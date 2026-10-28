# Copyright 2022 InstaDeep Ltd
#
# Licensed under the Creative Commons BY-NC-SA 4.0 License (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      https://creativecommons.org/licenses/by-nc-sa/4.0/
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Check released Flax artifacts against the fixture exported beside them.

The unit tests run on tiny random models, so they cannot catch a bad weight
conversion or a wrong vocabulary. This loads the real artifacts through the
package's loaders and replays ``roundtrip_fixture.npz``, the inputs and outputs
the exporting model produced. It needs the weights, so it is run by hand::

    python -m nucleotide_transformer_v3.check_release pretrained NTv3_pretrained
    python -m nucleotide_transformer_v3.check_release tracks NTv3_tracks
    python -m nucleotide_transformer_v3.check_release annotation NTv3_annot_mammals

The model argument is a released name or a directory holding the fixture.
"""

import argparse

import jax
import jax.numpy as jnp
import numpy as np
from flax import nnx

from nucleotide_transformer_v3.pretrained import (
    VARIABLES,
    count_parameters,
    load_annotation_model,
    load_pretrained_model,
    load_tracks_model,
    resolve_model_dir,
)
from nucleotide_transformer_v3.tokenizers import BaseTokenizer
from nucleotide_transformer_v3.with_tracks_heads import (
    DiscreteConditionedNTv3WithTracksHeads,
)

FIXTURE_FILE = "roundtrip_fixture.npz"

# The output the export compared, per family: categorical logits for the
# pretrained LM head and the annotation head, regression-valued tracks for the
# bigwig head.
COMPARED_OUTPUT = {
    "pretrained": "logits",
    "tracks": "bigwig_tracks_logits",
    "annotation": "bed_tracks_logits",
}
# Families whose output is decoded by argmax, so every position must agree.
CATEGORICAL = frozenset({"pretrained", "annotation"})

# Smallest change in the tracks, relative to their scale, that switching the
# species from human to mouse must produce. The released model measures 0.72;
# a conditioning that never reaches the head measures exactly 0.
MIN_SPECIES_EFFECT = 1e-3


def check(family: str, model: str, tolerance: float) -> dict[str, float]:
    """Replay the fixture and return the measured agreement.

    Raises:
        ValueError: If the outputs differ by more than ``tolerance`` of their
            scale, if a categorical output's argmax disagrees anywhere, or if
            the species conditioning does not reach the tracks head.
    """
    model_dir = resolve_model_dir(model)
    fixture = np.load(model_dir / FIXTURE_FILE)
    tokens = jnp.asarray(fixture["tokens"])
    key = COMPARED_OUTPUT[family]

    if family == "pretrained":
        loaded, _ = load_pretrained_model(model_dir, use_bf16_compute=False)
        got = np.asarray(loaded(tokens)[key])
        extra: dict[str, float] = {}
    elif family == "annotation":
        loaded, _ = load_annotation_model(model_dir, use_bf16_compute=False)
        # The released annotation models carry no conditions.
        got = np.asarray(loaded(tokens, [])[key])
        extra = {}
    else:
        loaded, _, (species,) = load_tracks_model(model_dir, use_bf16_compute=False)
        rows = tokens.shape[0]
        # The export ran with condition id 0 and every condition masked out,
        # so the fixture pins the unconditioned path only.
        got = np.asarray(
            loaded(
                tokens,
                [jnp.zeros((rows,), jnp.int32)],
                conditions_masks=[jnp.zeros((rows,), bool)],
            )[key]
        )
        extra = {"species_effect": _species_effect(loaded, species, tokens[:1])}
        # The fixture cannot see this, so it is asserted here or nowhere.
        if extra["species_effect"] < MIN_SPECIES_EFFECT:
            raise ValueError(
                "Species conditioning does not reach the bigwig head: human and "
                f"mouse differ by {extra['species_effect']:.3e} of the scale."
            )

    want = fixture[key]
    if want.shape != got.shape:
        raise ValueError(f"Shape mismatch: fixture {want.shape}, loaded {got.shape}.")
    delta = float(np.abs(want - got).max())
    scale = float(np.abs(want).max()) or 1.0
    result = {
        "leaves": float(len(jax.tree_util.tree_leaves(nnx.state(loaded, VARIABLES)))),
        "parameters": float(count_parameters(loaded)),
        "bitwise": float(np.array_equal(want, got)),
        "max_delta": delta,
        "relative": delta / scale,
        **extra,
    }
    if family in CATEGORICAL:
        result["argmax_agreement"] = float((want.argmax(-1) == got.argmax(-1)).mean())
        if result["argmax_agreement"] < 1.0:
            raise ValueError(f"argmax agreement {result['argmax_agreement']:.6%}")
    # Relative, because the fixture's probe is random sequence with a padded
    # tail: far out of distribution, it drives the pretrained logits to ~1e6,
    # where float32 noise alone exceeds any absolute bound worth setting.
    if delta / scale > tolerance:
        raise ValueError(
            f"max|delta| {delta:.3e} is {delta / scale:.3e} of the output scale, "
            f"above the tolerance {tolerance:.1e}."
        )
    return result


def _species_effect(
    model: DiscreteConditionedNTv3WithTracksHeads,
    species: BaseTokenizer,
    tokens: jax.Array,
) -> float:
    """Max relative change in the tracks between human and mouse conditioning.

    Zero would mean the species never reaches the bigwig head. The fixture
    cannot show that, since it masks the conditioning out.
    """
    ids = jnp.asarray([species.token_to_id("human"), species.token_to_id("mouse")])
    outputs = model(
        jnp.tile(tokens, (2, 1)),
        [ids],
        track_indices=jnp.arange(min(8, model.config.number_bigwig_tracks)),
    )
    tracks = np.asarray(outputs["bigwig_tracks_logits"])
    return float(np.abs(tracks[0] - tracks[1]).max() / (np.abs(tracks).max() or 1.0))


def main(argv: list[str] | None = None) -> int:
    """Run the check and print its measurements. Returns an exit status."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("family", choices=sorted(COMPARED_OUTPUT))
    parser.add_argument("model", help="Directory or released model name.")
    parser.add_argument(
        "--tolerance",
        type=float,
        default=1e-4,
        help="Largest max|delta| allowed, relative to the largest |output|.",
    )
    args = parser.parse_args(argv)

    result = check(args.family, args.model, args.tolerance)
    for name, value in result.items():
        print(f"{name}: {int(value) if float(value).is_integer() else f'{value:.3e}'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
