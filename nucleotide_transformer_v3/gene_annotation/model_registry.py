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

"""Named model configurations: weights that carry their own window geometry.

A checkpoint is not enough to annotate a genome with. It has to be paired with
the window length it was trained (or length-extended) for and with the crop the
head is scored at, and getting either wrong does not raise: it produces a
plausible-looking GFF3 that is simply wrong. Before this module every caller
retyped ``sequence_length=1048576 keep_target_center_fraction=0.75`` beside an
opaque ``s3://`` path, which made that silent failure one typo away.

A model config is that pairing, written down once and referred to by name:

    --model_cfg mammals

The layout follows `Tiberius <https://github.com/Gaius-Augustus/Tiberius>`_,
whose ``model_cfg/<name>.yaml`` is the interface this community already knows.
Files live in :data:`MODEL_CFG_DIR` and are addressed by bare basename, so
adding a clade means adding a YAML file and nothing else.

The geometry in a config is *not* recoverable from the checkpoint, which is the
whole reason it is written down here. The mammal checkpoint records
``keep_target_center_fraction: 1.0`` as its trained crop, while the validated
annotation geometry is ``0.75``, which the crop sweep measured rather than
assumed (see :data:`FIXED_KEEP_TARGET_CENTER_FRACTION`). Nothing here can reach
the trained value: :func:`apply_model_config` pins that constant on the
no-``model_cfg`` path too, and the export already substitutes ``0.75`` into its
config, keeping the trained crop under
``provenance.trained_keep_target_center_fraction``.

Which half of the geometry travels per entry, precisely, because the two are not
the same and saying "the geometry travels with the weights" overstated it:

* ``sequence_length`` **does**. It is 1,048,576 for all three registered
  clades, and it is the field an entry exists to carry: it genuinely varies per
  checkpoint, and a model annotated at the wrong window produces a worse GFF3
  rather than an error.
* ``keep_target_center_fraction`` **does not**. Production annotates at
  :data:`FIXED_KEEP_TARGET_CENTER_FRACTION` and only at that value, so the field
  is a **declaration cross-checked against that constant**, not a choice: an
  entry naming anything else is rejected by :func:`apply_model_config` rather
  than honoured. It is kept as a field rather than dropped so that ``--show_cfg``
  can still answer "what crop will this run at" without a reader having to go
  find a constant in the source, and so that editing a YAML to a different crop
  fails loudly instead of being silently ignored.
"""

import dataclasses
import os
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import yaml

from nucleotide_transformer_v3.gene_annotation.config import ProductionAnnotationConfig
from nucleotide_transformer_v3.logging import get_logger

logger = get_logger()

MODEL_CFG_DIR = Path(__file__).parent / "model_cfg"

# Set this to a directory of per-clade subdirectories to load weights from
# disk instead of the Hub. Optional: the shipped entries work without it.
_WEIGHTS_DIR_ENV = "NTV3_ANNOTATE_WEIGHTS_DIR"

# ${NAME} or $NAME, which is what os.path.expandvars honours.
_VARIABLE_PATTERN = re.compile(r"\$\{(\w+)\}|\$(\w+)")


class _Loader(yaml.SafeLoader):
    """SafeLoader without YAML 1.1's date and yes/no/on/off coercions.

    Matches the OmegaConf loader these files were written against. Without it
    an unquoted ``verified_on: 2026-09-22`` parses as a ``datetime.date`` and
    ``softmasking: yes`` as a bool, both of which the field-type check below
    rejects with a message about the wrong thing.

    ``yaml_implicit_resolvers`` is a class attribute that subclasses share by
    reference, so it is deep-copied before being edited. Editing it in place
    reaches ``yaml.safe_load`` and every other SafeLoader subclass in the
    process, OmegaConf's included.
    """

    yaml_implicit_resolvers = {
        first: list(resolvers)
        for first, resolvers in yaml.SafeLoader.yaml_implicit_resolvers.items()
    }


for _first, _resolvers in _Loader.yaml_implicit_resolvers.items():
    _Loader.yaml_implicit_resolvers[_first] = [
        (tag, regexp)
        for tag, regexp in _resolvers
        if tag not in ("tag:yaml.org,2002:timestamp", "tag:yaml.org,2002:bool")
    ]

# true/false stay boolean, as OmegaConf has them; yes/no/on/off do not.
_Loader.add_implicit_resolver(
    "tag:yaml.org,2002:bool",
    re.compile(r"^(?:true|True|TRUE|false|False|FALSE)$"),
    list("tTfF"),
)


# Fields every config file must set, and the type each is read back as. Declared
# rather than inferred so that a typo becomes an error naming the offending key
# instead of a field silently keeping its default.
_REQUIRED_FIELDS: dict[str, type] = {
    "description": str,
    "model_name": str,
    "sequence_length": int,
    "keep_target_center_fraction": float,
    "softmasking": bool,
    "species_conditioned": bool,
    "verified": bool,
}

# An entry names its weights ONE of two ways, and exactly one is required:
# Where the converted weights live. One field, because there is one way in:
# for a raw checkpoint that has not been (or will not be) released.
_WEIGHTS_FIELDS = ("weights_uri",)

_OPTIONAL_FIELDS: dict[str, type] = {
    "weights_uri": str,
    "verified_on": str,
    "species_vocabulary": list,
    "notes": str,
}


@dataclass(frozen=True)
class ModelConfig:
    """One named checkpoint and the geometry it must be run at.

    Attributes:
        name: The file's basename, which is what ``--model_cfg`` takes.
        description: One human-readable line naming the clade and the training
            run, shown by ``--list_cfg``.
        model_name: The released model these weights are an export of. Checked
            at load against the export's own ``provenance.source_model_name``,
            because ``weights_uri`` names a directory the user chose and
            nothing else distinguishes one clade's export from another's.
        weights_uri: A released model name, downloaded from the Hub, or a
            directory of converted artifacts. Loaded by
            :func:`~nucleotide_transformer_v3.pretrained.load_annotation_model`.
            This is how every shipped clade is addressed. Exactly one of this
            and ``weights_uri`` is set.
        weights_uri: Where a raw, unreleased checkpoint lives, for a clade not
            worth releasing yet. A URI rather than a path, so an ``s3://`` key
            and a ``gs://`` or ``https://`` URL are the same field: the resolver
            passes it straight through to ``weights_dir`` and the loader
            accepts any URI ``epath`` understands.
        sequence_length: Window length in bp fed to the model.
        keep_target_center_fraction: Fraction of each window the head is scored
            at. **Not a per-entry choice**: production runs only at
            :data:`FIXED_KEEP_TARGET_CENTER_FRACTION`, and an entry declaring
            anything else is rejected rather than honoured, so this is a
            declaration that must agree with that constant. Present so
            ``--show_cfg`` reports the crop and so a YAML edit fails loudly. It
            is not in general the checkpoint's own trained crop.
        softmasking: Whether the checkpoint expects a soft-masked vocabulary,
            i.e. whether lowercase repeat-masked bases carry their own token
            ids. Informational: the loader re-derives this from the checkpoint
            and cross-checks it against the model's ``alphabet_size``, so a
            wrong value here cannot corrupt a run. It tells a user whether to
            stage a soft-masked assembly, which nothing else does.
        species_conditioned: Whether the checkpoint conditions on species, and
            so whether ``species`` (or ``mask_species_condition``) is required.
        verified: Whether the ``weights_uri`` was confirmed to exist and to be
            an annotation checkpoint. False entries load with a loud warning.
        verified_on: ISO date the check was made, empty if never.
        species_vocabulary: The species names the checkpoint's own tokenizer
            knows, when they were recovered at registration time. Present so
            that ``--show_cfg`` can answer "what may I pass to ``--species``"
            without bucket access, which is the one question whose wrong answer
            is silent: an unknown name encodes to ``<unk>`` rather than failing.
            The checkpoint remains authoritative; this list is read back from
            the released model at registration time.
        notes: Anything a user should know before trusting the output.
    """

    name: str
    description: str
    sequence_length: int
    keep_target_center_fraction: float
    softmasking: bool
    species_conditioned: bool
    verified: bool
    model_name: str = ""
    weights_uri: str = ""
    verified_on: str = ""
    species_vocabulary: tuple[str, ...] = ()
    notes: str = ""


def list_model_config_names() -> list[str]:
    """Every registered configuration name, sorted.

    Returns:
        The basenames of the YAML files in :data:`MODEL_CFG_DIR`.
    """
    if not MODEL_CFG_DIR.is_dir():
        return []
    return sorted(path.stem for path in MODEL_CFG_DIR.glob("*.yaml"))


def _known_names_hint() -> str:
    """The 'known configurations' clause shared by every unknown-name error."""
    names = list_model_config_names()
    if not names:
        return f"No configurations are registered in {MODEL_CFG_DIR}."
    return f"Known configurations: {', '.join(names)}."


def load_model_config(name: str) -> ModelConfig:
    """Read one named configuration and validate it.

    Args:
        name: A bare basename, as ``--model_cfg`` takes it.

    Returns:
        The parsed configuration.

    Raises:
        ValueError: If the name is not a bare basename, is not registered, or
            names a file whose contents are not a valid configuration.
    """
    if not name:
        raise ValueError("model_cfg is empty.")
    if Path(name).name != name or name in {".", ".."}:
        raise ValueError(
            f"model_cfg={name!r} must be a bare configuration name, not a path. "
            f"{_known_names_hint()}"
        )

    path = MODEL_CFG_DIR / f"{name}.yaml"
    if not path.is_file():
        raise ValueError(f"Unknown model_cfg {name!r}. {_known_names_hint()}")

    loaded = yaml.load(path.read_text(), Loader=_Loader)
    if not isinstance(loaded, dict):
        raise ValueError(f"{path} must contain a YAML mapping, got {type(loaded)}.")
    return _build_entry(name, loaded, path)


def resolve_weights_uri(entry_name: str, uri: str, path: Path) -> str:
    """Where this entry's weights actually are, for this machine.

    The shipped entries name a released model, so ``pip install`` and one
    command go from FASTA to GFF3 with nothing staged: the loader downloads
    ``InstaDeepAI/<name>/jax_model/`` from the Hub and caches it. The repos are
    gated, so that needs a Hugging Face login with approved access.

    ``NTV3_ANNOTATE_WEIGHTS_DIR`` overrides it with a local directory holding
    one subdirectory per clade, named as the entry is. That is the route for an
    air-gapped host, a shared mirror, or weights you converted yourself; it
    takes precedence because someone who sets it means it.

    ``~`` and ``$VAR`` are still expanded in a hand-written ``weights_uri``,
    with the variable name read from the ORIGINAL string rather than recovered
    from the expanded one: recovering it produced "``${$MISSING_VAR/x}`` is not
    set" for a brace-less ``$VAR/x``, and let a bare ``$VAR`` through as a
    literal. An empty expansion counts as unset, since
    ``NTV3_ANNOTATE_WEIGHTS_DIR=""`` used to yield ``/mammals``.

    Args:
        entry_name: The registry entry's name, which is the subdirectory looked
            for under the override.
        uri: The entry's declared ``weights_uri``.
        path: The entry's YAML, named in any error.

    Returns:
        A released model name, or a local directory.

    Raises:
        ValueError: If a variable in ``uri`` is unset or empty.
    """
    override = os.environ.get(_WEIGHTS_DIR_ENV, "").strip()
    if override:
        return str(Path(os.path.expanduser(override)) / entry_name)

    for braced, bare in _VARIABLE_PATTERN.findall(uri):
        name = braced or bare
        if os.environ.get(name):
            continue
        state = "empty" if name in os.environ else "not set"
        raise ValueError(
            f"{path}: weights_uri is {uri!r} and ${name} is {state} in the "
            f"environment. Point {name} at the directory holding the converted "
            f"weights, e.g. export {name}=/data/ntv3_annotate"
        )
    return os.path.expandvars(os.path.expanduser(uri))


def _build_entry(name: str, raw: dict[Any, Any], path: Path) -> ModelConfig:
    """Validate one loaded mapping and turn it into a :class:`ModelConfig`.

    Args:
        name: The configuration name, i.e. the file's basename.
        raw: The mapping read from the YAML file.
        path: The file it came from, quoted in every error.

    Returns:
        The validated configuration.

    Raises:
        ValueError: If a required key is missing, an unknown key is present, a
            value has the wrong type, or the geometry is out of range.
    """
    keys = {str(key) for key in raw}
    missing = sorted(set(_REQUIRED_FIELDS) - keys)
    unknown = sorted(keys - set(_REQUIRED_FIELDS) - set(_OPTIONAL_FIELDS))
    if missing or unknown:
        # Reported together because the common cause is one misspelling, which
        # shows up as both at once: naming only the missing key would leave the
        # reader hunting for a value that is sitting in the file under a typo.
        problems = []
        if missing:
            problems.append(f"missing required key(s): {', '.join(missing)}")
        if unknown:
            problems.append(f"unknown key(s): {', '.join(unknown)}")
        allowed = ", ".join(sorted(set(_REQUIRED_FIELDS) | set(_OPTIONAL_FIELDS)))
        raise ValueError(
            f"{path} has {'; '.join(problems)}. Allowed keys: {allowed}. A "
            "misspelled key would otherwise silently drop the value it was "
            "meant to set."
        )

    values: dict[str, Any] = {}
    for field, expected in (_REQUIRED_FIELDS | _OPTIONAL_FIELDS).items():
        if field not in raw:
            continue
        value = raw[field]
        # bool is a subclass of int, so an int field must reject True explicitly.
        wrong_bool = expected is not bool and isinstance(value, bool)
        if expected is float and isinstance(value, int) and not isinstance(value, bool):
            value = float(value)
        if wrong_bool or not isinstance(value, expected):
            raise ValueError(
                f"{path}: {field} must be {expected.__name__}, got "
                f"{type(value).__name__} ({value!r})."
            )
        values[field] = value

    if values["sequence_length"] <= 0:
        raise ValueError(
            f"{path}: sequence_length must be positive, got "
            f"{values['sequence_length']}."
        )
    fraction = values["keep_target_center_fraction"]
    if not 0.0 < fraction <= 1.0:
        raise ValueError(
            f"{path}: keep_target_center_fraction must be in (0, 1], got {fraction}."
        )
    named = [field for field in _WEIGHTS_FIELDS if values.get(field)]
    if len(named) != 1:
        raise ValueError(
            f"{path} sets {named or 'neither'} of {list(_WEIGHTS_FIELDS)}; "
            "weights_uri is required. It points at a directory of converted "
            "model artifacts."
        )

    vocabulary = values.pop("species_vocabulary", [])
    if any(not isinstance(item, str) for item in vocabulary):
        raise ValueError(f"{path}: species_vocabulary must be a list of strings.")

    return ModelConfig(name=name, species_vocabulary=tuple(vocabulary), **values)


# The crop the pipeline annotates at, fixed. Not a default: the only value a
# production run may use.
#
# It is fixed rather than exposed because the parameter is easy to get wrong in
# a way that produces a worse annotation rather than an error. Above a
# checkpoint's trained crop the head is asked for window-edge positions it was
# never supervised on, which collapses accuracy (a model trained at 0.375 and
# evaluated at 1.0 measured Gene-F1 7.9). Below it, quality holds but the
# forward cost rises as 1/crop for nothing.
#
# 0.75 IS MEASURED. It was 0.375, which was inherited rather than chosen, and a
# whole-genome sweep on cow and beluga (1 Mb windows, one H100-80GB, our own
# evaluator) put accuracy flat across the range and falling only at 1.0:
#
#   crop    cow s/Mb  cow exon/gene   beluga s/Mb  beluga exon/gene
#   0.375   1.908     90.52 / 57.54   1.891        89.26 / 53.14
#   0.625   1.197     90.54 / 57.56   1.235        89.26 / 53.15
#   0.75    1.019     90.58 / 57.53   1.042        89.28 / 53.16
#   1.0     0.795     90.24 / 57.44   -            -
#
# 1.87x faster on cow, and the accuracy deltas are the size of run-to-run
# variation: 0.375 -> 0.75 moves Exon +0.06 and Gene -0.01, against 0.01 on both
# when the 0.625 arm was re-run on a second node. These are single runs per crop
# - the median-of-three in that harness applies to its one-chromosome phase, not
# to these - so that repeat is the whole basis for calling the deltas noise.
#
# The long-gene advantage that justifies a 1 Mb context survives at the
# resolution available: >200kb gene recall 41.47% against 41.82% at 0.375
# (paired McNemar p=0.50) against Tiberius' 32.01%. The minimum detectable
# effect there is 1.11 points (n=1,131, alpha=0.05, power=0.8), so the null
# excludes a loss of 1.11 points or worse and says nothing about smaller ones.
#
# 1.0 is excluded on measurement, on the POOLED two-genome tests rather than on
# any single raw p: >200kb at p=0.0002 and 20+ exon at p=0.0028, both under the
# 0.005 Bonferroni threshold for ten buckets per crop. Its raw cow p of 0.038 is
# one of thirty tests and Holm-adjusts to 0.114, so it establishes nothing alone.
#
# All three registered models are trained at 1.0, so 0.75 is a narrowing for
# each and legal. A checkpoint trained below 0.75 could not be registered at
# this crop, which is what ``_reject_registry_conditioning_drift``'s sibling
# checks catch at load rather than at annotate time.
#
# Every published number comes from this value. Making it settable invites a
# user to produce numbers that cannot be compared with ours, from a model that
# looks like it is working. See benchmarks/annotation_timing/README.md.
FIXED_KEEP_TARGET_CENTER_FRACTION = 0.75


def _reject_contradiction(
    *,
    name: str,
    field: str,
    given: Any,
    registered: Any,
    unset: Any,
) -> None:
    """Fail if the caller set ``field`` to something the registry disagrees with.

    Silently letting the explicit value win would put back exactly the failure
    this module exists to remove, and silently letting the registry win would
    ignore what the caller asked for. Neither is safe, so the disagreement is
    the error.

    Args:
        name: The configuration name, quoted in the message.
        field: The config field in question.
        given: What the caller passed.
        registered: What the named configuration says.
        unset: The value meaning "the caller did not pass this".

    Raises:
        ValueError: If ``given`` is set and differs from ``registered``.
    """
    if given == unset or given == registered:
        return
    raise ValueError(
        f"model_cfg={name!r} sets {field}={registered!r}, but {field}={given!r} "
        "was also passed. The registered value is the geometry this checkpoint "
        "was validated at, and a mismatched one produces a wrong annotation "
        "instead of an error, so the two cannot both be honoured. Drop "
        f"{field}. To annotate at a geometry this pipeline does not offer, "
        "drop model_cfg and pass weights_dir with an explicit "
        f"sequence_length; the crop stays {FIXED_KEEP_TARGET_CENTER_FRACTION} "
        "on every path, so that is not unlocked either."
    )


def _reject_a_chosen_crop(requested: float | None) -> None:
    """Refuse any crop other than the fixed one.

    Args:
        requested: What the caller asked for, or None for "did not ask".

    Raises:
        ValueError: If a crop other than the fixed value was requested.
    """
    if requested is None or requested == FIXED_KEEP_TARGET_CENTER_FRACTION:
        return
    raise ValueError(
        f"keep_target_center_fraction={requested} was requested, but production "
        f"annotates at {FIXED_KEEP_TARGET_CENTER_FRACTION} and only at that. "
        "Every published number comes from that crop, a larger one asks the head "
        "for window-edge positions it was never supervised on, and a smaller one "
        "costs forward compute for no gain. Drop the setting."
    )


def _reject_a_length_the_head_crops_differently(sequence_length: int) -> None:
    """Refuse a window where the tiler and the head disagree about the centre.

    Two expressions compute the same crop and are not equal in general:

    * ``compute_window_geometry`` (``gene_annotation/annotate_utils.py``) takes
      ``(L - int(L * k)) // 2``, and the tiler places every window with it.
    * ``crop_sequence_center`` (``nucleotide_transformer_v3.layers``) takes
      ``int((1 - k) * L) // 2``, and the head crops its output with it.

    At ``L=1048576, k=0.75`` both give 131,072, which is why every registered
    entry is unaffected. At ``L=999983`` the head keeps ``[124997, 874986)`` and
    the tiler assumes ``[124998, 874985)``: one bp of offset, two of length.
    Measured, they disagree for **375,000 of the 1,000,000** lengths in
    ``[100_000, 1_100_000)`` at ``k=0.75``, so this is the common case rather
    than a corner.

    The consequence depends on the decoder, and neither outcome is acceptable:
    ``hmm_constrained`` fails on a broadcast between the base codes and the
    states, while ``hmm`` and ``argmax`` shift every coordinate on the contig by
    a base and say nothing. This branch is the module's one surface that accepts
    an arbitrary window, and its whole purpose is that a wrong geometry must
    error rather than produce a plausible GFF3, so it rejects rather than
    rounding one way and hoping.

    The arithmetic is inherited from ``annotate/`` and is not changed here:
    reconciling the two expressions would change what every existing checkpoint
    predicts, which is not this pipeline's call to make.

    Args:
        sequence_length: The window the caller asked for.

    Raises:
        ValueError: When the two expressions disagree at this length.
    """
    from nucleotide_transformer_v3.gene_annotation.annotate_utils import (
        compute_window_geometry,
    )

    crop = FIXED_KEEP_TARGET_CENTER_FRACTION
    tiler_offset = compute_window_geometry(sequence_length, crop)
    head_offset = int((1.0 - crop) * sequence_length) // 2
    if tiler_offset == head_offset:
        return
    raise ValueError(
        f"sequence_length={sequence_length:,} cannot be annotated at crop {crop}: "
        f"the tiler places windows with a {tiler_offset:,} bp margin while the "
        f"model head crops {head_offset:,} bp, so the decoded block and the "
        f"predicted block differ by {abs(tiler_offset - head_offset) * 2} bp. "
        "That is a broadcast error under hmm_constrained and a silent "
        "coordinate shift under hmm/argmax. Use a length where "
        "(L - int(L * k)) // 2 == int((1 - k) * L) // 2 "
        f"(for example {_nearest_consistent_length(sequence_length):,}), or "
        "pass model_cfg=<name> to take a validated geometry from the registry."
    )


def _nearest_consistent_length(sequence_length: int) -> int:
    """Closest length at which the tiler and the head agree, for the error above.

    Searched outward rather than derived: the two expressions differ only in
    where ``int()`` truncates, so the consistent lengths are not a simple
    residue class and a wrong closed form would be worse than a short loop.
    """
    from nucleotide_transformer_v3.gene_annotation.annotate_utils import (
        compute_window_geometry,
    )

    crop = FIXED_KEEP_TARGET_CENTER_FRACTION
    for delta in range(1, 64):
        for candidate in (sequence_length - delta, sequence_length + delta):
            if (
                candidate > 0
                and compute_window_geometry(candidate, crop)
                == int((1.0 - crop) * candidate) // 2
            ):
                return candidate
    return sequence_length


def apply_model_config(
    config: ProductionAnnotationConfig,
) -> ProductionAnnotationConfig:
    """Resolve ``model_cfg`` into concrete weights and window geometry.

    The single entry point both surfaces go through: the Hydra script and the
    CLI each call this before anything else, so a name resolves the same way
    whichever one was used.

    With no ``model_cfg`` the config is returned untouched, except that a
    missing ``sequence_length`` is now an error rather than a default. There is
    no window length that is right for an unknown checkpoint, so a default one
    could only ever be a silent wrong answer.

    Args:
        config: The composed configuration, before any geometry is needed.

    Returns:
        ``config`` itself when no ``model_cfg`` was named, or a copy with
        ``weights_dir``, ``sequence_length`` and
        ``keep_target_center_fraction`` taken from the registry.

    Raises:
        ValueError: If the name is unknown, if an explicitly passed checkpoint
            or geometry contradicts the registered one, or if neither a
            ``model_cfg`` nor a full explicit geometry was given.
    """
    if not config.model_cfg:
        if not config.weights_dir:
            raise ValueError(
                "Pass model_cfg=<name> to annotate with a registered model, or "
                "weights_dir=<directory> with an explicit sequence_length. "
                f"{_known_names_hint()}"
            )
        if not config.sequence_length:
            raise ValueError(
                "sequence_length is required without a model_cfg: it "
                "must match the window the checkpoint was trained (or "
                "length-extended) for, and a mismatched window produces a wrong "
                "annotation rather than an error. Pass model_cfg=<name> instead "
                f"to take it from the registry. {_known_names_hint()}"
            )
        _reject_a_chosen_crop(config.keep_target_center_fraction)
        _reject_a_length_the_head_crops_differently(config.sequence_length)
        return dataclasses.replace(
            config,
            keep_target_center_fraction=FIXED_KEEP_TARGET_CENTER_FRACTION,
        )

    entry = load_model_config(config.model_cfg)
    if not entry.verified:
        logger.warning(
            f"model_cfg={entry.name!r} is marked unverified: its weights_uri "
            f"({entry.weights_uri}) has not been confirmed to exist and to be an "
            "annotation checkpoint. Treat its output accordingly."
        )

    _reject_contradiction(
        name=entry.name,
        field="weights_dir",
        given=config.weights_dir,
        # Resolved, because that is what the caller would have passed.
        registered=resolve_weights_uri(
            entry.name, entry.weights_uri, MODEL_CFG_DIR / f"{entry.name}.yaml"
        ),
        unset="",
    )
    _reject_contradiction(
        name=entry.name,
        field="sequence_length",
        given=config.sequence_length,
        registered=entry.sequence_length,
        unset=0,
    )
    # On THIS branch too, not only the explicit-checkpoint one. A YAML is the
    # other surface an arbitrary window arrives through, and the likelier of the
    # two: an entry at sequence_length=999983 used to be accepted happily while
    # the tiler placed windows at offset 124,998 and the head cropped at
    # 124,997. Under `hmm` or `argmax` every coordinate on the contig then
    # shifts by a base silently, because the slice stays in bounds either way.
    # Latent for the shipped clades only because all are 1,048,576, and powers
    # of two agree at this crop.
    _reject_a_length_the_head_crops_differently(entry.sequence_length)

    if entry.species_conditioned and not config.species:
        if not config.mask_species_condition:
            raise ValueError(
                f"model_cfg={entry.name!r} conditions on species. Pass "
                "species=<common name>, or mask_species_condition=true to "
                "annotate an organism this checkpoint never saw."
                + (
                    f" Known species: {', '.join(entry.species_vocabulary)}."
                    if entry.species_vocabulary
                    else ""
                )
            )
    elif not entry.species_conditioned and config.species:
        logger.warning(
            f"model_cfg={entry.name!r} is not species-conditioned, so "
            f"species={config.species!r} has no effect."
        )

    logger.info(
        f"model_cfg={entry.name!r}: {entry.weights_uri} at "
        f"{entry.sequence_length:,} bp, keep fraction "
        f"{entry.keep_target_center_fraction}."
    )
    _reject_a_chosen_crop(config.keep_target_center_fraction)
    if entry.keep_target_center_fraction != FIXED_KEEP_TARGET_CENTER_FRACTION:
        raise ValueError(
            f"model_cfg={entry.name!r} declares "
            f"keep_target_center_fraction={entry.keep_target_center_fraction}, but "
            f"production annotates only at {FIXED_KEEP_TARGET_CENTER_FRACTION}. The "
            "entry is registering a geometry the pipeline will not run, so it "
            "would be silently ignored. Fix the entry."
        )
    return dataclasses.replace(
        config,
        model_name=entry.model_name,
        weights_dir=resolve_weights_uri(
            entry.name, entry.weights_uri, MODEL_CFG_DIR / f"{entry.name}.yaml"
        ),
        sequence_length=entry.sequence_length,
        keep_target_center_fraction=FIXED_KEEP_TARGET_CENTER_FRACTION,
    )


def format_model_config_list() -> str:
    """The ``--list_cfg`` table: every registered name with its description.

    Returns:
        One line per configuration, plus a trailing hint, or a message saying
        the registry is empty.
    """
    names = list_model_config_names()
    if not names:
        return f"No model configurations are registered in {MODEL_CFG_DIR}."

    entries = [load_model_config(name) for name in names]
    width = max(len(entry.name) for entry in entries)
    lines = [
        f"{entry.name:<{width}}  {entry.description}"
        + ("" if entry.verified else "  [UNVERIFIED]")
        for entry in entries
    ]
    lines.append("")
    lines.append("Show one in full with --show_cfg <name>.")
    return "\n".join(lines)


def format_model_config(name: str) -> str:
    """The ``--show_cfg`` detail view for one configuration.

    Args:
        name: The configuration to show.

    Returns:
        Every field of the configuration, one per line.

    Raises:
        ValueError: If the name is not registered.
    """
    entry = load_model_config(name)
    rows: list[tuple[str, str]] = [
        ("name", entry.name),
        ("description", entry.description),
        ("weights_uri", entry.weights_uri),
        ("sequence_length", f"{entry.sequence_length:,}"),
        (
            "keep_target_center_fraction",
            f"{entry.keep_target_center_fraction} (fixed for all models)",
        ),
        ("softmasking", str(entry.softmasking)),
        ("species_conditioned", str(entry.species_conditioned)),
        (
            "verified",
            str(entry.verified)
            + (f" ({entry.verified_on})" if entry.verified_on else ""),
        ),
    ]
    if entry.species_vocabulary:
        rows.append(("species_vocabulary", ", ".join(entry.species_vocabulary)))
    if entry.notes:
        rows.append(("notes", entry.notes))

    width = max(len(key) for key, _ in rows)
    return "\n".join(f"{key:<{width}}  {value}" for key, value in rows)
