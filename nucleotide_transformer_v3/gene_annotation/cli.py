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

r"""A Tiberius-shaped command line over the production annotation pipeline.

`Tiberius <https://github.com/Gaius-Augustus/Tiberius>`_ is the interface this
community already uses, so the same three flags do the same three things here:

    ntv3-annotate --genome input.fasta --model_cfg mammals --out ntv3.gff3

Add ``--codingseq cds.fa --protseq prot.fa`` for the sequence outputs, and
``--resume`` to continue a killed run. ``--out`` chooses GFF3 or GTF by its
extension.

with ``--list_cfg`` and ``--show_cfg <name>`` to find out what may go after
``--model_cfg``.

This is a surface, not a second pipeline. Every flag is composed into the same
``ProductionAnnotationConfig`` ``annotate_genome.run`` takes, so the two entry
points cannot drift: there is one engine and one code path, reachable two ways.

Anything without a flag of its own is still reachable through ``--override``:

    ... --override decoder=hmm --override use_bf16_compute=false

Named flags are applied after ``--override``, so a named flag wins if both set
the same field.
"""

import argparse
import dataclasses
import logging
import sys
from types import UnionType
from typing import Any, Union, get_args, get_origin, get_type_hints

from nucleotide_transformer_v3.gene_annotation.config import ProductionAnnotationConfig
from nucleotide_transformer_v3.gene_annotation.model_registry import (
    format_model_config,
    format_model_config_list,
    list_model_config_names,
)


def build_parser() -> argparse.ArgumentParser:
    """The argument parser, with Tiberius' flag names where they line up.

    Returns:
        A parser whose namespace :func:`compose_config` understands.
    """
    parser = argparse.ArgumentParser(
        prog="annotate_genome",
        description=(
            "Annotate a genome FASTA with a trained NTv3 annotation head and "
            "write a GFF3."
        ),
        epilog=(
            "Any config field can be set with --override key=value; see "
            "ProductionAnnotationConfig for the full list."
        ),
    )
    parser.add_argument(
        "--genome",
        nargs="+",
        metavar="PATH",
        help=(
            "FASTA file(s), and/or directories of per-contig FASTAs, read in "
            "the order given. Every record across every input is annotated."
        ),
    )
    parser.add_argument(
        "--model_cfg",
        metavar="NAME",
        help=(
            "Registered model to annotate with. Brings its checkpoint and the "
            "window geometry it must be run at. See --list_cfg."
        ),
    )
    parser.add_argument("--out", metavar="PATH", help="GFF3 to write, gzipped if .gz.")
    parser.add_argument(
        "--species",
        metavar="NAME",
        help=(
            "Organism common name, for a species-conditioned model. No shipped "
            "model is one: all three clades ignore it, so leave it out. Kept "
            "for a future conditioned checkpoint, and warned about rather "
            "than silently accepted."
        ),
    )
    parser.add_argument(
        "--batch_size",
        type=int,
        metavar="N",
        help=(
            "Windows per forward pass, exactly that many. The output is "
            "invariant to it, so it trades only memory against speed. Leave it "
            "at 1 unless you have measured that more fits."
        ),
    )
    parser.add_argument(
        "--codingseq",
        metavar="PATH",
        help=(
            "Also write a FASTA of each gene's spliced CDS, gzipped if .gz. "
            "Tiberius' --codingseq."
        ),
    )
    parser.add_argument(
        "--protseq",
        metavar="PATH",
        help=(
            "Also write a FASTA of the protein each CDS translates to, gzipped "
            "if .gz. Tiberius' --protseq."
        ),
    )
    parser.add_argument(
        "--resume",
        action="store_true",
        default=None,
        help=(
            "Skip contigs already annotated in a previous run of the same "
            "input and settings, so a killed genome-scale job restarts where "
            "it stopped. A part produced from different sequence or different "
            "settings is re-annotated rather than reused."
        ),
    )
    parser.add_argument(
        "--seq_len",
        type=int,
        metavar="BP",
        help=(
            "Window length in bp. Supplied by --model_cfg; passing a value that "
            "disagrees with the registered one is an error, not an override."
        ),
    )
    parser.add_argument(
        "--override",
        "-O",
        action="append",
        default=[],
        metavar="KEY=VALUE",
        help="Set any config field: --override key=value. Repeatable.",
    )
    parser.add_argument(
        "--list_cfg",
        action="store_true",
        help="List the registered models and exit.",
    )
    parser.add_argument(
        "--show_cfg",
        metavar="NAME",
        help="Show one registered model in full and exit.",
    )
    return parser


_ANNOTATIONS: dict[str, Any] = get_type_hints(ProductionAnnotationConfig)


def _coerce(field: dataclasses.Field, raw: str) -> Any:
    """Turn an --override string into the field's declared type.

    The research pipeline got this from Hydra's override grammar, which typed
    and validated before the model loaded. Without it a mistyped value would
    survive and fail on the first batch instead, which on a genome-scale job is
    after the weights and the JIT compile have already been paid for.

    The type is read with ``typing.get_origin``/``get_args`` rather than by
    substring, so it stays right for whatever fields get added later: matching
    on text made ``list[str]`` win before ``int``, which would have coerced a
    future ``list[int]`` to a scalar and sent ``dict[str, int]`` down the int
    branch.
    """
    annotation = field.type
    if isinstance(annotation, str):
        annotation = _ANNOTATIONS.get(field.name, annotation)
    # Unwrap a union, and ONLY a union: get_args(list[str]) is (str,) too, so
    # unwrapping unconditionally turned every list field into its item type.
    optional = False
    target: Any = annotation
    if get_origin(annotation) in (Union, UnionType):
        members = [a for a in get_args(annotation) if a is not type(None)]
        optional = type(None) in get_args(annotation)
        target = members[0] if members else str

    # Hydra spelled "unset" as null, and every Optional field in this config
    # defaults to None, so there has to be a way back to it.
    if optional and raw.lower() in ("null", "none", "~"):
        return None

    if get_origin(target) is list:
        # Deliberately not split on commas: a path may contain one, and
        # `fasta_paths=/my,dir/a.fa` used to become ['/my', 'dir/a.fa'].
        # Repeat the override instead, which also needs no quoting.
        item_type = (get_args(target) or (str,))[0]
        return [_scalar(field.name, item_type, raw)]
    return _scalar(field.name, target, raw)


def _scalar(name: str, target: Any, raw: str) -> Any:
    """Parse one scalar of the declared type, or say what was expected."""
    if target is bool:
        lowered = raw.lower()
        if lowered not in ("true", "false", "1", "0"):
            raise ValueError(f"{name} expects true/false, got {raw!r}")
        return lowered in ("true", "1")
    if target in (int, float):
        try:
            return target(raw)
        except ValueError:
            article = "an" if target is int else "a"
            raise ValueError(
                f"{name} expects {article} {target.__name__}, got {raw!r}"
            ) from None
    return raw


def apply_overrides(
    config: ProductionAnnotationConfig, overrides: list[str]
) -> ProductionAnnotationConfig:
    """Apply ``key=value`` overrides, typed and validated against the dataclass.

    Args:
        config: The config to copy and update.
        overrides: Strings of the form ``field=value``.

    Returns:
        A new config with the overrides applied.

    Raises:
        ValueError: On a malformed override, an unknown field, or a value that
            does not parse as the field's type.
    """
    fields = {f.name: f for f in dataclasses.fields(config)}
    updates: dict[str, Any] = {}
    for override in overrides:
        key, separator, raw = override.partition("=")
        if not separator:
            raise ValueError(f"Override {override!r} is not of the form key=value.")
        key = key.strip()
        if key not in fields:
            raise ValueError(
                f"Unknown config field {key!r}. Known fields: "
                f"{', '.join(sorted(fields))}."
            )
        value = _coerce(fields[key], raw)
        if isinstance(value, list) and isinstance(updates.get(key), list):
            updates[key] = updates[key] + value
        else:
            updates[key] = value
    return dataclasses.replace(config, **updates)


def compose_config(args: argparse.Namespace) -> ProductionAnnotationConfig:
    """Build the config the flags describe.

    Overrides are applied first and named flags second, so a named flag wins
    when both set the same field -- the order the module docstring promises.

    Args:
        args: The parsed command line.

    Returns:
        The config, ready for ``run()``.
    """
    config = apply_overrides(ProductionAnnotationConfig(), list(args.override))
    return _apply_flags(config, args)


def _apply_flags(
    config: ProductionAnnotationConfig, args: argparse.Namespace
) -> ProductionAnnotationConfig:
    """Write the flags that were actually passed onto the config.

    A flag left off the command line is None and must not be written, or every
    unset flag would clobber the dataclass's own default with it.

    Args:
        config: The config to copy and update.
        args: The parsed command line.

    Returns:
        A new config with the passed flags applied.
    """
    flag_to_field: dict[str, str] = {
        "genome": "fasta_paths",
        "model_cfg": "model_cfg",
        "out": "output_path",
        "species": "species",
        "batch_size": "batch_size",
        "seq_len": "sequence_length",
        "codingseq": "codingseq_path",
        "protseq": "protseq_path",
        "resume": "resume",
    }
    updates: dict[str, Any] = {}
    for flag, field in flag_to_field.items():
        value: Any = getattr(args, flag)
        if value is None:
            continue
        updates[field] = list(value) if flag == "genome" else value
    return dataclasses.replace(config, **updates)


def main(argv: list[str] | None = None) -> int:
    """Run the CLI.

    Args:
        argv: Arguments to parse, defaulting to ``sys.argv[1:]``.

    Returns:
        A process exit status: 0 on success, 2 for a usage error.
    """
    parser = build_parser()
    args = parser.parse_args(argv)

    # The two query flags answer from the registry alone, so they stay fast and
    # work on a host with no accelerator and no bucket credentials.
    if args.list_cfg:
        print(format_model_config_list())
        return 0
    if args.show_cfg:
        try:
            print(format_model_config(args.show_cfg))
        except ValueError as error:
            parser.error(str(error))
        return 0

    if not args.genome:
        parser.error("--genome is required (or --list_cfg / --show_cfg).")
    if not args.out:
        parser.error("--out is required.")
    if not args.model_cfg and not any(
        override.strip().startswith("weights_dir=") for override in args.override
    ):
        parser.error(
            "--model_cfg is required, or --override weights_dir=<directory> "
            "with --seq_len for a model that is not registered. "
            f"Registered models: {', '.join(list_model_config_names()) or 'none'}."
        )

    try:
        config = compose_config(args)
    except ValueError as error:
        parser.error(str(error))

    # Nothing else installs a root handler, so a genome-scale run would
    # otherwise print nothing for its whole duration.
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")

    # Imported here rather than at module scope so that --list_cfg and
    # --show_cfg do not pay for JAX's import.
    from nucleotide_transformer_v3.gene_annotation.annotate_genome import run

    run(config)
    return 0


if __name__ == "__main__":
    sys.exit(main())
