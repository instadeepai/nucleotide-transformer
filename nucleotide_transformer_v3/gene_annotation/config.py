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

"""Config schema for the genome-annotation entry point."""

from dataclasses import dataclass, field


@dataclass
class ProductionAnnotationConfig:
    """Knobs for annotating one genome from FASTA to GFF3.

    Deliberately smaller than ``AnnotationConfig``: there is no MLflow run, no
    reference GFF, no metrics and no dataloader here, so none of the knobs that
    serve those exist. What remains is the model, the geometry, and the output.

    Attributes:
        fasta_paths: FASTA files and/or directories of per-contig FASTAs, read in
            the order given. Every record across every input is annotated.
        model_cfg: Name of a registered model configuration (see
            ``model_cfg/``), which supplies ``weights_dir`` together with
            the window geometry it must be run at. The recommended way in:
            geometry travelling with the weights is what stops a mismatched
            window from silently producing a wrong annotation. Empty means the
            checkpoint and geometry are given explicitly instead, and passing
            both a name and a value it disagrees with is an error rather than a
            silent override.
        model_name: The released model the weights are an export of, e.g.
            ``new_ntv3_700M_annotate_mammals_1mb``. Set from ``model_cfg``, and
            cross-checked against the export's own recorded provenance, so a
            directory holding the wrong clade is caught rather than annotated.
            It is also what the GFF3 provenance header records, since a
            directory path identifies nothing to whoever receives the file.
        weights_dir: Directory holding the converted model artifacts
            (``hyperparameters.json``, ``flax_model.msgpack`` and optionally
            ``tokenizer_config.json``). Set from ``model_cfg`` when one is
            named, or given directly for a model that is not registered.
        output_path: Annotation file to write; gzipped when it ends in ``.gz``.
            One file for the whole input, with gene ids unique across it. The
            format follows the extension, as Tiberius' does: ``.gff3`` / ``.gff``
            give GFF3 and ``.gtf`` gives GTF. Any other extension is rejected
            rather than guessed at.
        codingseq_path: Optional FASTA of every predicted gene's spliced CDS,
            one record per gene named by its GFF3 ``ID``. Tiberius'
            ``--codingseq``. Empty writes none.
        protseq_path: Optional FASTA of the protein each spliced CDS translates
            to under the standard genetic code, named the same way. Tiberius'
            ``--protseq``. Empty writes none.
        species: Organism common name for a species-conditioned checkpoint.
            Ignored (with a warning) by an unconditioned one.
        mask_species_condition: Drop the species conditioning instead of encoding
            ``species``. The explicit way to annotate an organism the checkpoint's
            species vocabulary does not know.
        decoder: Which decoder to run. ``hmm_constrained`` is the constrained
            Viterbi; ``hmm`` is unconstrained; ``argmax`` is per-position.
        sequence_length: Window length fed to the model (bp). Must match what the
            checkpoint was trained (or length-extended) for. 0 means unset, and
            is an error unless ``model_cfg`` supplies it: there is no window
            length that is right for an arbitrary checkpoint, so a default one
            could only ever be a silent wrong answer.
        keep_target_center_fraction: Not a setting. Production annotates at
            ``FIXED_KEEP_TARGET_CENTER_FRACTION`` (0.75) and nothing else, and
            asking for another value is an error rather than an override. The
            field remains because the geometry code downstream takes it as an
            argument; ``apply_model_config`` fills it in.
        backend: JAX backend ("gpu" or "cpu") selecting the device the whole
            annotation pass runs on: checkpoint load, model forward and decode.
            The first device of that backend is used. Empty accepts whatever
            backend JAX picks by default, which is what makes a CPU-only host
            work without configuration.
        batch_size: Windows per forward pass. Exactly that many; the output is
            the same at any batch size.
        use_bf16_compute: Rebuild the loaded model with bfloat16 compute dtypes.
            On by default: measured at ~1.6x the throughput of fp32 with no
            annotation-quality cost.
        apply_biological_filters: Drop genes whose spliced CDS has an internal
            in-frame stop codon. On by default, for parity with the documented
            recipe. The evaluator does *not* apply this filter, so it changes
            scored output and must be reported alongside any quoted metric.
            Needs per-base codes, so leaving it on makes the per-window
            ``base_codes`` computation unconditional.
        min_cds_length: Drop genes whose total CDS is shorter than this (0 off).
            Defaults to the recipe's 201. Idempotent with the evaluator's own
            ``min_pred_gene_cds_length=201`` (both count the same quantity: this
            counts non-intron positions in state space, the evaluator sums CDS
            feature lengths, and every non-intergenic non-intron class maps to
            type CDS), so it changes the shipped GFF3 rather than any reported
            metric. That is the point: nothing downstream of production would
            otherwise clean the artifact a user receives.
        decode_chunk_len: Viterbi block length for the blocked parallel decode.
        log_every_n_batches: Progress cadence within a contig.
        output_parts_dir: Where per-record GFF3 parts are staged before being
            merged into ``output_path``. Empty defaults to
            ``<output_path>.parts``. Each record is written as soon as it is
            finished, so a whole genome is never held in memory. Parts hold GFF3
            whatever the output format, so a GTF job resumes from a GFF3 job's
            parts and vice versa. The sequence outputs stage the same way, under
            ``<codingseq_path>.parts`` and ``<protseq_path>.parts``.
        resume: Reuse any record whose part file already exists and was produced
            from the same sequence under the same settings, so a killed
            whole-genome job restarts where it stopped rather than from the
            beginning.

            A part carries what it was made from (see ``PartIdentity``), and a
            part that disagrees with the record or the run in front of it is not
            finished work for this job: it is logged as stale and the record is
            annotated again. That covers an assembly edited under an unchanged
            record id and a run resumed with a different checkpoint, geometry,
            decoder, filters or species, either of which would otherwise merge
            features describing one input under a provenance header describing
            another. The output *format* is not part of that identity: parts
            always hold GFF3, so a GTF job still resumes from a GFF3 job's parts.

            The sequence check is a length plus a digest of a bounded sample, so
            it costs the same on a chromosome as on a scaffold. It sees every
            insertion, deletion or truncation, and any substitution inside the
            sampled blocks, which for a record up to 1 MiB is all of it. It does
            not see a length-preserving edit falling between two sampled blocks
            of a larger record.

            Gene ids are assigned from a counter threaded across records, and a
            resumed run seeds that counter past every id already staged in the
            parts directory, so the merged file cannot carry an ID twice. The
            numbering may therefore have gaps. Ids stay opaque labels, not
            identifiers whose stability we promise.
    """

    fasta_paths: list[str] = field(default_factory=list)
    model_cfg: str = ""
    model_name: str = ""
    weights_dir: str = ""
    output_path: str = ""
    codingseq_path: str = ""
    protseq_path: str = ""

    species: str | None = None
    mask_species_condition: bool = False

    decoder: str = "hmm_constrained"

    backend: str = ""

    sequence_length: int = 0
    keep_target_center_fraction: float | None = None
    batch_size: int = 1

    use_bf16_compute: bool = True

    apply_biological_filters: bool = True
    min_cds_length: int = 201
    decode_chunk_len: int = 512

    log_every_n_batches: int = 20

    output_parts_dir: str = ""
    resume: bool = False
