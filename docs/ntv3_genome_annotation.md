# 🧬 NTv3 Genome Annotation

Annotate a genome with a trained NTv3 annotation head and get a GFF3 or a GTF,
optionally with the spliced-CDS and protein FASTAs. One input, one command.

```bash
ntv3-annotate --genome /data/cow/genome.fa --model_cfg mammals --out /data/cow.gff3.gz
```

---

## 📑 Table of Contents

1. [⚠️ Hardware you will need](#️-hardware-you-will-need)
2. [📦 Models](#-models)
3. [🚀 Usage](#-usage)
4. [🧬 Inputs](#-inputs)
5. [📄 Output formats](#-output-formats)
6. [🔁 Resume and provenance](#-resume-and-provenance)
7. [🔬 Post-processing filters](#-post-processing-filters)
8. [📏 Why the window and the crop are not settings](#-why-the-window-and-the-crop-are-not-settings)
9. [🐍 Python API](#-python-api)
10. [⚖️ Known limitations](#️-known-limitations)

---

## ⚠️ Hardware you will need

Read this first; it decides whether this tool is usable for you.

| | |
| --- | --- |
| **GPU memory** | one 1 Mb window peaks at **25.3 GB**; allow **~32 GB** |
| **Practical minimum** | **A100-40GB or larger. A 24 GB card will not run it.** |
| Throughput | ~1.0 s/Mb, so a 2.7 Gb genome is roughly 45 minutes |
| `--batch_size` | 1 by default; 2 is the most that fits on an 80 GB card |
| CPU | supported and far slower: ~188 s/window, ~3,400 windows, 7.5 days |

Measured, not estimated. Both totals are per-window rates over the ~3,400
windows a 2.7 Gb genome tiles into at the 786,432 bp stride (a 1,048,576 bp
window less the 131,072 bp context margin each side). If you do not have a
data-centre GPU, [Tiberius](https://github.com/Gaius-Augustus/Tiberius)
recommends 8 GB and is the tool that will run.

### Install the GPU extra, or you will silently get the CPU

```bash
pip install 'nucleotide_transformer[gpu]'
```

PyPI's `jaxlib` carries no CUDA. Install without the extra on a GPU machine
and **nothing fails** -- JAX prints one line among a dozen backend probes and
falls back, so a 45-minute genome becomes a 7.5-day one and the only symptom is
that it seems slow. `ntv3-annotate` warns when it detects this, but the warning
scrolls past too, so install the extra.

Already installed the plain package? Check before starting a genome:

```bash
python -c "import jax; print(jax.devices())"   # want CudaDevice, not CpuDevice
```

## 📦 Models

Three clade models. `--model_cfg <name>` brings a model's weights **and the
window length it must be run at** with it.

| Name | Coverage |
| --- | --- |
| `mammals` | NTv3 703M, 1 Mb, 28 mammal assemblies |
| `plants` | NTv3 703M, 1 Mb, 53 plant assemblies |
| `insects` | NTv3 703M, 1 Mb, 93 insect assemblies |

```bash
ntv3-annotate --list_cfg            # what is registered
ntv3-annotate --show_cfg mammals    # everything one entry says
```

Both work with no accelerator and no weights downloaded.

### Where the weights come from

Nothing to stage. Each entry names a released model, and the loader downloads
`InstaDeepAI/<name>/jax_model/` from the Hub on first use and caches it — about
2.7 GiB per clade.

```bash
pip install .
export HF_TOKEN=...          # the repos are gated: request access on the Hub first
ntv3-annotate --genome genome.fa --model_cfg mammals --out out.gff3
```

| Clade | Repository |
| --- | --- |
| `mammals` | `InstaDeepAI/NTv3_annot_mammals` |
| `plants` | `InstaDeepAI/NTv3_annot_plants` |
| `insects` | `InstaDeepAI/NTv3_annot_insects` |

To load from disk instead — an air-gapped host, a shared mirror, or weights you
converted yourself — point `NTV3_ANNOTATE_WEIGHTS_DIR` at a directory holding
one subdirectory per clade, named as the clade is:

```bash
export NTV3_ANNOTATE_WEIGHTS_DIR=/data/ntv3_annotate   # holds mammals/ plants/ insects/
```

It takes precedence over the Hub when set. Each directory holds
`hyperparameters.json`, `flax_model.msgpack` and `tokenizer_config.json`.

The registry records which model each clade expects and checks it against the
weights' own `source_model_name` at load, so a mirror holding the wrong clade
is an error rather than a plausible, wrong GFF3.

## 🚀 Usage

```bash
# GFF3, gzipped
ntv3-annotate --genome genome.fa --model_cfg mammals --out out.gff3.gz

# GTF, plus the spliced CDS and its protein translation
ntv3-annotate --genome genome.fa --model_cfg plants --out out.gtf \
  --codingseq cds.fa --protseq prot.fa

# A model that is not registered: give the geometry explicitly
ntv3-annotate --genome genome.fa --out out.gff3 \
  --override weights_dir=/data/my_model --seq_len 1048576
```

Any config field is reachable with `--override key=value`, repeatable. Named
flags win when both set the same field.

## 🧬 Inputs

`--genome` takes files and directories interchangeably, in the order given. A
directory is globbed for `.fa` / `.fasta` / `.fna` and their `.gz` forms in any
case, and sorted by name — so a genome staged as one FASTA per contig needs no
concatenation.

Record ids are normalised into legal GFF3 seqids as they are read, since the
FASTA header goes straight into column 1. A record id appearing in two inputs
is an error rather than two overlapping feature sets under one seqid.

Remote paths are not supported; download the FASTA first.

## 📄 Output formats

The format follows `--out`'s extension. Anything else is rejected before the
model loads rather than guessed at.

| Extension | Format |
| --- | --- |
| `.gff3`, `.gff` | GFF3 |
| `.gtf` | GTF (Ensembl dialect) |

Either may be followed by `.gz`.

GFF3 rows are `gene`, `CDS` and `intron` only. There is no `mRNA` row, so
parsing a GFF3 for one finds nothing and raises nothing. GTF adds the
transcript level GFF3 leaves out: one per gene, named `<gene_id>.t1`, because
the decoder emits one state path per strand and there are no isoforms to name.
**CDS rows include the terminal stop codon**, which is AUGUSTUS' and Tiberius'
convention; Ensembl and GTF2.2 exclude it.

## 🔁 Resume and provenance

Each record's features are written to their own part file as soon as that
record finishes, and merged at the end, so a whole genome is never buffered in
memory and a job killed hours in keeps everything it finished.

```bash
ntv3-annotate ... --resume
```

A part records the length and a sampled digest of its sequence, and a digest of
the settings that decide what the features are. A resumed run reuses a part
only when both still match, so editing the assembly or changing the decoder
re-annotates the affected records rather than merging features that describe
one input under a header describing another.

The output carries `##ntv3-*` comments recording the weights, geometry, decoder
and filters. They are inert for tabix and IGV, and without them nobody
receiving the file can tell how it was produced.

## 🔬 Post-processing filters

Both are on by default and both change the file you receive.

| Setting | Default | What it drops |
| --- | --- | --- |
| `min_cds_length` | 201 | genes whose total CDS is shorter |
| `apply_biological_filters` | true | genes whose spliced CDS has an internal in-frame stop |

This matters when comparing numbers with anything else: **a metric computed on
filtered output is not comparable to one computed on unfiltered output.** Say
which you used. Set both off to reproduce raw decoder output.

## 📏 Why the window and the crop are not settings

Passing the wrong window length does not raise. It produces a plausible GFF3
that is wrong, which is why the length travels with the weights instead of
being retyped at each call site, and why `sequence_length` is required when you
bypass the registry.

The inference crop is fixed at **0.75** and is deliberately not exposed. Above a
model's trained crop the head is asked for window-edge positions it was never
supervised on, which collapses accuracy rather than failing; below it, quality
holds but each forward pass costs `1/crop` for no gain. Every published number
comes from 0.75.

## 🔍 Comparing two annotations

`check_equivalence` diffs two GFF3s and says whether the difference looks
structural. Differences that cluster at window seams are a stitching, cropping
or padding problem; differences that do not are numerical noise. It reports
each differing gene's distance to the nearest seam against what a uniform
distribution would predict, which is what makes the answer readable rather
than a matter of opinion.

```bash
python -m nucleotide_transformer_v3.gene_annotation.check_equivalence \
  a.gff3 b.gff3
```

Comment lines are ignored, so the provenance header is not a difference. Exits
non-zero when the two differ, so it can gate a job. Run both sides with
`use_bf16_compute=false`: bf16 on a GPU is not bit-reproducible.

## 🐍 Python API

```python
from nucleotide_transformer_v3.gene_annotation.annotate_genome import run
from nucleotide_transformer_v3.gene_annotation.config import (
    ProductionAnnotationConfig,
)

run(ProductionAnnotationConfig(
    fasta_paths=["/data/cow/genome.fa"],
    model_cfg="mammals",
    output_path="/data/cow.gff3.gz",
))
```

`run` logs through the `nucleotide_transformer_v3` logger, which attaches
only a `NullHandler`. Call `logging.basicConfig(level=logging.INFO)` to see
progress; the CLI does this for you.

## ⚖️ Known limitations

- **A resumed run may leave gaps in the gene numbering.** It seeds the counter
  past every id already staged rather than reusing one. Ids stay unique within
  a file, which is what GFF3 requires; they are opaque labels, not identifiers
  whose stability is promised.
- **A resumed run cannot see a length-preserving edit in the middle of a
  chromosome-sized record**, because the sequence check digests a bounded
  sample. Checksum the FASTA yourself if you need that.
- **No isoforms, no UTRs**, and no `start_codon` / `stop_codon` / `intron` rows
  in GTF.
- **Species conditioning is not supported.** None of the three released models
  takes a species input: every export sets `conditions_names` and
  `conditions_vocab_size` to `[]`. The mammal and plant checkpoints they were
  made from did carry a species vocabulary that nothing feeding the annotation
  head read, and the export records it under
  `provenance.dropped_conditioning`; the insect checkpoint never had one. So
  `--species` cannot change the output and is warned about.
