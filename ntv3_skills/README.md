# NTv3 skill

The NTv3 skill is a Claude Code skill, `/ntv3`, and the `ntv3`
command-line tool it drives. You describe a locus and the tracks you want in
one sentence. The skill asks one question at a time for anything missing,
shows one **Resolved** block (each value and the rule that chose it), and runs
only after you approve it. Each run writes one new directory under `results/`
and one self-contained report page.

| Ask | What it does |
|---|---|
| Predict | The tracks you name, at a gene or an exact region, for the species you name (human or mouse), from `NTv3_tracks`. Only the requested model heads are computed. |
| Saliency | Per base, the gradient of each track's mean over the crop, centred across A, C, G and T (their mean subtracted) and kept at the reference base. |
| Compare | One SNV: each track for the reference and the alternate allele, and the change (alt − ref). |
| Annotate | Gene models from `InstaDeepAI/NTv3_annot_<clade>` for a FASTA, or for the 1 Mb around a Predict run, drawn in the same browser. The forward runs on a GPU after approval. |
| Explain | The model, species, window, crop, bins, and model head behind each number, read from the run's files. |

## Requirements

- macOS or Linux, with network access to `api.genome.ucsc.edu`,
  `www.encodeproject.org`, `encode-public.s3.amazonaws.com`, and
  `huggingface.co`.
- `git`, Python 3.11, and [`uv`](https://docs.astral.sh/uv/).
- A Claude account that can use Claude Code.
- Disk space for the model weights in the Hugging Face cache. The first
  Predict downloads them.

Predict, Saliency, and Compare run on CPU or GPU at any length up to
1,048,576 bp; on CPU, 16,384 bp or more is slow and comes with a warning. A GPU
host also runs Annotate's forward.

## Install

### 1. Claude Code

Follow the install guide at <https://code.claude.com/docs>. With Node.js 18
or newer, one way is:

```bash
npm install -g @anthropic-ai/claude-code
claude --version
```

### 2. uv and Python 3.11

```bash
curl -LsSf https://astral.sh/uv/install.sh | sh
uv python install 3.11
```

### 3. This repository and its two environments

```bash
git clone https://github.com/instadeepai/nucleotide-transformer.git
cd nucleotide-transformer
"$(uv python find 3.11)" -m venv .venv
.venv/bin/pip install -e .           # laptop or CPU host
# .venv/bin/pip install -e '.[gpu]'  # GPU host: long windows and Annotate
uv sync --project ntv3_skills --python 3.11 --frozen
```

- `.venv` (repository root) runs the model: `jax` and
  `nucleotide_transformer_v3`. Pick one install line. Plain `-e .` gives CPU
  wheels, even on a GPU host, because PyPI's `jaxlib` is CPU-only; everything
  then runs on CPU. The `gpu` extra adds CUDA wheels.
- `ntv3_skills/.venv` runs the planner: `ntv3 plan`, `lookup`, `annotate`,
  `approve`, `doctor`, `critic`, `literature`, and `report`. It has no JAX.
  Do not use `uv run` after setup; it can re-resolve the root project.

Installing downloads no model weights.

### 4. Check the install

```bash
ntv3_skills/.venv/bin/ntv3 doctor
```

`doctor` checks `$NTV3_JAX_PYTHON` (when set), then `.venv/bin/python`. The
first one that imports `jax` and `nucleotide_transformer_v3` is the JAX
Python. If none passes, it prints the install command. When the JAX Python
runs on CPU, it also prints a `gpu` line: Annotate's forward refuses and long
Predict windows are slow, with the `gpu` install. It loads no
checkpoint and reads no `.env`.

To check whether ground truth can be compared:

```bash
.venv/bin/python -c "import pyBigWig; print(pyBigWig.remote, pyBigWig.numpy)"
```

`1 1` means predictions are compared with ENCODE signal. Anything else: see
Troubleshooting.

## Credentials

- **Claude Code**: sign in with your Claude account the first time you run
  `claude`.
- **UCSC and ENCODE**: public, no key.
- **Hugging Face**: `NTv3_tracks` downloads on the first Predict. The
  Resolved block's Weights line says whether a download will happen. The NTv3
  models are gated: request access once on any model's Hub page, then run
  `.venv/bin/huggingface-cli login` or set `HF_TOKEN`.

The tool reads no `.env` file and prints no credential.

## Models

- **Tracks**: `NTv3_tracks` (`InstaDeepAI/NTv3_tracks`), loaded with
  `nucleotide_transformer_v3.pretrained.load_tracks_model("NTv3_tracks")` and
  called as in `docs/nucleotide_transformer_v3.md`. It is
  species-conditioned, keeps lowercase (repeat) bases, takes any input length
  that is a multiple of 128, and predicts the central 37.5% of the window. For a
  4,096 bp window, predictions cover the central 1,536 bp; at its released
  1,048,576 bp, 393,216 bp. Predict always passes `track_indices`, so only the
  requested head columns are computed. No other model name loads.
- **Annotation**: one model per clade, `InstaDeepAI/NTv3_annot_mammals`,
  `_plants`, and `_insects` (`insecta` is accepted and recorded as
  `insects`), registered only at 1,048,576 bp windows with the central 75%
  kept. The decoder is `hmm_constrained` unless you name `argmax` or `hmm`.
  The released filters drop genes with less than 201 bp of CDS and genes with
  an internal in-frame stop. One window needs about 25 GB of GPU memory, so
  the forward runs on a GPU host.
- Paper: <https://www.biorxiv.org/content/10.64898/2025.12.22.695963>.

## CPU and GPU

- `ntv3 plan` picks the device from the JAX python `doctor` finds: `cuda`
  when it sees a GPU, else `cpu`. `--device cpu` or `--device cuda` overrides
  it, and the Window line of the Resolved block says which and why.
- Predict, Saliency, and Compare run at any multiple of 128 up to 1,048,576 bp
  on either device. Longer windows stop with a question.
- On CPU, a window of 16,384 bp or more adds a CPU warning line to the
  Resolved block (with a Saliency memory note above 131,072 bp), and predict
  logs it again. Measured on 16 cores: about 9 GB
  of RAM and 33 s per forward at 131,072 bp; Saliency needs roughly 0.25 GB of
  RAM per extra kb (16 GB at 32,768 bp). A GPU host runs both faster.
- On a GPU, Saliency above 131,072 bp adds a GPU warning line: on an 80 GB H100
  it ran at 131,072 bp and ran out of memory at 1,048,576 bp, where Predict
  and Compare fit.
- `--device cuda` on a machine where JAX sees no GPU stops before anything
  loads. It does not fall back to CPU.
- The `gpu` extra (`.venv/bin/pip install -e '.[gpu]'`) is for GPU Predict, the
  Annotate forward (`ntv3 annotate-forward`), and whole-genome `ntv3-annotate`
  on a GPU host (install step 3). See
  [`docs/ntv3_genome_annotation.md`](../docs/ntv3_genome_annotation.md).
  Do not install it on a laptop.

## Safety

- Commands write only to a new `results/<run>/` at the repository root; the
  CLI rejects other output paths. An existing run directory is never reused
  or overwritten.
- Nothing commits, pushes, uploads, or writes to S3 or GCS.
- The only weight downloads are `NTv3_tracks` on predict's first load and
  `NTv3_annot_<clade>` on the annotation forward's, and the Resolved block
  says whether one will happen before you approve. Both repos are gated, so a
  download needs a Hugging Face login with approved access.
- UCSC and ENCODE calls are read-only. Ground truth comes from public ENCODE
  HTTPS URLs.

## Your first ask

Start Claude Code in the repository root, so it finds the skill in
`.claude/skills/`:

```bash
claude
```

Accept the folder when asked, then type:

```text
/ntv3 HBE1 at chr11 5266144 5270240 on hg38 in K562. Predict K562 CAGE with a 4,096 bp window and the GENCODE and ClinVar lanes.
```

The skill replies with one Resolved block. It names:

- the assembly, locus, window, and crop;
- the model and whether its weights are cached;
- the model head and ENCODE experiment for each track;
- the browser lanes, the analyses, the critic, the run directory, and the JAX
  Python.

Reply `approve`, or strike a lane or change a line. The skill then runs
`plan`, `approve`, `lookup`, `predict`, and `report`, and publishes the
report page.

If your sentence leaves something out, the skill asks for it, one question at
a time, in this order: species, assembly, tissue, tracks, window length,
browser lanes. It never picks a species, tissue, track, window, or lane for
you. It also stops and asks when:

- a gene has two loci;
- your span is not the window length;
- an assay has no single ENCODE match.

To add Saliency, ask which bases drive a track. To add Compare, name one SNV
as `chrom position ref alt` (1-based).

## What a run writes

Each run is a new directory, `results/<assembly>_<chrom>_<start>_<end>_<analyses>/`,
analyses joined with `+` in plan order (for example
`hg38_chr7_55113274_55117370_predict+saliency`).

```text
results/<run>/
├── plan.yaml                  # the approved plan
├── resolved.txt, approval.json
├── locus.json                 # gene, transcript, strand, and TSS for a gene lookup
├── track_resolve.json         # each track: your words, ENCODE experiment, model head
├── gt/manifest.json           # public ENCODE file for each track
├── gt/observed.npy            # measured signal, when ground truth was compared
├── tables/track_agreement.csv # Pearson and Spearman per track, when compared
├── ucsc/*.json                # the browser lanes you approved
├── tracks/tracks.npy, accessions.json, metadata.json
├── interpret/<head>.npy|json  # Saliency, when asked
├── disrupt/                   # Compare, when asked
├── critic/audit_log.json      # KEEP / SUPPRESS for each claim
├── run_manifest.json          # model, device used, commit, input sequence hash, ground-truth status
├── report.json, report.md
├── ntv3_report.html           # the report page: sections and genome browser, no network needed
├── browser.html               # older static page
├── code/run_analysis.py       # replays the plan into a new directory
├── analysis.ipynb             # only when requested
└── ro-crate-metadata.json     # file hashes
```

The report page shows the lanes you approved, each predicted track, the
measured ENCODE signal when it was compared, the Saliency and alt − ref lanes
when they ran, a Compare table, and an Explain section. With the critic on,
`report` audits the run's claims and records KEEP or SUPPRESS for each.
Predictions are model output, not measurements.

The page embeds only what it draws (values at four significant digits) to stay
under the Claude artifact size limit.

An Annotate run writes `annotation_plan.json` and its Resolved block, then,
after the forward, `annotation.gff3` and `annotation_manifest.json` (model,
decoder, GPU, seconds, genes). `report` draws the GFF3 as the NTv3 annotation
lane (genes and CDS blocks, both strands); for a locus run it sits next to the
predicted tracks. A FASTA with several records shows the contig with the most
genes; `report --refresh --contig <seqid>` shows another. With the critic on,
`report` audits the gene claim against the GFF3 and the manifest.

## Command line without Claude Code

```bash
ntv3_skills/.venv/bin/ntv3 plan --resolve-only \
  --species human --assembly hg38 --gene HBE1 \
  --chrom chr11 --start 5266144 --end 5270240 \
  --tissue K562 --tracks CAGE --window-bp 4096 \
  --ucsc-tracks knownGene clinvarMain --analyses predict --want-report
# The same command without --resolve-only writes the plan, then:
ntv3_skills/.venv/bin/ntv3 approve --run-dir results/<run>
ntv3_skills/.venv/bin/ntv3 lookup --plan results/<run>/plan.yaml
PYTHONPATH=ntv3_skills .venv/bin/python -m ntv3_tools.cli predict --plan results/<run>/plan.yaml
ntv3_skills/.venv/bin/ntv3 report --run-dir results/<run>
```

- Leave out `--chrom`, `--start`, and `--end` to look the gene up on UCSC.
  The window is centred on its TSS. Distinct loci stop with a question.
- `--ucsc-tracks none` adds no lane. Leaving it out stops with the lanes
  question.
- `--analyses` takes `predict`, `saliency`, and `compare`. Compare needs
  `--variant CHROM POS REF ALT` (1-based); with no span, the window is centred
  on the variant.
- Annotate has its own command, for a FASTA or for the 1 Mb around a Predict
  run:

  ```bash
  ntv3_skills/.venv/bin/ntv3 annotate --resolve-only --clade mammals \
    --fasta <absolute-path> --out results/<run>   # or --locus-run results/<predict run>
  # The same command without --resolve-only writes the plan, then:
  ntv3_skills/.venv/bin/ntv3 approve --run-dir results/<run>
  PYTHONPATH=ntv3_skills .venv/bin/python -m ntv3_tools.cli annotate-forward \
    --plan results/<run>/annotation_plan.json
  ntv3_skills/.venv/bin/ntv3 report --run-dir results/<run>   # --refresh for a locus run
  ```

  `annotate-forward` refuses an unapproved or edited plan, and a JAX that sees
  no GPU, before any load. It runs the `nucleotide_transformer_v3.gene_annotation`
  pipeline, the one behind `ntv3-annotate`. A locus run's GFF3 is written in
  chromosome coordinates.
- Integrated Gradients, SmoothGrad, ISM, attention maps, interval knockout,
  sequence design, and motif scans are refused before any fetch or load.

`approve` is the only way to set `approved: true`; it records digests of the
plan and the block in `approval.json`. Predict refuses a plan that is not
approved, or that changed after approval. `approve --no-critic` turns the
critic off for that run.

To rank a Compare change against seeded random SNVs in the same window, run
`PYTHONPATH=ntv3_skills .venv/bin/python -m ntv3_tools.cli compare --plan … --chrom … --position … --ref … --alt … --null-n N --seed S`
on a plan without `compare`. Each null SNV is one more forward.

Replay a run into a new directory:

```bash
ntv3_skills/.venv/bin/python results/<run>/code/run_analysis.py \
  --research-python .venv/bin/python --output-run results/<run>-replay
```

## Track identifiers

One assay can have three identifiers, and the tool keeps them apart:

- a UCSC `track=` name for a browser lane;
- an ENCODE `ENCSR*` experiment for ground truth;
- a model head, one column of the model's track list (an `ENCSR*` or a
  FANTOM `CNhs*` id).

`K562 CAGE` is ground truth `ENCSR000CIO` and model head `CNhs11250_P`. A
track with no known model head stops with a question; no other head is
substituted.

## Troubleshooting

**"Ground truth not compared: this pyBigWig build cannot open remote files."**
The PyPI `pyBigWig` wheel is built without libcurl. Predict still saves the
tracks, records the reason in `run_manifest.json` (`ground_truth_note`), and
the page shows the tracks without a measured lane. To compare ground truth,
rebuild pyBigWig against libcurl (libcurl development headers are needed, for
example `libcurl4-openssl-dev` on Debian or Ubuntu):

```bash
.venv/bin/pip install --no-binary pyBigWig --no-build-isolation --force-reinstall pyBigWig
.venv/bin/python -c "import pyBigWig; print(pyBigWig.remote, pyBigWig.numpy)"   # expect 1 1
```

**Replay fails with `No module named 'numpy'`.**
Pass `--research-python` as the venv path that `ntv3 doctor` prints (for
example `.venv/bin/python`), not the interpreter that path links to.

**`JAX python: none`.**
The repository `.venv` is missing or cannot import the package. From the
repository root, run `"$(uv python find 3.11)" -m venv .venv && .venv/bin/pip install -e .`.

**`/ntv3` is not offered.**
Start `claude` from the repository root.

**`--device cuda … sees no CUDA device`.**
Plan again with `--device cpu`.

## Known gaps

- No raw-sequence Predict (regions come from UCSC), and one variant per run.
- Annotate's forward needs a GPU (about 25 GB per window); on a laptop it
  stops before any load.
- Saliency at 1,048,576 bp runs out of memory on an 80 GB GPU.

## Development

```bash
ntv3_skills/.venv/bin/python -m unittest discover -s ntv3_skills/tests -t ntv3_skills
.venv/bin/python -m unittest ntv3_skills/tests/test_predict_public_loader.py
ntv3_skills/.venv/bin/ruff check ntv3_skills/ntv3_tools ntv3_skills/tests
.venv/bin/python -c "import pyBigWig; print(pyBigWig.remote)"
```

The last line prints 1 when this install can compare ground truth, and 0 when
predict will record why it did not.

`tests/test_predict_public_loader.py` runs in the repository `.venv` with a
fake model, or a tiny randomly initialised model of the released class. No
weights are downloaded.

The report page shell is `ntv3_tools/report_app.html`. It is built from the
report app in `ntv3_skills/report/` (React, Vite, zod; the report schema is
`src/schema/report.ts`). Rebuild it after the app changes (Node 20+ and pnpm 9):

```bash
cd ntv3_skills/report && pnpm install --frozen-lockfile && pnpm test && cd ../..
PYTHONPATH=ntv3_skills ntv3_skills/.venv/bin/python -m ntv3_tools.report_app --app ntv3_skills/report
```

Commit the rebuilt `report_app.html` with the source change. To work on the
app against one run, start `pnpm dev` with `VITE_REPORT_URL` set to that run's
`report.json`.
