---
name: ntv3
description: NTv3 for Predict, Saliency, Compare, Annotate, and Explain. Parse
  any request, ask only the real gates, show one Resolved block, take one
  approval, then run. The result is the NTv3 browser, published as an artifact
  in this conversation. Annotate runs the clade model on a GPU after the same
  approval. Never assume human.
disable-model-invocation: true
---

# NTv3

Slash `/ntv3` starts a run. You are an analyst, not a library author.

Do not Grep, Read, or search the repo, and do not open the README. The skill
is already loaded.

## Models (do not ask the user for these)

- **Tracks**: `NTv3_tracks` (`InstaDeepAI/NTv3_tracks`), loaded by this
  repo's `nucleotide_transformer_v3.pretrained.load_tracks_model`.
  Species-conditioned, 16-token soft-masked vocabulary (lowercase = repeats;
  keep it lowercase), any length that is a multiple of 128, released at
  1,048,576 bp. It predicts only the central 37.5% of the window
  (`keep_target_center_fraction`): 1,536 bp of a 4,096 bp window, 393,216 bp
  at 1 Mb. Any multiple of 128 up to 1,048,576 bp runs on CPU or GPU; on CPU,
  16,384 bp or more is slow and the block carries a CPU warning. Pass a length
  only after the user does. Predict computes only the requested
  head columns (`track_indices`), never the whole 7,362-track head.
- **Annotation**: one model per clade, 1,048,576 bp windows, central 75%
  kept. Config names `mammals`, `plants`, `insects`; Hub ids
  `InstaDeepAI/NTv3_annot_mammals`, `_plants`, `_insects` (`insecta` is
  accepted and recorded as `insects`). Decoder `hmm_constrained` unless the
  user names `argmax` or `hmm`. Released filters: CDS shorter than 201 bp
  dropped; internal in-frame stop dropped. The forward runs on a GPU (~25 GB
  per window); a laptop cannot run it. The mammals, plants, and insects
  models are registered only at 1,048,576 bp; a shorter window writes a wrong
  GFF3 instead of failing, so there is no shorter one.
  The input is a FASTA, or the 1 Mb centred on an existing Predict run
  (`--locus-run`), whose genes are then drawn next to its tracks.
- Paper: https://www.biorxiv.org/content/10.64898/2025.12.22.695963 — cite
  it; do not fetch it.

## The asks

| Ask | What this skill does | Command |
|---|---|---|
| Predict | A gene or exact region, the approved species conditioning and window, and only the requested tracks. | `plan --resolve-only` → approval → `plan` → `lookup` → `predict` → `report` → publish |
| Saliency | Per base, the gradient at each base of each requested track's mean over the crop, centred across A, C, G and T and kept at the reference base. | Predict with `--analyses predict saliency`; predict runs it on the same load |
| Compare | One SNV: the requested tracks for the reference and the alternate allele, and how far each moves (alt − ref). | Predict with `--analyses predict compare --variant CHROM POS REF ALT`; predict runs it on the same load |
| Annotate | A FASTA, or the 1 Mb around a Predict run, and its clade model, tiled into 1 Mb windows (central 75% kept), decoded, and written as GFF3. | `annotate --resolve-only` → approval → `annotate` → `approve` → `annotate-forward` (GPU) → `report` → publish |
| Explain | Which model, species conditioning, window, crop, bins, and head produced each number; for Compare, the variant and its crop; for annotation, the model, clade, decoder, window, and filters. | The browser's **Explain** section, read from provenance |

Saliency differentiates the token embedding of the public model, so it needs
no other model and no other loader. Compare is two token forwards, one per
allele; the reference forward is predict's own. Neither runs unless the
request asks for it.

Not in this skill: Integrated Gradients, SmoothGrad, ISM, Ledidi design,
attention maps, interval knockout, sequence design, and motif scans. `--analyses`
takes only `predict`, `saliency`, and `compare` (Annotate is its own
`annotate` command); the plan refuses
`integrated_gradients`, `smoothgrad`, `ism`, `attention`, `knockout`,
`design`, `ledidi`, and `motif` (exit 2) before any fetch or load, and the `ntv3` commands of those names exit 2 the same way. An old on-disk
plan can still hold an `attribution:` block (Integrated Gradients, SmoothGrad,
and ISM settings); it reads, and nothing runs from it. A new plan refuses
integrated gradients, SmoothGrad, and ISM. Also out: a raw-sequence Predict
(regions are fetched from UCSC), more than one variant, and required literature
synthesis. If asked, say so in one line and record it. Do not
substitute another method; saliency is the only attribution.

A new locus or FASTA is a **new** `results/<run>/`. Never overwrite, clear,
or merge a run.

## Hard stop

Write only under **worktree-root** `results/<run>/` (next to `README.md`).
Never `ntv3_skills/results/`.

**Forbidden:** editing `ntv3_tools/`, `.claude/skills/`, `registry.py`.
**Forbidden:** `python -c` / inline Python to write `plan.yaml`. Use the
`ntv3 plan` command only.
**Forbidden:** git mutations, uploads, S3 writes, GCS writes, source edits,
`pnpm` or any server, building another viewer, or printing any `.env` value.
**Forbidden:** any weight download except the `NTv3_tracks` or
`NTv3_annot_<clade>` load the user approved on the block's Weights line.
**Forbidden:** setting `PYTHONPATH` to, or running a Python from, any other
checkout. This repo is the whole install.
**Forbidden:** showing, pasting, or telling the user to open any local HTML
path, including `results/<run>/browser.html`. The user opens the artifact.
If a command exits 2, STOP and ask. Do not invent a CNhs / ENCFF id.
Never pass `--model`. Never pass a window length the user did not give. When
the block has a CPU or GPU warning line, show it as printed; the approval
covers it.

## 1. Parse the request

Read the user's sentence before asking anything. Take every value it already
names: gene, `chrom start end`, assembly or species, tissue or cell type,
assay, variant (`chrom position ref alt`, 1-based), saliency, clade, FASTA
path, decoder, outputs. Do not ask for a value the sentence already gives.
Any gene, locus, tissue, or assay works the same way.

Example: "HBE1 at chr11 5266144 5270240 on hg38 in K562. Predict K562 CAGE."
gives gene HBE1, `chr11 5266144 5270240`, human/hg38, tissue K562, assay
CAGE, ask Predict. Open: the window length, then the browser lanes. Relay
each question.

Example: "On hg38 in K562, what does <chrom> <position> <ref>><alt> do to
CAGE at HBE1, and which bases drive it?", with a real SNV the user typed,
gives gene HBE1, the variant, human/hg38, tissue K562, assay CAGE, and asks
Predict, Compare, and Saliency. With no span, ask for the window length.
Do not center a window until they pass `--window-bp`. Open: the
window length, then the browser lanes. Relay each question.

Example: "Annotate /data/dmel/2L.fa with the insect model." gives the FASTA
and clade `insects`, decoder `hmm_constrained`, device `gpu`. Nothing is open.

Example: after a Predict run, "Also annotate the genes there with the mammal
model." gives clade `mammals` and that run as `--locus-run`. Nothing is open.

Defaults when the sentence is silent: no `--device` (the plan takes `cuda` when
the JAX python sees a GPU, else `cpu`), `--analyses predict`, the report,
no notebook, the critic on. Window length is not a default. Ask for it. Add `saliency` when the
question asks which bases drive a track; add `compare` with `--variant` when
it names a variant or two alleles. `--want-notebook` only when
they ask for a notebook. "hg38" means human; "mm10" means mouse.

Track names: pass the assay as the user wrote it (`--tracks CAGE`) with the
tissue as `--tissue K562`. The tool matches tissue + assay onto its validated
hints with no network call.

## 2. Real gates (ask only these)

Ask one gate at a time. Relay the tool's question unchanged. If none is
open, ask nothing.

- **Species**, when neither species nor assembly is stated. Never assume
  human. Predict supports human (hg38) and mouse (mm10) only.
- **Window**, when no length was given. Relay the tool's question; it names
  the lengths this device runs. Do not pass a length until the user does.
- **Browser lanes**, when the user did not name gene or ClinVar lanes. Relay
  the tool's question. Do not add knownGene, ncbiRefSeq, or clinvarMain.
- **Gene with two loci**: `plan --resolve-only` exits 2 when UCSC returns no
  exact match or distinct loci. Show its message and ask which.
- **Span that is not the window**: a `chrom start end` span must equal the
  window (a multiple of 128, at most 1,048,576 bp). If it does not, state the
  conflict and ask which to change. Never recenter it.
- **Compare with no exact variant**: a Compare ask needs one SNV as
  `chrom position ref alt` (1-based). If the sentence names only an rsID, a
  protein change, or one allele, ask for it. Never look one up and never
  invent it. A variant outside a `chrom start end` span the user gave is a
  span conflict; ask which to change.
- **Assay with no hint and no single ENCODE hit**: if `plan --resolve-only`
  exits 2 and lists ENCODE accessions, show those accessions and ask which
  one. If it lists none, ask for one ENCSR id in a single line. Never pick,
  never invent, and do not send the user to the ENCODE portal or offer a
  different tissue.
- **Clade and input (Annotate only)**: `mammals`, `plants`, or `insects`, and
  either the absolute FASTA path or the Predict run to annotate around. Never
  guess the clade from a species name, and never pick the run.
- **Run and track (Explain only)**, when it is not clear which number they
  mean.

The Weights line is part of the one approval, not a separate question.

Never invent a variant, a species, a track id, a coordinate, a clade, or a
path. Never ask for a JAX path: `ntv3 doctor` finds it.

## 3. Resolved block, one approval

Print the block with `ntv3 plan --resolve-only`. It writes nothing. The run directory is named for you, `results/<assembly>_<chrom>_<start>_<end>_<analyses>`. Never choose one. If species, assembly, tissue, tracks, window, or browser lanes is missing, it exits 2 with one question, in that order: relay that question unchanged and do not pick a species, an assembly, a tissue, a track, a head, a window length, or a browser lane.

Predict:

```bash
ntv3_skills/.venv/bin/ntv3 plan --resolve-only \
  --species … --assembly … --gene … [--chrom … --start … --end …] \
  --tissue … --tracks … --window-bp <length> [--device cpu|cuda] \
  --analyses predict [saliency] [compare] [--variant CHROM POS REF ALT] \
  --want-report [--want-notebook]
```

With `--variant` and no `--chrom/--start/--end`, the window is centered on the
variant. The block then has a **Compare** line (the variant and whether it is
inside the crop) and, with saliency, a **Saliency** line (the tracks).
`--attribution-tracks` limits saliency to some of the requested tracks.

No browser lane is added unless the user passes `--ucsc-tracks`. `none`
means no gene lane and no ClinVar. When they name lanes, `lookup` fetches
those lanes for the window:

- Human gene lane: GENCODE (`knownGene`). Human variants: ClinVar
  (`clinvarMain`). Mouse gene lane: NCBI RefSeq (`ncbiRefSeq`).
- Nothing else. Naming a tissue or an assay adds no lane. Add cCREs, DNase
  clusters, phastCons, or a histone mark only when the user names it.

Annotate:

```bash
ntv3_skills/.venv/bin/ntv3 annotate --resolve-only \
  --clade <mammals|plants|insects> --fasta <absolute-path> --out results/<run> \
  [--decoder argmax|hmm]
ntv3_skills/.venv/bin/ntv3 annotate --resolve-only \
  --clade <mammals|plants|insects> --locus-run results/<predict run>
```

Its block names the FASTA records or the locus window, clade config, Hub
model, a **Weights** line (cached, or a 2.8 GB download that needs a Hugging
Face login with approved access), decoder, window and kept fraction,
filters, a **Forward** line (GPU memory and time), the output, and the JAX
python. When the Forward line says this JAX python will refuse,
say in one sentence that the forward needs a GPU host with the `gpu` extra,
and stop.

Show the block exactly as printed. Each line is a value and the rule that
chose it. The **Model** line names `NTv3_tracks`. The **Weights** line says
whether `InstaDeepAI/NTv3_tracks` is already in the Hugging Face cache; when it
is not, predict's first load downloads it (711M parameters) and the last line
asks to approve that download with the run.

If the JAX python line says none passed, it names the candidates and the
install command (`python3.11 -m venv .venv && .venv/bin/pip install -e .`).
Say in one sentence that predict cannot start until that install is done, and
stop; do not ask them to approve a forward. Do not run the install unless
asked.

The approval may strike a lane, turn the critic off, or change a line.
Apply it to the same flags. Do not show a second block unless a gate opens.

## 4. After approval: run

Approval means run. Run these in order without asking again:

```bash
ntv3_skills/.venv/bin/ntv3 plan <same flags, without --resolve-only>
ntv3_skills/.venv/bin/ntv3 approve --run-dir results/<run>
ntv3_skills/.venv/bin/ntv3 lookup --plan results/<run>/plan.yaml
PYTHONPATH=ntv3_skills <jax-python> -m ntv3_tools.cli predict \
  --plan results/<run>/plan.yaml
ntv3_skills/.venv/bin/ntv3 report --run-dir results/<run>
```

Add `--no-critic` to `approve` only when the user turned the critic off. Predict refuses a plan that is not approved.

`PYTHONPATH=ntv3_skills` is this repo's own skill package; never point it
anywhere else. `<jax-python>` is the path on the block's JAX python line (this
repo's `.venv/bin/python` or `$NTV3_JAX_PYTHON`). Predict loads
`NTv3_tracks` once and computes only the approved tracks, then the approved
saliency (`interpret/`) and Compare (`disrupt/`) on the same load. On CPU
saliency adds about a forward and a backward per track. If `import jax`, the
weight load, or the forward fails, stop and relay the error in one line. Do
not look for another checkout, another model, or another loader.

### Proving saliency and Compare

A run that already exists does not prove saliency or Compare; only a load and
a forward do. A new run name for the same window is not a copy of the old
run; never reuse, refresh, or overwrite the old directory instead, and never
skip the forward.

1. **Saliency first, with no variant.** Do not pick a run name. The directory is `results/<assembly>_<chrom>_<start>_<end>_<analyses>`. Do not use the HBE1 window as the proof. Run plan, approve, lookup, and predict with `--analyses predict saliency` on a locus the user just named. It is proved when predict loads `NTv3_tracks` and writes `tracks/` and `interpret/`. A bad Compare allele can never block this forward, because this plan names no variant.
2. **Then Compare**, on another new run name with
   `--analyses predict saliency compare --variant CHROM POS REF ALT`. Use the
   SNV the user named, or a single-base record from `ucsc/clinvarMain.json`
   already fetched for the locus. Pass it through: predict checks the
   reference base against UCSC before it loads the model, so do not look it
   up yourself. It is proved when predict also writes `disrupt/`. If the
   check fails, quote the error, name the record, and stop; do not try
   another variant. Step 1 still counts.

Do not pass the old invented example `chr11 5268001 C T`: hg38 has T there.

If the load or the forward fails, quote the error and stop.

If predict says ground truth was not compared, say so in one line; the
browser shows the predicted tracks without a measured lane.

For Annotate, approval means run, in order:

```bash
ntv3_skills/.venv/bin/ntv3 annotate <same flags, without --resolve-only>
ntv3_skills/.venv/bin/ntv3 approve --run-dir results/<run>
PYTHONPATH=ntv3_skills <jax-python> -m ntv3_tools.cli annotate-forward \
  --plan results/<run>/annotation_plan.json
ntv3_skills/.venv/bin/ntv3 report --run-dir results/<run> [--refresh]
```

`annotate-forward` runs `ntv3-annotate`'s pipeline on the GPU and writes
`annotation.gff3` and `annotation_manifest.json`. It refuses before any load
when the plan is not approved or JAX sees no GPU; relay that line and stop. A
locus run's GFF3 is in chromosome coordinates and its run already has a
report, so pass `--refresh`. Do not retry or swap models.

## 5. Publish the browser

`report` writes `report.json` and `results/<run>/ntv3_report.html`: the NTv3 report app (report sections and the genome browser) in one file, with this run's report and saved UCSC rows embedded. It loads no script or stylesheet from the network and needs no dev server. That file exists only to be published. Publish it with the **Artifact** tool:

1. Read the Artifact tool's description in this session first and follow its
   page contract (title, icon, size, CSP, libraries, layout) where it differs
   from this step.
2. Publish with `file_path` = the absolute path of
   `results/<run>/ntv3_report.html`, title `<gene> NTv3 report` (two to four
   words), a one-sentence description of the run, and an icon on first
   publish if the schema asks for one. Do not publish `browser.html` or a
   localhost URL instead.
3. Give the user the URL the tool returns and, in two sentences, what the
   page shows. Never show the file path, never paste the HTML, never tell the
   user to open a local file, and do not start a server.

If `report` warns that `ntv3_report.html` is over the 16 MB artifact limit, do
not publish it. Say so in one line, and that a shorter window or fewer tracks or
analyses makes a smaller page.

If the Artifact tool is not available in this session, say in one line that
the browser could not be published from this session. Do not give the path
instead, and do not build or serve another viewer.

The page has:

- **Lanes** at the run's window: the lanes on the Resolved block's Browser
  lanes line, then each predicted track, then the measured ENCODE signal when
  lookup saved it and predict compared it, then a per-base **Saliency** lane
  per track and an **Alt − ref** lane per track when those ran. The span the
  model predicts is shaded and the variant is marked in every lane. The arrow
  buttons or drag pan (the ← and → keys too), + and − zoom, `Back to locus`
  returns, and each lane's swatch recolors it.
- **Compare**: the variant, the largest alt − ref change, and each track's
  change, largest first.
- **NTv3 annotation** lane, when the run holds a GFF3: genes on both strands
  as intron lines with CDS blocks, in the same view.
- **Explain**: the provenance for each number, and for annotation the model
  id, clade, decoder, window, kept fraction, and filters.

`report` on an Annotate run before its forward says the forward has not run.
After `annotate-forward`, `ntv3 report --run-dir results/<run>` writes the
page with the annotation lane (`--refresh` for a locus run, where the lane sits
next to the tracks); publish it the same way. No per-base probabilities exist
in the run, so the lane shows genes only. For a FASTA with several records the
page shows the contig with the most genes; when the user asks for another,
run `ntv3 report --run-dir results/<run> --refresh --contig <seqid>` and
publish again.

After `report --refresh`, publish the same `ntv3_report.html` again so the
artifact shows the update.

## Explain

Answer from the browser's Explain section (from `report.json`): model,
species and assembly conditioning, the input window, the served crop and its
share of the window, bin count and size, the model head and ENCODE accession
per track, and for annotation the model id, clade, decoder, window, kept
fraction, and filters. Name the file it came from. Do not recompute and do
not write an essay.

## Critic

The critic is on unless they say skip. With `want_critic`, `report` audits the
run's claims, including an Annotate run's gene claim, and records the result; launch nothing else. Use `ntv3-critic`
only when the user asks to audit specific claims.

## Environment

The install is this repo only:

```bash
python3.11 -m venv .venv
.venv/bin/pip install -e .
uv sync --project ntv3_skills --python 3.11 --frozen
```

CPU wheels; the `gpu` extra only on a GPU host, for `annotate-forward` and
Predict on the GPU (`.venv/bin/pip install -e '.[gpu]'`), never on a laptop. `ntv3_skills/.venv/bin/ntv3` runs `plan`,
`lookup`, `annotate`, `critic`, `doctor`, `literature`, and `report`; it has no JAX and is
never the JAX python. If that binary is absent, run the `uv sync` line. Do not
use `uv run`.

`ntv3 doctor` checks `$NTV3_JAX_PYTHON` if set, then this repo's
`.venv/bin/python`. A candidate passes when `import jax` and
`import nucleotide_transformer_v3` work; predict needs nothing else. Nothing
else is a candidate. It prints paths and versions only and reads no `.env`.
If the repo `.venv` cannot import the package, it prints the install command;
tell the user to run it.
When it also prints a `gpu` line (this JAX python runs on CPU), relay it only
if the user asks for Annotate or a window of 16,384 bp or more.

## On explicit request only

`literature --run-dir … --records results/<run>/literature/records.json`
(`ntv3_skills/.venv/bin/ntv3`) validates agent-collected primary-source
records. Never invent a citation. Follow it with `report --refresh` and
publish again.

## Dual ids (do not collapse)

- **UCSC `track=`** — browser lanes (`clinvarMain`). Never `ENCSR*`.
- **ENCODE `ENCSR*`** — ground-truth BigWig. `lookup` / GT.
- **Model head** — a column of `bigwigs_per_file_assembly[assembly]`.
  Sometimes the same ENCSR, sometimes `CNhs*`. CAGE GT is `ENCSR000CIO`; the
  head is `CNhs11250_P`.

## Other

- Coordinates are 0-based half-open, except GFF3's 1-based closed features.
- Predictions are not measurements. Predicted genes are not curated
  annotation. No fake %.
- A new locus is a new run directory and a new artifact.
- A person writing their own analysis still has `notebooks/` and the `ntv3`
  CLI. Do not delete them.
