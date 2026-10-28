---
name: ntv3-critic
description: LLM critic subagent. Use only when the user asks to audit
  specific claims in an NTv3 run (Predict or Annotate). Do not run in the
  same thread that wrote the report.
disable-model-invocation: true
---

# NTv3 Critic

## Use when

The user asks to audit specific claims. With `want_critic`, `report` already
audits the run's claims; this skill is not launched for that.

- Parent launches this as a **subagent**. Same thread will rubber-stamp.

## Subagent brief (parent pastes this)

You only see `results/<run>/`. You did not write the analysis. Read the artifacts.
Write `claims.json`. Then run the mechanical gate. KEEP or SUPPRESS only.
Do not invent evidence or a confidence %.

Evidence is the files the numbers came from: `report.json` (its `explain` and
`annotation` blocks), `tracks/`, `gt/`, `annotation_plan.json`,
`annotation_manifest.json`, and the run's GFF3. `ntv3_report.html` is only
the page the parent publishes, not evidence of its own; cite the file behind the number. Do not open a server or fetch the
network, never cite a file outside the run, and never show a local HTML path.
Every claim names at least one evidence path; a literature claim cites
`literature/citations.bib`.

Claims to SUPPRESS on sight: a prediction stated as a measurement, a
ground-truth agreement claim when predict did not compare ground truth,
predicted genes stated as curated annotation, an
annotation value (decoder, window, filters) that neither the GFF3 `##ntv3-`
header nor `annotation_plan.json` records, and any GFF3 or gene claim for an
Annotate run with no `annotation_manifest.json`.

## Steps (critic agent)

1. List files under the run dir. Missing file → SUPPRESS that claim.
2. Write `results/<run>/claims.json` as a list of
   `{claim_id, text, evidence_paths, citekey, kind}`.
3. Run:

```bash
ntv3_skills/.venv/bin/ntv3 critic --run-dir results/<run>
```

4. Open `results/<run>/critic/audit_log.json`. Do not paste it. One sentence:
   N KEEP / N SUPPRESS and the path.
5. Return to the parent. The parent then runs `ntv3 report --run-dir results/<run> --refresh`
   and publishes the rewritten browser page again as the artifact.
