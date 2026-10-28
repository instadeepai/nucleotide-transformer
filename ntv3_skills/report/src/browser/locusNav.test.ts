import assert from "node:assert/strict";
import { readFileSync } from "node:fs";
import { test } from "node:test";
import * as locusNav from "./locusNav.ts";
import {
  clampRange,
  locusAnnotationLabel,
  reportAnnotationTracks,
  parseCoord,
} from "./locusNav.ts";

test("parseCoord keeps the chromosome", () => {
  const parsed = parseCoord("chr2:60,400,000-60,620,000");
  assert.deepEqual(parsed, {
    chrom: "chr2",
    start: 60_400_000,
    end: 60_620_000,
  });
});

test("parseCoord normalizes chrX", () => {
  const parsed = parseCoord("CHRX:48700000-48900000");
  assert.equal(parsed?.chrom, "chrX");
});

test("parseCoord rejects a backwards range", () => {
  assert.equal(parseCoord("chr11:200-100"), null);
});

test("no hardcoded locus list or gene jump shortcut is exported", () => {
  assert.equal("AGARWAL_LOCI" in locusNav, false);
  assert.equal("activeAgarwalLocus" in locusNav, false);
  assert.equal("LOCUS_ANNOTATION_TRACKS" in locusNav, false);
});

test("clampRange keeps the span when it still fits the chromosome", () => {
  assert.deepEqual(clampRange(100, 200, 150), { start: 50, end: 150 });
});

test("clampRange fills a chromosome that is shorter than the window", () => {
  assert.deepEqual(clampRange(0, 400, 150), { start: 0, end: 150 });
});

test("reportAnnotationTracks keeps only the plan's lanes, once each", () => {
  assert.deepEqual(reportAnnotationTracks(["phastCons100way", "phastCons100way"]), [
    "phastCons100way",
  ]);
});

test("reportAnnotationTracks adds no gene or ClinVar lane for a plan with none", () => {
  assert.deepEqual(reportAnnotationTracks([]), []);
});

test("the browser draws no Agarwal tiles, K562 DNase lanes, or gene jump bar", () => {
  for (const file of ["../sections/BrowserSection.tsx", "./GenomeBrowser.tsx"]) {
    const source = readFileSync(new URL(file, import.meta.url), "utf8");
    for (const absent of [
      "agarwal",
      "K562_DNASE",
      "AGARWAL_LOCI",
      "gb-loci",
      "TileLane",
    ]) {
      assert.equal(source.includes(absent), false, `${file} still has ${absent}`);
    }
  }
});

test("locusAnnotationLabel keeps GENCODE V50 for hg38 knownGene only", () => {
  assert.equal(locusAnnotationLabel("knownGene", "hg38"), "GENCODE V50");
  assert.equal(locusAnnotationLabel("knownGene", "mm10"), "UCSC Genes");
  assert.equal(locusAnnotationLabel("cCREregistry", "hg38"), "ENCODE4 cCREs");
});
