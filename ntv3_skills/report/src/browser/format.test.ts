import assert from "node:assert/strict";
import { test } from "node:test";
import { fitLabel, formatScale, placeGeneLabels } from "./format.ts";

test("scale values keep two significant figures", () => {
  assert.equal(formatScale(0.003369), "0.0034");
  assert.equal(formatScale(5.7e-5), "5.7e-5");
  assert.equal(formatScale(1.77), "1.8");
  assert.equal(formatScale(12.34), "12.3");
  assert.equal(formatScale(250.4), "250");
  assert.equal(formatScale(0), "0");
});

test("a long lane label is shortened, a short one is kept", () => {
  assert.equal(fitLabel("NTv3 · CAGE", 170), "NTv3 · CAGE");
  const cut = fitLabel("Saliency · a very long track name that cannot fit", 62);
  assert.equal(cut.length, 10);
  assert.ok(cut.endsWith("…"));
});

test("gene names stay inside the plot and do not overlap", () => {
  const xs = placeGeneLabels(
    [
      { name: "HBE1", x0: 300, x1: 1000, span: 1_600 },
      { name: "HBG2", x0: 310, x1: 1000, span: 1_600 },
      { name: "OR51B5", x0: 100, x1: 200, span: 900 },
    ],
    180,
    1000,
  );
  assert.equal(xs[2], 204);
  assert.equal(xs[0], 300);
  assert.equal(xs[1], null);
});

test("a readthrough spanning the view does not hide a gene that starts in it", () => {
  // As on chr11 at HBE1: names collapse to one span each, all past the right edge.
  const xs = placeGeneLabels(
    [
      { name: "ENSG00000239920", x0: 184, x1: 992, span: 424_043 },
      { name: "HBE1", x0: 272, x1: 992, span: 237_308 },
      { name: "HBG2", x0: 184, x1: 992, span: 252_418 },
    ],
    184,
    992,
  );
  assert.deepEqual(xs, [null, 272, 184]);
});

test("when every gene spans the view, the shortest one is named", () => {
  const xs = placeGeneLabels(
    [
      { name: "HBG2", x0: 184, x1: 992, span: 252_418 },
      { name: "ENSG00000239920", x0: 184, x1: 992, span: 424_043 },
      { name: "HBE1", x0: 184, x1: 992, span: 237_308 },
    ],
    184,
    992,
  );
  assert.deepEqual(xs, [null, null, 184]);
});
