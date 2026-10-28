import assert from "node:assert/strict";
import { test } from "node:test";
import { explainLaneRows } from "./lanes.ts";

const TRACK = {
  label: "CAGE CNhs11250_P",
  modelHead: "CNhs11250_P",
  gtAccession: "ENCSR000CIO",
  gtUnit: "signal of unique reads",
};
const LETTERS = [{ label: "CAGE CNhs11250_P", modelHead: "CNhs11250_P" }];

test("each lane gets its unit: prediction, ground truth, saliency", () => {
  const rows = explainLaneRows(
    [TRACK],
    [{ modelHead: "CNhs11250_P", unit: "model logit", observed: [null, 1.5] }],
    LETTERS,
  );
  assert.deepEqual(
    rows.map((row) => [row.lane, row.unit]),
    [
      ["NTv3 · CAGE CNhs11250_P", "model logit"],
      ["Signal · CAGE CNhs11250_P", "signal of unique reads"],
      ["Saliency · CAGE CNhs11250_P", "contribution"],
    ],
  );
});

test("no stored ground truth, no signal row; an old report says signal", () => {
  const assays = [{ modelHead: "CNhs11250_P", unit: "model logit", observed: [] }];
  assert.equal(explainLaneRows([TRACK], assays, []).length, 1);
  const old = { label: "CAGE", modelHead: "CNhs11250_P", gtAccession: "ENCSR000CIO" };
  const rows = explainLaneRows(
    [old],
    [{ modelHead: "CNhs11250_P", unit: "model logit", observed: [0.2] }],
    [],
  );
  assert.equal(rows[1]?.unit, "signal");
});
