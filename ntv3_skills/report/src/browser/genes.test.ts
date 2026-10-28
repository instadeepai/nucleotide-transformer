import assert from "node:assert/strict";
import { test } from "node:test";
import { exonSegments } from "./genes.ts";
import { genesFromRows } from "./ucsc.ts";

// hg38 knownGene rows over chr11:5,266,144-5,270,240 (HBE1), trimmed to the
// fields the parser reads.
const HBE1_WINDOW = [
  { name: "ENST00000380252.7", geneName: "HBG2", strand: "-", chromStart: 5253187, chromEnd: 5505605, thickStart: 5253276, thickEnd: 5254441, blockSizes: "218,223,255,", chromStarts: "0,1104,252163,", tier: "basic,all" },
  { name: "ENST00000380259.7", geneName: "ENSG00000239920", strand: "-", chromStart: 5254391, chromEnd: 5678434, thickStart: 5645851, thickEnd: 5678434, blockSizes: "123,529,127,199,119,254,179,231,", chromStarts: "0,245,91381,336433,340840,372209,391293,423812,", tier: "canonical,basic,all" },
  { name: "ENST00000646569.1", geneName: "ENSG00000239920", strand: "-", chromStart: 5254879, chromEnd: 5328413, thickStart: 5254879, thickEnd: 5254879, blockSizes: "286,116,66,39,58,", chromStarts: "0,3129,3831,4619,73476,", tier: "all" },
  { name: "ENST00000396895.3", geneName: "HBE1", strand: "-", chromStart: 5268344, chromEnd: 5269945, thickStart: 5268468, thickEnd: 5269890, blockSizes: "253,223,147,", chromStarts: "0,1109,1454,", tier: "canonical,basic,all" },
  { name: "ENST00000380237.5", geneName: "HBE1", strand: "-", chromStart: 5268344, chromEnd: 5505604, thickStart: 5268468, thickEnd: 5269890, blockSizes: "253,223,358,43,36,", chromStarts: "0,1109,1454,13564,237224,", tier: "basic,all" },
];

test("a gene with a canonical transcript draws only that transcript", () => {
  const hbe1 = genesFromRows(HBE1_WINDOW).find((g) => g.name === "HBE1");

  assert.deepEqual(
    { start: hbe1?.start, end: hbe1?.end, cds: hbe1?.cds },
    { start: 5268344, end: 5269945, cds: { start: 5268468, end: 5269890 } },
  );
  assert.deepEqual(hbe1?.exons, [
    { start: 5268344, end: 5268597 },
    { start: 5269453, end: 5269676 },
    { start: 5269798, end: 5269945 },
  ]);
});

test("canonical selection is per gene; a gene without one keeps its rows", () => {
  const genes = genesFromRows(HBE1_WINDOW);

  assert.deepEqual(
    genes.map((g) => [g.name, g.start, g.end]),
    [
      ["HBG2", 5253187, 5505605],
      ["ENSG00000239920", 5254391, 5678434],
      ["HBE1", 5268344, 5269945],
    ],
  );
});

test("rows without a tier are merged by name, CDS spans joined", () => {
  const [gene] = genesFromRows([
    { name: "NM_1", name2: "Hbb-bs", strand: "-", txStart: 100, txEnd: 900, cdsStart: 150, cdsEnd: 800, exonStarts: "100,700,", exonEnds: "300,900," },
    { name: "NM_2", name2: "Hbb-bs", strand: "-", txStart: 100, txEnd: 1000, cdsStart: 150, cdsEnd: 950, exonStarts: "100,900,", exonEnds: "300,1000," },
  ]);

  assert.deepEqual(
    { start: gene.start, end: gene.end, cds: gene.cds, exons: gene.exons.length },
    { start: 100, end: 1000, cds: { start: 150, end: 950 }, exons: 4 },
  );
});

test("a non-coding transcript has an empty CDS span", () => {
  const [gene] = genesFromRows([HBE1_WINDOW[2]]);

  assert.deepEqual(gene.cds, { start: 5254879, end: 5254879 });
});

test("an exon splits into thin UTR and thick CDS parts", () => {
  assert.deepEqual(
    exonSegments({ start: 5268344, end: 5268597 }, { start: 5268468, end: 5269890 }),
    [
      { start: 5268344, end: 5268468, thick: false },
      { start: 5268468, end: 5268597, thick: true },
    ],
  );
  assert.deepEqual(
    exonSegments({ start: 5269798, end: 5269945 }, { start: 5268468, end: 5269890 }),
    [
      { start: 5269798, end: 5269890, thick: true },
      { start: 5269890, end: 5269945, thick: false },
    ],
  );
});

test("an exon of a non-coding gene is thin, and one without a CDS is thick", () => {
  assert.deepEqual(exonSegments({ start: 10, end: 20 }, { start: 5, end: 5 }), [
    { start: 10, end: 20, thick: false },
  ]);
  assert.deepEqual(exonSegments({ start: 10, end: 20 }, undefined), [
    { start: 10, end: 20, thick: true },
  ]);
});
