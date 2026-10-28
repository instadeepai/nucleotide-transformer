import assert from "node:assert/strict";
import { readFileSync } from "node:fs";
import { test } from "node:test";
import {
  AGARWAL_HS_CONTROLS,
  K562_DNASE_GT,
  K562_DNASE_GT_TRACKS,
  averageOrientationScores,
  indexTiles,
  stripOrientationSuffix,
  tilesInRange,
} from "./agarwalTiles.ts";

test("K562 DNase GT list covers the four training heads", () => {
  assert.equal(K562_DNASE_GT_TRACKS.length, 4);
  assert.deepEqual(
    K562_DNASE_GT_TRACKS.map((track) => track.headIndex),
    [1524, 1525, 1526, 1527],
  );
  assert.deepEqual(
    K562_DNASE_GT_TRACKS.map((track) => track.accession),
    ["ENCSR000EKN", "ENCSR000EKP", "ENCSR000EKQ", "ENCSR921NMD"],
  );
  assert.equal(new Set(K562_DNASE_GT_TRACKS.map((track) => track.url)).size, 4);
  assert.equal(K562_DNASE_GT.accession, "ENCSR921NMD");
  assert.match(K562_DNASE_GT.url, /ENCFF971AHO\.bigWig$/);
});

test("Agarwal HS controls rank HS2 strongest", () => {
  const ranked = [...AGARWAL_HS_CONTROLS].sort((a, b) => b.score - a.score);
  assert.equal(ranked[0]?.name, "HS2");
  assert.equal(ranked[0]?.score, 1.614);
});

test("pilot pack keeps HS2 and covers all seven loci", () => {
  const pack = JSON.parse(
    readFileSync(new URL("../data/agarwalK562PilotTiles.json", import.meta.url), "utf8"),
  ) as { tiles: { chrom: string; name: string; score: number }[] };
  const hs2 = pack.tiles.find((tile) => tile.name === "HS2");
  assert.equal(hs2?.score, 1.614);
  const chroms = new Set(pack.tiles.map((tile) => tile.chrom));
  for (const chrom of ["chr11", "chr16", "chrX", "chr8", "chr2", "chr20"]) {
    assert.ok(chroms.has(chrom), chrom);
  }
});

test("stripOrientationSuffix drops the reverse tag", () => {
  assert.equal(stripOrientationSuffix("peak14457_Reversed:"), "peak14457");
  assert.equal(stripOrientationSuffix("peak14457"), "peak14457");
});

test("averageOrientationScores means the two strands", () => {
  const tiles = averageOrientationScores([
    {
      chrom: "chr11",
      start: 5_280_660,
      end: 5_280_860,
      name: "peak14457",
      score: 0.508,
    },
    {
      chrom: "chr11",
      start: 5_280_660,
      end: 5_280_860,
      name: "peak14457_Reversed:",
      score: 0.257,
    },
  ]);
  assert.equal(tiles.length, 1);
  assert.equal(tiles[0]?.name, "peak14457");
  assert.equal(tiles[0]?.score, 0.3825);
});

test("tilesInRange keeps overlapping intervals only", () => {
  const index = indexTiles([
    {
      chrom: "chr16",
      start: 160_000,
      end: 160_200,
      name: "a",
      score: 1,
    },
    {
      chrom: "chr16",
      start: 400_000,
      end: 400_200,
      name: "b",
      score: 1,
    },
  ]);
  const visible = tilesInRange(index, "chr16", 100_000, 230_000);
  assert.deepEqual(
    visible.map((tile) => tile.name),
    ["a"],
  );
});
