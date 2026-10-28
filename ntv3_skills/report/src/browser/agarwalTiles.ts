import type { TileFeat } from "./types";

export type TileRow = {
  chrom: string;
  start: number;
  end: number;
  name: string;
  score: number;
};

/**
 * Agarwal K562 pilot / β-globin control oligos (ENCSR460LZI ENCFF006GZL).
 * HS1–HS4 log2 values match n_04 EXPERIMENTAL_DATA means. HS2 is the strongest.
 * HS5 is not in that pilot BED; the mean is from the n_04 replicate table.
 */
export const AGARWAL_HS_CONTROLS: TileFeat[] = [
  { chrom: "chr11", start: 5_275_890, end: 5_276_090, name: "HS1", score: 0.315 },
  { chrom: "chr11", start: 5_280_660, end: 5_280_860, name: "HS2", score: 1.614 },
  { chrom: "chr11", start: 5_284_690, end: 5_284_890, name: "HS3", score: 0.176 },
  { chrom: "chr11", start: 5_288_180, end: 5_288_380, name: "HS4", score: 0.452 },
  { chrom: "chr11", start: 5_291_330, end: 5_291_530, name: "HS5", score: -0.561 },
];

export type K562DnaseGtTrack = {
  id: string;
  headIndex: number;
  accession: string;
  fileAccession: string;
  url: string;
  label: string;
  origin: string;
};

/**
 * The four K562 DNase-seq heads on the 7362-track metadata (indices 1524–1527).
 * URLs are the ENCODE preferred_default GRCh38 read-depth signals for the
 * files listed in ntv3_7362_metadata.tsv. Do not invent ENCFF ids.
 *
 * Head 1527 trained on ENCFF009SCD (rep 2) + ENCFF971AHO (rep 1). The
 * browser shows ENCFF971AHO, the preferred_default file.
 */
export const K562_DNASE_GT_TRACKS: readonly K562DnaseGtTrack[] = [
  {
    id: "gt-dnase-1524",
    headIndex: 1524,
    accession: "ENCSR000EKN",
    fileAccession: "ENCFF899YDP",
    url: "https://encode-public.s3.amazonaws.com/2017/10/05/88725272-9c6d-4eb1-9ee5-167144b00b7f/ENCFF899YDP.bigWig",
    label: "K562 DNase · ENCSR000EKN",
    origin: "1524 · ENCFF899YDP · n_04 paper head",
  },
  {
    id: "gt-dnase-1525",
    headIndex: 1525,
    accession: "ENCSR000EKP",
    fileAccession: "ENCFF515UNC",
    url: "https://encode-public.s3.amazonaws.com/2017/10/04/6629b70e-1191-42ac-9ebb-30fdef2c8219/ENCFF515UNC.bigWig",
    label: "K562 DNase · ENCSR000EKP",
    origin: "1525 · ENCFF515UNC",
  },
  {
    id: "gt-dnase-1526",
    headIndex: 1526,
    accession: "ENCSR000EKQ",
    fileAccession: "ENCFF708UIS",
    url: "https://encode-public.s3.amazonaws.com/2017/10/05/a80607fe-8118-4fb9-9600-4836a0282009/ENCFF708UIS.bigWig",
    label: "K562 DNase · ENCSR000EKQ",
    origin: "1526 · ENCFF708UIS",
  },
  {
    id: "gt-dnase",
    headIndex: 1527,
    accession: "ENCSR921NMD",
    fileAccession: "ENCFF971AHO",
    url: "https://encode-public.s3.amazonaws.com/2018/01/18/0ff64bdc-4a8d-4975-a348-6091a3f18341/ENCFF971AHO.bigWig",
    label: "K562 DNase · ENCSR921NMD",
    origin: "1527 · ENCFF971AHO (rep 1; training also used ENCFF009SCD)",
  },
];

/** Skill-alias / 1527 representative. Prefer K562_DNASE_GT_TRACKS. */
export const K562_DNASE_GT = K562_DNASE_GT_TRACKS[3];

export function stripOrientationSuffix(name: string): string {
  return name.replace(/_Reversed:?$/u, "");
}

export function averageOrientationScores(rows: TileRow[]): TileFeat[] {
  const grouped = new Map<string, TileRow[]>();
  for (const row of rows) {
    const name = stripOrientationSuffix(row.name);
    const key = `${row.chrom}:${row.start}-${row.end}:${name}`;
    const bucket = grouped.get(key);
    if (bucket) {
      bucket.push(row);
    } else {
      grouped.set(key, [row]);
    }
  }
  const tiles: TileFeat[] = [];
  for (const [key, bucket] of grouped) {
    const [chrom, span, name] = splitTileKey(key);
    const [startText, endText] = span.split("-");
    tiles.push({
      chrom,
      start: Number(startText),
      end: Number(endText),
      name,
      score: bucket.reduce((sum, row) => sum + row.score, 0) / bucket.length,
    });
  }
  tiles.sort((a, b) => a.chrom.localeCompare(b.chrom) || a.start - b.start);
  return tiles;
}

export function indexTiles(tiles: TileFeat[]): Record<string, TileFeat[]> {
  const index: Record<string, TileFeat[]> = {};
  for (const tile of tiles) {
    const list = index[tile.chrom];
    if (list) {
      list.push(tile);
    } else {
      index[tile.chrom] = [tile];
    }
  }
  return index;
}

export function tilesInRange(
  index: Record<string, TileFeat[]>,
  chrom: string,
  start: number,
  end: number,
): TileFeat[] {
  const list = index[chrom] ?? [];
  return list.filter((tile) => tile.end > start && tile.start < end);
}

function splitTileKey(key: string): [string, string, string] {
  const chromEnd = key.indexOf(":");
  const spanEnd = key.indexOf(":", chromEnd + 1);
  return [key.slice(0, chromEnd), key.slice(chromEnd + 1, spanEnd), key.slice(spanEnd + 1)];
}
