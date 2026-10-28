import type { TrackGroup, TrackKind } from "../browser/types";
import type { FetchSpec } from "../browser/ucsc";
import { inferMode } from "../browser/ucsc";

export type CatalogLeaf = {
  id: string;
  label: string;
  type: string;
};

export type CatalogTrack = {
  id: string;
  label: string;
  long: string;
  type: string;
  group: string;
  visibility: string;
  nLeaves: number;
  leaves: CatalogLeaf[];
  priority: number;
};

export const UCSC_GROUP_ORDER = [
  "map",
  "genes",
  "phenDis",
  "varRep",
  "hprc",
  "rna",
  "expression",
  "singleCell",
  "regulation",
  "compGeno",
  "rep",
  "other",
] as const;

export const UCSC_GROUP_LABEL: Record<string, string> = {
  map: "Mapping & sequencing",
  genes: "Genes & predictions",
  phenDis: "Phenotypes & literature",
  varRep: "Variation",
  hprc: "Human pangenome (HPRC)",
  rna: "RNA & transcriptome",
  expression: "Expression",
  singleCell: "Single-cell RNA-seq",
  regulation: "Regulation",
  compGeno: "Comparative genomics",
  rep: "Repeats",
  other: "Other",
};

const GROUP_TO_SLOT: Record<string, TrackGroup> = {
  map: "Mapping",
  genes: "Annotations",
  phenDis: "Phenotypes",
  varRep: "Variation",
  hprc: "HPRC",
  rna: "RNA",
  expression: "Expression",
  singleCell: "Single cell",
  regulation: "Regulation",
  compGeno: "Conservation",
  rep: "Repeats",
  other: "Regulation",
};

const GROUP_RANK: Record<string, number> = {
  map: 4,
  genes: 14,
  phenDis: 24,
  varRep: 25,
  hprc: 26,
  rna: 15,
  expression: 36,
  singleCell: 37,
  regulation: 38,
  compGeno: 90,
  rep: 83,
  other: 39,
};

export function catalogToSpec(
  track: string,
  label: string,
  type: string,
): FetchSpec {
  const mode = inferMode(type, track);
  return {
    id: `ucsc:${track}`,
    track,
    label,
    mode,
    signed: /phylop/i.test(track),
    unit: mode === "signal" ? "value" : undefined,
  };
}

export function catalogSlotKind(type: string, track: string): TrackKind {
  return inferMode(type, track);
}

export function catalogGroup(group: string): TrackGroup {
  return GROUP_TO_SLOT[group] ?? "Regulation";
}

export function catalogRank(group: string): number {
  return GROUP_RANK[group] ?? 40;
}
