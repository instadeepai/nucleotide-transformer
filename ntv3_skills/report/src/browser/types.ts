export type Strand = "+" | "-";

export type GeneFeat = {
  name: string;
  start: number;
  end: number;
  strand: Strand;
  exons: { start: number; end: number }[];
  /** Coding span (thickStart..thickEnd); empty for a non-coding transcript. */
  cds?: { start: number; end: number };
};

export type CcreFeat = {
  name: string;
  start: number;
  end: number;
  klass: "PLS" | "pELS" | "dELS" | "DNase-H3K4me3" | "CTCF-only" | "CA";
};

export type IntervalFeat = {
  start: number;
  end: number;
  name: string;
  klass?: string;
};

export type TileFeat = {
  chrom: string;
  start: number;
  end: number;
  name: string;
  score: number;
};

export type SignalTrack = {
  id: string;
  name: string;
  unit: string;
  source: string;
  values: number[];
  signed?: boolean;
  spanStart?: number;
  spanEnd?: number;
  /** Genomic sequence aligned to spanStart..spanEnd (dynseq letters). */
  seq?: string;
};

export type ElementMark = {
  name: string;
  start: number;
  end: number;
};

export type LocusPack = {
  assembly: string;
  chrom: string;
  start: number;
  end: number;
  chromSize: number;
  predicted: SignalTrack[];
  saliency: SignalTrack[];
  /** Gene models from report.json, keyed by slot id (no UCSC fetch). */
  genes: Record<string, GeneFeat[]>;
  /** UCSC getData payloads saved with the run, keyed by track name. */
  saved: Record<string, unknown>;
  /** Ground truth stored in the report (prediction bins); never fetched. */
  observed: SignalTrack[];
};

export type TrackKind =
  | "genes"
  | "ccre"
  | "intervals"
  | "ticks"
  | "signal"
  | "gt"
  | "pred"
  | "saliency";

export type TrackGroup =
  | "Annotations"
  | "Variation"
  | "Regulation"
  | "MPRA"
  | "Ground truth"
  | "Predictions"
  | "Importance"
  | "Repeats"
  | "Conservation"
  | "Mapping"
  | "Phenotypes"
  | "HPRC"
  | "RNA"
  | "Expression"
  | "Single cell";

export type TrackSlot = {
  id: string;
  kind: TrackKind;
  label: string;
  group: TrackGroup;
  /** Shared key so GT DNase sits on top of NTv3 DNase. */
  assay?: string;
  /** UCSC /getData/track name. Live-filled when set. */
  ucscTrack?: string;
  defaultOn?: boolean;
  signed?: boolean;
  stackRank: number;
  origin: string;
  /** Shown in the row when it has no features, instead of an empty lane. */
  empty?: string;
};

export type QueryStat = {
  id: string;
  track: string;
  label: string;
  url: string;
  http: number | null;
  ok: boolean;
  empty: boolean;
  ms: number;
  bytes: number;
  items: number;
  itemsReturned?: number;
  mode: string;
  trackType: string;
  bigDataUrl?: string;
  error?: string;
  note?: string;
};

export type LiveBundle = {
  genes: Record<string, GeneFeat[]>;
  ccres: CcreFeat[];
  intervals: Record<string, IntervalFeat[]>;
  signals: Record<string, SignalTrack>;
  ok: Record<string, number>;
  errors: Record<string, string>;
  stats: QueryStat[];
  wallMs: number;
};
