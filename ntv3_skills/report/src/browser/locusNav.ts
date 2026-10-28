export type GenomicRange = {
  chrom: string;
  start: number;
  end: number;
};

const MIN_SPAN = 40;

export const LOCUS_ANNOTATION_LABELS: Record<string, string> = {
  knownGene: "GENCODE V50",
  clinvarMain: "ClinVar",
  cCREregistry: "ENCODE4 cCREs",
  wgEncodeRegDnaseClustered: "DNase Clusters",
  phastCons100way: "phastCons 100-way",
};

/** Label for a locus default track. GENCODE V50 is the hg38 knownGene build. */
export function locusAnnotationLabel(
  name: string,
  assembly = "hg38",
): string | undefined {
  if (name === "knownGene" && assembly !== "hg38") {
    return "UCSC Genes";
  }
  return LOCUS_ANNOTATION_LABELS[name];
}

export function parseCoord(text: string): GenomicRange | null {
  const match = text.replace(/,/g, "").match(/^(chr[\w.]+):(\d+)-(\d+)$/i);
  if (!match) {
    return null;
  }
  const start = Number(match[2]);
  const end = Number(match[3]);
  if (!Number.isFinite(start) || !Number.isFinite(end) || end <= start) {
    return null;
  }
  return { chrom: normalizeChrom(match[1]), start, end };
}

export function normalizeChrom(chrom: string): string {
  const trimmed = chrom.trim();
  if (!/^chr/i.test(trimmed)) {
    return trimmed;
  }
  const rest = trimmed.slice(3);
  if (/^(x|y|m|mt)$/i.test(rest)) {
    return `chr${rest.toUpperCase() === "MT" ? "M" : rest.toUpperCase()}`;
  }
  return `chr${rest}`;
}

export function formatCoord(range: GenomicRange): string {
  return `${range.chrom}:${range.start.toLocaleString()}-${range.end.toLocaleString()}`;
}

export function clampRange(
  start: number,
  end: number,
  chromSize: number,
  minSpan: number = MIN_SPAN,
): { start: number; end: number } {
  let rangeStart = start;
  let rangeEnd = end;
  if (rangeEnd - rangeStart < minSpan) {
    const mid = (rangeStart + rangeEnd) / 2;
    rangeStart = mid - minSpan / 2;
    rangeEnd = mid + minSpan / 2;
  }
  if (rangeEnd - rangeStart > chromSize) {
    return { start: 0, end: chromSize };
  }
  if (rangeStart < 0) {
    rangeEnd += -rangeStart;
    rangeStart = 0;
  }
  if (rangeEnd > chromSize) {
    rangeStart -= rangeEnd - chromSize;
    rangeEnd = chromSize;
  }
  rangeStart = Math.max(0, rangeStart);
  rangeEnd = Math.min(chromSize, rangeEnd);
  return { start: Math.round(rangeStart), end: Math.round(rangeEnd) };
}

/**
 * The UCSC lanes the approved plan named, once each. No lane is added for the
 * run: a plan that passed `--ucsc-tracks none` draws no genes and no ClinVar.
 */
export function reportAnnotationTracks(reportTracks: string[]): string[] {
  return [...new Set(reportTracks)];
}
