import type { GeneFeat } from "./types";

type Span = { start: number; end: number };

/** Split an exon into thick (CDS) and thin (UTR) parts; no CDS span draws thick. */
export function exonSegments(
  exon: Span,
  cds: GeneFeat["cds"],
): (Span & { thick: boolean })[] {
  if (!cds) return [{ ...exon, thick: true }];
  const parts = [
    { start: exon.start, end: Math.min(exon.end, cds.start), thick: false },
    {
      start: Math.max(exon.start, cds.start),
      end: Math.min(exon.end, cds.end),
      thick: true,
    },
    { start: Math.max(exon.start, cds.end), end: exon.end, thick: false },
  ];
  return parts.filter((part) => part.end > part.start);
}
