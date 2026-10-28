/**
 * Number and label formatting for the browser lanes. Pure functions, so the
 * node tests cover them without a DOM.
 */

/** A lane's scale value with two significant figures: 0.0034, 5.7e-5, 1.8, 12.3. */
export function formatScale(v: number): string {
  if (!Number.isFinite(v) || v === 0) return "0";
  const a = Math.abs(v);
  if (a >= 100) return v.toFixed(0);
  if (a >= 10) return v.toFixed(1);
  if (a >= 0.001) return String(Number(v.toPrecision(2)));
  return v.toExponential(1);
}

/** Rough width of 11 px semibold sans text, for fitting labels in the gutter. */
export const LABEL_CHAR_PX = 6.2;

/** Shorten a lane label with an ellipsis so it fits ``maxPx`` of gutter. */
export function fitLabel(label: string, maxPx: number, charPx = LABEL_CHAR_PX): string {
  const maxChars = Math.max(1, Math.floor(maxPx / charPx));
  if (label.length <= maxChars) return label;
  return `${label.slice(0, Math.max(1, maxChars - 1))}…`;
}

/** ``x0``/``x1``: visible ends in px; ``span``: genomic length in bp. */
export type GeneLabelInput = { name: string; x0: number; x1: number; span: number };

/**
 * Place gene names in the plot without overlap. ``x0``/``x1`` are the gene's
 * visible ends. A name goes after its gene when it fits; else above the
 * gene's visible start. Genes whose end is in view are placed first, then
 * genes whose start is in view, then genes spanning the whole view, shortest
 * gene first at a tie, so a long readthrough cannot hide a gene inside the
 * window. A name that would overlap
 * a placed one is dropped (it stays in the gene's hover title).
 */
export function placeGeneLabels(
  genes: GeneLabelInput[],
  left: number,
  right: number,
  charPx = 6.5,
): (number | null)[] {
  const plan = genes.map((gene) => {
    const width = gene.name.length * charPx;
    const after = gene.x1 + 4;
    if (after + width <= right) return { x: after, width, tier: 0 };
    const tier = gene.x0 > left ? 1 : 2;
    return { x: Math.max(left, Math.min(gene.x0, right - width)), width, tier };
  });
  const order = genes
    .map((_, i) => i)
    .sort(
      (a, b) =>
        plan[a]!.tier - plan[b]!.tier ||
        plan[a]!.x - plan[b]!.x ||
        genes[a]!.span - genes[b]!.span,
    );
  const out: (number | null)[] = genes.map(() => null);
  const placed: [number, number][] = [];
  for (const i of order) {
    const { x, width } = plan[i]!;
    if (placed.some(([a, b]) => x < b + 6 && a < x + width + 6)) continue;
    out[i] = x;
    placed.push([x, x + width]);
  }
  return out;
}
