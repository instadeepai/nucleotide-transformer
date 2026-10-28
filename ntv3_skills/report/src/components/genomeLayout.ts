/**
 * Shared genomic x-scale for Inference panels (annotations + tracks).
 * Keep GUT / PW / ticks identical so coords line up top-to-bottom.
 */
export const GENOME_FW = 920;
export const GENOME_GUT = 120;
export const GENOME_RIGHT = 16;
export const GENOME_PW = GENOME_FW - GENOME_GUT - GENOME_RIGHT;

export function fmtBp(n: number): string {
  return n.toLocaleString("en-US");
}

export function genomePx(
  bp: number,
  winStart: number,
  winEnd: number,
): number {
  const span = Math.max(1, winEnd - winStart);
  return GENOME_GUT + ((bp - winStart) / span) * GENOME_PW;
}

/** Nice major ticks shared by annotation rail and signal tracks. */
export function genomeTicks(winStart: number, winEnd: number): number[] {
  const span = Math.max(1, winEnd - winStart);
  const tickEvery = span > 20000 ? 10000 : 5000;
  const first = Math.ceil(winStart / tickEvery) * tickEvery;
  const ticks: number[] = [];
  for (let t = first; t < winEnd; t += tickEvery) ticks.push(t);
  return ticks;
}
