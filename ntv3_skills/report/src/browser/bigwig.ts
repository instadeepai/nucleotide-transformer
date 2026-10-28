import { BigWig } from "@gmod/bbi";
import type { Feature } from "@gmod/bbi";
import type { QueryStat, SignalTrack } from "./types";

export type BwSpec = {
  id: string;
  url: string;
  label: string;
  unit: string;
  signed?: boolean;
};

const readers = new Map<string, BigWig>();

function reader(url: string): BigWig {
  const hit = readers.get(url);
  if (hit) return hit;
  const bw = new BigWig({ url });
  readers.set(url, bw);
  return bw;
}

export async function sampleBwTrack(
  spec: BwSpec,
  chrom: string,
  start: number,
  end: number,
  nBins: number,
  signal?: AbortSignal,
): Promise<{ track: SignalTrack; stat: QueryStat }> {
  const url = spec.url;
  const t0 = performance.now();
  let bytes = 0;
  const basesPerSpan = (end - start) / nBins;
  const features = await reader(url).getFeatures(chrom, start, end, {
    signal,
    ...(basesPerSpan >= 2 ? { basesPerSpan } : {}),
    onProgress: (done: number, total: number) => {
      if (total > 0) bytes = done;
    },
  });
  const values = rasterize(features, start, end, nBins);
  const ms = performance.now() - t0;
  const track: SignalTrack = {
    id: spec.id,
    name: spec.label,
    unit: spec.unit,
    source: url.split("/").pop() ?? url,
    signed: spec.signed,
    values,
  };
  const stat: QueryStat = {
    id: spec.id,
    track: spec.id,
    label: spec.label,
    url,
    http: 206,
    ok: true,
    empty: features.length === 0,
    ms,
    bytes,
    items: features.length,
    mode: "signal",
    trackType: "bigWig",
    note: `${nBins} px bins`,
  };
  return { track, stat };
}

function rasterize(
  features: Feature[],
  start: number,
  end: number,
  nBins: number,
): number[] {
  const acc = Array.from({ length: nBins }, () => 0);
  const wt = Array.from({ length: nBins }, () => 0);
  const span = Math.max(1, end - start);
  for (const f of features) {
    const score = f.score;
    if (score == null || !Number.isFinite(score)) continue;
    const s = Math.max(f.start, start);
    const e = Math.min(f.end, end);
    if (e <= s) continue;
    const i0 = Math.max(0, Math.floor(((s - start) / span) * nBins));
    const i1 = Math.min(nBins, Math.ceil(((e - start) / span) * nBins));
    const lo = Math.min(i0, i1);
    const hi = Math.max(i0, i1, lo + 1);
    for (let i = lo; i < hi; i++) {
      acc[i] += score;
      wt[i] += 1;
    }
  }
  return acc.map((v, i) => (wt[i] ? v / wt[i] : 0));
}
