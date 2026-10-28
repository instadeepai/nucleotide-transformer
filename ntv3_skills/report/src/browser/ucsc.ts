import type {
  CcreFeat,
  GeneFeat,
  IntervalFeat,
  LiveBundle,
  QueryStat,
} from "./types";

const BASE = "https://api.genome.ucsc.edu";

type JsonObj = Record<string, unknown>;

export type FetchSpec = {
  id: string;
  track: string;
  label: string;
  mode: "genes" | "ccre" | "intervals" | "ticks" | "signal";
  signed?: boolean;
  unit?: string;
};

export function ucscTrackUrl(
  genome: string,
  track: string,
  chrom: string,
  start: number,
  end: number,
  maxItems = 50_000,
): string {
  const q = new URLSearchParams({
    genome,
    track,
    chrom,
    start: String(start),
    end: String(end),
    maxItemsOutput: String(maxItems),
  });
  return `${BASE}/getData/track?${q.toString()}`;
}

export async function fetchChromSize(
  genome: string,
  chrom: string,
  signal?: AbortSignal,
): Promise<number> {
  const response = await fetch(
    `${BASE}/list/chromosomes?${new URLSearchParams({ genome })}`,
    { signal },
  );
  if (!response.ok) {
    throw new Error(`UCSC chromosome request failed with HTTP ${response.status}.`);
  }
  const payload = (await response.json()) as JsonObj;
  const chromosomes = payload.chromosomes;
  if (
    typeof chromosomes !== "object" ||
    chromosomes === null ||
    Array.isArray(chromosomes)
  ) {
    throw new Error(`UCSC returned no size for ${chrom} in ${genome}.`);
  }
  const chromSize = (chromosomes as JsonObj)[chrom];
  if (typeof chromSize !== "number") {
    throw new Error(`UCSC returned no size for ${chrom} in ${genome}.`);
  }
  return chromSize;
}

export function inferMode(
  type: string,
  track: string,
): FetchSpec["mode"] {
  if (track === "cCREregistry") return "ccre";
  const t = type.toLowerCase();
  if (t.includes("genepred") || t === "psl" || t === "bigpsl") return "genes";
  if (t.includes("wig")) return "signal";
  if (t.includes("snp") || t.includes("vcf") || t === "gvf") return "ticks";
  return "intervals";
}

export function emptyBundle(): LiveBundle {
  return {
    genes: {},
    ccres: [],
    intervals: {},
    signals: {},
    ok: {},
    errors: {},
    stats: [],
    wallMs: 0,
  };
}

export async function fetchTracks(
  specs: FetchSpec[],
  genome: string,
  chrom: string,
  start: number,
  end: number,
): Promise<LiveBundle> {
  const bundle = emptyBundle();
  const t0 = performance.now();
  const results = await Promise.all(
    specs.map((spec) => fetchOneSpec(spec, genome, chrom, start, end)),
  );
  bundle.wallMs = performance.now() - t0;
  for (const { spec, json, stat } of results) {
    ingest(bundle, spec, json, stat);
  }
  return bundle;
}

export async function fetchOneSpec(
  spec: FetchSpec,
  genome: string,
  chrom: string,
  start: number,
  end: number,
): Promise<{ spec: FetchSpec; json: JsonObj | null; stat: QueryStat }> {
  const url = ucscTrackUrl(genome, spec.track, chrom, start, end);
  const t0 = performance.now();
  try {
    const res = await fetch(url);
    const buf = await res.arrayBuffer();
    const ms = performance.now() - t0;
    const bytes = buf.byteLength;
    if (!res.ok) {
      const msg = `UCSC ${res.status} ${res.statusText}`;
      return {
        spec,
        json: null,
        stat: baseStat(spec, url, res.status, ms, bytes, msg),
      };
    }
    const json = JSON.parse(new TextDecoder().decode(buf)) as JsonObj;
    const rows = rowsOf(json, spec.track);
    const trackType = String(json.trackType ?? spec.mode);
    const itemsReturned = num(json.itemsReturned) ?? rows.length;
    const bigDataUrl =
      typeof json.bigDataUrl === "string" ? json.bigDataUrl : undefined;
    return {
      spec,
      json,
      stat: {
        id: spec.id,
        track: spec.track,
        label: spec.label,
        url,
        http: res.status,
        ok: true,
        empty: rows.length === 0,
        ms,
        bytes,
        items: rows.length,
        itemsReturned,
        mode: spec.mode,
        trackType,
        bigDataUrl,
      },
    };
  } catch (err) {
    const ms = performance.now() - t0;
    const msg = err instanceof Error ? err.message : "request failed";
    return {
      spec,
      json: null,
      stat: baseStat(spec, url, null, ms, 0, msg),
    };
  }
}

/**
 * Serve specs from UCSC payloads saved with the run, through the same parser
 * a live response goes through. A payload is used only when it is for this
 * chromosome and covers start..end; every other spec is returned to fetch.
 */
export function savedBundle(
  specs: FetchSpec[],
  saved: Record<string, unknown>,
  chrom: string,
  start: number,
  end: number,
): { bundle: LiveBundle; rest: FetchSpec[]; covered: { start: number; end: number } } {
  const bundle = emptyBundle();
  const rest: FetchSpec[] = [];
  let covered = { start: -Infinity, end: Infinity };
  for (const spec of specs) {
    const json = saved[spec.track] as JsonObj | undefined;
    const from = num(json?.start);
    const to = num(json?.end);
    if (!json || json.chrom !== chrom || from == null || to == null || from > start || to < end) {
      rest.push(spec);
      continue;
    }
    const rows = rowsOf(json, spec.track);
    ingest(bundle, spec, json, {
      ...baseStat(spec, String(json._request_url ?? ""), 200, 0, 0, ""),
      ok: true,
      error: undefined,
      empty: rows.length === 0,
      items: rows.length,
      note: "saved with the run",
    });
    covered = { start: Math.max(covered.start, from), end: Math.min(covered.end, to) };
  }
  return { bundle, rest, covered };
}

export function mergeBundles(base: LiveBundle, extra: LiveBundle): LiveBundle {
  const statsById = new Map(base.stats.map((s) => [s.id, s]));
  for (const s of extra.stats) statsById.set(s.id, s);
  return {
    genes: { ...base.genes, ...extra.genes },
    ccres: extra.ccres.length ? extra.ccres : base.ccres,
    intervals: { ...base.intervals, ...extra.intervals },
    signals: { ...base.signals, ...extra.signals },
    ok: { ...base.ok, ...extra.ok },
    errors: { ...withoutKeys(base.errors, Object.keys(extra.ok)), ...extra.errors },
    stats: [...statsById.values()],
    wallMs: base.wallMs || extra.wallMs,
  };
}

function withoutKeys(
  rec: Record<string, string>,
  keys: string[],
): Record<string, string> {
  const out = { ...rec };
  for (const k of keys) delete out[k];
  return out;
}

function ingest(
  bundle: LiveBundle,
  spec: FetchSpec,
  json: JsonObj | null,
  stat: QueryStat,
): void {
  bundle.stats.push(stat);
  if (!stat.ok || json == null) {
    bundle.errors[spec.id] = stat.error ?? "request failed";
    return;
  }
  try {
    const rows = rowsOf(json, spec.track);
    if (spec.mode === "genes") {
      const genes = genesFromRows(rows);
      bundle.genes[spec.id] = genes;
      bundle.ok[spec.id] = genes.length;
      stat.items = genes.length;
    } else if (spec.mode === "ccre") {
      bundle.ccres = parseCcres(rows);
      bundle.ok[spec.id] = bundle.ccres.length;
      stat.items = bundle.ccres.length;
    } else if (spec.mode === "intervals" || spec.mode === "ticks") {
      const items = parseIntervals(rows);
      bundle.intervals[spec.id] = items;
      bundle.ok[spec.id] = items.length;
      stat.items = items.length;
    } else {
      const start = num(json.start) ?? 0;
      const end = num(json.end) ?? start + 1;
      const values = binWig(rows, start, end, 900);
      bundle.signals[spec.id] = {
        id: spec.id,
        name: spec.track,
        unit: spec.unit ?? "value",
        source: `UCSC ${spec.track}`,
        signed: spec.signed,
        values,
      };
      bundle.ok[spec.id] = rows.length;
    }
  } catch (err) {
    const msg = err instanceof Error ? err.message : "parse failed";
    bundle.errors[spec.id] = msg;
    stat.ok = false;
    stat.error = msg;
  }
}

function baseStat(
  spec: FetchSpec,
  url: string,
  http: number | null,
  ms: number,
  bytes: number,
  error: string,
): QueryStat {
  return {
    id: spec.id,
    track: spec.track,
    label: spec.label,
    url,
    http,
    ok: false,
    empty: false,
    ms,
    bytes,
    items: 0,
    mode: spec.mode,
    trackType: "",
    error,
  };
}

function rowsOf(json: JsonObj, track: string): JsonObj[] {
  const direct = json[track];
  if (Array.isArray(direct)) return direct as JsonObj[];
  for (const val of Object.values(json)) {
    if (Array.isArray(val) && val.length && typeof val[0] === "object") {
      return val as JsonObj[];
    }
  }
  return [];
}

function num(v: unknown): number | null {
  if (typeof v === "number" && Number.isFinite(v)) return v;
  if (typeof v === "string" && v !== "") {
    const n = Number(v);
    return Number.isFinite(n) ? n : null;
  }
  return null;
}

/** One gene per name: its canonical transcript when UCSC marks one, else all merged. */
export function genesFromRows(rows: JsonObj[]): GeneFeat[] {
  const canonical = new Set(
    rows.filter(isCanonical).map((row) => geneName(row)),
  );
  return collapseByName(
    parseGenes(
      rows.filter((row) => !canonical.has(geneName(row)) || isCanonical(row)),
    ),
  );
}

// hg38 knownGene tier, e.g. "canonical,basic,all".
function isCanonical(row: JsonObj): boolean {
  return String(row.tier || "").split(",").includes("canonical");
}

function geneName(row: JsonObj): string {
  return String(row.geneName || row.name2 || row.name || "gene");
}

function parseGenes(rows: JsonObj[]): GeneFeat[] {
  const out: GeneFeat[] = [];
  for (const row of rows) {
    const start =
      num(row.chromStart) ?? num(row.txStart) ?? num(row.genoStart);
    const end = num(row.chromEnd) ?? num(row.txEnd) ?? num(row.genoEnd);
    if (start == null || end == null) continue;
    const strand = row.strand === "-" ? "-" : "+";
    const cdsStart = num(row.thickStart) ?? num(row.cdsStart);
    const cdsEnd = num(row.thickEnd) ?? num(row.cdsEnd);
    out.push({
      name: geneName(row),
      start,
      end,
      strand,
      exons: exonsFromRow(row, start),
      ...(cdsStart != null && cdsEnd != null
        ? { cds: { start: cdsStart, end: Math.max(cdsStart, cdsEnd) } }
        : {}),
    });
  }
  return out;
}

function exonsFromRow(row: JsonObj, chromStart: number): GeneFeat["exons"] {
  const exonStarts = parseCsvInts(row.exonStarts);
  const exonEnds = parseCsvInts(row.exonEnds);
  if (exonStarts.length && exonEnds.length) {
    return exonStarts.map((s, i) => ({ start: s, end: exonEnds[i] ?? s }));
  }
  const sizes = parseCsvInts(row.blockSizes);
  const rel = parseCsvInts(row.chromStarts);
  if (sizes.length && rel.length) {
    return rel.map((off, i) => ({
      start: chromStart + off,
      end: chromStart + off + (sizes[i] ?? 0),
    }));
  }
  return [];
}

function parseCcres(rows: JsonObj[]): CcreFeat[] {
  const out: CcreFeat[] = [];
  for (const row of rows) {
    const start = num(row.chromStart);
    const end = num(row.chromEnd);
    if (start == null || end == null) continue;
    const raw = String(row.cCRE_class || row.ccre || row.encodeLabel || "CA");
    const klass = CCRE_KLASS[raw] ?? "CA";
    out.push({
      start,
      end,
      name: String(row.name || ""),
      klass,
    });
  }
  return out;
}

const CCRE_KLASS: Record<string, CcreFeat["klass"]> = {
  PLS: "PLS",
  pELS: "pELS",
  dELS: "dELS",
  "DNase-H3K4me3": "DNase-H3K4me3",
  "CA-H3K4me3": "DNase-H3K4me3",
  "CTCF-only": "CTCF-only",
  "CA-CTCF": "CTCF-only",
  CA: "CA",
};

function parseIntervals(rows: JsonObj[]): IntervalFeat[] {
  const out: IntervalFeat[] = [];
  for (const row of rows) {
    const start =
      num(row.chromStart) ?? num(row.start) ?? num(row.genoStart);
    const end = num(row.chromEnd) ?? num(row.end) ?? num(row.genoEnd);
    if (start == null || end == null) continue;
    const name = String(
      row.repName || row.name || row.geneName || row.description || "",
    );
    out.push({
      start,
      end,
      name,
      klass: String(row.repClass || row.cCRE_class || ""),
    });
  }
  return out;
}

function binWig(
  rows: JsonObj[],
  winStart: number,
  winEnd: number,
  nBins: number,
): number[] {
  const acc = Array.from({ length: nBins }, () => 0);
  const cnt = Array.from({ length: nBins }, () => 0);
  const span = Math.max(1, winEnd - winStart);
  for (const row of rows) {
    const s = num(row.start) ?? num(row.chromStart);
    const e = num(row.end) ?? num(row.chromEnd);
    const v = num(row.value);
    if (s == null || e == null || v == null) continue;
    const mid = (s + e) / 2;
    const i = Math.max(
      0,
      Math.min(nBins - 1, Math.floor(((mid - winStart) / span) * nBins)),
    );
    acc[i] += v;
    cnt[i] += 1;
  }
  return acc.map((v, i) => (cnt[i] ? v / cnt[i] : 0));
}

function collapseByName(genes: GeneFeat[]): GeneFeat[] {
  const byName = new Map<string, GeneFeat>();
  for (const g of genes) {
    const prev = byName.get(g.name);
    if (!prev) {
      byName.set(g.name, g);
      continue;
    }
    byName.set(g.name, {
      ...prev,
      start: Math.min(prev.start, g.start),
      end: Math.max(prev.end, g.end),
      exons: [...prev.exons, ...g.exons],
      cds: joinCds(prev.cds, g.cds),
    });
  }
  return [...byName.values()];
}

/** The span covering both coding spans; an empty (non-coding) span adds nothing. */
function joinCds(a: GeneFeat["cds"], b: GeneFeat["cds"]): GeneFeat["cds"] {
  if (!a || a.end <= a.start) return b ?? a;
  if (!b || b.end <= b.start) return a;
  return { start: Math.min(a.start, b.start), end: Math.max(a.end, b.end) };
}

function parseCsvInts(v: unknown): number[] {
  if (typeof v !== "string" || !v) return [];
  return v
    .split(",")
    .map((x) => Number.parseInt(x, 10))
    .filter((n) => Number.isFinite(n));
}
