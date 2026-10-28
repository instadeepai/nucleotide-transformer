import {
  lazy,
  Suspense,
  useEffect,
  useMemo,
  useRef,
  useState,
  type MouseEvent as ReactMouseEvent,
  type PointerEvent as ReactPointerEvent,
} from "react";
import {
  catalogGroup,
  catalogRank,
  catalogSlotKind,
  catalogToSpec,
} from "../data/ucscGroups";
import { sampleBwTrack, type BwSpec } from "./bigwig";
import { ColorPalette } from "./ColorPicker";
import {
  formatCoord,
  parseCoord,
  type GenomicRange,
} from "./locusNav";
import { exonSegments } from "./genes";
import type {
  CcreFeat,
  GeneFeat,
  IntervalFeat,
  LiveBundle,
  LocusPack,
  SignalTrack,
  TrackGroup,
  TrackKind,
  TrackSlot,
} from "./types";
import {
  emptyBundle,
  fetchChromSize,
  fetchTracks,
  mergeBundles,
  savedBundle,
  type FetchSpec,
} from "./ucsc";
import { fitLabel, formatScale, placeGeneLabels } from "./format";
import { useViewport } from "./useViewport";

const CatalogPanel = lazy(async () => {
  const m = await import("./CatalogPanel");
  return { default: m.CatalogPanel };
});

// A wider gutter so lane labels such as "Saliency · CAGE CNhs11250_P" fit.
const GUT = 184;
const RIGHT = 12;
const PW = 808;
const FW = GUT + PW + RIGHT;
const CCRE: Record<CcreFeat["klass"], string> = {
  PLS: "#e11d48",
  pELS: "#f59e0b",
  dELS: "#eab308",
  "DNase-H3K4me3": "#fda4af",
  "CTCF-only": "#0ea5e9",
  CA: "#b0b0b0",
};

const NUC_COLOR: Record<string, string> = {
  A: "#00b050",
  C: "#3973fc",
  G: "#fdb932",
  T: "#ff554d",
};

const PRED_BLUE = "#1f63c4";
const GT_GRAY = "#3d4653";

function isColorable(slot: TrackSlot): boolean {
  return (
    slot.kind === "gt" ||
    slot.kind === "pred" ||
    slot.kind === "signal" ||
    slot.kind === "saliency"
  );
}

function defaultColor(slot: TrackSlot): string {
  if (slot.kind === "pred") return PRED_BLUE;
  if (slot.kind === "gt") return GT_GRAY;
  return GT_GRAY;
}

function rowHeight(slot: TrackSlot): number {
  if (slot.kind === "genes") return 44;
  if (slot.kind === "ccre" || slot.kind === "intervals" || slot.kind === "ticks") {
    return 22;
  }
  if (slot.kind === "saliency") return 56;
  return 40;
}

type Populate = {
  annotations: boolean;
  gt: boolean;
  pred: boolean;
  saliency: boolean;
};

type Props = {
  pack: LocusPack;
  populate: Populate;
  slots: TrackSlot[];
  bwTracks: Record<string, BwSpec>;
};

function defaultOnMap(slots: TrackSlot[]): Record<string, boolean> {
  return Object.fromEntries(slots.map((s) => [s.id, s.defaultOn !== false]));
}

function defaultColorMap(slots: TrackSlot[]): Record<string, string> {
  return Object.fromEntries(slots.map((s) => [s.id, defaultColor(s)]));
}

export function GenomeBrowser({ pack, populate, slots, bwTracks }: Props) {
  const slotList = slots;
  const bwMap = bwTracks;
  const storageKey = useMemo(
    () => browserStorageKey(pack),
    [pack.assembly, pack.chrom, pack.start, pack.end],
  );
  const initialExtraSlots = useMemo(
    () => readPersistedExtraSlots(storageKey),
    [storageKey],
  );
  const bwKey = Object.entries(bwMap)
    .map(([id, spec]) => `${id}:${spec.url}`)
    .sort()
    .join("|");
  const [viewChrom, setViewChrom] = useState(pack.chrom);
  const [viewChromSize, setViewChromSize] = useState(pack.chromSize);
  const { viewStart, viewEnd, setView, shiftFrac, zoomAt } = useViewport(
    pack.start,
    pack.end,
    viewChromSize,
  );
  const [on, setOn] = useState<Record<string, boolean>>(() =>
    defaultOnMap([...slotList, ...initialExtraSlots]),
  );
  const [catalogQ, setCatalogQ] = useState("");
  const [coordText, setCoordText] = useState(
    formatCoord({ chrom: pack.chrom, start: pack.start, end: pack.end }),
  );
  const [live, setLive] = useState<LiveBundle | null>(null);
  const [liveState, setLiveState] = useState<"idle" | "loading" | "ok" | "error">(
    "idle",
  );
  const [liveError, setLiveError] = useState<string | null>(null);
  const [extraSlots, setExtraSlots] =
    useState<TrackSlot[]>(initialExtraSlots);
  const [pendingTrack, setPendingTrack] = useState<string | null>(null);
  const [bwSignals, setBwSignals] = useState<Record<string, SignalTrack>>({});
  const [bwErrors, setBwErrors] = useState<Record<string, string>>({});
  const [colors, setColors] = useState<Record<string, string>>(() =>
    defaultColorMap([...slotList, ...initialExtraSlots]),
  );
  const [colorEdit, setColorEdit] = useState<{
    id: string;
    label: string;
    top: number;
    left: number;
  } | null>(null);
  const annRange = useRef({
    chrom: pack.chrom,
    start: pack.start,
    end: pack.end,
  });
  const [jumpError, setJumpError] = useState<string | null>(null);
  const savedRef = useRef<{ chrom: string; bundle: LiveBundle } | null>(null);
  const allSlots = useMemo(
    () => [...slotList, ...extraSlots],
    [slotList, extraSlots],
  );
  const drag = useRef<{
    x: number;
    start: number;
    end: number;
    moved: boolean;
  } | null>(null);

  const span = Math.max(1, viewEnd - viewStart);
  const px = (bp: number) => GUT + ((bp - viewStart) / span) * PW;
  const bpAt = (clientX: number, el: Element) => {
    const r = el.getBoundingClientRect();
    const x = clientX - r.left;
    const t = (x - GUT * (r.width / FW)) / (PW * (r.width / FW));
    return viewStart + t * span;
  };

  useEffect(() => {
    persistExtraSlots(storageKey, extraSlots);
  }, [storageKey, extraSlots]);

  useEffect(() => {
    setOn((prev) => {
      const next = { ...prev };
      let changed = false;
      for (const slot of slotList) {
        if (!(slot.id in next)) {
          next[slot.id] = slot.defaultOn !== false;
          changed = true;
        }
      }
      return changed ? next : prev;
    });
  }, [slotList]);

  useEffect(() => {
    if (viewChrom === pack.chrom) {
      setViewChromSize(pack.chromSize);
    }
  }, [pack.chrom, pack.chromSize, viewChrom]);

  const goToRange = (range: GenomicRange) => {
    setJumpError(null);
    setCoordText(formatCoord(range));
    if (range.chrom === viewChrom) {
      setView(range.start, range.end, viewChromSize);
      return;
    }
    setViewChrom(range.chrom);
    setView(range.start, range.end, 300_000_000);
    setLive(null);
    annRange.current = { chrom: range.chrom, start: -1, end: -1 };
    void fetchChromSize(pack.assembly, range.chrom)
      .then((size) => {
        setViewChromSize(size);
        setView(range.start, range.end, size);
      })
      .catch((err: unknown) => {
        setJumpError(
          err instanceof Error ? err.message : `No size for ${range.chrom}.`,
        );
      });
  };

  useEffect(() => {
    let cancelled = false;
    setLiveState("loading");
    const specs = slotList
      .filter((slot) => slot.ucscTrack)
      .map(fetchSpecForSlot);
    const pad = Math.max(1, viewEnd - viewStart);
    const start = Math.max(0, viewStart - pad);
    const end = Math.min(viewChromSize, viewEnd + pad);
    // Tracks saved with the run draw from the report; only the rest are fetched.
    const stored = savedBundle(specs, pack.saved, viewChrom, viewStart, viewEnd);
    savedRef.current = { chrom: viewChrom, bundle: stored.bundle };
    const fetched = stored.rest.length
      ? fetchTracks(stored.rest, pack.assembly, viewChrom, start, end)
      : Promise.resolve(emptyBundle());
    fetched
      .then((bundle) => {
        if (cancelled) return;
        annRange.current = {
          chrom: viewChrom,
          start: Math.max(start, stored.covered.start),
          end: Math.min(end, stored.covered.end),
        };
        setLive(mergeBundles(stored.bundle, bundle));
        setLiveState("ok");
      })
      .catch((err: unknown) => {
        if (cancelled) return;
        setLiveState("error");
        setLiveError(err instanceof Error ? err.message : "UCSC fetch failed");
      });
    return () => {
      cancelled = true;
    };
  }, [pack.assembly, viewChrom, slotList]);

  useEffect(() => {
    const ids = Object.entries(on)
      .filter(([id, isOn]) => isOn && id in bwMap)
      .map(([id]) => id)
      .filter((id) => !(id.startsWith("gt-") && !populate.gt));
    if (!ids.length) return;
    const ac = new AbortController();
    const span = Math.max(1, viewEnd - viewStart);
    const fetchStart = Math.max(0, viewStart - span);
    const fetchEnd = Math.min(viewChromSize, viewEnd + span);
    const nBins = Math.min(
      2400,
      Math.max(PW, Math.round(PW * ((fetchEnd - fetchStart) / span))),
    );
    const handle = window.setTimeout(() => {
      void (async () => {
        const next: Record<string, SignalTrack> = {};
        const nextErrors: Record<string, string> = {};
        await Promise.all(
          ids.map(async (id) => {
            const spec = bwMap[id];
            if (!spec) return;
            try {
              const { track } = await sampleBwTrack(
                spec,
                viewChrom,
                fetchStart,
                fetchEnd,
                nBins,
                ac.signal,
              );
              track.spanStart = fetchStart;
              track.spanEnd = fetchEnd;
              next[id] = track;
            } catch (error: unknown) {
              nextErrors[id] =
                error instanceof Error ? error.message : "BigWig load failed";
            }
          }),
        );
        if (ac.signal.aborted) return;
        setBwSignals((prev) => ({ ...prev, ...next }));
        setBwErrors((prev) => {
          const updated = { ...prev };
          for (const id of ids) delete updated[id];
          return { ...updated, ...nextErrors };
        });
      })();
    }, 80);
    return () => {
      ac.abort();
      window.clearTimeout(handle);
    };
  }, [viewStart, viewEnd, on, populate.gt, viewChrom, viewChromSize, bwKey]);

  useEffect(() => {
    if (
      annRange.current.chrom === viewChrom &&
      viewStart >= annRange.current.start &&
      viewEnd <= annRange.current.end
    ) {
      return;
    }
    const span = viewEnd - viewStart;
    const start = Math.max(0, viewStart - span);
    const end = Math.min(viewChromSize, viewEnd + span);
    const handle = window.setTimeout(() => {
      const specs = allSlots
        .filter((slot) => on[slot.id] && slot.ucscTrack)
        .map(fetchSpecForSlot);
      if (!specs.length) return;
      fetchTracks(specs, pack.assembly, viewChrom, start, end)
        .then((bundle) => {
          annRange.current = { chrom: viewChrom, start, end };
          // A failed live fetch keeps the rows saved with the run on their chrom.
          const saved =
            savedRef.current?.chrom === viewChrom ? savedRef.current.bundle : emptyBundle();
          setLive(mergeBundles(saved, bundle));
        })
        .catch(() => undefined);
    }, 200);
    return () => window.clearTimeout(handle);
  }, [
    viewStart,
    viewEnd,
    allSlots,
    on,
    pack.assembly,
    viewChrom,
    viewChromSize,
  ]);

  const visible = orderStack(slotList, on, populate, extraSlots);

  const onPointerDown = (e: ReactPointerEvent<SVGSVGElement>) => {
    const bounds = e.currentTarget.getBoundingClientRect();
    const gutterWidth = GUT * (bounds.width / FW);
    if (e.clientX - bounds.left < gutterWidth) return;
    e.currentTarget.setPointerCapture(e.pointerId);
    drag.current = { x: e.clientX, start: viewStart, end: viewEnd, moved: false };
  };
  const onPointerMove = (e: ReactPointerEvent<SVGSVGElement>) => {
    if (!drag.current) return;
    const dx = e.clientX - drag.current.x;
    if (!drag.current.moved && Math.abs(dx) < 4) return;
    drag.current.moved = true;
    const r = e.currentTarget.getBoundingClientRect();
    const dbp =
      (-dx / (PW * (r.width / FW))) * (drag.current.end - drag.current.start);
    setView(drag.current.start + dbp, drag.current.end + dbp);
  };
  const onPointerUp = () => {
    drag.current = null;
  };
  const onDoubleClick = (e: ReactMouseEvent<SVGSVGElement>) => {
    const el = e.currentTarget;
    const r = el.getBoundingClientRect();
    const gutPx = GUT * (r.width / FW);
    if (e.clientX - r.left < gutPx) return;
    const zoomOut = e.shiftKey || e.altKey;
    zoomAt(zoomOut ? 1 / 3 : 3, bpAt(e.clientX, el));
  };

  const goCoord = () => {
    const parsed = parseCoord(coordText);
    if (!parsed) return;
    goToRange(parsed);
  };

  const resetLocus = () => {
    goToRange({ chrom: pack.chrom, start: pack.start, end: pack.end });
  };

  const annotationWindow = () => {
    const span = Math.max(1, viewEnd - viewStart);
    return {
      start: Math.max(0, viewStart - span),
      end: Math.min(viewChromSize, viewEnd + span),
    };
  };

  const reloadUcsc = () => {
    setLiveState("loading");
    const specs = allSlots
      .filter((slot) => on[slot.id] && slot.ucscTrack)
      .map(fetchSpecForSlot);
    const window = annotationWindow();
    fetchTracks(specs, pack.assembly, viewChrom, window.start, window.end)
      .then((bundle) => {
        setLive(bundle);
        setLiveState("ok");
        setLiveError(null);
      })
      .catch((err: unknown) => {
        setLiveState("error");
        setLiveError(err instanceof Error ? err.message : "UCSC fetch failed");
      });
  };

  const queryCatalog = (
    track: string,
    label: string,
    type: string,
    group: string,
  ) => {
    const existing =
      slotList.find((s) => s.ucscTrack === track) ??
      extraSlots.find((s) => s.ucscTrack === track);
    const spec: FetchSpec = existing
      ? {
          id: existing.id,
          track,
          label: existing.label,
          mode: kindToMode(existing.kind),
          signed: existing.signed,
        }
      : catalogToSpec(track, label, type);
    if (!existing) {
      const slot: TrackSlot = {
        id: spec.id,
        kind: catalogSlotKind(type, track),
        label,
        group: catalogGroup(group),
        ucscTrack: track,
        defaultOn: true,
        signed: /phylop/i.test(track),
        stackRank: catalogRank(group),
        origin: `UCSC ${track}`,
      };
      setExtraSlots((prev) => [...prev, slot]);
      setColors((prev) => ({ ...prev, [slot.id]: defaultColor(slot) }));
    }
    setOn((prev) => ({ ...prev, [spec.id]: true }));
    setPendingTrack(track);
    const window = annotationWindow();
    fetchTracks([spec], pack.assembly, viewChrom, window.start, window.end)
      .then((extra) => {
        setLive((prev) => (prev ? mergeBundles(prev, extra) : extra));
        setLiveState("ok");
      })
      .catch((err: unknown) => {
        setLiveError(err instanceof Error ? err.message : "UCSC fetch failed");
      })
      .finally(() => setPendingTrack(null));
  };

  const hideCatalog = (track: string) => {
    const builtin = slotList.find((s) => s.ucscTrack === track);
    if (builtin) {
      setOn((prev) => ({ ...prev, [builtin.id]: false }));
      return;
    }
    const extra = extraSlots.find((s) => s.ucscTrack === track);
    if (extra) removeExtra(extra.id);
  };

  const removeExtra = (id: string) => {
    setExtraSlots((prev) => prev.filter((s) => s.id !== id));
    setOn((prev) => {
      const next = { ...prev };
      delete next[id];
      return next;
    });
  };

  const loadedTracks = new Set(
    allSlots
      .filter(
        (s) =>
          s.ucscTrack &&
          on[s.id] &&
          live &&
          (s.id in live.ok || s.id in live.errors),
      )
      .map((s) => s.ucscTrack as string),
  );

  return (
    <div className="gb">
      <div className="gb-toolbar">
        <div className="gb-nav">
          <div className="gb-cluster">
            <span className="gb-kicker">Move</span>
            <button type="button" onClick={() => shiftFrac(-0.95)}>
              {"<<<"}
            </button>
            <button type="button" onClick={() => shiftFrac(-0.5)}>
              {"<<"}
            </button>
            <button type="button" onClick={() => shiftFrac(-0.1)}>
              {"<"}
            </button>
            <button type="button" onClick={() => shiftFrac(0.1)}>
              {">"}
            </button>
            <button type="button" onClick={() => shiftFrac(0.5)}>
              {">>"}
            </button>
            <button type="button" onClick={() => shiftFrac(0.95)}>
              {">>>"}
            </button>
          </div>
          <div className="gb-cluster">
            <span className="gb-kicker">Zoom in</span>
            <button type="button" onClick={() => zoomAt(1.5, (viewStart + viewEnd) / 2)}>
              1.5x
            </button>
            <button type="button" onClick={() => zoomAt(3, (viewStart + viewEnd) / 2)}>
              3x
            </button>
            <button type="button" onClick={() => zoomAt(10, (viewStart + viewEnd) / 2)}>
              10x
            </button>
          </div>
          <div className="gb-cluster">
            <span className="gb-kicker">Zoom out</span>
            <button type="button" onClick={() => zoomAt(1 / 1.5, (viewStart + viewEnd) / 2)}>
              1.5x
            </button>
            <button type="button" onClick={() => zoomAt(1 / 3, (viewStart + viewEnd) / 2)}>
              3x
            </button>
            <button type="button" onClick={resetLocus}>
              back to locus
            </button>
          </div>
        </div>
        <form
          className="gb-coord"
          onSubmit={(e) => {
            e.preventDefault();
            goCoord();
          }}
        >
          <input
            value={coordText}
            onChange={(e) => setCoordText(e.target.value)}
            aria-label="Genomic range"
          />
          <button type="submit">go</button>
        </form>
      </div>

      <div className="gb-pos">
        {viewChrom}:{viewStart.toLocaleString()}-{viewEnd.toLocaleString()}{" "}
        <span>{span.toLocaleString()} bp</span>
        {viewChrom !== pack.chrom ? (
          <span className="gb-pos-note"> · GT + UCSC only</span>
        ) : null}
        {jumpError ? <span className="err"> · {jumpError}</span> : null}
      </div>

      <Ideogram
        chrom={viewChrom}
        chromSize={viewChromSize}
        viewStart={viewStart}
        viewEnd={viewEnd}
      />

      <svg
        className="gb-stack"
        viewBox={`0 0 ${FW} ${stackHeight(visible)}`}
        role="img"
        aria-label="Aligned genome tracks"
        onPointerDown={onPointerDown}
        onPointerMove={onPointerMove}
        onPointerUp={onPointerUp}
        onPointerCancel={onPointerUp}
        onDoubleClick={onDoubleClick}
      >
        <TrackStack
          pack={pack}
          viewChrom={viewChrom}
          live={live}
          bwSignals={bwSignals}
          bwErrors={bwErrors}
          viewStart={viewStart}
          viewEnd={viewEnd}
          px={px}
          slots={visible}
          populate={populate}
          onColorSlot={(slot, top, left) =>
            setColorEdit({ id: slot.id, label: slot.label, top, left })
          }
          colors={colors}
        />
        </svg>
      {colorEdit ? (
        <ColorPalette
          label={colorEdit.label}
          value={colors[colorEdit.id] ?? GT_GRAY}
          top={colorEdit.top}
          left={colorEdit.left}
          onChange={(hex) =>
            setColors((prev) => ({ ...prev, [colorEdit.id]: hex }))
          }
          onClose={() => setColorEdit(null)}
        />
      ) : null}

      <div className="gb-picker">
        <div className="gb-picker-head">
          <span>Visible tracks</span>
          <button type="button" onClick={reloadUcsc} disabled={liveState === "loading"}>
            {liveState === "loading" ? "UCSC…" : "Reload UCSC tracks"}
          </button>
          {liveState === "ok" && live ? (
            <span className="ok">
              live {Object.keys(live.ok).length}/{Object.keys(live.ok).length + Object.keys(live.errors).length} tracks
            </span>
          ) : null}
          {liveState === "error" ? <span className="err">{liveError}</span> : null}
        </div>
        {(
          [
            "Mapping",
            "Annotations",
            "Phenotypes",
            "Variation",
            "HPRC",
            "RNA",
            "Expression",
            "Single cell",
            "Regulation",
            "MPRA",
            "Ground truth",
            "Predictions",
            "Importance",
            "Repeats",
            "Conservation",
          ] as TrackGroup[]
        ).map((group) => {
          const slots = allSlots.filter((s) => s.group === group);
          if (!slots.length) return null;
          return (
            <div key={group} className="gb-group">
              <div className="gb-group-lab">{group}</div>
              {slots.map((s) => {
                const extra = extraSlots.some((e) => e.id === s.id);
                return (
                  <div key={s.id} className="gb-track-opt">
                    <label>
                      <input
                        type="checkbox"
                        checked={Boolean(on[s.id])}
                        onChange={() => {
                          if (extra && on[s.id]) {
                            removeExtra(s.id);
                            return;
                          }
                          const next = !on[s.id];
                          setOn((prev) => ({ ...prev, [s.id]: next }));
                          if (
                            next &&
                            s.ucscTrack &&
                            !(live && s.id in live.ok)
                          ) {
                            const spec = fetchSpecForSlot(s);
                            fetchTracks(
                              [spec],
                              pack.assembly,
                              viewChrom,
                              annRange.current.start,
                              annRange.current.end,
                            ).then((bundle) => {
                              setLive((prev) =>
                                prev ? mergeBundles(prev, bundle) : bundle,
                              );
                            }).catch(() => undefined);
                          }
                        }}
                      />
                      {s.label}
                    </label>
                    {extra ? (
                      <button
                        type="button"
                        className="gb-remove"
                        aria-label={`Remove ${s.label}`}
                        onClick={() => removeExtra(s.id)}
                      >
                        ×
                      </button>
                    ) : null}
                  </div>
                );
              })}
            </div>
          );
        })}
      </div>

      <Suspense fallback={null}>
        <CatalogPanel
          genome={pack.assembly}
          loaded={loadedTracks}
          pending={pendingTrack}
          filterQ={catalogQ}
          onFilterChange={setCatalogQ}
          onQuery={queryCatalog}
          onHide={hideCatalog}
        />
      </Suspense>
    </div>
  );
}

function browserStorageKey(pack: LocusPack): string {
  return [
    "ntv3-browser-tracks",
    pack.assembly,
    pack.chrom,
    pack.start,
    pack.end,
  ].join(":");
}

function readPersistedExtraSlots(storageKey: string): TrackSlot[] {
  try {
    const raw = window.localStorage.getItem(storageKey);
    if (!raw) return [];
    const parsed: unknown = JSON.parse(raw);
    if (!Array.isArray(parsed)) return [];
    return parsed.filter(isPersistedTrackSlot);
  } catch {
    return [];
  }
}

function persistExtraSlots(storageKey: string, slots: TrackSlot[]): void {
  try {
    window.localStorage.setItem(storageKey, JSON.stringify(slots));
  } catch {
    // The report still works when browser storage is unavailable.
  }
}

function isPersistedTrackSlot(value: unknown): value is TrackSlot {
  if (!value || typeof value !== "object") return false;
  const slot = value as Partial<TrackSlot>;
  return (
    typeof slot.id === "string" &&
    typeof slot.kind === "string" &&
    typeof slot.label === "string" &&
    typeof slot.group === "string" &&
    typeof slot.ucscTrack === "string" &&
    typeof slot.stackRank === "number" &&
    typeof slot.origin === "string"
  );
}

function kindToMode(kind: TrackKind): FetchSpec["mode"] {
  if (
    kind === "genes" ||
    kind === "ccre" ||
    kind === "ticks" ||
    kind === "signal"
  ) {
    return kind;
  }
  return "intervals";
}

function fetchSpecForSlot(slot: TrackSlot): FetchSpec {
  return {
    id: slot.id,
    track: slot.ucscTrack ?? slot.id,
    label: slot.label,
    mode: kindToMode(slot.kind),
    signed: slot.signed,
  };
}

/**
 * Stack order: annotations → GT/pred pairs (interleaved when pred is on
 * and populated) → saliency → phyloP last.
 * Pred rows are omitted until they exist; then they sit under their GT.
 */
function orderStack(
  slots: TrackSlot[],
  on: Record<string, boolean>,
  populate: Populate,
  extraSlots: TrackSlot[],
): TrackSlot[] {
  const checked = [...slots, ...extraSlots].filter((s) => on[s.id]);
  const assays: string[] = [];
  for (const s of checked) {
    if (s.assay && !assays.includes(s.assay)) assays.push(s.assay);
  }
  const pairIds = new Set<string>();
  const pairs: TrackSlot[] = [];
  for (const assay of assays) {
    const gt = checked.find((s) => s.kind === "gt" && s.assay === assay);
    const pred = checked.find((s) => s.kind === "pred" && s.assay === assay);
    if (gt) {
      pairs.push(gt);
      pairIds.add(gt.id);
    }
    if (pred && populate.pred) {
      pairs.push(pred);
      pairIds.add(pred.id);
    }
  }
  const rest = checked.filter((s) => !pairIds.has(s.id) && s.kind !== "pred");
  rest.sort((a, b) => a.stackRank - b.stackRank);
  const before = rest.filter((s) => s.stackRank < 50);
  const after = rest.filter((s) => s.stackRank >= 50);
  return [...before, ...pairs, ...after];
}

function stackHeight(slots: TrackSlot[]): number {
  let h = 36;
  for (const s of slots) h += rowHeight(s);
  return h + 28;
}

function rowReady(
  slot: TrackSlot,
  pack: LocusPack,
  populate: Populate,
  live: LiveBundle | null,
  bwSignals: Record<string, SignalTrack>,
): boolean {
  if (slot.kind === "gt") {
    if (!populate.gt) return false;
    return slot.id in bwSignals || pack.observed.some((track) => track.id === slot.id);
  }
  if (slot.ucscTrack) return live != null && slot.id in live.ok;
  if (slot.kind === "pred") return populate.pred;
  if (slot.kind === "saliency") return populate.saliency;
  return populate.annotations;
}

function TrackStack({
  pack,
  viewChrom,
  live,
  bwSignals,
  bwErrors,
  viewStart,
  viewEnd,
  px,
  slots,
  populate,
  onColorSlot,
  colors,
}: {
  pack: LocusPack;
  viewChrom: string;
  live: LiveBundle | null;
  bwSignals: Record<string, SignalTrack>;
  bwErrors: Record<string, string>;
  viewStart: number;
  viewEnd: number;
  px: (bp: number) => number;
  slots: TrackSlot[];
  populate: Populate;
  onColorSlot: (slot: TrackSlot, top: number, left: number) => void;
  colors: Record<string, string>;
}) {
  let y = 8;
  const rows = [];
  for (const slot of slots) {
    const ready = rowReady(slot, pack, populate, live, bwSignals);
    const error = bwErrors[slot.id] ?? live?.errors[slot.id];
    const h = rowHeight(slot);
    rows.push(
      <g key={slot.id}>
        <g
          className={isColorable(slot) ? "gb-lab-hit recolor" : "gb-lab-hit"}
          data-track-id={slot.id}
          onPointerDown={(e) => {
            e.stopPropagation();
          }}
          onClick={(e) => {
            e.stopPropagation();
            if (isColorable(slot)) {
              const box = (
                e.currentTarget as SVGGElement
              ).getBoundingClientRect();
              onColorSlot(slot, box.bottom + 4, box.left);
            }
          }}
        >
          <title>
            {isColorable(slot) ? `Color ${slot.label}` : slot.label}
          </title>
          <rect x={0} y={y} width={GUT - 2} height={h} fill="transparent" />
          <text
            x={GUT - 10}
            y={y + (slot.kind === "genes" ? 18 : 13)}
            textAnchor="end"
            className="gb-lab"
            style={
              isColorable(slot)
                ? { fill: colors[slot.id] ?? defaultColor(slot) }
                : undefined
            }
          >
            {fitLabel(slot.label, GUT - 14)}
          </text>
        </g>
        {ready ? (
          <RowBody
            slot={slot}
            pack={pack}
            live={live}
            bwSignals={bwSignals}
            y={y}
            h={h}
            viewStart={viewStart}
            viewEnd={viewEnd}
            px={px}
            fill={colors[slot.id] ?? defaultColor(slot)}
          />
        ) : (
          <EmptyRow
            y={y}
            h={h}
            message={error ? `failed · ${error}` : `waiting · ${slot.origin}`}
          />
        )}
      </g>,
    );
    y += h + 4;
  }
  rows.push(
    <Axis
      key="axis"
      chrom={viewChrom}
      viewStart={viewStart}
      viewEnd={viewEnd}
      y={y + 4}
      px={px}
    />,
  );
  return <g>{rows}</g>;
}

function RowBody({
  slot,
  pack,
  live,
  bwSignals,
  y,
  h,
  viewStart,
  viewEnd,
  px,
  fill,
}: {
  slot: TrackSlot;
  pack: LocusPack;
  live: LiveBundle | null;
  bwSignals: Record<string, SignalTrack>;
  y: number;
  h: number;
  viewStart: number;
  viewEnd: number;
  px: (bp: number) => number;
  fill: string;
}) {
  if (slot.kind === "genes") {
    const genes = live?.genes[slot.id] ?? pack.genes[slot.id] ?? [];
    if (!genes.length && slot.empty) return <EmptyRow y={y} h={h} message={slot.empty} />;
    return <GeneLane genes={genes} y={y} px={px} viewStart={viewStart} viewEnd={viewEnd} />;
  }
  if (slot.kind === "ccre") {
    return (
      <CcreLane
        ccres={live?.ccres ?? []}
        y={y}
        px={px}
        viewStart={viewStart}
        viewEnd={viewEnd}
      />
    );
  }
  if (slot.kind === "intervals") {
    return (
      <IntervalLane
        items={live?.intervals[slot.id] ?? []}
        y={y}
        px={px}
        viewStart={viewStart}
        viewEnd={viewEnd}
        fill="#5b6b7c"
      />
    );
  }
  if (slot.kind === "ticks") {
    return (
      <TickLane
        items={live?.intervals[slot.id] ?? []}
        y={y}
        px={px}
        viewStart={viewStart}
        viewEnd={viewEnd}
        color={slot.id === "clinvar" ? "#a81f26" : "#1e4b8c"}
      />
    );
  }
  const track = resolveSignal(slot, pack, live, bwSignals);
  if (!track) return null;
  if (slot.kind === "saliency") {
    return (
      <DynseqLane
        track={track}
        y={y}
        h={h}
        viewStart={viewStart}
        viewEnd={viewEnd}
        winStart={track.spanStart ?? pack.start}
        winEnd={track.spanEnd ?? pack.end}
        fill={fill}
      />
    );
  }
  return (
    <SignalLane
      track={track}
      y={y}
      h={h}
      viewStart={viewStart}
      viewEnd={viewEnd}
      winStart={track.spanStart ?? pack.start}
      winEnd={track.spanEnd ?? pack.end}
      signed={Boolean(track.signed || slot.signed)}
      fill={fill}
    />
  );
}

function resolveSignal(
  slot: TrackSlot,
  pack: LocusPack,
  live: LiveBundle | null,
  bwSignals: Record<string, SignalTrack>,
): SignalTrack | null {
  if (bwSignals[slot.id]) return bwSignals[slot.id] ?? null;
  if (live?.signals[slot.id]) return live.signals[slot.id] ?? null;
  if (slot.ucscTrack) return null;
  if (slot.kind === "gt") {
    return pack.observed.find((track) => track.id === slot.id) ?? null;
  }
  const saliency = pack.saliency.find((track) => track.id === slot.id);
  if (saliency) return saliency;
  const pr = pack.predicted.find((t) => t.id === slot.id);
  return pr ?? null;
}

function IntervalLane({
  items,
  y,
  px,
  viewStart,
  viewEnd,
  fill,
}: {
  items: IntervalFeat[];
  y: number;
  px: (bp: number) => number;
  viewStart: number;
  viewEnd: number;
  fill: string;
}) {
  return (
    <g>
      {items
        .filter((c) => c.end > viewStart && c.start < viewEnd)
        .map((c, i) => (
          <rect
            key={`${c.start}-${c.name}-${i}`}
            x={px(Math.max(c.start, viewStart))}
            y={y + 4}
            width={Math.max(
              1,
              px(Math.min(c.end, viewEnd)) - px(Math.max(c.start, viewStart)),
            )}
            height={12}
            fill={fill}
            opacity={0.85}
          >
            <title>{c.name}</title>
          </rect>
        ))}
    </g>
  );
}

function TickLane({
  items,
  y,
  px,
  viewStart,
  viewEnd,
  color,
}: {
  items: IntervalFeat[];
  y: number;
  px: (bp: number) => number;
  viewStart: number;
  viewEnd: number;
  color: string;
}) {
  return (
    <g>
      {items
        .filter((c) => c.end > viewStart && c.start < viewEnd)
        .slice(0, 800)
        .map((c, i) => {
          const x = px((c.start + c.end) / 2);
          return (
            <line
              key={`${c.start}-${c.name}-${i}`}
              x1={x}
              x2={x}
              y1={y + 3}
              y2={y + 17}
              stroke={color}
              strokeWidth={1}
            >
              <title>{c.name}</title>
            </line>
          );
        })}
    </g>
  );
}

function GeneLane({
  genes,
  y,
  px,
  viewStart,
  viewEnd,
}: {
  genes: GeneFeat[];
  y: number;
  px: (bp: number) => number;
  viewStart: number;
  viewEnd: number;
}) {
  const visible = genes.filter((g) => g.end > viewStart && g.start < viewEnd);
  const labelX = placeGeneLabels(
    visible.map((g) => ({
      name: g.name,
      x0: px(Math.max(g.start, viewStart)),
      x1: px(Math.min(g.end, viewEnd)),
      span: g.end - g.start,
    })),
    GUT,
    GUT + PW,
  );
  return (
    <g>
      {visible.map((g, gi) => {
        const x0 = px(Math.max(g.start, viewStart));
        const x1 = px(Math.min(g.end, viewEnd));
        const lx = labelX[gi];
        const midY = y + 16;
        const dir = g.strand === "+" ? 1 : -1;
        return (
          <g key={`${g.name}-${g.start}`}>
            <title>{g.name}</title>
            <line x1={x0} y1={midY} x2={x1} y2={midY} stroke="#1e4b8c" strokeWidth={1} />
            {chevrons(x0, x1, midY, dir)}
            {g.exons
              .flatMap((e) => exonSegments(e, g.cds))
              .filter((e) => e.end > viewStart && e.start < viewEnd)
              .map((e, i) => (
                <rect
                  key={`${g.name}-${e.start}-${e.end}-${i}`}
                  x={px(Math.max(e.start, viewStart))}
                  y={midY - (e.thick ? 5 : 2.5)}
                  width={Math.max(
                    1,
                    px(Math.min(e.end, viewEnd)) - px(Math.max(e.start, viewStart)),
                  )}
                  height={e.thick ? 10 : 5}
                  fill="#1e4b8c"
                />
              ))}
            {lx == null ? null : (
              <text x={lx} y={midY - 8} className="gb-gene" fontStyle="italic">
                {g.name}
              </text>
            )}
          </g>
        );
      })}
    </g>
  );
}

function RangeLabel({ y, h, text }: { y: number; h: number; text: string }) {
  return (
    <text x={GUT - 10} y={y + h - 2} textAnchor="end" className="gb-range">
      {text}
    </text>
  );
}

function chevrons(x0: number, x1: number, y: number, dir: number) {
  const nodes = [];
  const left = Math.min(x0, x1);
  const right = Math.max(x0, x1);
  for (let x = left + 10; x < right - 8; x += 22) {
    const tip = x + 4 * dir;
    nodes.push(
      <path
        key={x}
        d={`M${x - 3 * dir},${y - 3.5} L${tip},${y} L${x - 3 * dir},${y + 3.5}`}
        fill="none"
        stroke="#1e4b8c"
        strokeWidth={1}
      />,
    );
  }
  return nodes;
}

function CcreLane({
  ccres,
  y,
  px,
  viewStart,
  viewEnd,
}: {
  ccres: CcreFeat[];
  y: number;
  px: (bp: number) => number;
  viewStart: number;
  viewEnd: number;
}) {
  return (
    <g>
      {ccres
        .filter((c) => c.end > viewStart && c.start < viewEnd)
        .map((c, i) => (
          <rect
            key={`${c.name}-${c.start}-${i}`}
            x={px(Math.max(c.start, viewStart))}
            y={y + 4}
            width={Math.max(
              1,
              px(Math.min(c.end, viewEnd)) - px(Math.max(c.start, viewStart)),
            )}
            height={12}
            fill={CCRE[c.klass]}
          >
            <title>{`${c.name} · ${c.klass}`}</title>
          </rect>
        ))}
    </g>
  );
}

function SignalLane({
  track,
  y,
  h,
  viewStart,
  viewEnd,
  winStart,
  winEnd,
  signed,
  fill,
}: {
  track: SignalTrack;
  y: number;
  h: number;
  viewStart: number;
  viewEnd: number;
  winStart: number;
  winEnd: number;
  signed: boolean;
  fill: string;
}) {
  const stableMax = useMemo(() => absMax(track.values), [track.values]);
  const { dFill, dStroke, maxAbs } = useMemo(
    () =>
      pathForWindow(
        track.values,
        winStart,
        winEnd,
        viewStart,
        viewEnd,
        y,
        h,
        signed,
        signed ? stableMax : undefined,
      ),
    [track.values, winStart, winEnd, viewStart, viewEnd, y, h, signed, stableMax],
  );
  const outside = viewEnd <= winStart || viewStart >= winEnd;
  const zeroY = signed ? y + h / 2 : y + h - 6;
  if (outside && track.id.startsWith("pr-")) {
    return (
      <g>
        <line x1={GUT} y1={zeroY} x2={GUT + PW} y2={zeroY} stroke="#e8ebf0" />
        <text x={GUT + 8} y={y + h / 2 + 3} className="gb-empty">
          outside NTv3 predicted span
        </text>
        <RangeLabel y={y} h={h} text="—" />
      </g>
    );
  }
  return (
    <g>
      <line x1={GUT} y1={zeroY} x2={GUT + PW} y2={zeroY} stroke="#e8ebf0" />
      <path d={dFill} fill={fill} fillOpacity={signed ? 0.14 : 0.18} />
      <path d={dStroke} fill="none" stroke={fill} strokeWidth={0.9} />
      <RangeLabel
        y={y}
        h={h}
        text={
          signed
            ? `±${formatScale(maxAbs)}${track.unit ? ` ${track.unit}` : ""}`
            : `0–${formatScale(maxAbs)}${track.unit ? ` ${track.unit}` : ""}`
        }
      />
    </g>
  );
}

function DynseqLane({
  track,
  y,
  h,
  viewStart,
  viewEnd,
  winStart,
  winEnd,
  fill,
}: {
  track: SignalTrack;
  y: number;
  h: number;
  viewStart: number;
  viewEnd: number;
  winStart: number;
  winEnd: number;
  fill: string;
}) {
  const span = Math.max(1, viewEnd - viewStart);
  const pxPerBp = PW / span;
  const seq = track.seq;
  const lo = Math.max(winStart, Math.floor(viewStart));
  const hi = Math.min(winEnd, Math.ceil(viewEnd));
  const showLetters = Boolean(seq) && pxPerBp >= 4 && hi - lo <= 4000;
  if (!showLetters || !seq) {
    return (
      <SignalLane
        track={track}
        y={y}
        h={h}
        viewStart={viewStart}
        viewEnd={viewEnd}
        winStart={winStart}
        winEnd={winEnd}
        signed
        fill={fill}
      />
    );
  }
  const winSpan = Math.max(1, winEnd - winStart);
  const n = track.values.length;
  const baseY = y + h / 2;
  const usable = (h - 12) / 2;
  const i0 = Math.max(0, Math.floor(((lo - winStart) / winSpan) * n));
  const i1 = Math.min(n, Math.ceil(((hi - winStart) / winSpan) * n));
  const maxAbs = Math.max(
    1e-9,
    ...track.values.slice(i0, i1).map((v) => Math.abs(v)),
  );
  const glyph = 10;
  const letters = [];
  for (let bp = lo; bp < hi; bp++) {
    const idx = bp - winStart;
    if (idx < 0 || idx >= seq.length) continue;
    const base = seq[idx]?.toUpperCase() ?? "N";
    const vi = Math.min(n - 1, Math.max(0, Math.floor((idx / winSpan) * n)));
    const score = track.values[vi] ?? 0;
    const letterH = Math.max(2, (Math.abs(score) / maxAbs) * usable);
    const x = GUT + ((bp + 0.5 - viewStart) / span) * PW;
    const color = NUC_COLOR[base] ?? "#6b7482";
    const sx = pxPerBp / (glyph * 0.72);
    const sy = (score >= 0 ? 1 : -1) * (letterH / glyph);
    letters.push(
      <g key={bp} transform={`translate(${x} ${baseY}) scale(${sx},${sy})`}>
        <text
          textAnchor="middle"
          dominantBaseline="ideographic"
          fontSize={glyph}
          fontFamily="Arial, sans-serif"
          fontWeight={700}
          fill={color}
        >
          {base}
        </text>
      </g>,
    );
  }
  return (
    <g>
      <line x1={GUT} y1={baseY} x2={GUT + PW} y2={baseY} stroke="#e8ebf0" />
      {letters}
      <RangeLabel y={y} h={h} text={`±${formatScale(maxAbs)}`} />
    </g>
  );
}

function absMax(values: number[]): number {
  let m = 1e-6;
  for (const v of values) {
    const a = Math.abs(v);
    if (a > m) m = a;
  }
  return m;
}

function pathForWindow(
  values: number[],
  winStart: number,
  winEnd: number,
  viewStart: number,
  viewEnd: number,
  y: number,
  h: number,
  signed: boolean,
  scaleMax?: number,
): { dFill: string; dStroke: string; maxAbs: number } {
  const n = values.length;
  const winSpan = Math.max(1, winEnd - winStart);
  const viewSpan = Math.max(1, viewEnd - viewStart);
  const pad = 3;
  const baseY = y + h - 6;
  const usable = h - 12;
  const cols = PW;
  const colVals = Array.from({ length: cols }, () => 0);
  for (let c = 0; c < cols; c++) {
    const bp0 = viewStart + (c / cols) * viewSpan;
    const bp1 = viewStart + ((c + 1) / cols) * viewSpan;
    const j0 = Math.max(0, Math.floor(((bp0 - winStart) / winSpan) * n));
    const j1 = Math.min(n, Math.max(j0 + 1, Math.ceil(((bp1 - winStart) / winSpan) * n)));
    if (signed) {
      let pick = 0;
      let pickAbs = 0;
      for (let j = j0; j < j1; j++) {
        const v = values[j] ?? 0;
        const a = Math.abs(v);
        if (a > pickAbs) {
          pickAbs = a;
          pick = v;
        }
      }
      colVals[c] = pick;
    } else {
      let pick = 0;
      for (let j = j0; j < j1; j++) {
        const v = values[j] ?? 0;
        if (v > pick) pick = v;
      }
      colVals[c] = pick;
    }
  }
  let windowMax = 1e-6;
  for (const v of colVals) {
    const a = Math.abs(v);
    if (a > windowMax) windowMax = a;
  }
  const maxAbs = scaleMax ?? windowMax;
  const zeroY = signed ? baseY - usable / 2 : baseY;
  const yOf = (v: number) => {
    if (signed) return zeroY - (v / maxAbs) * (usable / 2 - pad);
    return baseY - pad - (v / maxAbs) * (usable - pad);
  };
  let dFill = "";
  let dStroke = "";
  for (let c = 0; c < cols; c++) {
    const x0 = GUT + c;
    const x1 = GUT + c + 1;
    const yy = yOf(colVals[c] ?? 0);
    dFill += `M${x0},${zeroY} L${x0},${yy} L${x1},${yy} L${x1},${zeroY} Z`;
    dStroke += dStroke ? ` L${x0},${yy}` : `M${x0},${yy}`;
    dStroke += ` L${x1},${yy}`;
  }
  return { dFill, dStroke, maxAbs };
}

function EmptyRow({ y, h, message }: { y: number; h: number; message: string }) {
  return (
    <g>
      <rect
        x={GUT}
        y={y + 2}
        width={PW}
        height={h - 6}
        fill="none"
        stroke="#cdd4dc"
        strokeDasharray="4 3"
      />
      <text x={GUT + 8} y={y + h / 2 + 3} className="gb-empty">
        {message}
      </text>
    </g>
  );
}

function Axis({
  chrom,
  viewStart,
  viewEnd,
  y,
  px,
}: {
  chrom: string;
  viewStart: number;
  viewEnd: number;
  y: number;
  px: (bp: number) => number;
}) {
  const span = viewEnd - viewStart;
  const step = span > 20000 ? 10000 : span > 5000 ? 2000 : span > 1000 ? 500 : 100;
  const first = Math.ceil(viewStart / step) * step;
  const ticks = [];
  for (let t = first; t < viewEnd; t += step) ticks.push(t);
  // 9 px mono is about 5.4 px a character. A tick label that would run into
  // the chromosome name at the left keeps its tick and drops its text.
  const chromRight = GUT + chrom.length * 5.4 + 6;
  return (
    <g>
      <line x1={GUT} y1={y} x2={GUT + PW} y2={y} stroke="#3d4653" strokeWidth={0.8} />
      {ticks.map((t) => {
        const label = `${(t / 1000).toFixed(span > 5000 ? 0 : 1)} kb`;
        return (
          <g key={t}>
            <line x1={px(t)} y1={y} x2={px(t)} y2={y + 5} stroke="#3d4653" strokeWidth={0.8} />
            {px(t) - (label.length * 5.4) / 2 < chromRight ? null : (
              <text x={px(t)} y={y + 16} textAnchor="middle" className="gb-tick">
                {label}
              </text>
            )}
          </g>
        );
      })}
      <text x={GUT} y={y + 16} className="gb-tick">
        {chrom}
      </text>
    </g>
  );
}

function Ideogram({
  chrom,
  chromSize,
  viewStart,
  viewEnd,
}: {
  chrom: string;
  chromSize: number;
  viewStart: number;
  viewEnd: number;
}) {
  const w = 640;
  const x = (bp: number) => (bp / chromSize) * w;
  return (
    <svg className="gb-ideo" viewBox={`0 0 780 28`} aria-hidden="true">
      <text x="0" y="18" className="gb-lab">
        {chrom}
      </text>
      <rect x="80" y="10" width={w} height="8" rx="4" fill="#d8dee6" />
      <rect x={80 + x(50_000_000)} y="10" width={x(8_000_000)} height="8" fill="#9aa3ad" />
      <rect
        x={80 + x(viewStart)}
        y="6"
        width={Math.max(2, x(viewEnd) - x(viewStart))}
        height="16"
        fill="none"
        stroke="#c41e3a"
        strokeWidth="1.4"
      />
    </svg>
  );
}
