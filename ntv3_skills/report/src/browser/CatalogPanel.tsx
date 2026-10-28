import { useMemo, useState } from "react";
import type { CatalogTrack } from "../data/ucscGroups";
import { UCSC_GROUP_LABEL, UCSC_GROUP_ORDER } from "../data/ucscGroups";
import { UCSC_CATALOG } from "../data/ucscCatalog";

type Props = {
  genome: string;
  loaded: Set<string>;
  pending: string | null;
  filterQ: string;
  onFilterChange: (q: string) => void;
  onQuery: (track: string, label: string, type: string, group: string) => void;
  onHide: (track: string) => void;
};

export function CatalogPanel({
  genome,
  loaded,
  pending,
  filterQ,
  onFilterChange,
  onQuery,
  onHide,
}: Props) {
  const q = filterQ;
  const setQ = onFilterChange;
  const [open, setOpen] = useState<Record<string, boolean>>({});
  const [openComp, setOpenComp] = useState<Record<string, boolean>>({});
  const [custom, setCustom] = useState("");

  const query = q.trim().toLowerCase();
  const hasCatalog = genome === UCSC_CATALOG.genome;
  const tracks = useMemo(() => {
    if (!hasCatalog) return [];
    if (!query) return UCSC_CATALOG.tracks;
    return UCSC_CATALOG.tracks.filter((t) => {
      if (matches(t.label, query) || matches(t.id, query) || matches(t.long, query)) {
        return true;
      }
      return t.leaves.some(
        (leaf) => matches(leaf.label, query) || matches(leaf.id, query),
      );
    });
  }, [hasCatalog, query]);

  const byGroup = useMemo(() => {
    const map = new Map<string, CatalogTrack[]>();
    for (const t of tracks) {
      const list = map.get(t.group) ?? [];
      list.push(t);
      map.set(t.group, list);
    }
    return map;
  }, [tracks]);

  const submitCustom = () => {
    const name = custom.trim();
    if (!name) return;
    onQuery(name, name, "", "other");
  };

  return (
    <section className="gb-catalog">
      <div className="gb-picker-head">
        <span>UCSC catalog</span>
      </div>
      <form
        className="gb-catalog-query"
        onSubmit={(e) => {
          e.preventDefault();
          submitCustom();
        }}
      >
        <input
          value={custom}
          onChange={(e) => setCustom(e.target.value)}
          placeholder="any UCSC track= name (e.g. gnomadVariantsV4.1)"
          aria-label="UCSC track name"
        />
        <button type="submit" disabled={!custom.trim() || pending != null}>
          Fetch track
        </button>
      </form>
      <input
        className="gb-catalog-filter"
        value={q}
        onChange={(e) => setQ(e.target.value)}
        placeholder="filter catalog…"
        aria-label="Filter UCSC catalog"
      />
      {!hasCatalog ? (
        <p className="muted">
          The browsable catalog snapshot is for {UCSC_CATALOG.genome}. Enter a
          documented {genome} track name above.
        </p>
      ) : null}
      {UCSC_GROUP_ORDER.map((group) => {
        const items = byGroup.get(group) ?? [];
        if (!items.length) return null;
        const shown = open[group] || Boolean(query);
        return (
          <div key={group} className="gb-cat-group">
            <button
              type="button"
              className="gb-cat-toggle"
              onClick={() => setOpen((p) => ({ ...p, [group]: !p[group] }))}
            >
              {shown ? "▾" : "▸"} {UCSC_GROUP_LABEL[group] ?? group}
              <span className="muted">{items.length}</span>
            </button>
            {shown
              ? items.map((t) => (
                  <CatalogRow
                    key={t.id}
                    track={t}
                    loaded={loaded}
                    pending={pending}
                    open={Boolean(openComp[t.id] || query)}
                    onOpen={() =>
                      setOpenComp((p) => ({ ...p, [t.id]: !p[t.id] }))
                    }
                    onQuery={onQuery}
                    onHide={onHide}
                  />
                ))
              : null}
          </div>
        );
      })}
    </section>
  );
}

function CatalogRow({
  track,
  loaded,
  pending,
  open,
  onOpen,
  onQuery,
  onHide,
}: {
  track: CatalogTrack;
  loaded: Set<string>;
  pending: string | null;
  open: boolean;
  onOpen: () => void;
  onQuery: (track: string, label: string, type: string, group: string) => void;
  onHide: (track: string) => void;
}) {
  const isOn = loaded.has(track.id);
  const busy = pending === track.id;
  return (
    <div className="gb-cat-row">
      <label>
        <input
          type="checkbox"
          checked={isOn}
          disabled={busy}
          onChange={() => {
            if (isOn) onHide(track.id);
            else onQuery(track.id, track.label, track.type, track.group);
          }}
        />
        {track.label}
      </label>
      <code>{track.id}</code>
      {track.nLeaves > 0 ? (
        <button type="button" className="linkish" onClick={onOpen}>
          {track.nLeaves.toLocaleString()} leaves
        </button>
      ) : (
        <span className="chip">{track.type || "track"}</span>
      )}
      {open && track.nLeaves > 0 ? (
        <div className="gb-cat-leaves">
          {track.leaves.map((leaf) => (
            <label key={leaf.id}>
              <input
                type="checkbox"
                checked={loaded.has(leaf.id)}
                disabled={pending === leaf.id}
                onChange={() => {
                  if (loaded.has(leaf.id)) onHide(leaf.id);
                  else onQuery(leaf.id, leaf.label, leaf.type, track.group);
                }}
              />
              {leaf.label}
              <code>{leaf.id}</code>
            </label>
          ))}
          {track.nLeaves > track.leaves.length ? (
            <p className="muted">
              Showing {track.leaves.length} of {track.nLeaves.toLocaleString()}.
              Paste a leaf name in Fetch track for the rest.
            </p>
          ) : null}
        </div>
      ) : null}
    </div>
  );
}

function matches(s: string, q: string): boolean {
  return s.toLowerCase().includes(q);
}
