import { explainLaneRows } from "../browser/lanes";
import { Section } from "../components/Section";
import type { Report } from "../schema/report";

function fmt(n: number): string {
  return n.toLocaleString("en-US");
}

function span(s: { chrom: string; start: number; end: number; bp: number }) {
  return `${s.chrom}:${fmt(s.start)}-${fmt(s.end)} · ${fmt(s.bp)} bp`;
}

type Props = { report: Report };

/** Which model, window, crop, and head produced each number, from provenance. */
export function ExplainSection({ report }: Props) {
  const explain = report.explain;
  if (!explain) return null;
  const cropPercent = Math.round((explain.crop.bp / explain.window.bp) * 100);
  const laneRows = explainLaneRows(
    explain.tracks,
    report.tracks.assays,
    explain.letters,
  );

  return (
    <Section id="g-explain" title="Explain">
      <p className="interp interp--muted">
        Read from {explain.sources.join(", ")}. Nothing here is recomputed.
      </p>
      <dl className="meta-dl">
        <div>
          <dt>Model</dt>
          <dd>
            <code>{explain.model}</code>
            {explain.manifestModel && explain.manifestModel !== explain.model ? (
              <>
                {" "}
                (manifest: <code>{explain.manifestModel}</code>)
              </>
            ) : null}
          </dd>
        </div>
        <div>
          <dt>Species</dt>
          <dd>
            {explain.species} · {explain.assembly} conditioning
          </dd>
        </div>
        <div>
          <dt>Window</dt>
          <dd>
            <code>{span(explain.window)}</code> read by the model
          </dd>
        </div>
        <div>
          <dt>Crop</dt>
          <dd>
            <code>{span(explain.crop)}</code> predicted (central {cropPercent}%)
          </dd>
        </div>
        <div>
          <dt>Bins</dt>
          <dd>
            {fmt(explain.nBins)} × {fmt(explain.binSize)} bp
          </dd>
        </div>
        <div>
          <dt>Device</dt>
          <dd>{explain.device}</dd>
        </div>
      </dl>
      <table className="ledger">
        <thead>
          <tr>
            <th>Lane</th>
            <th>Model head</th>
            <th>Unit</th>
            <th>Ground truth</th>
          </tr>
        </thead>
        <tbody>
          {laneRows.map((row) => (
            <tr key={row.lane}>
              <td>{row.lane}</td>
              <td>{row.modelHead ? <code>{row.modelHead}</code> : "—"}</td>
              <td>{row.unit}</td>
              <td>{row.groundTruth ? <code>{row.groundTruth}</code> : "—"}</td>
            </tr>
          ))}
        </tbody>
      </table>
      {explain.compare ? (
        <p className="interp">
          <span className="summary-kicker">Compare</span>
          <code>{explain.compare.label}</code> used{" "}
          <code>{explain.compare.checkpoint}</code> over{" "}
          <code>{span(explain.compare.window)}</code>; its deltas cover{" "}
          <code>{span(explain.compare.crop)}</code> at{" "}
          {fmt(explain.compare.binSize)} bp bins.
        </p>
      ) : null}
      {explain.letters.map((row) => (
        <p className="interp" key={row.modelHead}>
          <span className="summary-kicker">Letters</span>
          Saliency for {row.label} (<code>{row.modelHead}</code>) over{" "}
          <code>
            {explain.window.chrom}:{fmt(row.start)}-{fmt(row.end)}
          </code>
          . Letters appear in the browser when zoomed in. Saliency is not a causal
          effect.
        </p>
      ))}
    </Section>
  );
}
