import { Section } from "../components/Section";
import type { Report } from "../schema/report";

type Props = { report: Report };

/** Optional approved Ledidi model-output design result. */
export function DesignSection({ report }: Props) {
  const design = report.design;
  if (!design) return null;
  const interval = design.editableInterval;
  const selected = design.selectedCandidate;

  return (
    <Section id="g-design" title="Sequence design">
      <p className="interp">{design.caveat}</p>
      <p>
        Ledidi optimized <strong>{design.trackLabel}</strong> over{" "}
        <strong>
          {interval.chrom}:{interval.start.toLocaleString()}–
          {interval.end.toLocaleString()}
        </strong>
        . The requested model-output delta was {design.targetDelta.toFixed(3)}. The
        objective averages output bins {design.scoredBins.start.toLocaleString()}–
        {(design.scoredBins.end - 1).toLocaleString()} around peak bin{" "}
        {design.scoredBins.peak.toLocaleString()}.
      </p>
      <table className="ledger">
        <thead>
          <tr>
            <th>Baseline</th>
            <th>Desired</th>
            <th>Selected</th>
            <th>Error</th>
            <th>Edits</th>
          </tr>
        </thead>
        <tbody>
          <tr>
            <td>{design.baselineTarget.toFixed(3)}</td>
            <td>{design.desiredTarget.toFixed(3)}</td>
            <td>{selected.achievedTarget.toFixed(3)}</td>
            <td>{selected.absoluteTargetError.toFixed(3)}</td>
            <td>{selected.editCount.toLocaleString()}</td>
          </tr>
        </tbody>
      </table>
      {selected.edits.length ? (
        <p>
          Selected edits:{" "}
          {selected.edits
            .map((edit) => `${edit.position.toLocaleString()} ${edit.ref}>${edit.alt}`)
            .join(", ")}
        </p>
      ) : (
        <p>No edit was selected by this bounded optimization run.</p>
      )}
      <p className="interp">
        Candidate {selected.rank} of {design.candidateCount} · seed{" "}
        {design.seed.toLocaleString()} · <code>{design.modelHead}</code>
      </p>
      <figure className="saliency-logo">
        <img
          src={design.figure.png}
          alt="Ledidi model-output trajectory and before-versus-after track"
        />
        <figcaption>
          Optimization trajectory and selected-track profile before versus after ·{" "}
          <a href={design.figure.pdf}>vector PDF</a>
        </figcaption>
      </figure>
    </Section>
  );
}
