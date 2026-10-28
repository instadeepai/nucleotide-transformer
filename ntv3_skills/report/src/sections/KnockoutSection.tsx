import { Section } from "../components/Section";
import type { Report } from "../schema/report";

type Props = { report: Report };

/** Optional approved interval-shuffle comparison. */
export function KnockoutSection({ report }: Props) {
  const knockout = report.knockout;
  if (!knockout) return null;
  const interval = knockout.interval;

  return (
    <Section id="g-knockout" title="Interval shuffle">
      <p className="interp">{knockout.caveat}</p>
      <p>
        Approved interval{" "}
        <strong>
          {interval.chrom}:{interval.start.toLocaleString()}–
          {interval.end.toLocaleString()}
        </strong>{" "}
        · seed {knockout.seed.toLocaleString()}
        {knockout.sourceEvidence ? ` · evidence ${knockout.sourceEvidence}` : ""}
      </p>
      <table className="ledger">
        <thead>
          <tr>
            <th>Model head</th>
            <th>Signed max shuffle effect</th>
            <th>Mean shuffle effect</th>
          </tr>
        </thead>
        <tbody>
          {knockout.perTrack.map((row) => (
            <tr key={row.modelHead}>
              <td>
                <code>{row.modelHead}</code>
              </td>
              <td>{row.signedMaxEffect.toFixed(3)}</td>
              <td>{row.meanEffect.toFixed(3)}</td>
            </tr>
          ))}
        </tbody>
      </table>
      <figure className="saliency-logo">
        <img src={knockout.figure.png} alt="Wild-type and shuffled track comparison" />
        <figcaption>
          Wild type versus shuffled sequence ·{" "}
          <a href={knockout.figure.pdf}>vector PDF</a>
        </figcaption>
      </figure>
    </Section>
  );
}
