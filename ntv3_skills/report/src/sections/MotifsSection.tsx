import { Section } from "../components/Section";
import type { Report } from "../schema/report";

type Props = { report: Report };

/** Optional explicit JASPAR motif occurrence evidence. */
export function MotifsSection({ report }: Props) {
  const motifs = report.motifs;
  if (!motifs) return null;

  return (
    <Section id="g-motifs" title="Motif occurrences">
      <p className="interp">{motifs.caveat}</p>
      <p>
        {motifs.scanner} {motifs.scannerVersion} found{" "}
        <strong>{motifs.hits.length.toLocaleString()}</strong> occurrence
        {motifs.hits.length === 1 ? "" : "s"} at p ≤ {motifs.threshold} for{" "}
        {motifs.attributionMethod} on <strong>{motifs.trackLabel}</strong>.
      </p>
      {motifs.sources.map((source) => (
        <p key={source.matrixId}>
          <a href={source.url}>
            {source.matrixId} · {source.name}
          </a>
        </p>
      ))}
      {motifs.hits.length ? (
        <table className="ledger">
          <thead>
            <tr>
              <th>Motif</th>
              <th>Interval</th>
              <th>Strand</th>
              <th>p-value</th>
              <th>|attribution| sum</th>
            </tr>
          </thead>
          <tbody>
            {motifs.hits.map((hit, index) => (
              <tr key={`${hit.matrixId}-${hit.start}-${hit.strand}-${index}`}>
                <td>{hit.matrixName}</td>
                <td>
                  {report.locus.chrom}:{hit.start.toLocaleString()}–
                  {hit.end.toLocaleString()}
                </td>
                <td>{hit.strand}</td>
                <td>{hit.pValue.toExponential(2)}</td>
                <td>{hit.attributionAbsoluteSum.toExponential(2)}</td>
              </tr>
            ))}
          </tbody>
        </table>
      ) : (
        <p>No motif occurrences passed the approved threshold.</p>
      )}
    </Section>
  );
}
