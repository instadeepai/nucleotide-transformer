import { Section } from "../components/Section";
import type { Report } from "../schema/report";

type Props = { report: Report };

function correlation(value: number | null): string {
  return value == null ? "undefined" : value.toFixed(3);
}

export function AgreementSection({ report }: Props) {
  return (
    <Section id="g-agreement" title="Prediction agreement">
      <p className="interp">
        Pearson and Spearman use finite paired model bins. They compare shape, not
        absolute scale.
      </p>
      <table className="ledger">
        <thead>
          <tr>
            <th>Track</th>
            <th>Model / GT</th>
            <th>Valid bins</th>
            <th>Pearson</th>
            <th>Spearman</th>
          </tr>
        </thead>
        <tbody>
          {report.trackAgreement.map((row) => (
            <tr key={`${row.modelHead}-${row.gtAccession}`}>
              <td>{row.trackQuery}</td>
              <td>
                <code>{row.modelHead}</code>
                <br />
                <code>
                  {row.gtAccession} / {row.fileAccession}
                </code>
              </td>
              <td>{row.nValid.toLocaleString()}</td>
              <td title={row.caveat}>{correlation(row.pearson)}</td>
              <td title={row.caveat}>{correlation(row.spearman)}</td>
            </tr>
          ))}
        </tbody>
      </table>
    </Section>
  );
}
