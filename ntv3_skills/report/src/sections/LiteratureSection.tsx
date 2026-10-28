import { Section } from "../components/Section";
import type { Report } from "../schema/report";

type Props = { report: Report };

/** Optional validated primary-source literature records. */
export function LiteratureSection({ report }: Props) {
  const literature = report.literature;
  if (!literature) return null;

  return (
    <Section id="g-literature" title="Literature">
      <p className="interp">
        Primary-source records inspected for this run. Search queries and retrieval
        dates are retained for reproducibility.
      </p>
      <table className="ledger">
        <thead>
          <tr>
            <th>Source</th>
            <th>Relevance</th>
            <th>Retrieved</th>
          </tr>
        </thead>
        <tbody>
          {literature.records.map((record) => (
            <tr key={record.citekey}>
              <td>
                <a href={record.sourceUrl}>{record.title}</a>
                <br />
                {record.authors.join(", ")} ({record.year})
                <br />
                <code>{record.citekey}</code>
                {record.doi ? (
                  <>
                    {" · "}
                    <a href={`https://doi.org/${record.doi}`}>{record.doi}</a>
                  </>
                ) : null}
              </td>
              <td>{record.relevance}</td>
              <td>{record.retrievalDate}</td>
            </tr>
          ))}
        </tbody>
      </table>
      <p>
        <a href={literature.bib}>Download BibTeX</a>
      </p>
    </Section>
  );
}
