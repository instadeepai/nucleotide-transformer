import { useState } from "react";
import { Section } from "../components/Section";
import type { Report } from "../schema/report";

type Props = { report: Report };

export function ClaimsSection({ report }: Props) {
  const [showSuppressed, setShowSuppressed] = useState(false);
  const nSup = report.claims.filter((c) => !c.keep).length;
  const rows = showSuppressed
    ? report.claims
    : report.claims.filter((c) => c.keep);

  return (
    <Section id="g-claims" title="Claims">
      {nSup > 0 ? (
        <div className="figtools">
          <button
            type="button"
            className="tgl"
            aria-pressed={showSuppressed}
            onClick={() => setShowSuppressed((v) => !v)}
          >
            {showSuppressed
              ? `hide suppressed (${nSup})`
              : `show suppressed (${nSup})`}
          </button>
        </div>
      ) : null}
      <table className="ledger">
        <thead>
          <tr>
            <th className="label">ID</th>
            <th className="label">Claim</th>
            <th className="label">Support</th>
            <th className="label">Status</th>
          </tr>
        </thead>
        <tbody>
          {rows.map((c) => (
            <tr key={c.id} className={c.keep ? undefined : "is-sup"}>
              <td className="label">{c.id}</td>
              <td className="label">{c.text}</td>
              <td>
                {c.support}
                {c.command ? <span className="cmd">$ {c.command}</span> : null}
              </td>
              <td className="label">
                <span className={`flag ${c.keep ? "keep" : "sup"}`}>
                  {c.keep ? "retained" : "withheld"}
                </span>
              </td>
            </tr>
          ))}
        </tbody>
      </table>
    </Section>
  );
}
