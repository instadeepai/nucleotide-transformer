import { Section } from "../components/Section";
import type { Report } from "../schema/report";

type Props = { report: Report };

/** Compare ask: one approved SNV, alternate against reference allele. */
export function DisruptionSection({ report }: Props) {
  const disruption = report.disruption;
  if (!disruption) return null;

  return (
    <Section id="g-disruption" title="Compare: single SNV">
      <p className="interp">{disruption.caveat}</p>
      <p>
        <strong>{disruption.variant.label}</strong> produced a signed effect of{" "}
        <strong>{disruption.signedEffect.toFixed(3)}</strong>.{" "}
        {disruption.percentile == null ? (
          "No null SNVs were run, so the effect is not ranked."
        ) : (
          <>
            Its absolute effect ranked at the{" "}
            <strong>{disruption.percentile.toFixed(1)}th percentile</strong> against{" "}
            {disruption.nullN.toLocaleString()} seeded null SNVs.
          </>
        )}
      </p>
      <table className="ledger">
        <thead>
          <tr>
            <th>Track</th>
            <th>Model head</th>
            <th>Signed effect</th>
            <th>Null percentile</th>
          </tr>
        </thead>
        <tbody>
          {disruption.perTrack.map((row) => (
            <tr key={row.modelHead}>
              <td>{row.track}</td>
              <td>
                <code>{row.modelHead}</code>
              </td>
              <td>{row.signedEffect.toFixed(3)}</td>
              <td>{row.percentile == null ? "—" : row.percentile.toFixed(1)}</td>
            </tr>
          ))}
        </tbody>
      </table>
      <p className="interp">
        Seed {disruption.seed.toLocaleString()} · {disruption.method}
      </p>
      <figure className="saliency-logo">
        <img
          src={disruption.figure.png}
          alt="Per-bin alternate-versus-reference log2 fold change"
        />
        <figcaption>
          Paper-style per-bin log2(Alt/Ref), centered on the SNV ·{" "}
          <a href={disruption.figure.pdf}>vector PDF</a>
        </figcaption>
      </figure>
    </Section>
  );
}
