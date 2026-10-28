import { Section } from "../components/Section";
import type { Report } from "../schema/report";

type Props = { report: Report };

/** Optional saved last-layer mean attention pyramid. */
export function AttentionSection({ report }: Props) {
  const attention = report.attention;
  if (!attention) return null;

  return (
    <Section id="g-attention" title="Attention">
      <p className="interp">{attention.caveat}</p>
      <figure className="saliency-logo">
        <img
          src={attention.png}
          alt={`Last-layer mean attention pyramid for ${report.locus.chrom}:${report.locus.start}-${report.locus.end}`}
        />
        <figcaption>
          Layer {attention.layer.toLocaleString()} ·{" "}
          {attention.headCount.toLocaleString()} heads ·{" "}
          {attention.binSize.toLocaleString()} bp/bin ·{" "}
          <a href={attention.pdf}>vector PDF</a>
        </figcaption>
      </figure>
    </Section>
  );
}
