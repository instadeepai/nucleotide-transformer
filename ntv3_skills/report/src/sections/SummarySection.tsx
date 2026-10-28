import { Section } from "../components/Section";
import type { Report } from "../schema/report";

function fmt(n: number): string {
  return n.toLocaleString("en-US");
}

type Props = { report: Report };

export function SummarySection({ report }: Props) {
  const trackSpan = `${report.tracks.chrom}:${fmt(report.tracks.start)}-${fmt(
    report.tracks.end,
  )}`;
  const predictionBp = report.tracks.end - report.tracks.start;
  const predictionPercent = Math.round((predictionBp / report.windowBp) * 100);

  return (
    <Section id="g-summary" title="Summary">
      <div className="summary-stack">
        <p className="interp">
          <span className="summary-kicker">Prediction</span>
          NTv3 generated {report.tracks.assays.length.toLocaleString()} functional
          track{report.tracks.assays.length === 1 ? "" : "s"} over{" "}
          <code>{trackSpan}</code>. This is the centered {predictionPercent}% (
          {predictionBp.toLocaleString()} bp) of the requested{" "}
          {report.windowBp.toLocaleString()} bp window. These values are model
          outputs, not experimental measurements. Explain lists the model, window,
          crop, and head behind each one.
        </p>
        {report.disruption ? (
          <p className="interp">
            <span className="summary-kicker">Compare</span>
            {report.disruption.variant.label} against the reference allele. The
            tracks that move most are listed under Compare.
          </p>
        ) : null}
        {report.gtSources.length > 0 ? (
          <p className="interp">
            <span className="summary-kicker">Ground truth</span>
            Computed prediction-vs-ground-truth comparison for{" "}
            {report.trackAgreement.length.toLocaleString()} public ENCODE track
            {report.trackAgreement.length === 1 ? "" : "s"} over the exact model
            bins. Model logits and experimental signal differ; correlations compare
            shape.
          </p>
        ) : null}
        {report.saliencyTracks.length > 0 ? (
          <p className="interp">
            <span className="summary-kicker">Saliency</span>
            Signed saliency for{" "}
            {report.saliencyTracks.length.toLocaleString()} track
            {report.saliencyTracks.length === 1 ? "" : "s"} is a row in the genome
            browser; zoom in to read the bases. Saliency is not a causal effect.
          </p>
        ) : null}
      </div>
    </Section>
  );
}
