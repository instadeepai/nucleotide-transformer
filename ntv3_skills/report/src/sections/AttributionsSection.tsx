import { Section } from "../components/Section";
import type { Report } from "../schema/report";

type Props = { report: Report };

const METHOD_LABELS: Record<string, string> = {
  saliency: "Saliency",
  integrated_gradients: "Integrated Gradients",
  smoothgrad: "SmoothGrad (Saliency)",
  ism: "ISM substitution deltas",
};

/** Method-aware static logos backed by full-locus arrays on disk. */
export function AttributionsSection({ report }: Props) {
  if (report.attributionTracks.length === 0) {
    return (
      <Section id="g-attrib" title="Attributions">
        {report.saliencyTracks.map((track) => (
          <div className="mod" key={track.name}>
            <div className="mod-lab">
              {track.name} · saliency · signed reference contribution
            </div>
            <div className="mod-pane">
              <figure className="saliency-logo">
                <img src={track.logo.png} alt={`${track.name} saliency logo`} />
                <figcaption>
                  Display crop only · full-span compatibility array is saved ·{" "}
                  <a href={track.logo.pdf}>vector PDF</a>
                </figcaption>
              </figure>
            </div>
          </div>
        ))}
      </Section>
    );
  }
  return (
    <Section id="g-attrib" title="Attributions">
      {report.attributionTracks.map((track) => (
        <div className="mod" key={`${track.method}-${track.modelHead}`}>
          <div className="mod-lab">
            {track.name} · {track.methodLabel ?? METHOD_LABELS[track.method]} ·{" "}
            {track.referenceSemantic.replaceAll("_", " ")}
          </div>
          <div className="mod-pane">
            <figure className="saliency-logo">
              <img
                src={track.logo.png}
                alt={`${track.name} ${track.method} reference-display logo`}
              />
              <figcaption>
                {report.locus.chrom}:{track.logo.start.toLocaleString()}–
                {track.logo.end.toLocaleString()} · display crop only · full-span{" "}
                <code>{track.hypotheticalPath}</code> and{" "}
                <code>{track.referencePath}</code> are saved ·{" "}
                <a href={track.logo.pdf}>vector PDF</a>
              </figcaption>
            </figure>
            {track.caveats.map((caveat) => (
              <p className="interp" key={caveat}>
                {caveat}
              </p>
            ))}
          </div>
        </div>
      ))}
    </Section>
  );
}
