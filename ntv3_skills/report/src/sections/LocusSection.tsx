import { Section } from "../components/Section";
import type { Report } from "../schema/report";

function fmt(n: number): string {
  return n.toLocaleString("en-US");
}

type Props = { report: Report };

export function LocusSection({ report }: Props) {
  const { locus } = report;
  const windowLabel = `${locus.chrom}:${fmt(locus.start)}-${fmt(locus.end)}`;
  const tssLabel =
    locus.tss != null ? `${fmt(locus.tss)} (${locus.strand ?? "?"})` : null;

  return (
    <Section id="g-locus" title="Locus">
      <dl className="meta-dl">
        <div>
          <dt>Window</dt>
          <dd>
            <code>{windowLabel}</code>
          </dd>
        </div>
        {tssLabel ? (
          <div>
            <dt>TSS</dt>
            <dd>
              <code>{tssLabel}</code>
            </dd>
          </div>
        ) : null}
        {locus.elements.length > 0 ? (
          <div>
            <dt>Elements</dt>
            <dd>
              <span className="pill-row">
                {locus.elements.map((name) => (
                  <span className="pill" key={name}>
                    {name}
                  </span>
                ))}
              </span>
            </dd>
          </div>
        ) : null}
        {locus.source?.url ? (
          <div>
            <dt>Source</dt>
            <dd>
              <a href={locus.source.url} target="_blank" rel="noreferrer">
                UCSC {locus.source.track}
              </a>
            </dd>
          </div>
        ) : null}
        {locus.focusStart !== locus.start || locus.focusEnd !== locus.end ? (
          <div>
            <dt>Focus</dt>
            <dd>
              <code>
                {locus.chrom}:{fmt(locus.focusStart)}-{fmt(locus.focusEnd)}
              </code>
            </dd>
          </div>
        ) : null}
        <div>
          <dt>Assays</dt>
          <dd>
            <span className="pill-row">
              {locus.assays.map((name) => (
                <span className="pill" key={name}>
                  {name}
                </span>
              ))}
            </span>
          </dd>
        </div>
      </dl>
    </Section>
  );
}
