import { Section } from "../components/Section";
import type { Report } from "../schema/report";

type Props = { report: Report };

/** Trace without strikethrough on skipped. Methods replaced by run parameters. */
export function ProvenanceSection({ report }: Props) {
  const params = Object.entries(report.runParameters);

  return (
    <Section id="g-trace" title="Provenance">
      <div className="block">
        <p className="subhead">Trace</p>
        <div className="trace">
          {report.trace.map((step, i) => (
            <div
              key={`${step.name}-${i}`}
              className={`span${step.status === "skipped" ? " skip" : ""}`}
            >
              <span className="n">{String(i + 1).padStart(2, "0")}</span>
              <span className="name">{step.name}</span>
              <span className="art">{step.artifact ?? "—"}</span>
              <span className={`st ${step.status === "ok" ? "ok" : "skip"}`}>
                {step.status === "ok" ? "ok" : "skipped"}
              </span>
            </div>
          ))}
        </div>
      </div>

      {report.critic ? (
        <div className="block">
          <p className="subhead">Critic</p>
          <p className="critic-line">
            {report.critic.enabled
              ? `The critic ran: ${report.critic.kept} kept, ${report.critic.suppressed} suppressed.`
              : "The critic was off for this run; it was turned off at approval."}
          </p>
        </div>
      ) : null}

      <div className="block">
        <p className="subhead">Run parameters</p>
        <pre className="run-params" aria-label="Run parameters">
          <span className="c"># run parameters</span>
          {"\n{\n"}
          {params.map(([k, v], i) => {
            const isNum = typeof v === "number";
            const val = isNum ? String(v) : `"${v}"`;
            const comma = i < params.length - 1 ? "," : "";
            return (
              <span key={k}>
                {"  "}
                <span className="k">"{k}"</span>:{" "}
                <span className={isNum ? "v-num" : "v"}>{val}</span>
                {comma}
                {"\n"}
              </span>
            );
          })}
          {"}"}
        </pre>
      </div>
    </Section>
  );
}
