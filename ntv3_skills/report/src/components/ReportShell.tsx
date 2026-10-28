import { useEffect, useMemo, useState, type ReactNode } from "react";
import type { Report } from "../schema/report";

const NAV = [
  { id: "g-locus", label: "Locus" },
  { id: "g-summary", label: "Summary" },
  { id: "g-agreement", label: "Agreement" },
  { id: "g-browser", label: "Browser" },
  { id: "g-disruption", label: "Compare" },
  { id: "g-explain", label: "Explain" },
  { id: "g-literature", label: "Literature" },
  { id: "g-trace", label: "Provenance" },
] as const;

/** Only sections with results are listed. Nothing is shown as pending. */
function visibleNav(report: Report) {
  const present: Record<(typeof NAV)[number]["id"], boolean> = {
    "g-locus": true,
    "g-summary": true,
    "g-agreement": report.trackAgreement.length > 0,
    "g-browser": true,
    "g-disruption": report.disruption != null,
    "g-explain": report.explain != null,
    "g-literature": report.literature != null,
    "g-trace": true,
  };
  return NAV.filter((item) => present[item.id]);
}

type ReportShellProps = {
  report: Report;
  children: ReactNode;
};

export function ReportShell({ report, children }: ReportShellProps) {
  const [current, setCurrent] = useState<string>("g-locus");
  const nav = useMemo(() => visibleNav(report), [report]);

  useEffect(() => {
    document.title = `${report.gene} in ${report.cell}`;
  }, [report.gene, report.cell]);

  useEffect(() => {
    const sections = nav.map((item) => document.getElementById(item.id)).filter(
      (el): el is HTMLElement => el != null,
    );
    if (!sections.length) return;
    const updateCurrent = () => {
      const anchor = window.innerHeight * 0.25;
      let active = sections[0];
      for (const section of sections) {
        if (section.getBoundingClientRect().top <= anchor) active = section;
      }
      setCurrent(active.id);
    };
    updateCurrent();
    window.addEventListener("scroll", updateCurrent, { passive: true });
    window.addEventListener("resize", updateCurrent);
    return () => {
      window.removeEventListener("scroll", updateCurrent);
      window.removeEventListener("resize", updateCurrent);
    };
  }, [nav]);

  return (
    <div className="wrap">
      <header className="rep">
        <h1>
          {report.gene} in {report.cell}
        </h1>
        <button type="button" onClick={() => window.location.reload()}>
          Refresh report
        </button>
      </header>

      <nav className="obs-nav" aria-label="Report sections">
        {nav.map((item) => (
          <a
            key={item.id}
            href={`#${item.id}`}
            aria-current={current === item.id ? "true" : undefined}
            onClick={(event) => {
              event.preventDefault();
              document.getElementById(item.id)?.scrollIntoView({ block: "start" });
              window.history.replaceState(null, "", `#${item.id}`);
              setCurrent(item.id);
            }}
          >
            {item.label}
          </a>
        ))}
      </nav>

      <div className="stream">{children}</div>
    </div>
  );
}
