import { useEffect, useState } from "react";
import { ReportShell } from "./components/ReportShell";
import { loadRun } from "./loadRun";
import type { Report } from "./schema/report";
import { AgreementSection } from "./sections/AgreementSection";
import { BrowserSection } from "./sections/BrowserSection";
import { DisruptionSection } from "./sections/DisruptionSection";
import { ExplainSection } from "./sections/ExplainSection";
import { LocusSection } from "./sections/LocusSection";
import { LiteratureSection } from "./sections/LiteratureSection";
import { ProvenanceSection } from "./sections/ProvenanceSection";
import { SummarySection } from "./sections/SummarySection";

export function App() {
  const [report, setReport] = useState<Report | null>(null);
  const [error, setError] = useState<string | null>(null);

  useEffect(() => {
    void loadRun()
      .then(setReport)
      .catch((err: unknown) => {
        setError(err instanceof Error ? err.message : "Failed to load report.json");
      });
  }, []);

  if (error) {
    return (
      <main style={{ padding: "2rem", fontFamily: "sans-serif" }}>
        <h1>No run loaded</h1>
        <p>{error}</p>
        <p>
          After <code>ntv3 report</code>, start the app with{" "}
          <code>VITE_REPORT_URL</code> pointing at that run&apos;s{" "}
          <code>report.json</code>.
        </p>
      </main>
    );
  }

  if (!report) {
    return (
      <main style={{ padding: "2rem", fontFamily: "sans-serif" }}>
        <p>Loading this run…</p>
      </main>
    );
  }

  return (
    <ReportShell report={report}>
      <LocusSection report={report} />
      <SummarySection report={report} />
      {report.trackAgreement.length > 0 ? (
        <AgreementSection report={report} />
      ) : null}
      <BrowserSection report={report} />
      {report.disruption ? <DisruptionSection report={report} /> : null}
      <ExplainSection report={report} />
      {report.literature ? <LiteratureSection report={report} /> : null}
      <ProvenanceSection report={report} />
    </ReportShell>
  );
}
