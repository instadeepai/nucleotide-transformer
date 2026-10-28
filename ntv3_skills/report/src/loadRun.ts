import { ReportSchema, type Report } from "./schema/report";

/** Id of the JSON script tag `ntv3 report` fills when it writes one HTML file. */
const EMBEDDED_REPORT_ID = "ntv3-report";

/**
 * Load the run: the report embedded in this page when `ntv3 report` wrote it,
 * else results/<run>/report.json from VITE_REPORT_URL (the dev server).
 * Does not fall back to a hardcoded example run.
 */
export async function loadRun(): Promise<Report> {
  const embedded = document.getElementById(EMBEDDED_REPORT_ID)?.textContent?.trim();
  if (embedded && embedded !== "null") {
    return withBase(ReportSchema.parse(JSON.parse(embedded)), document.baseURI);
  }
  const url = import.meta.env.VITE_REPORT_URL ?? "/report.json";
  const response = await fetch(url, { cache: "no-store" });
  if (!response.ok) {
    throw new Error(
      `No report at ${url} (${response.status}). Point VITE_REPORT_URL at this run's report.json.`,
    );
  }
  return withBase(ReportSchema.parse(await response.json()), response.url);
}

/** Resolve the report's relative figure links against where the run lives. */
function withBase(report: Report, base: string): Report {
  return {
    ...report,
    saliencyTracks: report.saliencyTracks.map((track) => ({
      ...track,
      logo: {
        ...track.logo,
        png: new URL(track.logo.png, base).href,
        pdf: new URL(track.logo.pdf, base).href,
      },
    })),
    attributionTracks: report.attributionTracks.map((track) => ({
      ...track,
      logo: {
        ...track.logo,
        png: new URL(track.logo.png, base).href,
        pdf: new URL(track.logo.pdf, base).href,
      },
    })),
    attention: report.attention
      ? {
          ...report.attention,
          png: new URL(report.attention.png, base).href,
          pdf: new URL(report.attention.pdf, base).href,
        }
      : null,
    disruption: report.disruption
      ? {
          ...report.disruption,
          figure: {
            png: new URL(report.disruption.figure.png, base).href,
            pdf: new URL(report.disruption.figure.pdf, base).href,
          },
        }
      : null,
    knockout: report.knockout
      ? {
          ...report.knockout,
          figure: {
            png: new URL(report.knockout.figure.png, base).href,
            pdf: new URL(report.knockout.figure.pdf, base).href,
          },
        }
      : null,
    design: report.design
      ? {
          ...report.design,
          figure: {
            png: new URL(report.design.figure.png, base).href,
            pdf: new URL(report.design.figure.pdf, base).href,
          },
        }
      : null,
    literature: report.literature
      ? {
          ...report.literature,
          bib: new URL(report.literature.bib, base).href,
        }
      : null,
  };
}
