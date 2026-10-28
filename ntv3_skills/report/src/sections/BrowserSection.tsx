import { lazy, Suspense, useEffect, useState } from "react";
import type { BwSpec } from "../browser/bigwig";
import {
  locusAnnotationLabel,
  reportAnnotationTracks,
} from "../browser/locusNav";
import {
  hasStoredValues,
  predictionLabel,
  saliencyLabel,
  signalLabel,
} from "../browser/lanes";
import type { LocusPack, TrackSlot } from "../browser/types";
import { fetchChromSize } from "../browser/ucsc";
import { Section } from "../components/Section";
import { findCatalogTrack } from "../data/ucscCatalog";
import {
  catalogGroup,
  catalogRank,
  catalogSlotKind,
} from "../data/ucscGroups";
import type { Report } from "../schema/report";
import "../styles/browser.css";

const GenomeBrowser = lazy(async () => {
  const module = await import("../browser/GenomeBrowser");
  return { default: module.GenomeBrowser };
});

const POPULATE = {
  annotations: true,
  gt: true,
  pred: true,
  saliency: true,
};

type Props = { report: Report };

/**
 * Aligned genome browser as a report section. Reads this run only.
 */
export function BrowserSection({ report }: Props) {
  const [chromSize, setChromSize] = useState(report.locus.end);
  useEffect(() => {
    const abortController = new AbortController();
    setChromSize(report.locus.end);
    fetchChromSize(report.assembly, report.locus.chrom, abortController.signal)
      .then(setChromSize)
      .catch(() => undefined);
    return () => abortController.abort();
  }, [report.assembly, report.locus.chrom, report.locus.end]);
  const pack = packFromReport(report, chromSize);
  const { slots, bwTracks } = slotsFromReport(report);
  const crop = `${report.tracks.chrom}:${report.tracks.start.toLocaleString()}-${report.tracks.end.toLocaleString()}`;

  return (
    <Section id="g-browser" title="Genome browser">
      <p className="lede">
        NTv3 predictions cover <code>{crop}</code>. Pan past that span with the
        arrows: UCSC annotations stay visible.{" "}
        <code>back to locus</code> returns. Genes and ClinVar are fetched live
        for the window; add any other track from the UCSC catalog. Click a
        prediction&apos;s label to recolor it.{" "}
        {report.saliencyTracks.length > 0
          ? "Zoom in until bases appear to read the saliency letters."
          : "Saliency letters were not requested for this run."}
      </p>
      <Suspense fallback={<p className="muted">Loading genome browser…</p>}>
        <GenomeBrowser
          pack={pack}
          populate={POPULATE}
          slots={slots}
          bwTracks={bwTracks}
        />
      </Suspense>
    </Section>
  );
}

function packFromReport(report: Report, chromSize: number): LocusPack {
  const predStart = report.tracks.start;
  const predEnd = report.tracks.end;
  return {
    assembly: report.assembly,
    chrom: report.locus.chrom,
    start: report.locus.start,
    end: report.locus.end,
    chromSize,
    predicted: report.tracks.assays.map((assay, i) => ({
      id: `pr-${i}`,
      name: predictionLabel(assay.name),
      // The scale shows the number only; "model logit" is said in the Summary.
      unit: "",
      source: assay.source,
      values: assay.predicted,
      spanStart: predStart,
      spanEnd: predEnd,
    })),
    saliency: report.saliencyTracks.map((track, i) => ({
      id: `sal-${i}`,
      name: saliencyLabel(track.name),
      unit: "contribution",
      source: `${report.modelShort} · ${track.modelHead}`,
      signed: true,
      values: track.values,
      spanStart: track.start,
      spanEnd: track.end,
      seq: track.seq,
    })),
    genes: {
      [ANNOTATION_SLOT]: (report.annotation?.genes ?? []).filter(
        (gene) => gene.chrom === report.locus.chrom,
      ),
    },
    saved: report.ucscData,
    observed: report.tracks.assays.flatMap((assay, i) =>
      hasStoredValues(assay.observed)
        ? [
            {
              id: `gt-${i}`,
              name: signalLabel(assay.name),
              unit: "signal",
              source: "ground truth stored in the report",
              values: assay.observed.map((value) => value ?? Number.NaN),
              spanStart: predStart,
              spanEnd: predEnd,
            },
          ]
        : [],
    ),
  };
}

const ANNOTATION_SLOT = "ntv3-annotation";

/** The NTv3 gene-annotation row, for a run with an annotation output. */
function annotationSlot(report: Report): TrackSlot {
  return {
    id: ANNOTATION_SLOT,
    kind: "genes",
    label: "NTv3 annotation",
    group: "Annotations",
    defaultOn: true,
    stackRank: 15,
    origin: report.annotation?.model ?? "annotation_plan.json only",
    empty:
      "No NTv3 annotation in this run: the annotation models run only at " +
      "1,048,576 bp, and this laptop does not run them.",
  };
}

function slotsFromReport(report: Report): {
  slots: TrackSlot[];
  bwTracks: Record<string, BwSpec>;
} {
  const catalog: TrackSlot[] = reportAnnotationTracks(report.ucscTracks).map(
    (name) => {
      const entry = findCatalogTrack(name, report.assembly);
      const type = entry?.type ?? "";
      const group = entry?.group ?? "other";
      return {
        id: `ucsc:${name}`,
        kind: catalogSlotKind(type, name),
        label: locusAnnotationLabel(name, report.assembly) ?? entry?.label ?? name,
        group: catalogGroup(group),
        ucscTrack: name,
        defaultOn: true,
        signed: /phylop/i.test(name),
        stackRank: catalogRank(group),
        origin: `UCSC ${name}`,
      };
    },
  );
  // No annotation output, no rail: an empty dashed box says nothing useful.
  if (report.annotation) catalog.push(annotationSlot(report));
  const bwTracks: Record<string, BwSpec> = {};
  const assaySlots: TrackSlot[] = [];
  report.tracks.assays.forEach((assay, i) => {
    const key = slug(assay.name);
    const predId = `pr-${i}`;
    assaySlots.push({
      id: predId,
      kind: "pred",
      label: predictionLabel(assay.name),
      group: "Predictions",
      assay: key,
      defaultOn: true,
      stackRank: 51 + i * 2,
      origin: assay.source,
    });
    // A ground-truth lane only when its values are stored in the report; its
    // BigWig URL is never fetched, so an offline page shows no failed lane.
    if (hasStoredValues(assay.observed)) {
      const gt = report.gtSources.find(
        (source) => source.modelHead === assay.modelHead,
      );
      assaySlots.push({
        id: `gt-${i}`,
        kind: "gt",
        label: signalLabel(gt?.name || assay.name),
        group: "Ground truth",
        assay: key,
        defaultOn: true,
        stackRank: 50 + i * 2,
        origin: gt?.accession ?? "ground truth stored in the report",
      });
    }
  });
  report.saliencyTracks.forEach((track, i) => {
    assaySlots.push({
      id: `sal-${i}`,
      kind: "saliency",
      label: saliencyLabel(track.name),
      group: "Importance",
      defaultOn: true,
      stackRank: 70 + i,
      origin: `${report.modelShort} · ${track.modelHead}`,
    });
  });
  return { slots: [...catalog, ...assaySlots], bwTracks };
}

function slug(name: string): string {
  return name.toLowerCase().replace(/[^a-z0-9]+/g, "-").replace(/^-|-$/g, "") || "assay";
}
