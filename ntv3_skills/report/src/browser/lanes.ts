/**
 * Lane names and units, shared by the genome browser and the Explain table so
 * both say the same thing. Pure functions, covered by node tests.
 */

/** Lane labels: what drew the lane, then the track. */
export function predictionLabel(track: string): string {
  return `NTv3 · ${track}`;
}

export function signalLabel(track: string): string {
  return `Signal · ${track}`;
}

export function saliencyLabel(track: string): string {
  return `Saliency · ${track}`;
}

/** True when a lane has at least one finite stored value. */
export function hasStoredValues(values: (number | null)[]): boolean {
  return values.some((value) => value != null && Number.isFinite(value));
}

export type ExplainTrack = {
  label: string;
  modelHead: string;
  gtAccession: string;
  gtUnit?: string;
};
export type ExplainAssay = { modelHead: string; unit: string; observed: (number | null)[] };
export type ExplainLetters = { label: string; modelHead: string };
export type LaneRow = {
  lane: string;
  modelHead: string | null;
  unit: string;
  groundTruth: string | null;
};

/**
 * One Explain row per lane the browser draws, in track order: the prediction,
 * its ground-truth signal when values are stored, then its saliency.
 *
 * Units: the prediction's own unit from ``tracks.assays``; the ENCODE output
 * type of the BigWig (``gtUnit``), else "signal"; "contribution" for saliency.
 */
export function explainLaneRows(
  tracks: ExplainTrack[],
  assays: ExplainAssay[],
  letters: ExplainLetters[],
): LaneRow[] {
  const rows: LaneRow[] = [];
  for (const track of tracks) {
    const assay = assays.find((row) => row.modelHead === track.modelHead);
    rows.push({
      lane: predictionLabel(track.label),
      modelHead: track.modelHead,
      unit: assay?.unit || "—",
      groundTruth: track.gtAccession || null,
    });
    if (track.gtAccession && assay && hasStoredValues(assay.observed)) {
      rows.push({
        lane: signalLabel(track.label),
        modelHead: null,
        unit: track.gtUnit || "signal",
        groundTruth: track.gtAccession,
      });
    }
    for (const row of letters.filter((l) => l.modelHead === track.modelHead)) {
      rows.push({
        lane: saliencyLabel(row.label),
        modelHead: row.modelHead,
        unit: "contribution",
        groundTruth: null,
      });
    }
  }
  return rows;
}
