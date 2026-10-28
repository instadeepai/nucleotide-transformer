import { z } from "zod";

export const BaseWeightsSchema = z.tuple([
  z.number(),
  z.number(),
  z.number(),
  z.number(),
]);

const AttributionCalloutSchema = z.object({
  start: z.number(),
  end: z.number(),
  label: z.string(),
  pwm: z.array(BaseWeightsSchema).optional(),
});

export const AttributionTrackSchema = z.object({
  chrom: z.string(),
  start: z.number(),
  end: z.number(),
  originBp: z.number(),
  scores: z.array(BaseWeightsSchema),
  callouts: z.array(AttributionCalloutSchema).default([]),
});

const TrackAssaySchema = z.object({
  name: z.string(),
  modelHead: z.string(),
  unit: z.string(),
  range: z.tuple([z.number(), z.number()]),
  source: z.string(),
  observed: z.array(z.number().nullable()),
  predicted: z.array(z.number()),
});

const TracksSchema = z.object({
  chrom: z.string(),
  start: z.number(),
  end: z.number(),
  assays: z.array(TrackAssaySchema),
});

const ClaimSchema = z.object({
  id: z.string(),
  text: z.string(),
  support: z.string(),
  command: z.string().nullable(),
  keep: z.boolean(),
});

const TraceStepSchema = z.object({
  name: z.string(),
  artifact: z.string().nullable(),
  status: z.enum(["ok", "skipped"]),
});

const AnalysisStageSchema = z.object({
  id: z.enum([
    "attention",
    "disruption",
    "shuffle",
    "design",
    "motifs",
    "literature",
    "claims",
  ]),
  label: z.string(),
  status: z.enum(["complete", "pending", "not_requested"]),
  message: z.string(),
});

const SpanSchema = z.object({
  chrom: z.string(),
  start: z.number().int().nonnegative(),
  end: z.number().int().positive(),
  bp: z.number().int().positive(),
});

/** Which model, window, crop, and head produced each number. From provenance. */
const ExplainSchema = z.object({
  model: z.string(),
  manifestModel: z.string().nullable().optional(),
  species: z.string(),
  assembly: z.string(),
  device: z.string(),
  window: SpanSchema,
  crop: SpanSchema,
  binSize: z.number().int().positive(),
  nBins: z.number().int().positive(),
  tracks: z.array(
    z.object({
      label: z.string(),
      modelHead: z.string(),
      gtAccession: z.string(),
      // ENCODE output type of the ground-truth BigWig; absent in older reports.
      gtUnit: z.string().optional(),
    }),
  ),
  letters: z
    .array(
      z.object({
        label: z.string(),
        modelHead: z.string(),
        start: z.number().int().nonnegative(),
        end: z.number().int().positive(),
      }),
    )
    .default([]),
  compare: z
    .object({
      label: z.string(),
      checkpoint: z.string(),
      window: SpanSchema,
      crop: SpanSchema,
      binSize: z.number().int().positive(),
    })
    .nullable()
    .default(null),
  sources: z.array(z.string()),
});

/**
 * NTv3 gene models from a GFF3 in the run. Null when the run has none: the
 * annotation models run only at 1,048,576 bp, and a laptop does not run them.
 */
const AnnotationSchema = z.object({
  model: z.string().nullable().optional(),
  genes: z.array(
    z.object({
      chrom: z.string(),
      name: z.string(),
      start: z.number().int().nonnegative(),
      end: z.number().int().positive(),
      strand: z.enum(["+", "-"]),
      exons: z.array(z.object({ start: z.number().int(), end: z.number().int() })),
    }),
  ),
});

/** What the critic did for this run, as `ntv3 report` recorded it. */
const CriticSchema = z.union([
  z.object({
    enabled: z.literal(true),
    log: z.string(),
    claims: z.string(),
    kept: z.number().int().nonnegative(),
    suppressed: z.number().int().nonnegative(),
  }),
  z.object({ enabled: z.literal(false), note: z.string() }),
]);

export const ReportSchema = z.object({
  tool: z.string(),
  commit: z.string(),
  gene: z.string(),
  cell: z.string(),
  assembly: z.string(),
  modelShort: z.string(),
  windowBp: z.number().positive(),
  locus: z.object({
    chrom: z.string(),
    start: z.number().int().nonnegative(),
    end: z.number().int().positive(),
    tss: z.number().int().nullable(),
    strand: z.enum(["+", "-"]).nullable(),
    focusStart: z.number().int().nonnegative(),
    focusEnd: z.number().int().positive(),
    elements: z.array(z.string()),
    source: z
      .object({
        track: z.string(),
        url: z.string(),
      })
      .nullable()
      .optional(),
    assays: z.array(z.string()),
  }),
  tracks: TracksSchema,
  claims: z.array(ClaimSchema),
  trace: z.array(TraceStepSchema),
  runParameters: z.record(z.string(), z.union([z.string(), z.number()])),
  gtSources: z.array(
    z.object({
      name: z.string(),
      accession: z.string(),
      fileAccession: z.string(),
      modelHead: z.string(),
      url: z.string(),
      origin: z.string(),
    }),
  ),
  trackAgreement: z.array(
    z.object({
      trackQuery: z.string(),
      modelHead: z.string(),
      gtAccession: z.string(),
      fileAccession: z.string(),
      nValid: z.number().int().nonnegative(),
      pearson: z.number().nullable(),
      spearman: z.number().nullable(),
      caveat: z.string(),
    }),
  ),
  saliencyTracks: z.array(
    z.object({
      name: z.string(),
      modelHead: z.string(),
      values: z.array(z.number()),
      start: z.number().int().nonnegative(),
      end: z.number().int().positive(),
      seq: z.string(),
      logo: z.object({
        selection: z.enum(["strongest_peak", "explicit_center"]),
        centerBp: z.number().int().nonnegative(),
        start: z.number().int().nonnegative(),
        end: z.number().int().positive(),
        png: z.string(),
        pdf: z.string(),
      }),
    }),
  ),
  attributionTracks: z
    .array(
      z.object({
        name: z.string(),
        modelHead: z.string(),
        method: z.enum([
          "saliency",
          "integrated_gradients",
          "smoothgrad",
          "ism",
        ]),
        methodLabel: z.string().optional(),
        methodClass: z.string(),
        start: z.number().int().nonnegative(),
        end: z.number().int().positive(),
        hypotheticalSemantic: z.string(),
        referenceSemantic: z.string(),
        hypotheticalPath: z.string(),
        referencePath: z.string(),
        caveats: z.array(z.string()),
        logo: z.object({
          selection: z.enum(["strongest_peak", "explicit_center"]),
          centerBp: z.number().int().nonnegative(),
          start: z.number().int().nonnegative(),
          end: z.number().int().positive(),
          png: z.string(),
          pdf: z.string(),
          mode: z.literal("reference"),
        }),
      }),
    )
    .default([]),
  attention: z
    .object({
      layer: z.number().int().positive(),
      headCount: z.number().int().positive(),
      matrixShape: z.tuple([
        z.number().int().positive(),
        z.number().int().positive(),
      ]),
      binSize: z.number().int().positive(),
      checkpoint: z.string(),
      source: z.string(),
      caveat: z.string(),
      png: z.string(),
      pdf: z.string(),
    })
    .nullable()
    .optional(),
  disruption: z
    .object({
      variant: z.object({
        chrom: z.string(),
        position: z.number().int().positive(),
        ref: z.string().length(1),
        alt: z.string().length(1),
        label: z.string(),
      }),
      signedEffect: z.number(),
      absoluteEffect: z.number().nonnegative(),
      // Null when no null SNVs were run (nullN 0), the default Compare.
      percentile: z.number().min(0).max(100).nullable(),
      nullN: z.number().int().nonnegative(),
      seed: z.number().int(),
      perTrack: z.array(
        z.object({
          track: z.string(),
          modelHead: z.string(),
          signedEffect: z.number(),
          percentile: z.number().min(0).max(100).nullable(),
        }),
      ),
      method: z.string(),
      figure: z.object({ png: z.string(), pdf: z.string() }),
      caveat: z.string(),
    })
    .nullable()
    .default(null),
  knockout: z
    .object({
      mechanism: z.string(),
      interval: z.object({
        chrom: z.string(),
        start: z.number().int().nonnegative(),
        end: z.number().int().positive(),
        coordinateSystem: z.string(),
      }),
      seed: z.number().int(),
      sourceEvidence: z.string(),
      perTrack: z.array(
        z.object({
          modelHead: z.string(),
          signedMaxEffect: z.number(),
          meanEffect: z.number(),
        }),
      ),
      figure: z.object({ png: z.string(), pdf: z.string() }),
      caveat: z.string(),
    })
    .nullable()
    .default(null),
  design: z
    .object({
      wording: z.string(),
      modelHead: z.string(),
      trackLabel: z.string(),
      baselineTarget: z.number(),
      targetDelta: z.number(),
      desiredTarget: z.number(),
      selectedCandidate: z.object({
        rank: z.number().int().positive(),
        achievedTarget: z.number(),
        absoluteTargetError: z.number().nonnegative(),
        editCount: z.number().int().nonnegative(),
        edits: z.array(
          z.object({
            position: z.number().int().positive(),
            ref: z.string().length(1),
            alt: z.string().length(1),
          }),
        ),
      }),
      candidateCount: z.number().int().positive(),
      seed: z.number().int(),
      editableInterval: z.object({
        chrom: z.string(),
        start: z.number().int().nonnegative(),
        end: z.number().int().positive(),
        coordinateSystem: z.string(),
      }),
      scoredBins: z.object({
        peak: z.number().int().nonnegative(),
        start: z.number().int().nonnegative(),
        end: z.number().int().positive(),
        smoothing: z.number().int().positive(),
        searchRadius: z.number().int().positive(),
      }),
      figure: z.object({ png: z.string(), pdf: z.string() }),
      caveat: z.string(),
    })
    .nullable()
    .default(null),
  motifs: z
    .object({
      scanner: z.string(),
      scannerVersion: z.string(),
      threshold: z.number().positive(),
      attributionMethod: z.string(),
      modelHead: z.string(),
      trackLabel: z.string(),
      sources: z.array(
        z.object({
          matrixId: z.string(),
          name: z.string(),
          url: z.string(),
        }),
      ),
      hits: z.array(
        z.object({
          matrixId: z.string(),
          matrixName: z.string(),
          start: z.number().int().nonnegative(),
          end: z.number().int().positive(),
          strand: z.string(),
          score: z.number(),
          pValue: z.number().nonnegative(),
          attributionSignedSum: z.number(),
          attributionAbsoluteSum: z.number().nonnegative(),
        }),
      ),
      caveat: z.string(),
    })
    .nullable()
    .default(null),
  literature: z
    .object({
      records: z.array(
        z.object({
          citekey: z.string(),
          title: z.string(),
          authors: z.array(z.string()),
          year: z.number().int().positive(),
          sourceUrl: z.string(),
          doi: z.string().nullable(),
          retrievalDate: z.string(),
          searchQuery: z.string(),
          relevance: z.string(),
        }),
      ),
      bib: z.string(),
    })
    .nullable()
    .default(null),
  annotation: AnnotationSchema.nullable().default(null),
  ucscTracks: z.array(z.string()),
  /**
   * Raw UCSC getData payloads the run's lookup saved, keyed by track. The
   * browser draws these without a network call while the view stays inside
   * the window they cover.
   */
  ucscData: z.record(z.string(), z.unknown()).default({}),
  /** Null only for a report written before the critic was recorded. */
  critic: CriticSchema.nullable().default(null),
  explain: ExplainSchema.nullable().default(null),
  analysisStages: z.array(AnalysisStageSchema).default([]),
});

export type Report = z.infer<typeof ReportSchema>;
export type BaseWeights = z.infer<typeof BaseWeightsSchema>;
export type AttributionTrack = z.infer<typeof AttributionTrackSchema>;
export type AnalysisStage = z.infer<typeof AnalysisStageSchema>;
export type Explain = z.infer<typeof ExplainSchema>;
