import type { ReactNode } from "react";
import type {
  AttributionTrack as AttributionTrackData,
  BaseWeights,
} from "../schema/report";
import { NUCLEOTIDE_COLORS } from "../data/colors";
import { GENOME_FW } from "./genomeLayout";

const BASES = ["A", "C", "G", "T"] as const;
const GLYPH = 72;

/**
 * Local layout — tighter left gut than Inference (no lane-name column).
 * Keep full GENOME_FW so the SVG still spans the mod pane.
 */
const GUT = 48;
const RIGHT = 12;
const PW = GENOME_FW - GUT - RIGHT;
const STACK_H = 64;
const MID = 70;
const LABEL_GAP = 16;
const H = MID + STACK_H + LABEL_GAP + 4;

type Props = {
  track: AttributionTrackData;
};

function fmtRelKb(bp: number, origin: number): string {
  const kb = (bp - origin) / 1000;
  const sign = kb >= 0 ? "+" : "−";
  return `${sign}${Math.abs(kb).toFixed(1)} kb`;
}

function signedParts(col: BaseWeights): {
  pos: { b: (typeof BASES)[number]; v: number }[];
  neg: { b: (typeof BASES)[number]; v: number }[];
  extent: number;
} {
  const pos = BASES.map((b, bi) => ({ b, v: Math.max(0, col[bi] ?? 0) }))
    .filter((e) => e.v > 0)
    .sort((a, b) => a.v - b.v);
  const neg = BASES.map((b, bi) => ({ b, v: Math.max(0, -(col[bi] ?? 0)) }))
    .filter((e) => e.v > 0)
    .sort((a, b) => a.v - b.v);
  const sumPos = pos.reduce((s, e) => s + e.v, 0);
  const sumNeg = neg.reduce((s, e) => s + e.v, 0);
  return { pos, neg, extent: Math.max(sumPos, sumNeg) };
}

function stackSigned(
  col: BaseWeights,
  cx: number,
  maxExtent: number,
  sx: number,
): ReactNode[] {
  const { pos, neg } = signedParts(col);
  const nodes: ReactNode[] = [];

  let cursorUp = MID;
  for (const e of pos) {
    const h = (e.v / maxExtent) * STACK_H;
    if (h < 0.3) continue;
    nodes.push(
      <g
        key={`p-${e.b}`}
        transform={`translate(${cx}, ${cursorUp}) scale(${sx.toFixed(4)}, ${(h / GLYPH).toFixed(4)})`}
      >
        <text
          x={0}
          y={0}
          fontSize={GLYPH}
          fontWeight={700}
          fill={NUCLEOTIDE_COLORS[e.b]}
          textAnchor="middle"
          fontFamily="var(--mono)"
          style={{ letterSpacing: "-0.04em" }}
        >
          {e.b}
        </text>
      </g>,
    );
    cursorUp -= h;
  }

  let cursorDn = MID;
  for (const e of neg) {
    const h = (e.v / maxExtent) * STACK_H;
    if (h < 0.3) continue;
    const origin = cursorDn + h;
    nodes.push(
      <g
        key={`n-${e.b}`}
        transform={`translate(${cx}, ${origin}) scale(${sx.toFixed(4)}, ${(h / GLYPH).toFixed(4)})`}
      >
        <text
          x={0}
          y={0}
          fontSize={GLYPH}
          fontWeight={700}
          fill={NUCLEOTIDE_COLORS[e.b]}
          textAnchor="middle"
          fontFamily="var(--mono)"
          style={{ letterSpacing: "-0.04em" }}
        >
          {e.b}
        </text>
      </g>,
    );
    cursorDn += h;
  }
  return nodes;
}

/**
 * Per-position contribution logo (signed, mean-centered).
 * Compact left gut — only room for y-label + ticks.
 */
export function AttributionTrack({ track }: Props) {
  const { start, end, scores, originBp, chrom } = track;
  const n = Math.max(1, scores.length);
  const colW = PW / n;
  const maxExtent = Math.max(
    0.05,
    ...scores.map((c) => signedParts(c).extent),
  );
  const sx = Math.min(1.15, Math.max(0.45, (colW * 0.98) / 40));
  const yTop = MID - STACK_H;
  const yBot = MID + STACK_H;
  const tickLabelX = GUT - 6;
  const yTitleX = 12;

  return (
    <div className="attr-track">
      <svg
        viewBox={`0 0 ${GENOME_FW} ${H}`}
        preserveAspectRatio="xMinYMin meet"
        role="img"
        aria-label={`Contribution logo ${chrom}:${start}-${end}`}
      >
        <line
          x1={GUT}
          y1={yTop}
          x2={GUT}
          y2={yBot}
          stroke="#3d4653"
          strokeWidth={0.9}
        />
        {[
          { y: yTop, lab: maxExtent.toFixed(1) },
          { y: MID, lab: "0" },
          { y: yBot, lab: (-maxExtent).toFixed(1) },
        ].map((t) => (
          <g key={t.lab}>
            <line
              x1={GUT - 4}
              y1={t.y}
              x2={GUT}
              y2={t.y}
              stroke="#3d4653"
            />
            <text
              x={tickLabelX}
              y={t.y + 3}
              fontSize={9}
              fill="#6b7482"
              textAnchor="end"
              fontFamily="var(--mono)"
            >
              {t.lab}
            </text>
          </g>
        ))}
        <text
          x={yTitleX}
          y={MID}
          fontSize={11}
          fill="#3d4653"
          fontFamily="var(--sans)"
          transform={`rotate(-90 ${yTitleX} ${MID})`}
          textAnchor="middle"
        >
          Attribution
        </text>

        {scores.map((col, i) => {
          const cx = GUT + i * colW + colW / 2;
          return <g key={i}>{stackSigned(col, cx, maxExtent, sx)}</g>;
        })}

        <line
          x1={GUT}
          y1={MID}
          x2={GUT + PW}
          y2={MID}
          stroke="#3d4653"
          strokeWidth={0.9}
        />
        <text
          x={GUT}
          y={yBot + 14}
          fontSize={9}
          fill="#6b7482"
          textAnchor="start"
          fontFamily="var(--mono)"
        >
          {fmtRelKb(start, originBp)}
        </text>
        <text
          x={GUT + PW}
          y={yBot + 14}
          fontSize={9}
          fill="#6b7482"
          textAnchor="end"
          fontFamily="var(--mono)"
        >
          {fmtRelKb(end, originBp)}
        </text>
      </svg>
      <p className="attr-track__src">
        Mean-centered contribution scores (L×4) · {chrom}:
        {start.toLocaleString("en-US")}–{end.toLocaleString("en-US")} · relative
        to {originBp.toLocaleString("en-US")}.
      </p>
    </div>
  );
}
