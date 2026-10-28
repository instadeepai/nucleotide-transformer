import { useCallback, useState } from "react";
import { clampRange } from "./locusNav";

export function useViewport(
  initialStart: number,
  initialEnd: number,
  chromSize: number,
) {
  const [viewStart, setViewStart] = useState(initialStart);
  const [viewEnd, setViewEnd] = useState(initialEnd);

  const setView = useCallback(
    (start: number, end: number, size: number = chromSize) => {
      const next = clampRange(start, end, size);
      setViewStart(next.start);
      setViewEnd(next.end);
    },
    [chromSize],
  );

  const shiftFrac = useCallback(
    (frac: number) => {
      const span = viewEnd - viewStart;
      setView(viewStart + frac * span, viewEnd + frac * span);
    },
    [viewStart, viewEnd, setView],
  );

  const zoomAt = useCallback(
    (factor: number, anchorBp: number) => {
      const span = viewEnd - viewStart;
      const newSpan = span / factor;
      const t = (anchorBp - viewStart) / span;
      setView(anchorBp - t * newSpan, anchorBp + (1 - t) * newSpan);
    },
    [viewStart, viewEnd, setView],
  );

  const reset = useCallback(() => {
    setViewStart(initialStart);
    setViewEnd(initialEnd);
  }, [initialStart, initialEnd]);

  return { viewStart, viewEnd, setView, shiftFrac, zoomAt, reset };
}
