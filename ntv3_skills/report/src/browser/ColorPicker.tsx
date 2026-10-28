import { useEffect, useRef, useState } from "react";

const PALETTE = [
  "#1f63c4",
  "#3d4653",
  "#c45c1f",
  "#12704f",
  "#a81f26",
  "#6b4c9a",
  "#0e7c7b",
  "#b45309",
  "#3973fc",
  "#6b7482",
];

export function ColorPalette({
  label,
  value,
  top,
  left,
  onChange,
  onClose,
}: {
  label: string;
  value: string;
  top: number;
  left: number;
  onChange: (hex: string) => void;
  onClose: () => void;
}) {
  const [hex, setHex] = useState(value);
  const root = useRef<HTMLDivElement>(null);

  useEffect(() => {
    setHex(value);
  }, [value]);

  useEffect(() => {
    const onDoc = (e: MouseEvent) => {
      if (!root.current?.contains(e.target as Node)) onClose();
    };
    const onKey = (e: KeyboardEvent) => {
      if (e.key === "Escape") onClose();
    };
    document.addEventListener("mousedown", onDoc);
    document.addEventListener("keydown", onKey);
    return () => {
      document.removeEventListener("mousedown", onDoc);
      document.removeEventListener("keydown", onKey);
    };
  }, [onClose]);

  const applyHex = () => {
    const next = normalizeHex(hex);
    if (next) onChange(next);
    else setHex(value);
  };

  return (
    <div
      ref={root}
      className="gb-palette"
      role="dialog"
      aria-label={`${label} color`}
      style={{ top, left }}
    >
      <div className="gb-palette-grid">
        {PALETTE.map((c) => (
          <button
            key={c}
            type="button"
            className={c.toLowerCase() === value.toLowerCase() ? "on" : undefined}
            style={{ background: c }}
            aria-label={c}
            onClick={() => {
              onChange(c);
              setHex(c);
              onClose();
            }}
          />
        ))}
      </div>
      <input
        className="gb-hex"
        value={hex}
        spellCheck={false}
        placeholder="#1f63c4"
        aria-label={`${label} hex color`}
        onChange={(e) => setHex(e.target.value)}
        onBlur={applyHex}
        onKeyDown={(e) => {
          if (e.key === "Enter") {
            applyHex();
            onClose();
          }
        }}
      />
    </div>
  );
}

function normalizeHex(raw: string): string | null {
  let s = raw.trim().replace(/^#/, "");
  if (/^[0-9a-fA-F]{3}$/.test(s)) {
    s = [...s].map((ch) => ch + ch).join("");
  }
  if (!/^[0-9a-fA-F]{6}$/.test(s)) return null;
  return `#${s.toLowerCase()}`;
}
