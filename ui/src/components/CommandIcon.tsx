/**
 * The command glyphs, drawn on the same currentColor convention as the unit glyphs.
 *
 * The text label stays beside the icon, so the accessible name of every ribbon control is still
 * its word; the SVG is `aria-hidden`.
 */

export type CommandIconKind =
  | "solve"
  | "undo"
  | "redo"
  | "report"
  | "new"
  | "open"
  | "demo"
  | "examples"
  | "save"
  | "calculate"
  | "calculate-all";

export function CommandIcon({ kind }: { kind: CommandIconKind }) {
  return (
    <svg className="command-icon" viewBox="0 0 16 16" aria-hidden="true" focusable="false">
      {kind === "solve" ? <path className="fill" d="M4 3 L12 8 L4 13 Z" /> : null}
      {kind === "undo" ? (
        <>
          <path d="M6 4 L2 8 L6 12" />
          <path d="M2 8 H9 A4 4 0 0 1 9 15" />
        </>
      ) : null}
      {kind === "redo" ? (
        <>
          <path d="M10 4 L14 8 L10 12" />
          <path d="M14 8 H7 A4 4 0 0 0 7 15" />
        </>
      ) : null}
      {kind === "report" ? (
        <>
          <rect x="3" y="2" width="10" height="12" rx="1" />
          <path d="M5.5 5 H10.5 M5.5 8 H10.5 M5.5 11 H8.5" />
        </>
      ) : null}
      {kind === "new" ? <path d="M8 3 V13 M3 8 H13" /> : null}
      {kind === "open" ? <path d="M2 5 H7 L8.5 7 H14 V13 H2 Z" /> : null}
      {kind === "demo" ? (
        <>
          <rect x="2" y="3" width="12" height="10" rx="1" />
          <path d="M6 5 L10 8 L6 11 Z" className="fill" />
        </>
      ) : null}
      {kind === "examples" ? (
        <path d="M3 4 H13 M3 8 H13 M3 12 H13" />
      ) : null}
      {kind === "save" ? (
        <>
          <path d="M3 3 H11 L13 5 V13 H3 Z" />
          <path d="M5 3 V7 H10 V3 M5 13 V9 H11 V13" />
        </>
      ) : null}
      {kind === "calculate" ? <path d="M4 3 H12 L8 8 L12 13 H4" /> : null}
      {kind === "calculate-all" ? (
        <>
          <path d="M3 3 H9 L6 8 L9 13 H3" />
          <path d="M10 6 L11.5 8 L14 4" />
        </>
      ) : null}
    </svg>
  );
}
