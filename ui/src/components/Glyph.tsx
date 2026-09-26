/**
 * The drawings: one symbol per unit operation, on the 72 x 48 box the caller hangs the ports on.
 *
 * **A flowsheet is read at a glance, and a pump drawn as a rectangle is not read at all.** What
 * makes a canvas legible without opening a node is that the symbol states what the machine does: a
 * valve is a bowtie because it throttles, a column carries trays because it separates, and the two
 * machines are mirrored because one squeezes a gas and the other lets it expand.
 *
 * **No drawing names a colour.** The outline, the internals, and the one filled mark that says
 * which way a fluid moves are the three classes the stylesheet owns - `body`, `detail` and
 * `accent` - because the theme owns every colour in the editor, and a glyph that set its own stroke
 * would be the one element that ignored it.
 *
 * **A nozzle reaches the edge of the box.** The caller attaches the handles to the box's edges
 * rather than to points this file would have to agree with, so a vessel whose ports leave from the
 * top reaches `y = 0` and one whose ports leave to the side reaches `x = 0` or `x = 72`.
 */

import type { GlyphKind, Sense } from "../state/glyphs";

export function Glyph({
  kind,
  sense,
  selected,
  flagged,
}: {
  kind: GlyphKind;
  /** Which way a machine's moving part points, where it has one. */
  sense: Sense;
  selected: boolean;
  flagged: boolean; // a diagnostic names this node
}): React.ReactElement {
  return (
    <svg
      className={`glyph${selected ? " selected" : ""}${flagged ? " flagged" : ""}`}
      viewBox="0 0 72 48"
      width={72}
      height={48}
      aria-hidden="true"
      focusable="false"
    >
      {drawing(kind, sense)}
    </svg>
  );
}

/** The shapes for one kind, laid out on the box the caller hangs its ports on. */
function drawing(kind: GlyphKind, sense: Sense): React.ReactElement {
  switch (kind) {
    case "pump":
      return (
        <>
          <line className="body" x1={0} y1={24} x2={23} y2={24} />
          <line className="body" x1={49} y1={24} x2={72} y2={24} />
          <circle className="body" cx={36} cy={24} r={13} />
          <polygon className="accent" points="35,18 47,24 35,30" />
        </>
      );

    case "machine": {
      // Mirrored, because the narrow end is where the pressure is made: a compressor squeezes into
      // its outlet, an expander opens out of it.
      const squeezes = sense !== "out";
      return (
        <>
          <line className="body" x1={0} y1={24} x2={12} y2={24} />
          <line className="body" x1={58} y1={24} x2={72} y2={24} />
          {/* The flange face at the outlet end, which both machines discharge through. */}
          <line className="body" x1={58} y1={6} x2={58} y2={42} />
          <polygon
            className="body"
            points={squeezes ? "12,8 58,16 58,32 12,40" : "12,16 58,8 58,40 12,32"}
          />
        </>
      );
    }

    case "valve":
      return (
        <>
          <line className="body" x1={0} y1={24} x2={20} y2={24} />
          <line className="body" x1={52} y1={24} x2={72} y2={24} />
          <polygon className="body" points="20,12 36,24 20,36" />
          <polygon className="body" points="52,12 36,24 52,36" />
        </>
      );

    case "duty": {
      // Mirrored: a heater's arrow rises from below, a cooler's falls from above.
      const heats = sense !== "out";
      return (
        <>
          <line className="body" x1={0} y1={24} x2={14} y2={24} />
          <line className="body" x1={58} y1={24} x2={72} y2={24} />
          <rect className="body" x={14} y={8} width={44} height={32} />
          <polygon
            className="accent"
            points={heats ? "30,20 42,20 36,10" : "30,28 42,28 36,38"}
          />
          <polygon
            className="accent"
            points={heats ? "34,20 38,20 38,38 34,38" : "34,10 38,10 38,28 34,28"}
          />
        </>
      );
    }

    case "pipe":
      return (
        <>
          <line className="body" x1={0} y1={17} x2={72} y2={17} />
          <line className="body" x1={0} y1={31} x2={72} y2={31} />
          <line className="detail" x1={36} y1={17} x2={36} y2={31} />
        </>
      );

    case "filter":
      return (
        <>
          <line className="body" x1={0} y1={24} x2={14} y2={24} />
          <line className="body" x1={58} y1={24} x2={72} y2={24} />
          <rect className="body" x={14} y={8} width={44} height={32} />
          <line className="detail" x1={16} y1={38} x2={26} y2={10} />
          <line className="detail" x1={26} y1={38} x2={36} y2={10} />
          <line className="detail" x1={36} y1={38} x2={46} y2={10} />
          <line className="detail" x1={46} y1={38} x2={56} y2={10} />
        </>
      );

    case "mixer":
      return <polygon className="body" points="0,4 0,44 72,24" />;

    case "splitter":
      return <polygon className="body" points="72,4 72,44 0,24" />;

    case "ejector":
      return (
        <polygon className="body" points="0,6 26,20 40,20 72,8 72,40 40,28 26,28 0,42" />
      );

    case "reactor":
      return (
        <>
          <line className="body" x1={0} y1={24} x2={22} y2={24} />
          <line className="body" x1={50} y1={24} x2={72} y2={24} />
          <rect className="body" x={22} y={2} width={28} height={44} rx={6} />
          <line className="detail" x1={28} y1={24} x2={44} y2={24} />
          <polyline className="detail" points="31,21 28,24 31,27" />
          <polyline className="detail" points="41,21 44,24 41,27" />
        </>
      );

    case "stirred":
      return (
        <>
          <line className="body" x1={0} y1={24} x2={20} y2={24} />
          <line className="body" x1={52} y1={24} x2={72} y2={24} />
          <rect className="body" x={20} y={6} width={32} height={40} rx={5} />
          <line className="detail" x1={36} y1={6} x2={36} y2={36} />
          <line className="detail" x1={27} y1={36} x2={45} y2={36} />
        </>
      );

    case "tube":
      return (
        <>
          <rect className="body" x={0} y={14} width={72} height={20} rx={10} />
          <line className="detail" x1={14} y1={24} x2={58} y2={24} />
          <polyline className="detail" points="52,20 58,24 52,28" />
        </>
      );

    case "separator":
      return (
        <>
          <line className="body" x1={0} y1={24} x2={20} y2={24} />
          <line className="body" x1={36} y1={0} x2={36} y2={4} />
          <line className="body" x1={36} y1={44} x2={36} y2={48} />
          <rect className="body" x={20} y={4} width={32} height={40} rx={6} />
        </>
      );

    case "scrubber":
      return (
        <>
          <line className="body" x1={0} y1={24} x2={20} y2={24} />
          <line className="body" x1={36} y1={0} x2={36} y2={4} />
          <line className="body" x1={36} y1={44} x2={36} y2={48} />
          <rect className="body" x={20} y={4} width={32} height={40} rx={6} />
          <line className="detail" x1={24} y1={12} x2={48} y2={12} />
          <line className="detail" x1={24} y1={16} x2={48} y2={16} />
        </>
      );

    case "column":
      return (
        <>
          <line className="body" x1={0} y1={24} x2={22} y2={24} />
          <line className="body" x1={36} y1={0} x2={36} y2={2} />
          <line className="body" x1={36} y1={46} x2={36} y2={48} />
          <rect className="body" x={22} y={2} width={28} height={44} rx={4} />
          <line className="detail" x1={24} y1={9} x2={48} y2={9} />
          <line className="detail" x1={24} y1={14} x2={48} y2={14} />
          <line className="detail" x1={24} y1={19} x2={48} y2={19} />
          <line className="detail" x1={24} y1={24} x2={48} y2={24} />
          <line className="detail" x1={24} y1={29} x2={48} y2={29} />
          <line className="detail" x1={24} y1={34} x2={48} y2={34} />
          <line className="detail" x1={24} y1={39} x2={48} y2={39} />
        </>
      );

    case "packed":
      return (
        <>
          <line className="body" x1={0} y1={24} x2={22} y2={24} />
          <line className="body" x1={36} y1={0} x2={36} y2={2} />
          <line className="body" x1={36} y1={46} x2={36} y2={48} />
          <rect className="body" x={22} y={2} width={28} height={44} rx={4} />
          {/* The internals are the difference from a tray column, not the outline. */}
          <line className="detail" x1={24} y1={10} x2={32} y2={18} />
          <line className="detail" x1={32} y1={10} x2={40} y2={18} />
          <line className="detail" x1={40} y1={10} x2={48} y2={18} />
          <line className="detail" x1={32} y1={10} x2={24} y2={18} />
          <line className="detail" x1={40} y1={10} x2={32} y2={18} />
          <line className="detail" x1={48} y1={10} x2={40} y2={18} />
          <line className="detail" x1={24} y1={22} x2={32} y2={30} />
          <line className="detail" x1={32} y1={22} x2={40} y2={30} />
          <line className="detail" x1={40} y1={22} x2={48} y2={30} />
          <line className="detail" x1={32} y1={22} x2={24} y2={30} />
          <line className="detail" x1={40} y1={22} x2={32} y2={30} />
          <line className="detail" x1={48} y1={22} x2={40} y2={30} />
          <line className="detail" x1={24} y1={34} x2={32} y2={42} />
          <line className="detail" x1={32} y1={34} x2={40} y2={42} />
          <line className="detail" x1={40} y1={34} x2={48} y2={42} />
          <line className="detail" x1={32} y1={34} x2={24} y2={42} />
          <line className="detail" x1={40} y1={34} x2={32} y2={42} />
          <line className="detail" x1={48} y1={34} x2={40} y2={42} />
        </>
      );

    case "exchanger":
      return (
        <>
          <line className="body" x1={0} y1={24} x2={17} y2={24} />
          <line className="body" x1={59} y1={24} x2={72} y2={24} />
          <circle className="body" cx={26} cy={24} r={9} />
          <circle className="body" cx={50} cy={24} r={9} />
          <line className="body" x1={26} y1={15} x2={50} y2={15} />
          <line className="body" x1={26} y1={33} x2={50} y2={33} />
        </>
      );

    case "flare":
      return (
        <>
          <line className="body" x1={0} y1={44} x2={31} y2={44} />
          <line className="body" x1={36} y1={0} x2={36} y2={4} />
          <polygon className="body" points="28,48 31,44 32,12 40,12 41,44 44,48" />
          {/* The flame sits on the stack's mouth, and the flue gas leaves through its tip. */}
          <polygon className="accent" points="36,4 40,12 32,12" />
        </>
      );

    case "manifold":
      return (
        <>
          <rect className="body" x={0} y={18} width={72} height={8} />
          <line className="body" x1={18} y1={26} x2={18} y2={38} />
          <line className="body" x1={36} y1={26} x2={36} y2={38} />
          <line className="body" x1={54} y1={26} x2={54} y2={38} />
        </>
      );

    case "tank":
      return (
        <>
          <line className="body" x1={0} y1={24} x2={18} y2={24} />
          <line className="body" x1={36} y1={0} x2={36} y2={3} />
          <line className="body" x1={36} y1={45} x2={36} y2={48} />
          <rect className="body" x={18} y={8} width={36} height={32} />
          {/* Flat caps: the vessel stands on its axis, which is the side its ports leave from. */}
          <ellipse className="body" cx={36} cy={8} rx={18} ry={5} />
          <ellipse className="body" cx={36} cy={40} rx={18} ry={5} />
        </>
      );

    case "fallback":
      return <rect className="body" x={6} y={4} width={60} height={40} rx={8} />;
  }
}
