/**
 * What each unit operation is drawn as, and which side of it a port leaves from.
 *
 * **Two tables, and both are keyed on the leaf of `unit_ops.<leaf>`.** The leaf is unique across
 * all twenty-nine entries, which the *id* is not useful for: every id is `unit_ops.<thing>`, so
 * splitting it takes the thing for a family. The family is on the wire now (`Form.family`, from the
 * palette directory) and this is deliberately not that - a pump and a compressor are both
 * `two_port` and are drawn as different machines.
 *
 * **A drawing is a claim about the machinery**, so the two tables are held to the catalogue in both
 * directions by `test/glyphs.test.ts`: a new palette entry with no drawing fails, and a drawing for
 * an entry that no longer exists fails too.
 */

/** What a unit operation is drawn as. */
export type GlyphKind =
  | "pump"
  | "machine"
  | "valve"
  | "duty"
  | "pipe"
  | "filter"
  | "mixer"
  | "splitter"
  | "ejector"
  | "reactor"
  | "stirred"
  | "tube"
  | "separator"
  | "scrubber"
  | "column"
  | "packed"
  | "exchanger"
  | "flare"
  | "manifold"
  | "tank"
  | "fallback";

/** The direction a drawing's moving part points, where it has one. */
export type Sense = "in" | "out" | "none";

/**
 * Every palette entry, by the leaf of its id.
 *
 * `pump` is a circle with a discharge triangle and a compressor is a trapezoid, because that is
 * what they are; `column` carries tray lines and `packed` a hatch, because the difference between a
 * tray column and a packed one is the internals and not the outline.
 */
export const GLYPHS: Readonly<Record<string, GlyphKind>> = {
  absorption_column: "column",
  component_splitter: "splitter",
  compressor: "machine",
  cooler: "duty",
  distillation_column: "column",
  ejector: "ejector",
  expander: "machine",
  filter: "filter",
  flare: "flare",
  gibbs_reactor: "reactor",
  gas_scrubber: "scrubber",
  heat_exchanger: "exchanger",
  heater: "duty",
  manifold: "manifold",
  mixer: "mixer",
  packed_column: "packed",
  pipe: "pipe",
  plug_flow_reactor: "tube",
  pump: "pump",
  rate_based_packed_column: "packed",
  separator: "separator",
  shortcut_distillation_column: "column",
  simple_absorber: "column",
  splitter: "splitter",
  stirred_tank_reactor: "stirred",
  stripping_column: "column",
  tank: "tank",
  three_phase_separator: "separator",
  throttling_valve: "valve",
};

/** Which entries are drawn as a machine that moves a fluid one way rather than the other. */
export const SENSE: Readonly<Record<string, Sense>> = {
  compressor: "in",
  expander: "out",
  cooler: "out",
  heater: "in",
};

/** The leaf of a palette id, which is what both tables are keyed on. */
export function leafOf(unit: string | undefined): string {
  return (unit ?? "").replace(/^unit_ops\./, "");
}

/** The drawing for a palette id, or the fallback box for one nothing here knows. */
export function glyphFor(unit: string | undefined): GlyphKind {
  return GLYPHS[leafOf(unit)] ?? "fallback";
}

/** The sense of a drawing, where it has one. */
export function senseFor(unit: string | undefined): Sense {
  return SENSE[leafOf(unit)] ?? "none";
}

/** Which side of a drawing a port leaves from. */
export type Side = "top" | "bottom" | "left" | "right";

/**
 * The ports a drawing does not put where its *direction* would.
 *
 * **The fallback is the declaration's own sense of direction**: an inlet is on the left and an
 * outlet is on the right, which is right for most machinery and wrong for anything a fluid passes
 * *through* vertically. A column's distillate leaves the top and its bottoms the bottom; a
 * separator's vapour leaves the top; a stack's flue gas leaves upward. Those are the entries here,
 * and nothing else is listed - so an override cannot quietly become the rule.
 */
export const PORT_SIDES: Readonly<Record<string, Readonly<Record<string, Side>>>> = {
  distillation_column: { distillate: "top", bottoms: "bottom" },
  shortcut_distillation_column: { distillate: "top", bottoms: "bottom" },
  packed_column: { distillate: "top", bottoms: "bottom" },
  rate_based_packed_column: { gas_out: "top", liquid_out: "bottom" },
  absorption_column: { gas_out: "top", liquid_out: "bottom" },
  stripping_column: { overhead_gas: "top", lean_liquid: "bottom" },
  simple_absorber: { gas: "top", liquid: "bottom" },
  gas_scrubber: { gas: "top", liquid: "bottom" },
  three_phase_separator: { vapour: "top", heavy_liquid: "bottom" },
  separator: { vapour: "top", liquid: "bottom" },
  tank: { gas: "top", liquid: "bottom" },
  flare: { product: "top" },
  plug_flow_reactor: { product: "right" },
};

/** Whether a side is a vertical one, which is what the drawing's own height has to accommodate. */
export function isVertical(side: Side): boolean {
  return side === "top" || side === "bottom";
}

/**
 * Where one declared port sits, given its name and the direction the declaration gives it.
 *
 * The override first, then the declaration: an inlet on the left, an outlet on the right.
 */
export function sideOf(unit: string, port: string, direction: "in" | "out"): Side {
  const declared = PORT_SIDES[leafOf(unit)]?.[port];
  if (declared !== undefined) {
    return declared;
  }
  return direction === "in" ? "left" : "right";
}
