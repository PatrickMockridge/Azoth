/**
 * The unit a number is read in, and the one multiplication that gets it there.
 *
 * **The front end carries no conversion table.** A quantity arrives as `{magnitude_si, unit}`, and
 * the two facts needed to show it another way - which dimension that unit measures, and what the
 * unit a reader wants is worth - are both in the catalogue, computed by the library from the same
 * conversion a calculation runs. So the arithmetic here is a division and nothing else, and a unit
 * the vocabulary gains next year appears in the switcher without a line changing.
 *
 * **A dimension no set names is shown in the unit the library computed it in.** That is the
 * fallback rather than a blank, and it is deliberate: most of the vocabulary has no engineering
 * alternative - a `Pa*m**6/mol**2` attraction parameter is read in its own unit or in none - so the
 * sets name the dimensions a person actually switches.
 *
 * **Nothing is on the wire about which set is in force.** It is a reader's choice, like the theme,
 * and it is remembered the same way: one `localStorage` key and no document, because a flowsheet
 * opened by somebody else has no opinion about who is looking at it.
 */

import type { Catalogue } from "../wire/types";

/** One declared unit, as the catalogue carries it. */
export interface UnitDecl {
  dimension: string | null;
  factor: number | null;
  /**
   * The affine constant, in the same base unit `factor` is measured from - and zero for every
   * unit that is a scale, which is all but a display unit.
   *
   * **`si = (value + offset) * factor` is uom's convention and this is its inverse**:
   * `value = si / factor - offset`. One expression covers both kinds, so nothing here branches on
   * whether a unit happens to be shifted.
   */
  offset: number | null;
}

/** One named set: a unit per dimension, and what to call it. */
export interface UnitSetDecl {
  id: string;
  name: string;
  units: Record<string, string>;
}

/** The library's units and sets, and which set is in force. */
export interface Units {
  /** Every declared unit, keyed by the string a quantity carries. */
  declared: Record<string, UnitDecl>;
  sets: UnitSetDecl[];
  /** The id of the set in force, which may name one the catalogue no longer carries. */
  active: string;
}

/** Where the choice is remembered. One key, and nothing else in the app writes one. */
const KEY = "azoth.units";

/** The set a person last chose, or none. */
export function readUnitSet(): string | null {
  try {
    return window.localStorage.getItem(KEY);
  } catch {
    // A browser can refuse storage (a private window, a policy). The library's first set is still
    // the one shown, and a choice that cannot be remembered is not a fault worth a message.
    return null;
  }
}

/** Remember the set a person chose. */
export function rememberUnitSet(id: string): void {
  try {
    window.localStorage.setItem(KEY, id);
  } catch {
    // As above: the choice is in force either way, and the store is only how it survives a reload.
  }
}

/**
 * The library's units, with the set in force resolved against them.
 *
 * A remembered id the catalogue does not carry - a vocabulary that lost a set, a store written by
 * an older build - falls back to the first set rather than to no conversion, because a reader who
 * asked for field units and gets SI is better served than one who gets a number with no unit.
 */
export function unitsOf(catalogue: Catalogue | null, wanted: string | null): Units {
  const sets = catalogue?.unit_sets ?? [];
  const first = sets[0]?.id ?? "";
  const named = sets.some((set) => set.id === wanted) ? (wanted as string) : first;
  return { declared: catalogue?.units ?? {}, sets, active: named };
}

/**
 * The unit one quantity is read in, and the factor its SI magnitude is divided by.
 *
 * `null` factors and unknown units return the quantity's own unit with a factor of one, which
 * makes every absence the same absence: what the library wrote, unchanged.
 */
export function displayOf(
  units: Units,
  unit: string,
): { unit: string; factor: number; offset: number } {
  const same = { unit, factor: 1, offset: 0 };
  const declared = units.declared[unit];
  if (declared?.dimension == null) {
    return same;
  }
  const set = units.sets.find((candidate) => candidate.id === units.active);
  const wanted = set?.units[declared.dimension];
  if (wanted === undefined) {
    return same;
  }
  const target = units.declared[wanted];
  const factor = target?.factor;
  // A factor of zero would be a division by zero rather than a unit, and a negative one is not a
  // unit at all - neither is a conversion the library declares, so both mean "leave it alone".
  if (factor === undefined || factor === null || !(factor > 0)) {
    return same;
  }
  // **The offset is the *target* unit's**, because it is the unit being shown in that decides
  // where its zero sits - a degree Celsius is 273.15 above a kelvin, whatever the quantity
  // arrived in.
  return { unit: wanted, factor, offset: target?.offset ?? 0 };
}

/**
 * The number to draw for a magnitude the library wrote in SI, in the unit it is read in.
 *
 * **One formula for both kinds of unit**: a scale's offset is zero, so `si / factor - offset`
 * is the division it always was where a set says nothing about scale, and where a set says
 * `°C` it is the affine inverse as well.
 */
export function inDisplayUnit(units: Units, magnitudeSi: number, unit: string): number {
  const shown = displayOf(units, unit);
  return magnitudeSi / shown.factor - shown.offset;
}
