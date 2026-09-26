/**
 * A form parameter, as a control and a value.
 *
 * **The kind comes from the library and this file only reads it.** `Param` in the palette carries
 * a unit and no type, so the kind is the model's own input declaration compiled into the form —
 * which means the widget a field gets and the value the checker will accept are decided in one
 * place and not two. A kind the library cannot name arrives as `"unknown"`, and a field that
 * cannot choose a widget says so rather than guessing at a number box.
 */

import { displayOf, inDisplayUnit, toDisplayUnit, type Units } from "../state/units";
import type { Catalogue, FormParameter, FormRange, Kind } from "./types";

/** What a field renders as. */
export type Control = "number" | "list" | "switch" | "select" | "text" | "unknown";

/** The control a kind gets. `components` and `matrix` are declared and never a palette parameter. */
export function controlFor(kind: Kind): Control {
  switch (kind) {
    case "quantity":
      return "number";
    case "vector":
      return "list";
    case "boolean":
      return "switch";
    case "enum":
      return "select";
    case "string":
      return "text";
    default:
      return "unknown";
  }
}

/** A field's bounds, for the input's own `min`/`max` and for what to say on a violation. */
export interface Bounds {
  min: number | null;
  max: number | null;
  /**
   * Whether each endpoint is excluded — **each its own**, because the models' bounds are.
   * An efficiency is `(0, 1]`: one is the reversible limit and legal, and a single flag for
   * both ends would mark the capture's own efficiency-of-one row as a violation.
   */
  minExclusive: boolean;
  maxExclusive: boolean;
  /** Why the bound exists, which is the model's own sentence. */
  rationale: string | null;
  /** What a violation means: an error refuses the run, a warning does not. */
  severity: "error" | "warning" | null;
}

/**
 * The bound a parameter's range gives, or none.
 *
 * **The library's ranges are a list and usually one entry long.** A parameter can carry more than
 * one — a forbidden band inside an allowed interval is two — and a field draws the first, because
 * an input has one `min` and one `max`. The rationale of the one it drew travels with them, so a
 * violation is explained by the bound that caused it.
 */
export function boundsOf(parameter: FormParameter): Bounds {
  const range: FormRange | undefined = parameter.range[0];
  if (range === undefined) {
    return {
      min: null,
      max: null,
      minExclusive: false,
      maxExclusive: false,
      rationale: null,
      severity: null,
    };
  }
  return {
    min: range.min,
    max: range.max,
    minExclusive: !range.min_inclusive,
    maxExclusive: !range.max_inclusive,
    rationale: range.rationale,
    severity: range.severity,
  };
}

/**
 * The value a command should carry for what a field holds.
 *
 * **The kind shapes it, so the document gets the value its declaration describes** rather than
 * whatever a text input happened to contain: a vector is split on commas, a switch is a boolean,
 * and a quantity is a number or, if it is not one, the text as typed — which the checker then
 * refuses with `parameter_kind` rather than the run refusing it later.
 */
export function commandValue(kind: Kind, raw: string): unknown {
  switch (controlFor(kind)) {
    case "number": {
      const number = Number(raw);
      return raw.trim() !== "" && Number.isFinite(number) ? number : raw;
    }
    case "list":
      return raw
        .split(",")
        .map((entry) => Number(entry.trim()))
        .filter((entry) => Number.isFinite(entry));
    case "switch":
      return raw === "true";
    default:
      return raw;
  }
}

/**
 * The acceleration names the tool schema publishes, which is the class's own list.
 *
 * **Read rather than spelled.** The three names live in one place — `recycle::ACCELERATION_NAMES`
 * — and the schema carries them; a widget that typed them again would be the copy that drifts,
 * and the one nobody parses. The schema says `value` is a number *or* one of these, so this reads
 * the second arm; an absent schema is an empty list and the caller says nothing.
 */
export function accelerationNames(catalogue: Catalogue | null): string[] {
  const tool = catalogue?.tools?.find((entry) => entry.name === "set_recycle");
  const schema = tool?.input_schema as
    | { properties?: { value?: { oneOf?: { enum?: unknown }[] } } }
    | undefined;
  const values = schema?.properties?.value?.oneOf?.[1]?.enum;
  return Array.isArray(values) ? values.filter((value): value is string => typeof value === "string") : [];
}

/** What a field shows for a value the document already holds. */
export function fieldText(value: unknown): string {
  if (value === undefined || value === null) {
    return "";
  }
  if (Array.isArray(value)) {
    return value.join(", ");
  }
  return String(value);
}

/**
 * What a field shows for a value the document holds, in the unit set in force.
 *
 * **A document's unit and a reader's are two different things, and a field is where they meet**:
 * the value is stored in the unit the spec declares and shown in the one the set names. So the
 * two conversions live here, beside each other, rather than in a widget - and a caller with no
 * set gets exactly what the document says, which is what keeps every test that has no catalogue
 * reading the shipped units.
 *
 * A vector is converted entry by entry: a tray-temperature profile is one `K` column, and a
 * field that converted the first entry and not the rest would be showing two units in one box.
 */
export function fieldTextIn(
  value: unknown,
  declared: string | null,
  units: Units | null,
): string {
  if (declared === null || units === null) {
    return fieldText(value);
  }
  if (typeof value === "number") {
    return shownNumber(inDisplayUnit(units, value, declared), value);
  }
  if (Array.isArray(value)) {
    return value
      .map((entry) =>
        typeof entry === "number" ? shownNumber(inDisplayUnit(units, entry, declared), entry) : entry,
      )
      .join(", ");
  }
  return fieldText(value);
}

/**
 * A converted number as a box shows it.
 *
 * **Rounded only where the conversion did something.** An input is controlled, so its text is
 * rewritten on every keystroke: a field that always rounded would fight somebody typing a tenth
 * digit, which is exactly why the SI set - whose conversion is the identity - shows `String(value)`
 * and has always behaved. A converted value has no such history to preserve, and the alternative
 * is a box reading `366.4833333333333` where the reader wants a temperature.
 */
function shownNumber(converted: number, original: number): string {
  return converted === original ? String(converted) : String(Number(converted.toPrecision(6)));
}

/**
 * The value a command should carry for what a field holds, converted back to the document's unit.
 *
 * [`commandValue`] shapes it by kind first, so a field cannot write a document the checker will
 * refuse for its type; this converts whatever number came out of that into the declared unit -
 * which is the whole point, because the document stores what its spec declares and a reader types
 * what the set names.
 */
export function commandValueIn(
  kind: Kind,
  raw: string,
  declared: string | null,
  units: Units | null,
): unknown {
  const shaped = commandValue(kind, raw);
  if (declared === null || units === null) {
    return shaped;
  }
  if (typeof shaped === "number") {
    return toDisplayUnit(units, shaped, declared);
  }
  if (Array.isArray(shaped)) {
    return shaped.map((entry) =>
      typeof entry === "number" ? toDisplayUnit(units, entry, declared) : entry,
    );
  }
  return shaped;
}

/**
 * A scalar with a unit, shortened for a label: 2.00e+6 Pa, 477.96 K.
 *
 * **This is the one place a quantity becomes text**, which is what makes the unit set a parameter
 * rather than a second formatter: a caller with a set passes it and gets the number and the unit
 * a reader asked for, and a caller with none - a bound's hint under a parameter, a value the
 * reader cannot change - passes none and gets the library's own. Nothing here converts; the
 * division is `state/units.ts`'s and the factor is the catalogue's.
 */
export function formatQuantity(
  magnitude: number,
  unit: string,
  digits = 4,
  units?: Units,
): string {
  const shown = units === undefined ? { unit, factor: 1, offset: 0 } : displayOf(units, unit);
  const rounded = Number((magnitude / shown.factor - shown.offset).toPrecision(digits));
  return `${rounded} ${shown.unit}`;
}

/**
 * Whether a value is outside a bound, and how to say so.
 *
 * **The library's own comparison, not a second opinion about the range.** A field marks a value
 * the run would refuse; it does not decide that the value is wrong.
 */
export function violates(bounds: Bounds, value: unknown): boolean {
  if (bounds.severity === null || typeof value !== "number") {
    return false;
  }
  const { min, max, minExclusive, maxExclusive } = bounds;
  if (min !== null && (minExclusive ? value <= min : value < min)) {
    return true;
  }
  if (max !== null && (maxExclusive ? value >= max : value > max)) {
    return true;
  }
  return false;
}
