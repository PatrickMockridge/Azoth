/**
 * A form parameter, as a control and a value.
 *
 * **The kind comes from the library and this file only reads it.** `Param` in the palette carries
 * a unit and no type, so the kind is the model's own input declaration compiled into the form —
 * which means the widget a field gets and the value the checker will accept are decided in one
 * place and not two. A kind the library cannot name arrives as `"unknown"`, and a field that
 * cannot choose a widget says so rather than guessing at a number box.
 */

import type { FormParameter, FormRange, Kind } from "./types";

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

/** A scalar with a unit, shortened for a label: 2.00e+6 Pa, 477.96 K. */
export function formatQuantity(magnitude: number, unit: string, digits = 4): string {
  const rounded = Number(magnitude.toPrecision(digits));
  return `${rounded} ${unit}`;
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
