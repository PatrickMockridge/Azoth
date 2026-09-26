/**
 * What a unit operation published, read without knowing which unit operation it is.
 *
 * **The shape is the library's and not this front end's**, which is why nothing here enumerates it:
 * a result is whatever that operation's model declares, and the Rust half serialises it under the
 * field names the Python half publishes under. So the reader classifies *values* rather than looking
 * up units operation by operation - a number, a quantity, a vector of either - and a field added to
 * a model next year appears without a line here.
 *
 * **Nothing is computed.** A scalar is the number the kernel reached, a series is the vector it
 * filled in, and a warning is the library's own sentence. The one thing this module does decide is
 * *how many digits* a magnitude shows, which is a formatting choice and not arithmetic.
 */

import type { TabId } from "./layout";
import type { Envelope, Quantity, UnitResult } from "../wire/types";

/** One scalar a result carries: a number, a name, a switch, or a quantity with its unit. */
export interface Scalar {
  key: string;
  /** What to draw, already formatted. */
  text: string;
  /** The unit, where the field carries one. */
  unit: string | null;
  /** Whether the cell is a number, which right-aligns it. */
  numeric: boolean;
}

/** One vector a result carries: a profile, a set of tray values, a history. */
export interface Series {
  key: string;
  values: number[];
  unit: string | null;
}

/** One warning, as the wire writes it. */
export interface ResultWarning {
  code: string;
  message: string;
  field: string | null;
}

/** The result one instance published, or `undefined` where it published none. */
export function resultOf(envelope: Envelope, instance: string): UnitResult | undefined {
  return envelope.session?.results[instance];
}

/**
 * The scalars of a result, in the model's own order.
 *
 * A `null` is a row with a dash rather than a missing row: a skipped solve has no residual, and a
 * reader comparing two runs needs to see that it is absent rather than not look for it.
 */
export function scalarsOf(result: UnitResult): Scalar[] {
  const out: Scalar[] = [];
  for (const [key, value] of Object.entries(result)) {
    if (key === "warnings") {
      continue;
    }
    if (typeof value === "number") {
      out.push({ key, text: value.toFixed(6), unit: null, numeric: true });
    } else if (typeof value === "boolean") {
      out.push({ key, text: String(value), unit: null, numeric: false });
    } else if (typeof value === "string") {
      out.push({ key, text: value, unit: null, numeric: false });
    } else if (value === null) {
      out.push({ key, text: "—", unit: null, numeric: false });
    } else if (isQuantity(value)) {
      out.push({
        key,
        text: value.magnitude_si.toFixed(6),
        unit: value.unit,
        numeric: true,
      });
    }
  }
  return out;
}

/** The vectors of a result - the profiles, tray columns and histories - with their units. */
export function seriesOf(result: UnitResult): Series[] {
  const out: Series[] = [];
  for (const [key, value] of Object.entries(result)) {
    if (isComposition(key) || !Array.isArray(value) || value.length === 0) {
      continue;
    }
    if (value.every((entry) => typeof entry === "number")) {
      out.push({ key, values: value as number[], unit: null });
    } else if (value.every(isQuantity)) {
      const quantities = value as Quantity[];
      out.push({
        key,
        values: quantities.map((entry) => entry.magnitude_si),
        unit: quantities[0]?.unit ?? null,
      });
    }
  }
  return out;
}

/**
 * Whether a field is a composition rather than a profile.
 *
 * **The library names every one of them `*_z`** - the record's own `z`, a port's `outlet_z`, a
 * column's `distillate_z` - and their index is a *substance* rather than a stage, so a sparkline of
 * one is a picture of nothing. A heater's result carries `outlet_z`, which is why this rule exists:
 * without it every result with a composition offered a stage table.
 */
function isComposition(key: string): boolean {
  return key === "z" || key.endsWith("_z");
}

/** The warnings a result carries, in the library's own words. */
export function warningsOf(result: UnitResult): ResultWarning[] {
  const warnings = result["warnings"];
  if (!Array.isArray(warnings)) {
    return [];
  }
  return warnings.flatMap((warning) => {
    if (typeof warning !== "object" || warning === null) {
      return [];
    }
    const { code, message, field } = warning as Partial<ResultWarning>;
    return [
      {
        code: typeof code === "string" ? code : "warning",
        message: typeof message === "string" ? message : "",
        field: typeof field === "string" ? field : null,
      },
    ];
  });
}

/** Whether a value is the codec's `{magnitude_si, unit}` shape. */
function isQuantity(value: unknown): value is Quantity {
  return (
    typeof value === "object" &&
    value !== null &&
    "magnitude_si" in value &&
    typeof (value as Quantity).magnitude_si === "number"
  );
}

/**
 * The sheets a result adds to its own window, in the order they are shown.
 *
 * **A sheet is offered iff it has something to draw**, which is the rule that makes every absence
 * fall out rather than be asked about: a document that has not run has no result and so no result
 * sheets, and a unit operation that published only scalars has no stage table. The profile comes
 * first because it is the shape of the answer and the scalars are its summary.
 */
export function resultTabs(result: UnitResult | undefined): TabId[] {
  if (result === undefined) {
    return [];
  }
  const tabs: TabId[] = [];
  if (seriesOf(result).length > 0) {
    tabs.push("stages");
  }
  if (scalarsOf(result).length > 0) {
    tabs.push("results");
  }
  return tabs;
}
