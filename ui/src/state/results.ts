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

import { displayOf, unitLabelOf, type Units } from "./units";
import type { TabId } from "./layout";
import { humanLabel, explanationOf } from "./labels";
import { fieldOf, isQuantity } from "../wire/decode";
import type { CalcOutput, Envelope, Quantity, UnitResult } from "../wire/types";

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
export function scalarsOf(result: UnitResult, units?: Units): Scalar[] {
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
      // The same division `formatQuantity` does, kept here as a number and a unit because a
      // result sheet draws them in two cells.
      const read = shown(units, value.unit);
      out.push({
        key,
        text: (value.magnitude_si / read.factor).toFixed(6),
        unit: read.unit,
        numeric: true,
      });
    }
  }
  return out;
}

/**
 * One result field as a person reads it, and as the record spells it.
 *
 * **Two names and one value.** `key` is the model's own field name and goes in the
 * machine-readable appendix; `label` is the human summary's. The value is formatted but never
 * recomputed, and a field the run did not publish is not a row at all.
 */
export interface OutputRow {
  key: string;
  label: string;
  text: string;
  unit: string | null;
  numeric: boolean;
  /** The spec's sentence about the field, often the equation behind it. */
  description: string | null;
  /** The declared kind, or `null` for a field the form does not declare. */
  kind: string | null;
  /** A published record list's row count, for the fields that carry rows rather than a value. */
  records: number | null;
}

/**
 * The published fields of a result, in the model's declared order.
 *
 * **Only what exists.** A `null`, an `undefined` or an empty array is not a row - hx1 declares
 * 76 fields and publishes 8, and the 64 electrical and mechanical fields that were never
 * configured are absent rather than drawn as dashes. `false` and `0` are values and stay.
 *
 * The declared order is the spec's, which is the order the model's author put the record in;
 * a field the form does not declare (an older catalogue, a newer kernel) is appended in the
 * result's own order rather than dropped. Compositions are the streams' business and numeric
 * vectors are the Stages sheet's, so neither is a scalar row here.
 */
export function outputRowsOf(
  result: UnitResult,
  units?: Units,
  outputs?: readonly CalcOutput[],
): OutputRow[] {
  const declared = new Map((outputs ?? []).map((output) => [output.name, output]));
  const order: string[] = [];
  for (const output of outputs ?? []) {
    if (output.name !== "warnings" && output.name in result) {
      order.push(output.name);
    }
  }
  for (const key of Object.keys(result)) {
    if (key === "warnings" || declared.has(key)) {
      continue;
    }
    order.push(key);
  }

  const rows: OutputRow[] = [];
  for (const key of order) {
    if (isComposition(key)) {
      continue;
    }
    const row = outputRow(key, result[key], declared.get(key), units);
    if (row !== null) {
      rows.push(row);
    }
  }
  return rows;
}

function outputRow(
  key: string,
  value: unknown,
  output: CalcOutput | undefined,
  units: Units | undefined,
): OutputRow | null {
  if (value === null || value === undefined) {
    return null;
  }
  const base = {
    key,
    label: humanLabel(key),
    description: explanationOf(output?.description),
    kind: output?.kind ?? null,
    records: null,
  };
  if (typeof value === "number") {
    const read = readUnit(units, output?.unit ?? null);
    return {
      ...base,
      text: formatNumber(value / read.factor),
      unit: read.label,
      numeric: true,
    };
  }
  if (typeof value === "boolean") {
    // The record says `true`; a person reads yes. The appendix keeps the record's spelling.
    return { ...base, text: value ? "yes" : "no", unit: null, numeric: false };
  }
  if (typeof value === "string") {
    return { ...base, text: value, unit: null, numeric: false };
  }
  if (isQuantity(value)) {
    const read = readUnit(units, value.unit);
    return {
      ...base,
      text: formatNumber(value.magnitude_si / read.factor),
      unit: read.label,
      numeric: true,
    };
  }
  if (!Array.isArray(value)) {
    return null;
  }
  if (value.length === 0) {
    return null;
  }
  // A numeric or quantity vector is a profile, which the Stages sheet owns. A vector of names
  // and a record list are not profiles: the first is a set a reader can see, the second is a
  // table the machine appendix carries.
  if (value.every((entry) => typeof entry === "number" || isQuantity(entry))) {
    return null;
  }
  const records = value.every((entry) => typeof entry === "object");
  return {
    ...base,
    text: records
      ? `${String(value.length)} ${value.length === 1 ? "record" : "records"}`
      : value.map((entry) => (typeof entry === "string" ? entry : "…")).join(", "),
    unit: null,
    numeric: false,
    records: records ? value.length : null,
  };
}

/** A declared unit as the reader's unit set shows it, with the factor a magnitude divides by. */
function readUnit(
  units: Units | undefined,
  unit: string | null,
): { label: string | null; factor: number } {
  if (unit === null) {
    return { label: null, factor: 1 };
  }
  const read = units === undefined ? { unit, factor: 1 } : displayOf(units, unit);
  return { label: unitLabelOf(units ?? null, read.unit), factor: read.factor };
}

/**
 * A number as a result sheet shows it: six significant figures, no trailing zeros.
 *
 * Not rounded into a claim: the value's own `f64` is what the library wrote, and the digits are
 * a formatting choice. A non-finite number is shown as itself rather than hidden, because the
 * library uses `NaN` for a state it could not reach.
 */
export function formatNumber(value: number): string {
  if (!Number.isFinite(value)) {
    return String(value);
  }
  return String(Number(value.toPrecision(6)));
}

/** The vectors of a result - the profiles, tray columns and histories - with their units. */
export function seriesOf(result: UnitResult, units?: Units): Series[] {
  const out: Series[] = [];
  for (const [key, value] of Object.entries(result)) {
    if (isComposition(key) || !Array.isArray(value) || value.length === 0) {
      continue;
    }
    if (value.every((entry) => typeof entry === "number")) {
      out.push({ key, values: value as number[], unit: null });
    } else if (allQuantities(value)) {
      const quantities = value;
      // **One factor for the column**, taken from the first entry: a profile is one dimension, and
      // a vector whose entries disagreed about theirs would be a defect the library cannot write.
      const first = quantities[0];
      const read = first === undefined ? { unit: null, factor: 1 } : shown(units, first.unit);
      out.push({
        key,
        values: quantities.map((entry) => entry.magnitude_si / read.factor),
        unit: read.unit,
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
    const code = fieldOf(warning, "code");
    const message = fieldOf(warning, "message");
    const field = fieldOf(warning, "field");
    return [
      {
        code: typeof code === "string" ? code : "warning",
        message: typeof message === "string" ? message : "",
        field: typeof field === "string" ? field : null,
      },
    ];
  });
}

/**
 * One quantity as it is read: its unit and the factor it is divided by.
 *
 * **A caller with no set gets the library's own unit and a factor of one**, which is what makes the
 * unit set a reading of a result rather than a requirement on one: a caller that has no catalogue
 * yet - the first render of either door - sees exactly what the run wrote.
 */
function shown(units: Units | undefined, unit: string): { unit: string; factor: number } {
  return units === undefined ? { unit, factor: 1 } : displayOf(units, unit);
}

/** Whether a value is the codec's `{magnitude_si, unit}` shape. */
/** `Array.prototype.every` keeps its predicate's narrowing, which the built-in cannot express. */
function allQuantities(values: unknown[]): values is Quantity[] {
  return values.every(isQuantity);
}

/**
 * The sheets a result adds to its own window, in the order they are shown.
 *
 * **A sheet is offered iff it has something to draw**, which is the rule that makes every absence
 * fall out rather than be asked about: a document that has not run has no result and so no result
 * sheets, and a unit operation that published only scalars has no stage table. The profile comes
 * first because it is the shape of the answer and the scalars are its summary.
 */
export function resultTabs(
  result: UnitResult | undefined,
  outputs?: readonly CalcOutput[],
): TabId[] {
  if (result === undefined) {
    return [];
  }
  const tabs: TabId[] = [];
  if (seriesOf(result).length > 0) {
    tabs.push("stages");
  }
  // A result whose every field is absent has nothing to draw, so it offers no Results sheet.
  if (outputRowsOf(result, undefined, outputs).length > 0) {
    tabs.push("results");
  }
  return tabs;
}
