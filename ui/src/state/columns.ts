/**
 * What a table's columns are, and where the substance axis comes from.
 *
 * **A solved stream's `z` is a bare array**, positional and unnamed: the library's record is
 * `n, z, P, T, h` and the names live on the *feed* that declared the fluid. So the axis has to be
 * recovered from the document, and the rule for that is this module's first content - with the one
 * case it must refuse to guess about.
 *
 * **A column is a declaration and not a table's private business**, which is what lets a quantity
 * added to the record years from now appear with a label rather than a blank: [`streamColumns`]
 * names the ones it knows and then draws any *unclaimed* key that is a quantity, labelled from its
 * own key. Nothing here computes a value a kernel computed.
 *
 * **A column's unit is the unit a reader will see**, which is the unit set's and not the record's:
 * the cell is divided by the same factor the header names, so the two cannot disagree about which
 * unit a column is in.
 */

import { displayOf, type Units } from "./units";
import type { Catalogue, Envelope, Quantity, StreamRecord } from "../wire/types";

/** The substance axis a composition table is laid out over. */
export interface Axis {
  /** One entry per column. */
  names: readonly string[];
  /**
   * Whether the names are **positions rather than names**.
   *
   * True when the feeds disagree about the fluid - two lists of the same length with different
   * names, or lists of different lengths - because then no feed's own list describes every stream
   * and a header drawn from one of them would be a claim the document does not make.
   */
  positional: boolean;
  /** Why, when it is positional; `null` otherwise. */
  note: string | null;
}

/**
 * The substances every stream of this document is made of, read from the feeds that declare them.
 *
 * **A feed is the only place the names are stated.** A product has no list of its own and a solved
 * stream's `z` is a bare array, so there is nothing else to read - and where the feeds disagree,
 * `positional` says so rather than silently preferring the first one.
 */
export function componentAxis(envelope: Envelope): Axis {
  const declared = envelope.flowsheet.graph.nodes
    .filter((node) => node.role === "input")
    .map((node) => node.data.input?.components ?? [])
    .filter((components) => components.length > 0);

  const first = declared[0];
  if (first === undefined) {
    return {
      names: [],
      positional: true,
      note: "no feed this document declares names its substances, so a composition has no names to show",
    };
  }
  if (declared.length === 1) {
    return { names: first, positional: false, note: null };
  }

  const agrees = declared.every(
    (list) => list.length === first.length && list.every((name, i) => name === first[i]),
  );
  if (agrees) {
    return { names: first, positional: false, note: null };
  }

  const longest = declared.reduce((a, b) => (b.length > a.length ? b : a));
  const lists = declared.map((list) => list.join(" / ")).join("; and ");
  return {
    names: longest,
    positional: true,
    note: `the feeds disagree about the fluid (${lists}), so these columns are positions and not names`,
  };
}

/** How a fractional column's cell reads, to four decimals - enough for a mole fraction. */
export function fraction(value: number): string {
  return value.toFixed(4);
}

/** One column of the workbook, and how its cells read. */
export interface Column {
  /** The record's own key, which is what a cell is looked up by. */
  key: string;
  label: string;
  /** The unit, where the column carries one; `null` for a bare number or a fraction. */
  unit: string | null;
  /** How many decimal places a cell shows. */
  digits: number;
}

/**
 * The columns of a stream record, in reading order.
 *
 * **The declared list first, then anything it did not claim.** The five the palette declares, then
 * the three the run adds, and then - for a record that carries something this file has never heard
 * of - a column labelled from the key itself. That last part is the difference between a front end
 * that shows what arrived and one that shows what it was written to expect.
 */
export function streamColumns(record: StreamRecord | undefined, units?: Units): Column[] {
  const declared: Column[] = [
    { key: "n", label: "flow", unit: "mol/s", digits: 4 },
    { key: "mass_flow", label: "mass flow", unit: "kg/s", digits: 5 },
    { key: "molar_mass", label: "M", unit: "kg/mol", digits: 5 },
    // **Four significant figures, and one would do for Pa alone.** A pressure in pascals is
    // around 1e5, so one figure reads as the same number; the same quantity in psi is around
    // 70, where one figure is 70 against 72.52. The unit set is what makes the coarse
    // reading wrong, and the column carries one precision for every unit it may be read in.
    { key: "P", label: "P", unit: "Pa", digits: 4 },
    { key: "T", label: "T", unit: "K", digits: 3 },
    { key: "h", label: "h", unit: "J/mol", digits: 3 },
    { key: "vapour_fraction", label: "VF", unit: null, digits: 4 },
  ];
  const shown = (unit: string | null) =>
    unit === null || units === undefined ? unit : displayOf(units, unit).unit;
  const reading = declared.map((column) => ({ ...column, unit: shown(column.unit) }));
  if (record === undefined) {
    return reading;
  }
  const claimed = new Set([...declared.map((column) => column.key), "z"]);
  const extra = Object.entries(record)
    .filter(([key, value]) => !claimed.has(key) && isQuantity(value))
    .map(([key, value]) => ({
      key,
      label: key,
      unit: shown((value as Quantity).unit),
      digits: 4,
    }));
  return [...reading, ...extra];
}

/** Whether a record's value is the codec's `{magnitude_si, unit}` shape. */
function isQuantity(value: unknown): boolean {
  return (
    typeof value === "object" &&
    value !== null &&
    "magnitude_si" in value &&
    "unit" in value
  );
}

/** One row of a workbook: what the row *is*, and the record it read, where a run reached it. */
export interface Row {
  /** The stream's path, which is also what selects it on the canvas. */
  readonly path: string;
  /** The node a click should select. */
  readonly nodeId: string;
  /** The unit op or boundary that produces it, for the row's second column. */
  readonly from: string;
  readonly record: StreamRecord | undefined;
}

/**
 * Every row the workbook draws, **from the endpoints the document declares and in its own order**.
 *
 * **Only the *producing* side gets a row, and that is what keeps the table worth reading.** A unit
 * op's `outlet` handle and the next one's `inlet` are the same physical stream - the envelope keys
 * its value by whichever endpoint produced it - so a table with both would carry every stream twice
 * and every second copy would be blank. An inlet is wiring; the value lives on the producer.
 *
 * **This is the rule that makes it a workbook rather than a run log.** The paths are what the
 * document declares, so the grid has the same shape before and after Solve and a stream the run did
 * not produce is a row of dashes rather than a missing row.
 *
 * **A stream the run produced that no port names is appended**, which is the safety net for a
 * `many` outlet whose projected handles are shorter than what the run returned. The shipped demo
 * has none: its tear is *not* one, because the envelope carries a tear's four residuals and not its
 * stream - that stream is the loop's accelerated copy of the edge's `from`, which is a row already,
 * and the residuals are the Tears sheet.
 */
export function streamRows(envelope: Envelope): Row[] {
  const streams = envelope.session?.streams ?? {};
  const rows: Row[] = [];
  const seen = new Set<string>();

  for (const node of envelope.flowsheet.graph.nodes) {
    const paths =
      node.role === "instance"
        ? (node.data.ports?.outlets ?? []).flatMap((port) => port.handles)
        : // A boundary is an endpoint in its own right: a feed produces its own record and a
          // product is where a stream is keyed by the name the document gave it.
          [node.data.name];
    for (const path of paths) {
      if (seen.has(path)) {
        continue;
      }
      seen.add(path);
      rows.push({ path, nodeId: node.id, from: node.data.name, record: streams[path] });
    }
  }

  for (const path of Object.keys(streams)) {
    if (!seen.has(path)) {
      seen.add(path);
      rows.push({ path, nodeId: "", from: "—", record: streams[path] });
    }
  }
  return rows;
}

/**
 * Every substance's molar mass, by name, in kg/mol - the databank's own numbers.
 *
 * A name the databank does not carry is absent from the map and a name it carries without a mass
 * maps to `null`; both mean the same thing downstream, which is that a mass cannot be formed.
 */
export type MolarMasses = Record<string, number | null>;

/** The catalogue's substances, keyed by name. */
export function molarMasses(catalogue: Catalogue | null): MolarMasses {
  const out: MolarMasses = {};
  for (const component of catalogue?.components ?? []) {
    out[component.name] = component.molar_mass;
  }
  return out;
}

/** One row of the composition sheet. */
export interface CompositionRow {
  readonly component: string;
  readonly path: string;
  readonly nodeId: string;
  readonly fraction: number;
  /** `n · z`, or `null` where the run reached no flow for this stream. */
  readonly flow: number | null;
  /**
   * `z_i · M_i / Σ z_j · M_j`, or `null` where any substance's mass is unknown.
   *
   * **`null` rather than a partial fraction**, because a fraction is a ratio to a total: a
   * denominator missing one term is not a smaller answer, it is a different one.
   */
  readonly mass_fraction: number | null;
  /** `n · z_i · M_i`, in kg/s, or `null` where *this* substance's mass is unknown. */
  readonly mass_flow: number | null;
}

/**
 * Every stream's substances, one row per stream per substance.
 *
 * **The one derived cell is converted in the derived unit.** `n · z` is a molar flow, so it is read
 * in the set's molar-flow unit and not in the record's - the multiplication happens in SI and the
 * division is the same one `formatQuantity` would have done, which is what keeps the sheet's `n · z`
 * column in the same unit as the workbook's `flow` column.
 */
export function compositionRows(
  envelope: Envelope,
  units?: Units,
  masses: MolarMasses = {},
): CompositionRow[] {
  const axis = componentAxis(envelope);
  const perSecond = units === undefined ? 1 : displayOf(units, "mol/s").factor;
  const rows: CompositionRow[] = [];
  for (const row of streamRows(envelope)) {
    const z = row.record?.z;
    if (z === undefined) {
      continue;
    }
    const n = row.record?.n.magnitude_si;
    // The weights this stream's substances have, in its own axis order - so a substance the
    // databank does not carry is `null` here and takes both mass columns with it.
    const weights = z.map((_, index) => masses[axis.names[index] ?? ""] ?? null);
    const total = weights.every((mass) => mass !== null)
      ? z.reduce((sum, value, index) => sum + value * (weights[index] as number), 0)
      : null;
    z.forEach((value, index) => {
      const mass = weights[index];
      rows.push({
        component: axis.names[index] ?? `z[${index}]`,
        path: row.path,
        nodeId: row.nodeId,
        fraction: value,
        flow: n === undefined ? null : (value * n) / perSecond,
        mass_fraction: total === null || total === 0 ? null : (value * (mass as number)) / total,
        mass_flow: mass === null || mass === undefined || n === undefined ? null : value * n * mass,
      });
    });
  }
  return rows;
}
