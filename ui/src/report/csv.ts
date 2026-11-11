/**
 * The workbook, as files a spreadsheet can open.
 *
 * **RFC 4180 and one guard the format does not have.** A cell that begins with `=`, `+`, `-` or
 * `@` is a formula to a spreadsheet, and a substance name or a message is data; the guard prefixes
 * those cells with an apostrophe so an exported case cannot execute what a document contains.
 *
 * **The units are in the header, not repeated in every cell.** A column is one dimension, and a
 * reader sorting the sheet wants numbers; the header is where the unit is stated once.
 */

import { compositionRows, molarMasses, streamColumns, streamRows } from "../state/columns";
import { displayOf, inDisplayUnit, type Units } from "../state/units";
import { fieldOf, isQuantity } from "../wire/decode";
import type { Catalogue, Envelope } from "../wire/types";

/** One cell, quoted only where the format requires it, and guarded against formula injection. */
export function csvCell(value: string): string {
  const guarded = /^[=+\-@]/.test(value) ? `'${value}` : value;
  if (/[",\r\n]/.test(guarded)) {
    return `"${guarded.replaceAll('"', '""')}"`;
  }
  return guarded;
}

/** A table as RFC 4180 text, with a trailing newline. */
export function tableToCsv(headers: string[], rows: string[][]): string {
  return [headers, ...rows]
    .map((row) => row.map(csvCell).join(","))
    .join("\r\n")
    .concat("\r\n");
}

/** A cell of a stream record in the unit set, without its unit repeated. */
export function streamCellText(
  record: unknown,
  key: string,
  digits: number,
  units: Units,
): string {
  const value = fieldOf(record, key);
  if (value === undefined || value === null) {
    return "";
  }
  if (isQuantity(value)) {
    const shown = inDisplayUnit(units, value.magnitude_si, value.unit);
    return String(Number(shown.toPrecision(digits)));
  }
  if (typeof value === "number") {
    return String(Number(value.toPrecision(digits)));
  }
  return String(value);
}

/** The streams sheet: every producing endpoint, in the workbook's order and columns. */
export function streamsCsv(envelope: Envelope, units: Units): string {
  const rows = streamRows(envelope);
  const columns = streamColumns(rows[0]?.record, units);
  return tableToCsv(
    [
      "stream",
      ...columns.map((column) =>
        column.unit === null ? column.label : `${column.label} [${column.unit}]`,
      ),
    ],
    rows.map((row) => [
      row.path,
      ...columns.map((column) => streamCellText(row.record, column.key, column.digits, units)),
    ]),
  );
}

/** The composition sheet: one row per stream per substance. */
export function compositionCsv(
  envelope: Envelope,
  catalogue: Catalogue | null,
  units: Units,
): string {
  const flow = displayOf(units, "mol/s").unit;
  const mass = displayOf(units, "kg/s");
  const rows = compositionRows(envelope, units, molarMasses(catalogue));
  return tableToCsv(
    ["component", "stream", "mole fraction", `n·z [${flow}]`, "mass fraction", `mass flow [${mass.unit}]`],
    rows.map((row) => [
      row.component,
      row.path,
      row.fraction.toFixed(4),
      row.flow === null ? "" : String(Number(row.flow.toPrecision(5))),
      row.mass_fraction === null ? "" : row.mass_fraction.toFixed(4),
      row.mass_flow === null ? "" : String(Number((row.mass_flow / mass.factor).toPrecision(5))),
    ]),
  );
}

/** The tears sheet: the loops and the residuals the run last measured. */
export function tearsCsv(envelope: Envelope): string {
  const rows = envelope.flowsheet.graph.edges
    .filter((edge) => edge.data.kind === "recycle")
    .map((edge) => {
      const tear = envelope.session?.tears.find((candidate) => candidate.stream === edge.data.path);
      const at = (value: number | undefined) => (value === undefined ? "" : value.toExponential(3));
      return [
        edge.data.path,
        edge.data.from,
        String(tear?.iterations ?? 0),
        String(tear?.solved ?? false),
        String(tear?.active ?? false),
        at(tear?.residuals?.flow),
        at(tear?.residuals?.composition),
        at(tear?.residuals?.temperature),
        at(tear?.residuals?.pressure),
      ];
    });
  return tableToCsv(
    ["stream", "from", "passes", "solved", "active", "flow", "composition", "temperature", "pressure"],
    rows,
  );
}
