import { useState } from "react";

import { componentAxis, compositionRows, fraction, streamColumns, streamRows } from "../../state/columns";
import { displayOf, type Units } from "../../state/units";
import { formatQuantity } from "../../wire/field";
import type { Envelope, Quantity, StreamRecord } from "../../wire/types";
import { DockTabs } from "../property/DockTabs";
import { WorkbookGrid } from "./WorkbookGrid";

/** The workbook's sheets. */
type Sheet = "streams" | "composition" | "tears";

/** What each sheet is called on its tab. */
const LABELS: Record<Sheet, string> = {
  streams: "Streams",
  composition: "Composition",
  tears: "Tears",
};

/**
 * Every stream the document wires, down the page, with what the run reached across it.
 *
 * **This is the HYSYS workbook, and the one thing that makes it a workbook is where its rows come
 * from.** They are the *document's* port paths - see `streamRows` - so the grid has the same shape
 * before and after Solve, and a stream the run did not reach is a row of dashes rather than a row
 * that is not there. A run log would be the same table with the rows it happened to compute.
 *
 * **Nothing here computes a number.** A cell is the record's own value formatted, or a dash. The
 * one derived column in the composition sheet is `n · z`, two fields of one record multiplied, and
 * it says so.
 */
export function Workbook({
  envelope,
  units,
  onSelect,
}: {
  envelope: Envelope;
  /** The unit a reader wants the values in; a catalogue that has not loaded converts nothing. */
  units: Units;
  onSelect: (id: string | null) => void;
}) {
  const [sheet, setSheet] = useState<Sheet>("streams");
  const rows = streamRows(envelope);
  const axis = componentAxis(envelope);
  const stale = envelope.dirty;

  return (
    <div className="workbook">
      <DockTabs
        label="workbook"
        active={sheet}
        tabs={[
          { id: "streams", label: LABELS.streams, count: rows.length },
          { id: "composition", label: LABELS.composition },
          { id: "tears", label: LABELS.tears, count: envelope.session?.tears.length ?? 0 },
        ]}
        onTab={(id) => setSheet(id as Sheet)}
      />
      {sheet === "streams" ? (
        <WorkbookGrid
          columns={streamColumns(rows[0]?.record, units)}
          stale={stale}
          empty="this document wires no streams"
          onSelect={(id) => {
            const row = rows.find((candidate) => candidate.path === id);
            onSelect(row?.nodeId === undefined || row.nodeId === "" ? null : row.nodeId);
          }}
          rows={rows.map((row) => ({
            id: row.path,
            head: row.path,
            cells: streamColumns(row.record, units).map((column) =>
              cell(row.record, column.key, units),
            ),
          }))}
        />
      ) : null}
      {sheet === "composition" ? (
        <>
          {axis.note === null ? null : <p className="note">{axis.note}</p>}
          <WorkbookGrid
            columns={[
              { key: "component", label: "component" },
              { key: "stream", label: "stream" },
              { key: "fraction", label: "mole fraction", numeric: true },
              // **The label names the unit the cells are in**, which is the set's: the sheet
              // converts `n · z` with the same factor the workbook converts `n` with.
              { key: "flow", label: `n · z ${flowUnit(units)}`, numeric: true },
            ]}
            stale={stale}
            empty="no stream in this document carries a composition yet"
            onSelect={(id) => {
              const row = compositionRows(envelope, units).find(
                (candidate) => `${candidate.path}:${candidate.component}` === id,
              );
              onSelect(row?.nodeId === undefined || row.nodeId === "" ? null : row.nodeId);
            }}
            rows={compositionRows(envelope, units).map((row) => ({
              id: `${row.path}:${row.component}`,
              head: row.component,
              cells: [
                row.path,
                fraction(row.fraction),
                row.flow === null ? "—" : fraction(row.flow),
              ],
            }))}
          />
        </>
      ) : null}
      {sheet === "tears" ? (
        <WorkbookGrid
          columns={[
            { key: "stream", label: "tear" },
            { key: "from", label: "from" },
            { key: "iterations", label: "passes", numeric: true },
            { key: "solved", label: "solved" },
            { key: "active", label: "active" },
            { key: "flow", label: "flow", numeric: true },
            { key: "composition", label: "z", numeric: true },
            { key: "temperature", label: "T", numeric: true },
            { key: "pressure", label: "P", numeric: true },
          ]}
          stale={stale}
          empty="this document declares no tear"
          rows={envelope.flowsheet.graph.edges
            .filter((edge) => edge.data.kind === "recycle")
            .map((edge) => {
              const tear = envelope.session?.tears.find(
                (candidate) => candidate.stream === edge.data.path,
              );
              const at = (value: number | undefined) =>
                value === undefined ? "—" : value.toExponential(2);
              return {
                id: edge.data.path,
                head: edge.data.path,
                cells: [
                  edge.data.from,
                  tear?.iterations ?? 0,
                  String(tear?.solved ?? false),
                  String(tear?.active ?? false),
                  at(tear?.residuals?.flow),
                  at(tear?.residuals?.composition),
                  at(tear?.residuals?.temperature),
                  at(tear?.residuals?.pressure),
                ],
              };
            })}
        />
      ) : null}
    </div>
  );
}

/** The unit the composition sheet's derived column is read in. */
function flowUnit(units: Units): string {
  return displayOf(units, "mol/s").unit;
}

/** One cell: the record's value, or a dash where the run left none. */
function cell(record: StreamRecord | undefined, key: string, units: Units): string {
  const value = (record as Record<string, unknown> | undefined)?.[key];
  if (value === undefined || value === null) {
    return "—";
  }
  if (typeof value === "number") {
    return value.toFixed(4);
  }
  const quantity = value as Quantity;
  if (typeof quantity.magnitude_si === "number") {
    return formatQuantity(quantity.magnitude_si, quantity.unit, 4, units);
  }
  return String(value);
}
