/**
 * The case report, assembled from the envelope and nothing else.
 *
 * **A report is a reading, not a second run.** Every number here is a value the last Solve
 * reached, formatted in the unit set in force; the provenance is the library's own record; the
 * caveats and messages are its own sentences. Nothing is recomputed and nothing is rounded into
 * a claim.
 *
 * **The report is a list of tables.** A caller (the printable view, a future export) draws them
 * without knowing which one it has, so a field added to a result or a stream appears in the
 * report the same way it appears in the workbook.
 */

import {
  compositionRows,
  molarMasses,
  streamColumns,
  streamRows,
} from "../state/columns";
import { displayUnitOf, inDisplayUnit, type Units } from "../state/units";
import { scalarsOf, seriesOf, warningsOf } from "../state/results";
import { provenanceOf, shortHash } from "../state/provenance";
import { fieldTextIn } from "../wire/field";
import type { Catalogue, Envelope } from "../wire/types";
import { streamCellText } from "./csv";

/** One table the report prints. */
export interface ReportTable {
  title: string;
  headers: string[];
  rows: string[][];
}

/** The case, as the printable view and its tests read it. */
export interface CaseReport {
  title: string;
  id: string;
  generated: string;
  status: string;
  library: string;
  version: string;
  unitSet: string;
  tables: ReportTable[];
  diagnostics: ReportTable;
  summary: {
    unitOps: number;
    streams: number;
    warnings: number;
    errors: number;
  };
}

/** The status word, from the same facts the title strip reads. */
function statusWord(envelope: Envelope): string {
  if (!envelope.ok) {
    return "refused";
  }
  if (envelope.run_error !== null) {
    return "failed";
  }
  return envelope.dirty ? "stale" : "solved";
}

/**
 * Assemble the report.
 *
 * `generated` is passed in rather than read from the clock here, so a test can hold the report to
 * one document instead of to the day it ran.
 */
export function buildCaseReport(
  envelope: Envelope,
  catalogue: Catalogue | null,
  units: Units,
  generated: string,
): CaseReport {
  const nodes = envelope.flowsheet.graph.nodes;
  const unitOps = nodes.filter((node) => node.role === "instance");
  const feeds = nodes.filter((node) => node.role === "input");
  const set = units.sets.find((candidate) => candidate.id === units.active);

  const feedRows = feeds.flatMap((node) => {
    const record = node.data.input;
    if (record === undefined) {
      return [];
    }
    return [
      [
        node.data.name,
        record.components.join(", "),
        fieldTextIn(record.n, "mol/s", units),
        fieldTextIn(record.P, "Pa", units),
        fieldTextIn(record.T, "K", units),
      ],
    ];
  });

  const parameterRows = unitOps.flatMap((node) => {
    const form = catalogue?.unit_ops.find((entry) => entry.id === node.data.unit);
    if (form === undefined) {
      return [];
    }
    return form.parameters.map((parameter) => {
      const value = node.data.parameters?.[parameter.name];
      return [
        node.data.name,
        parameter.name,
        value === undefined ? "—" : fieldTextIn(value, parameter.unit, units),
        displayUnitOf(units, parameter.unit) ?? "",
      ];
    });
  });

  const rows = streamRows(envelope);
  const columns = streamColumns(rows[0]?.record, units);
  const streamRowsReport = rows.map((row) => [
    row.path,
    ...columns.map((column) => streamCellText(row.record, column.key, column.digits, units)),
  ]);

  const masses = molarMasses(catalogue);
  const composition = compositionRows(envelope, units, masses);
  const compositionTable = composition.map((row) => [
    row.component,
    row.path,
    row.fraction.toFixed(4),
    row.flow === null ? "" : String(Number(row.flow.toPrecision(5))),
    row.mass_fraction === null ? "—" : row.mass_fraction.toFixed(4),
    row.mass_flow === null
      ? "—"
      : String(Number(inDisplayUnit(units, row.mass_flow, "kg/s").toPrecision(5))),
  ]);

  const resultRows: string[][] = [];
  let warningCount = 0;
  for (const node of unitOps) {
    const result = envelope.session?.results[node.data.name];
    if (result === undefined) {
      continue;
    }
    for (const scalar of scalarsOf(result, units)) {
      resultRows.push([node.data.name, scalar.key, scalar.text, scalar.unit ?? ""]);
    }
    for (const series of seriesOf(result, units)) {
      resultRows.push([
        node.data.name,
        series.key,
        series.values.map((value) => String(Number(value.toPrecision(5)))).join(", "),
        series.unit ?? "",
      ]);
    }
    for (const warning of warningsOf(result)) {
      warningCount += 1;
      resultRows.push([node.data.name, warning.code, warning.message, warning.field ?? ""]);
    }
  }

  const provenanceRows = unitOps.map((node) => {
    const entry = provenanceOf(envelope, node.data.name);
    return [
      node.data.name,
      entry?.unit ?? node.data.unit ?? "—",
      entry?.model ?? "—",
      entry?.source ?? "—",
      entry?.verification ?? "—",
      shortHash(entry?.spec_sha256 ?? null) ?? "—",
      shortHash(entry?.rust_sha256 ?? null) ?? "—",
    ];
  });

  const tearRows = envelope.flowsheet.graph.edges
    .filter((edge) => edge.data.kind === "recycle")
    .map((edge) => {
      const tear = envelope.session?.tears.find((candidate) => candidate.stream === edge.data.path);
      return [
        edge.data.path,
        edge.data.from,
        String(tear?.iterations ?? 0),
        String(tear?.solved ?? false),
        String(tear?.active ?? false),
      ];
    });

  const errors = envelope.diagnostics.filter((diagnostic) => diagnostic.severity === "error");
  const diagnosticRows = envelope.diagnostics.map((diagnostic) => [
    diagnostic.severity,
    diagnostic.code,
    diagnostic.path === "" ? diagnostic.section : `${diagnostic.section}/${diagnostic.path}`,
    diagnostic.message,
  ]);

  return {
    title: envelope.flowsheet.name,
    id: envelope.flowsheet.id,
    generated,
    status: statusWord(envelope),
    library: envelope.provenance.library,
    version: envelope.provenance.version,
    unitSet: set?.name ?? units.active,
    tables: [
      {
        title: "Feeds",
        headers: [
          "feed",
          "components",
          `n [${displayUnitOf(units, "mol/s") ?? "mol/s"}]`,
          `P [${displayUnitOf(units, "Pa") ?? "Pa"}]`,
          `T [${displayUnitOf(units, "K") ?? "K"}]`,
        ],
        rows: feedRows,
      },
      {
        title: "Unit operations",
        headers: ["unit operation", "parameter", "value", "unit"],
        rows: parameterRows,
      },
      {
        title: "Streams",
        headers: [
          "stream",
          ...columns.map((column) =>
            column.unit === null ? column.label : `${column.label} [${column.unit}]`,
          ),
        ],
        rows: streamRowsReport,
      },
      {
        title: "Composition",
        headers: [
          "component",
          "stream",
          "mole fraction",
          `n·z [${displayUnitOf(units, "mol/s") ?? "mol/s"}]`,
          "mass fraction",
          `mass flow [${displayUnitOf(units, "kg/s") ?? "kg/s"}]`,
        ],
        rows: compositionTable,
      },
      {
        title: "Results and caveats",
        headers: ["unit operation", "field", "value", "unit or field"],
        rows: resultRows,
      },
      {
        title: "Recycles",
        headers: ["stream", "from", "passes", "solved", "active"],
        rows: tearRows,
      },
      {
        title: "Source and checks",
        headers: ["unit operation", "entry", "model", "source", "verification", "spec", "kernel"],
        rows: provenanceRows,
      },
    ],
    diagnostics: {
      title: "Messages",
      headers: ["severity", "code", "where", "message"],
      rows: diagnosticRows,
    },
    summary: {
      unitOps: unitOps.length,
      streams: Object.keys(envelope.session?.streams ?? {}).length,
      warnings: warningCount,
      errors: errors.length,
    },
  };
}
