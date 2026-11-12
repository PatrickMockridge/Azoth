/**
 * The case report, assembled from the envelope and nothing else.
 *
 * **A report is a reading, not a second run.** Every number here is a value the last Solve
 * reached, formatted in the unit set in force; the provenance is the library's own record; the
 * caveats and messages are its own sentences. Nothing is recomputed and nothing is rounded into
 * a claim.
 *
 * **Two layers, and the reader chooses.** The `cover` is for a person: the drawing, the streams
 * and a results table in plain labels, display units and the library's own sentence about how
 * each number is obtained - often the equation behind it. The `appendix` is the machine record:
 * the model's field names, the library's units, the design parameters that were set, the
 * profiles, and the provenance (model, raw verification status, hashes, data authority). A
 * reader can check the cover by hand; an agent can reproduce and compare the appendix.
 */

import {
  compositionRows,
  molarMasses,
  streamColumns,
  streamRows,
} from "../state/columns";
import { displayUnitOf, inDisplayUnit, type Units } from "../state/units";
import { formatNumber, outputRowsOf, seriesOf, warningsOf } from "../state/results";
import { provenanceOf, shortHash } from "../state/provenance";
import { fieldTextIn } from "../wire/field";
import type { Catalogue, Envelope } from "../wire/types";
import { streamCellText } from "./csv";
import { statusWord } from "../state/status";

/** One table the report prints. */
export interface ReportTable {
  title: string;
  headers: string[];
  rows: string[][];
  /** A sentence under the table, where the table needs one to be read correctly. */
  note?: string;
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
  /** What a person reads: results in plain labels, streams and caveats. */
  cover: ReportTable[];
  /** The machine record: field names, library units, parameters, profiles and provenance. */
  appendix: ReportTable[];
  diagnostics: ReportTable;
  summary: {
    unitOps: number;
    streams: number;
    warnings: number;
    errors: number;
  };
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

  const formOf = (unit: string | undefined) =>
    unit === undefined ? undefined : catalogue?.unit_ops.find((entry) => entry.id === unit);

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

  // **Only the parameters that were set.** The editor's Design sheet draws every declared field
  // so one can be filled in; a report is a record of what the run was given, and an unset field
  // is not part of that record. A dash would look like a value nobody supplied.
  const parameterRows = unitOps.flatMap((node) => {
    const form = formOf(node.data.unit);
    if (form === undefined) {
      return [];
    }
    return form.parameters.flatMap((parameter) => {
      const value = node.data.parameters?.[parameter.name];
      if (value === undefined || value === null || value === "") {
        return [];
      }
      return [
        [
          node.data.name,
          parameter.name,
          fieldTextIn(value, parameter.unit, units),
          displayUnitOf(units, parameter.unit) ?? "",
        ],
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

  const humanResultRows: string[][] = [];
  const machineResultRows: string[][] = [];
  const profileRows: string[][] = [];
  const caveatRows: string[][] = [];
  let warningCount = 0;
  for (const node of unitOps) {
    const result = envelope.session?.results[node.data.name];
    if (result === undefined) {
      continue;
    }
    const outputs = formOf(node.data.unit)?.outputs ?? [];
    for (const row of outputRowsOf(result, units, outputs)) {
      humanResultRows.push([
        node.data.name,
        row.label,
        row.text,
        row.unit ?? "",
        row.description ?? "",
      ]);
    }
    for (const row of outputRowsOf(result, undefined, outputs)) {
      machineResultRows.push([
        node.data.name,
        row.key,
        row.text,
        row.unit ?? "",
        row.description ?? "",
      ]);
    }
    for (const series of seriesOf(result)) {
      profileRows.push([
        node.data.name,
        series.key,
        series.values.map((value) => formatNumber(value)).join(", "),
        series.unit ?? "",
      ]);
    }
    for (const warning of warningsOf(result)) {
      warningCount += 1;
      caveatRows.push([node.data.name, warning.code, warning.field ?? "", warning.message]);
    }
  }

  const provenanceRows = unitOps.map((node) => {
    const entry = provenanceOf(envelope, node.data.name);
    return [
      node.data.name,
      entry?.model ?? "—",
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

  const cover: ReportTable[] = [
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
  ];
  if (humanResultRows.length > 0) {
    cover.push({
      title: "Results",
      headers: ["unit operation", "output", "value", "unit", "how it is calculated"],
      rows: humanResultRows,
    });
  }
  cover.push(
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
  );
  if (caveatRows.length > 0) {
    cover.push({
      title: "Caveats",
      headers: ["unit operation", "code", "field", "message"],
      rows: caveatRows,
    });
  }

  const appendix: ReportTable[] = [
    {
      title: "Design parameters",
      headers: ["unit operation", "parameter", "value", "unit"],
      rows: parameterRows,
      note: "Only parameters the run was given are listed; an unset field uses the model's declared default.",
    },
  ];
  if (machineResultRows.length > 0) {
    appendix.push({
      title: "Result fields",
      headers: ["unit operation", "field", "value", "unit", "description"],
      rows: machineResultRows,
      note: "The model's own field names, in the library's units, for comparison and reproduction.",
    });
  }
  if (profileRows.length > 0) {
    appendix.push({
      title: "Profiles",
      headers: ["unit operation", "field", "values", "unit"],
      rows: profileRows,
    });
  }
  if (tearRows.length > 0) {
    appendix.push({
      title: "Recycles",
      headers: ["stream", "from", "passes", "solved", "active"],
      rows: tearRows,
    });
  }
  appendix.push({
    title: "Model checks",
    headers: ["unit operation", "model", "verification", "spec", "kernel"],
    rows: provenanceRows,
    note: "Verification describes how far the model itself is checked; it is not the status of the run.",
  });
  appendix.push({
    title: "Data authority",
    headers: ["keycard"],
    rows: [[envelope.provenance.keycard]],
  });

  return {
    title: envelope.flowsheet.name,
    id: envelope.flowsheet.id,
    generated,
    status: statusWord(envelope),
    library: envelope.provenance.library,
    version: envelope.provenance.version,
    unitSet: set?.name ?? units.active,
    cover,
    appendix,
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
