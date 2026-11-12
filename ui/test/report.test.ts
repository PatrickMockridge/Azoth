/**
 * The report is a reading: every table comes from the envelope, and a refused document still
 * reports why rather than pretending to have results.
 *
 * **Two layers are the point.** The cover is for a person - plain labels, display units and the
 * library's sentence - and the appendix keeps the model's field names and the provenance. The
 * tests below hold the line between them: a null is not a row anywhere, and the raw spelling of a
 * status belongs in the appendix rather than beside a number.
 */

import { readFileSync } from "node:fs";
import { resolve } from "node:path";
import { describe, expect, it } from "vitest";

import { buildCaseReport, type CaseReport } from "../src/report/model";
import { unitsOf } from "../src/state/units";
import { parseCatalogue, parseEnvelope } from "../src/wire/decode";

const catalogue = parseCatalogue(readFileSync(resolve("test/fixtures/catalogue.json"), "utf8"));
const units = unitsOf(catalogue, "si");
const solved = parseEnvelope(readFileSync(resolve("test/fixtures/envelope.json"), "utf8"));
const refused = parseEnvelope(readFileSync(resolve("test/fixtures/broken.json"), "utf8"));

const tables = (report: CaseReport) => [...report.cover, ...report.appendix];
const table = (report: CaseReport, title: string) =>
  tables(report).find((entry) => entry.title === title);

describe("the case report", () => {
  it("names the case, its status, its units and its library", () => {
    const report = buildCaseReport(solved, catalogue, units, "now");
    expect(report.title).toBe("Feed, pump, heat, separate, recycle");
    expect(report.id).toBe("flowsheets.demo");
    expect(report.status).toBe("solved");
    expect(report.unitSet).toBe("SI");
    expect(report.version.length).toBeGreaterThan(0);
    expect(report.summary.unitOps).toBe(4);
    expect(report.summary.streams).toBeGreaterThan(0);
  });

  it("carries the feeds, the streams, the results and the source", () => {
    const report = buildCaseReport(solved, catalogue, units, "now");
    expect(table(report, "Feeds")?.rows[0]?.[0]).toBe("feed_1");
    expect(table(report, "Streams")?.rows.map((row) => row[0])).toContain("hx1.outlet");

    // **The cover's results are plain labels and no empty row.** hx1 publishes eight fields;
    // the 64 electrical and mechanical fields it never configured are absent, not dashes.
    const results = table(report, "Results")?.rows ?? [];
    expect(results.length).toBeGreaterThan(0);
    expect(results.every((row) => row[1] !== "—" && row[2] !== "—")).toBe(true);
    expect(results.some((row) => row[0] === "hx1" && row[1] === "Outlet duty")).toBe(true);
    expect(results.some((row) => row[1] === "Max utilization")).toBe(true);
    // The library's sentence travels with the value.
    expect(results.find((row) => row[1] === "Outlet duty")?.[4]).toContain("inlet_n");

    // The appendix is the record: the model's own names, in the library's units.
    const fields = table(report, "Result fields")?.rows ?? [];
    expect(fields.some((row) => row[0] === "hx1" && row[1] === "outlet_duty")).toBe(true);
    // Only the parameter the demo set.
    const parameters = table(report, "Design parameters")?.rows ?? [];
    expect(parameters.filter((row) => row[0] === "hx1").map((row) => row[1])).toEqual([
      "outlet_temperature",
    ]);
    const checks = table(report, "Model checks")?.rows ?? [];
    expect(checks.some((row) => row[1] === "process.heater")).toBe(true);
    // The raw status is the machine record's, and it stays spelled the library's way.
    expect(checks.some((row) => row[2] === "partially_verified")).toBe(true);
    expect(table(report, "Data authority")?.rows[0]?.[0]).toContain("standard keycard");
    expect(report.diagnostics.rows).toEqual([]);
  });

  it("splits caveats out of the results, and omits tables with nothing in them", () => {
    const report = buildCaseReport(solved, catalogue, units, "now");
    // The demo's heater raised no warning, so there is no caveats table at all.
    expect(table(report, "Caveats")).toBeUndefined();
    // A table with no rows is not printed by the view, so the model does not invent an empty one.
    expect(table(report, "Results")?.rows.length).toBeGreaterThan(0);
  });

  it("reports a refused document's diagnostics and no results", () => {
    const report = buildCaseReport(refused, catalogue, units, "now");
    expect(report.status).toBe("incomplete");
    expect(report.diagnostics.rows.length).toBeGreaterThan(0);
    expect(report.summary.errors).toBeGreaterThan(0);
    expect(table(report, "Results")).toBeUndefined();
    expect(table(report, "Result fields")).toBeUndefined();
  });
});
