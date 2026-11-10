/**
 * The report is a reading: every table comes from the envelope, and a refused document still
 * reports why rather than pretending to have results.
 */

import { readFileSync } from "node:fs";
import { resolve } from "node:path";
import { describe, expect, it } from "vitest";

import { buildCaseReport } from "../src/report/model";
import { unitsOf } from "../src/state/units";
import { parseCatalogue, parseEnvelope } from "../src/wire/decode";

const catalogue = parseCatalogue(readFileSync(resolve("test/fixtures/catalogue.json"), "utf8"));
const units = unitsOf(catalogue, "si");
const solved = parseEnvelope(readFileSync(resolve("test/fixtures/envelope.json"), "utf8"));
const refused = parseEnvelope(readFileSync(resolve("test/fixtures/broken.json"), "utf8"));

const table = (report: ReturnType<typeof buildCaseReport>, title: string) =>
  report.tables.find((entry) => entry.title === title);

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
    const results = table(report, "Results and caveats")?.rows ?? [];
    expect(results.some((row) => row[0] === "hx1" && row[1] === "outlet_duty")).toBe(true);
    const source = table(report, "Source and checks")?.rows ?? [];
    expect(source.some((row) => row[2] === "process.heater")).toBe(true);
    expect(report.diagnostics.rows).toEqual([]);
  });

  it("reports a refused document's diagnostics and no results", () => {
    const report = buildCaseReport(refused, catalogue, units, "now");
    expect(report.status).toBe("refused");
    expect(report.diagnostics.rows.length).toBeGreaterThan(0);
    expect(report.summary.errors).toBeGreaterThan(0);
    expect(table(report, "Results and caveats")?.rows ?? []).toEqual([]);
  });
});
