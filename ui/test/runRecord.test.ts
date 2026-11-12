/**
 * The run record: the machine-readable artefact an agent reproduces a run from.
 *
 * **It is the envelope, not a summary of it.** The test holds the record to the fixture's own
 * document, results, provenance and diagnostics, and to the absence of any display formatting:
 * a value is the library's magnitude and unit, not the converted, rounded number the cover
 * shows. A record that carried the cover's numbers would be a second answer.
 */

import { readFileSync } from "node:fs";
import { resolve } from "node:path";
import { describe, expect, it } from "vitest";

import { buildRunRecord } from "../src/report/runRecord";
import { unitsOf } from "../src/state/units";
import { parseCatalogue, parseEnvelope } from "../src/wire/decode";

const catalogue = parseCatalogue(readFileSync(resolve("test/fixtures/catalogue.json"), "utf8"));
const envelope = parseEnvelope(readFileSync(resolve("test/fixtures/envelope.json"), "utf8"));
const units = unitsOf(catalogue, "si");

describe("the run record", () => {
  it("carries the document, the raw session and the provenance", () => {
    const record = buildRunRecord(envelope, units, "now");
    expect(record.format).toBe("azoth.run-record/1");
    expect(record.case.id).toBe("flowsheets.demo");
    // The document is present so the run can be repeated rather than only read.
    expect(record.case.document).toContain("[[instances]]");
    // **The raw results**, exactly as the library wrote them: a quantity keeps its magnitude and
    // unit, and the absent fields stay absent rather than becoming display dashes.
    expect(record.session?.results.hx1?.outlet_duty).toEqual({
      magnitude_si: -4668.8484576717265,
      unit: "W",
    });
    expect(record.session?.results.hx1?.vfd_efficiency).toBeNull();
    expect(record.status).toBe("solved");
    expect(record.ok).toBe(true);
    expect(record.dirty).toBe(false);
    expect(record.execution_order).toBe("insertion");
    expect(record.library.name).toBe("azoth");
    expect(record.library.version).not.toBe("");
    expect(record.keycard).toContain("standard keycard");
    expect(record.units.active).toBe("si");
    expect(record.provenance.units.map((entry) => entry.instance)).toEqual([
      "mix1",
      "p1",
      "hx1",
      "sep1",
    ]);
    expect(record.diagnostics).toEqual([]);
    expect(record.run_error).toBeNull();
  });

  it("carries no display formatting", () => {
    const record = buildRunRecord(envelope, units, "now");
    // The library's outlet_n is a bare 1.0; the cover renders it as "1" with a unit from the
    // declaration, and the record does not carry that label or that unit.
    expect(record.session?.results.hx1?.outlet_n).toBe(1);
    const text = JSON.stringify(record);
    expect(text).not.toContain("Outlet duty");
    expect(text).not.toContain("how it is calculated");
  });
});
