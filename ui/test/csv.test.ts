/**
 * The workbook as a file: the quoting, the injection guard and the units in the header.
 */

import { readFileSync } from "node:fs";
import { resolve } from "node:path";
import { describe, expect, it } from "vitest";

import { compositionCsv, csvCell, streamsCsv, tableToCsv, tearsCsv } from "../src/report/csv";
import { unitsOf } from "../src/state/units";
import { parseCatalogue, parseEnvelope } from "../src/wire/decode";

const envelope = parseEnvelope(readFileSync(resolve("test/fixtures/envelope.json"), "utf8"));
const catalogue = parseCatalogue(readFileSync(resolve("test/fixtures/catalogue.json"), "utf8"));
const units = unitsOf(catalogue, "si");

describe("CSV", () => {
  it("quotes only what the format requires, and guards a formula", () => {
    expect(csvCell("methane")).toBe("methane");
    expect(csvCell("a,b")).toBe('"a,b"');
    expect(csvCell('say "hi"')).toBe('"say ""hi"""');
    expect(csvCell("line\nbreak")).toBe('"line\nbreak"');
    expect(csvCell("=1+1")).toBe("'=1+1");
    expect(csvCell("+cmd")).toBe("'+cmd");
    expect(csvCell("-2")).toBe("'-2");
    expect(csvCell("@here")).toBe("'@here");
    expect(tableToCsv(["a", "b"], [["1", "2"]])).toBe("a,b\r\n1,2\r\n");
  });

  it("exports the streams with their units in the header", () => {
    const text = streamsCsv(envelope, units);
    const [header, first] = text.split("\r\n");
    expect(header).toContain("flow [mol/s]");
    expect(header).toContain("T [K]");
    expect(header?.startsWith("stream,")).toBe(true);
    expect(first).toBeDefined();
    expect(text).toContain("feed_1");
    expect(text).toContain("hx1.outlet");
  });

  it("exports the composition one substance per stream, and the tears", () => {
    const composition = compositionCsv(envelope, catalogue, units);
    expect(composition.split("\r\n")[0]).toContain("mole fraction");
    expect(composition).toContain("methane");
    expect(composition).toContain("feed_1");

    const tears = tearsCsv(envelope);
    expect(tears).toContain("recycle_1");
    expect(tears.split("\r\n")[0]).toContain("passes");
  });
});
