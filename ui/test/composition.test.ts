/**
 * The feed draft: fractions, normalisation and the one conversion a feed needs.
 *
 * The conversion is not a table in the front end — it is the catalogue's molar masses, the same
 * ones the workbook's mass column reads. These tests hold the arithmetic to that source rather
 * than to a constant written here.
 */

import { readFileSync } from "node:fs";
import { resolve } from "node:path";
import { describe, expect, it } from "vitest";

import {
  draftFromInput,
  emptyDraft,
  moleFractions,
  normalise,
  readyToSend,
  recordOf,
  sumOf,
  type FeedDraft,
} from "../src/state/composition";
import { unitsOf } from "../src/state/units";
import { parseCatalogue } from "../src/wire/decode";

const catalogue = parseCatalogue(readFileSync(resolve("test/fixtures/catalogue.json"), "utf8"));
const si = unitsOf(catalogue, "si");
const field = unitsOf(catalogue, "field");

function draft(components: [string, string][]): FeedDraft {
  return {
    ...emptyDraft("feed_1", si),
    components: components.map(([name, fraction]) => ({ name, fraction })),
  };
}

describe("a feed draft", () => {
  it("sums and normalises only when asked", () => {
    const uneven = draft([
      ["methane", "1"],
      ["ethane", "1"],
    ]);
    expect(sumOf(uneven)).toBe(2);
    expect(readyToSend(uneven)).toBe(false);

    const normal = normalise(uneven);
    expect(normal).not.toBeNull();
    expect(sumOf(normal as FeedDraft)).toBeCloseTo(1, 12);
    expect(readyToSend(normal as FeedDraft)).toBe(true);
  });

  it("converts mass fractions with the catalogue's molar masses", () => {
    const masses = draft([
      ["methane", "0.5"],
      ["ethane", "0.5"],
    ]);
    const weight = (name: string) =>
      catalogue.components?.find((entry) => entry.name === name)?.molar_mass ?? 0;
    expect(weight("methane")).toBeGreaterThan(0);
    expect(weight("ethane")).toBeGreaterThan(0);
    const moles = moleFractions({ ...masses, mode: "mass" }, catalogue);
    expect(moles).not.toBeNull();
    const methane = (0.5 / weight("methane")) / (0.5 / weight("methane") + 0.5 / weight("ethane"));
    expect((moles as number[])[0]).toBeCloseTo(methane, 12);
    expect((moles as number[]).reduce((sum, value) => sum + value, 0)).toBeCloseTo(1, 12);
  });

  it("shows a record in the unit set and stores the document's unit", () => {
    const record = {
      components: ["methane"],
      n: 1,
      z: [1],
      P: 500000,
      T: 366.4833333333333,
    };
    const shown = draftFromInput("feed_1", record, field);
    expect(Number(shown.T)).toBeCloseTo(200, 3);
    const prepared = recordOf(shown, catalogue, field);
    expect(prepared).not.toBeNull();
    expect(prepared?.T).toBeCloseTo(366.4833, 3);
    // The field shows six significant figures, so the stored value is within a fifth of a pascal:
    // a reader's display precision, not a conversion error.
    expect(prepared?.P).toBeCloseTo(500000, 0);
    expect(prepared?.z).toEqual([1]);
  });

  it("opens a new feed at the physical defaults in the set in force", () => {
    const fresh = {
      ...emptyDraft("feed_1", field),
      components: [{ name: "methane", fraction: "1" }],
    };
    // 5e5 Pa is 72.52 psi and 300 K is 80.33 °F; the boxes say so.
    expect(Number(fresh.P)).toBeCloseTo(72.52, 1);
    expect(Number(fresh.T)).toBeCloseTo(80.33, 1);
    const prepared = recordOf(fresh, catalogue, field);
    expect(prepared).not.toBeNull();
    expect(prepared?.P).toBeCloseTo(500000, 0);
    expect(prepared?.T).toBeCloseTo(300, 2);
  });

  it("returns the draft in SI when the unit set is SI", () => {
    const prepared = recordOf(draft([["methane", "1"]]), catalogue, si);
    expect(prepared).not.toBeNull();
    expect(prepared?.n).toBe(1);
    expect(prepared?.P).toBe(500000);
  });
});
