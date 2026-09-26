/**
 * Where a composition table's substance axis comes from.
 *
 * **A solved stream's `z` is a bare array**: the record is `n, z, P, T, h` and the names live on
 * the *feed* that declared the fluid, so the axis has to be recovered from the document - and the
 * case that matters is the one where the document does not agree with itself, where the honest
 * answer is a position rather than a name.
 */

import { describe, expect, it } from "vitest";

import { componentAxis, compositionRows, fraction, streamColumns, streamRows } from "../src/state/columns";
import { unitsOf } from "../src/state/units";
import catalogue from "./fixtures/catalogue.json";
import fixture from "./fixtures/envelope.json";
import type { Catalogue, Envelope, GraphNode } from "../src/wire/types";

const envelope = fixture as unknown as Envelope;
const field = unitsOf(catalogue as unknown as Catalogue, "field");

/** A document whose feeds declare the given substance lists. */
function withFeeds(lists: readonly (readonly string[])[]): Envelope {
  const nodes: GraphNode[] = lists.map((components, index) => ({
    id: `input:feed_${index}`,
    type: "stream",
    role: "input",
    position: { x: 0, y: 0 },
    data: {
      name: `feed_${index}`,
      input: { components: [...components], n: 1, z: components.map(() => 1 / components.length), P: 1e5, T: 300 },
    },
  }));
  return {
    ...envelope,
    flowsheet: { ...envelope.flowsheet, graph: { ...envelope.flowsheet.graph, nodes } },
  };
}

it("takes the substances from the feed that declares them", () => {
  const axis = componentAxis(envelope);
  expect(axis.positional).toBe(false);
  expect(axis.note).toBeNull();
  expect(axis.names.length).toBeGreaterThan(0);
  // The shipped document's fluid, from its own feed rather than from a fixture in this file.
  expect(axis.names).toEqual(["methane", "n-butane"]);
});

it("keeps the names where two feeds declare the same fluid", () => {
  const axis = componentAxis(withFeeds([["methane", "n-butane"], ["methane", "n-butane"]]));
  expect(axis.positional).toBe(false);
  expect(axis.names).toEqual(["methane", "n-butane"]);
});

it("marks the axis positional where the feeds disagree, rather than picking one", () => {
  // Same length, different names: a header drawn from either feed is a claim the document does not
  // make, which is the case a front end must not guess at.
  const axis = componentAxis(withFeeds([["methane", "n-butane"], ["propane", "ethane"]]));
  expect(axis.positional).toBe(true);
  expect(axis.note).toContain("disagree");
  expect(axis.names).toEqual(["methane", "n-butane"]);

  // And different lengths, where one feed's list cannot even be a renaming of the other's.
  const uneven = componentAxis(withFeeds([["methane"], ["methane", "n-butane", "propane"]]));
  expect(uneven.positional).toBe(true);
  expect(uneven.names).toEqual(["methane", "n-butane", "propane"]);
});

it("says so where no feed declares a fluid", () => {
  const axis = componentAxis(withFeeds([]));
  expect(axis.names).toEqual([]);
  expect(axis.note).toContain("no feed");
});

it("writes a fraction to four places", () => {
  expect(fraction(0.123456)).toBe("0.1235");
  expect(fraction(0)).toBe("0.0000");
});

describe("the workbook's rows", () => {
  it("comes from the document's ports, in the graph's own order", () => {
    const rows = streamRows(envelope);
    // **The rule that makes it a workbook and not a run log.** These are the ports the *document*
    // declares, so the grid has this shape before a run as well as after one.
    expect(rows.map((row) => row.path)).toEqual([
      // **One row per stream and not one per port**: a producing outlet and the inlet it feeds are
      // the same stream, and only the producer's endpoint carries the value.
      "mix1.product",
      "p1.outlet",
      "hx1.outlet",
      "sep1.vapour",
      "sep1.liquid",
      "feed_1",
      "vapour_product",
    ]);
    // Every row is a path, the object it belongs to, and the record the run reached.
    expect(rows[0]?.from).toBe("mix1");
    expect(rows.every((row) => row.record !== undefined)).toBe(true);
  });

  it("has the same rows with no run at all, and the values are what is missing", () => {
    const unrun: Envelope = { ...envelope, session: null, dirty: true };
    const before = streamRows(unrun);
    const after = streamRows(envelope);
    expect(before.map((row) => row.path)).toEqual(after.map((row) => row.path));
    expect(before.every((row) => row.record === undefined)).toBe(true);
    expect(after.some((row) => row.record !== undefined)).toBe(true);
  });

  it("draws a column for a quantity it has never heard of", () => {
    // **The forward-compatibility claim, tested rather than asserted.** A field added to the record
    // after this file was written appears with its own key as its label - which is the difference
    // between a table that shows what arrived and one that shows what it expected.
    const record = {
      ...(envelope.session?.streams["sep1.vapour"] as object),
      entropy: { magnitude_si: -26.9, unit: "J/(mol*K)" },
    } as never;
    const columns = streamColumns(record);
    expect(columns.map((column) => column.key)).toContain("entropy");
    expect(columns.find((column) => column.key === "entropy")?.unit).toBe("J/(mol*K)");
    // And a composition is never a column of its own: it has a sheet.
    expect(columns.map((column) => column.key)).not.toContain("z");
  });

  it("reads every stream's substances against the feed's axis", () => {
    const rows = compositionRows(envelope);
    const vapour = rows.filter((row) => row.path === "sep1.vapour");
    expect(vapour.map((row) => row.component)).toEqual(["methane", "n-butane"]);
    // `n · z` is two fields of the record multiplied, and the flow is the run's own.
    const flow = envelope.session?.streams["sep1.vapour"]?.n.magnitude_si ?? 0;
    expect(vapour[0]?.flow).toBeCloseTo(flow * (vapour[0]?.fraction ?? 0), 10);

    // **And the derived cell is read in the set's unit.** `n · z` is a molar flow, so the field
    // set is the one that changes it - and it changes by the factor of `kmol/h`, not of `mol/s`.
    const kmolPerHour = 1000 / 3600;
    const inField = compositionRows(envelope, field).filter((row) => row.path === "sep1.vapour");
    expect(inField[0]?.flow).toBeCloseTo((vapour[0]?.flow ?? 0) / kmolPerHour, 10);
  });

  it("renames a column to the unit it is read in, and leaves the rest alone", () => {
    // **A column's unit is the unit its cells are in**, which is what keeps the header and the
    // numbers from disagreeing. The four that change are the four the field set names; `M` and
    // `VF` are not among them, because no set names their dimensions.
    const inSi = streamColumns(envelope.session?.streams["sep1.vapour"], unitsOf(null, null));
    const inField = streamColumns(envelope.session?.streams["sep1.vapour"], field);
    expect(inSi.map((column) => `${column.label} ${String(column.unit)}`)).toEqual([
      "flow mol/s",
      "mass flow kg/s",
      "M kg/mol",
      "P Pa",
      "T K",
      "h J/mol",
      "VF null",
    ]);
    expect(inField.map((column) => `${column.label} ${String(column.unit)}`)).toEqual([
      "flow kmol/h",
      "mass flow lb/h",
      "M kg/mol",
      "P psi",
      "T K",
      "h kJ/mol",
      "VF null",
    ]);
  });
});
