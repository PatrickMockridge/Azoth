/**
 * The drawing tables, held to the palette in both directions.
 *
 * **A glyph is a claim about machinery**, and the two ways that claim can rot are both silent: a
 * palette entry added with no drawing is a node that quietly becomes a grey box, and a drawing left
 * behind by a renamed port is a handle that quietly moves to the wrong side of the vessel. Neither
 * is visible to a test that only counts nodes, so the tables are read against the catalogue.
 *
 * The catalogue is the *captured* one (`test/fixtures/catalogue.json`, `azoth forms --tools`), which
 * `fixtures.test.ts` holds to the live output - so this is the shipped palette and not a fixture
 * invented here.
 */

import { expect, it } from "vitest";

import {
  GLYPHS,
  PORT_SIDES,
  glyphFor,
  leafOf,
  senseFor,
  sideOf,
} from "../src/state/glyphs";
import fixture from "./fixtures/catalogue.json";

const catalogue = fixture as unknown as {
  unit_ops: { id: string; ports: { name: string; direction: "in" | "out" }[] }[];
};

it("draws every entry the palette carries", () => {
  for (const form of catalogue.unit_ops) {
    const leaf = leafOf(form.id);
    expect(GLYPHS[leaf], `${form.id} has no drawing`).toBeDefined();
    // And it is a drawing rather than the fallback, which is what an omission would look like.
    expect(glyphFor(form.id), `${form.id} falls back`).not.toBe("fallback");
  }
});

it("draws nothing the palette does not carry", () => {
  const carried = new Set(catalogue.unit_ops.map((form) => leafOf(form.id)));
  for (const leaf of Object.keys(GLYPHS)) {
    expect(carried.has(leaf), `\`${leaf}\` is drawn and is not a palette entry`).toBe(true);
  }
});

it("names only ports the entries declare", () => {
  // **The failure this catches**: `stripping_column`'s overhead is `overhead_gas` and not
  // `gas_out`, and a table that said otherwise would move a handle to a side the drawing does not
  // have - rendered as a stream leaving a vessel's middle, with no test noticing.
  const ports = new Map(
    catalogue.unit_ops.map((form) => [leafOf(form.id), new Set(form.ports.map((p) => p.name))]),
  );
  for (const [leaf, sides] of Object.entries(PORT_SIDES)) {
    const declared = ports.get(leaf);
    expect(declared, `\`${leaf}\` has port sides and is not a palette entry`).toBeDefined();
    for (const port of Object.keys(sides)) {
      expect(declared?.has(port), `\`${leaf}\` declares no port called \`${port}\``).toBe(true);
    }
  }
});

it("places every declared port, on the side its direction implies where nothing overrides it", () => {
  for (const form of catalogue.unit_ops) {
    for (const port of form.ports) {
      const side = sideOf(form.id, port.name, port.direction);
      expect(["top", "bottom", "left", "right"]).toContain(side);
      // An inlet is never on the right and an outlet never on the left *unless* the table says so:
      // which is what makes the override table the only place a vertical port can come from.
      const overridden = PORT_SIDES[leafOf(form.id)]?.[port.name] !== undefined;
      if (!overridden) {
        expect(side).toBe(port.direction === "in" ? "left" : "right");
      }
    }
  }
});

it("gives an unknown entry the fallback and no sense", () => {
  expect(glyphFor("unit_ops.nonesuch")).toBe("fallback");
  expect(glyphFor(undefined)).toBe("fallback");
  expect(senseFor("unit_ops.nonesuch")).toBe("none");
  // The two machines differ in the direction their moving part points, which is how a compressor
  // and an expander are told apart at a glance.
  expect(senseFor("unit_ops.compressor")).toBe("in");
  expect(senseFor("unit_ops.expander")).toBe("out");
});
