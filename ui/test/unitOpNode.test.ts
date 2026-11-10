/**
 * The handles on a node, held to the catalogue and to the side they share.
 *
 * **The defect this exists for was invisible to every existing test.** Two single-handle ports on
 * the same side both sat at 50%, so the heat exchanger's two inlets and two outlets were drawn as
 * one of each; the catalogue and the executor were correct. A test that only counted nodes could
 * not see it.
 */

import { readFileSync } from "node:fs";
import { resolve } from "node:path";
import { describe, expect, it } from "vitest";

import { portSlots } from "../src/components/UnitOpNode";
import { parseCatalogue } from "../src/wire/decode";
import type { Form, GraphNode, Handle } from "../src/wire/types";

const catalogue = parseCatalogue(readFileSync(resolve("test/fixtures/catalogue.json"), "utf8"));

/** A graph node for a catalogue form, with `handlesPerPort` handles on every port. */
function nodeFor(form: Form, handlesPerPort = 1): GraphNode {
  const inlets: Handle[] = [];
  const outlets: Handle[] = [];
  for (const port of form.ports) {
    const handles = Array.from({ length: handlesPerPort }, (_, index) =>
      handlesPerPort === 1 ? port.name : `${port.name}[${String(index)}]`,
    );
    (port.direction === "in" ? inlets : outlets).push({
      name: port.name,
      multiplicity: port.multiplicity,
      handles,
    });
  }
  return {
    id: `instance:${form.id}`,
    type: "unit_op",
    role: "instance",
    position: { x: 0, y: 0 },
    data: {
      name: form.id,
      unit: form.id,
      unit_name: form.name,
      parameters: {},
      ports: { inlets, outlets },
    },
  };
}

describe("the port slots", () => {
  it("gives every handle on a side its own place", () => {
    for (const form of catalogue.unit_ops) {
      const slots = portSlots(form, nodeFor(form));
      const places = new Set(slots.map((slot) => `${slot.side}:${slot.at}`));
      expect(places.size, `${form.id} draws two handles on top of each other`).toBe(slots.length);
    }
  });

  it("draws the heat exchanger's four connections", () => {
    const form = catalogue.unit_ops.find((entry) => entry.id === "unit_ops.heat_exchanger");
    expect(form).toBeDefined();
    const slots = portSlots(form, nodeFor(form as Form));
    expect(slots).toHaveLength(4);
    expect(slots.filter((slot) => slot.side === "left").map((slot) => slot.at)).toEqual([1 / 3, 2 / 3]);
    expect(slots.filter((slot) => slot.side === "right").map((slot) => slot.at)).toEqual([1 / 3, 2 / 3]);
    expect(new Set(slots.map((slot) => slot.port))).toEqual(
      new Set(["hot_in", "hot_out", "cold_in", "cold_out"]),
    );
  });

  it("spreads a many outlet's streams along its side", () => {
    const form = catalogue.unit_ops.find((entry) => entry.id === "unit_ops.splitter");
    expect(form).toBeDefined();
    const slots = portSlots(form, nodeFor(form as Form, 3));
    expect(slots.filter((slot) => slot.side === "right").map((slot) => slot.at)).toEqual([
      0.25, 0.5, 0.75,
    ]);
  });
});
