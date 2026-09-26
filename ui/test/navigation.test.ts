/**
 * The navigator's rows.
 *
 * **What this file holds is the one thing the module decides** rather than reads: the order of the
 * families, which is the catalogue's so that this pane and the palette agree, and the two row kinds
 * that are not "a thing in the middle of the flowsheet" — the streams and the boundary.
 *
 * The demo's own capture is the fixture, so the rows below are the objects the shipped document
 * has; the rest of the cases build the smallest envelope that can carry the case being made.
 */

import { describe, expect, it } from "vitest";

import { familyLabel, familyOfForm, navigationOf } from "../src/state/navigation";
import catalogue from "./fixtures/catalogue.json";
import fixture from "./fixtures/envelope.json";
import type { Catalogue, Envelope, GraphEdge, GraphNode } from "../src/wire/types";

const demo = fixture as unknown as Envelope;
const forms = catalogue as unknown as Catalogue;

/** The shipped catalogue declares families this document does not use, which is the point of one
 * of the cases below; nothing here needs the rest of what a catalogue carries. */
const withFamilies = (families: (string | null)[]): Catalogue =>
  ({
    unit_ops: families.map((family, index) => ({
      id: `unit_ops.f${index}`,
      name: `f${index}`,
      source: null,
      model: null,
      family,
      runnable: true,
      refusal: null,
      ports: [],
      parameters: [],
      unmodelled_ranges: [],
    })),
  }) as Catalogue;

const node = (role: GraphNode["role"], name: string, unit?: string): GraphNode => ({
  id: role === "instance" ? `instance:${name}` : `${role}:${name}`,
  type: role === "instance" ? "unit_op" : "stream",
  role,
  position: { x: 0, y: 0 },
  data: unit === undefined ? { name } : { name, unit },
});

const envelopeOf = (nodes: GraphNode[], edges: GraphEdge[] = []): Envelope =>
  ({
    ok: true,
    dirty: false,
    execution_order: "insertion",
    flowsheet: {
      id: "f",
      name: "f",
      document: "",
      graph: { id: "f", name: "f", nodes, edges },
    },
    diagnostics: [],
    paths: [],
    session: null,
    run_error: null,
  }) as Envelope;

/** The labels of the groups, and of the rows under each, as a reader sees them. */
const shape = (groups: ReturnType<typeof navigationOf>) =>
  groups.map((group) => [group.label, group.rows.map((row) => row.label)] as const);

describe("the navigator", () => {
  it("lists the document's objects, families first and the boundary last", () => {
    // **The families come in the catalogue's order**, not the graph's: the document places `mix1`
    // before `p1` and `p1` before `hx1`, so first-appearance order would put `two port` third.
    // One order for families across both panes is the reason it does not.
    expect(shape(navigationOf(demo, forms))).toEqual([
      ["mixer", ["mix1"]],
      ["separator", ["sep1"]],
      ["two port", ["p1", "hx1"]],
      [
        "Streams",
        ["feed_1", "mix1.product", "p1.outlet", "hx1.outlet", "sep1.vapour", "recycle_1"],
      ],
      ["Feeds", ["feed_1"]],
      ["Products", ["vapour_product"]],
    ]);
  });

  it("selects by the id the canvas selects by, which is the whole point of it", () => {
    const rows = navigationOf(demo, forms).flatMap((group) => group.rows);
    const nodes = new Set(demo.flowsheet.graph.nodes.map((entry) => entry.id));
    const edges = new Set(demo.flowsheet.graph.edges.map((entry) => entry.id));

    // Every row is one of the two things a selection can be, and every object has a row.
    expect(rows.every((row) => nodes.has(row.id) || edges.has(row.id))).toBe(true);
    expect(rows.filter((row) => nodes.has(row.id))).toHaveLength(nodes.size);
    expect(rows.filter((row) => edges.has(row.id))).toHaveLength(edges.size);
    // A stream's row is the edge, not the endpoint: `mix1.product` is `e1` here.
    expect(rows.find((row) => row.label === "mix1.product")?.id).toBe("e1");
  });

  it("says a tear is a recycle and says nothing about a connection", () => {
    const streams = navigationOf(demo, forms).find((group) => group.label === "Streams")?.rows ?? [];
    expect(streams.map((row) => [row.label, row.tag])).toEqual([
      ["feed_1", null],
      ["mix1.product", null],
      ["p1.outlet", null],
      ["hx1.outlet", null],
      ["sep1.vapour", null],
      ["recycle_1", "recycle"],
    ]);
  });

  it("files a unit op the catalogue does not carry under unfiled, and shows no empty family", () => {
    const nodes = [node("instance", "u1", "unit_ops.f0"), node("instance", "u2", "unit_ops.gone")];
    // `f0` is a `column` and `u1` is one, `f1` is a `utility` and nothing here is, and `f2` states
    // no family at all - so the group drawn is where the *nodes* are, not where the catalogue is.
    expect(shape(navigationOf(envelopeOf(nodes), withFamilies(["column", "utility", null])))).toEqual(
      [
        ["column", ["u1"]],
        ["unfiled", ["u2"]],
      ],
    );

    // A catalogue that carries no family at all puts every unit op in the one visible group.
    expect(shape(navigationOf(envelopeOf(nodes), withFamilies([null])))).toEqual([
      ["unfiled", ["u1", "u2"]],
    ]);
  });

  it("lists nothing for a document with nothing in it", () => {
    expect(navigationOf(envelopeOf([]), forms)).toEqual([]);
  });

  it("reads a form's family with one word for the wire that stopped carrying it", () => {
    expect(familyOfForm({ family: "column" } as never)).toBe("column");
    expect(familyOfForm({ family: null } as never)).toBe("unfiled");
  });

  it("spells a family the way the palette beside it spells it", () => {
    // The key is the spec's directory and the label is two words, which is why they are two
    // functions rather than one that would have to be undone to group by.
    expect(familyLabel("two_port")).toBe("two port");
    expect(familyLabel("column")).toBe("column");
  });
});
