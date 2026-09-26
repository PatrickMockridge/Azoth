/**
 * The property view's tab rule.
 *
 * **This is the one piece of the view that is a decision rather than a drawing**, and its failure
 * mode is the reason it is a function: the tab remembered for the kind you were in is not a tab the
 * object you just selected offers — you were reading a column's stages and selected a pump — and a
 * component that had to notice that inline is a component that eventually does not.
 */

import { expect, it } from "vitest";

import { kindOf, tabFor, tabsFor } from "../src/state/layout";
import type { GraphEdge, GraphNode } from "../src/wire/types";

const node = (role: GraphNode["role"]): GraphNode => ({
  id: `${role}:x`,
  type: role === "instance" ? "unit_op" : "stream",
  role,
  position: { x: 0, y: 0 },
  data: { name: "x" },
});

const edge = (kind: "connection" | "recycle"): GraphEdge => ({
  id: "e1",
  source: "instance:a",
  target: "instance:b",
  sourceHandle: "a.outlet",
  targetHandle: "b.feed",
  data: { kind, from: "a.outlet", to: "b.feed", path: kind === "recycle" ? "r1" : "a.outlet" },
});

it("reads a selection as one of the four kinds", () => {
  expect(kindOf(node("instance"), null)).toBe("instance");
  expect(kindOf(node("input"), null)).toBe("feed");
  expect(kindOf(node("product"), null)).toBe("product");
  expect(kindOf(null, edge("connection"))).toBe("edge");
  expect(kindOf(null, null)).toBeNull();
});

it("offers a tear a convergence sheet and a plain connection none", () => {
  // **The seven settings live on `[[recycles]]`**, so a connection's window offering Convergence
  // would be a tab whose only possible content is "there is nothing here".
  expect(tabsFor("edge", edge("recycle"))).toContain("convergence");
  expect(tabsFor("edge", edge("connection"))).not.toContain("convergence");
});

it("shows a unit operation's design first and a stream's conditions first", () => {
  expect(tabsFor("instance", null)[0]).toBe("design");
  expect(tabsFor("feed", null)[0]).toBe("conditions");
  expect(tabsFor("product", null)[0]).toBe("conditions");
  expect(tabsFor("edge", edge("connection"))[0]).toBe("connections");
});

it("keeps the tab remembered for the kind, and falls back where it is gone", () => {
  const available = tabsFor("instance", null);
  // Remembered and still offered.
  expect(tabFor("instance", { instance: "worksheet" }, available)).toBe("worksheet");
  // **Remembered and gone**: a column's stages on a pump, whose window has no such sheet.
  expect(tabFor("instance", { instance: "profiles" }, available)).toBe("design");
  // Nothing remembered for this kind, though something is for another.
  expect(tabFor("feed", { instance: "worksheet" }, tabsFor("feed", null))).toBe("conditions");
});
