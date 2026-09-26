/**
 * Where a composition table's substance axis comes from.
 *
 * **A solved stream's `z` is a bare array**: the record is `n, z, P, T, h` and the names live on
 * the *feed* that declared the fluid, so the axis has to be recovered from the document - and the
 * case that matters is the one where the document does not agree with itself, where the honest
 * answer is a position rather than a name.
 */

import { expect, it } from "vitest";

import { componentAxis, fraction } from "../src/state/columns";
import fixture from "./fixtures/envelope.json";
import type { Envelope, GraphNode } from "../src/wire/types";

const envelope = fixture as unknown as Envelope;

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
