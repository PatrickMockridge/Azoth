/**
 * A diagnostic's target, as the node a click selects.
 *
 * The interesting cases are the ones that select nothing: a target is derived from the checker's
 * own variant fields, and where it names no node — an endpoint nobody declared, a palette entry, a
 * fact about the document — a front-end that selected something arbitrary would be pointing a user
 * at the wrong thing to look responsive.
 */

import { describe, expect, it } from "vitest";

import envelopeJson from "./fixtures/envelope.json";
import { selectedNode, targetNodeId } from "../src/state/selection";
import type { Envelope, Target } from "../src/wire/types";

// **Imported, not read.** `vite` resolves a JSON import and TypeScript infers its shape, so the
// fixture needs no `@types/node`.
const envelope = envelopeJson as unknown as Envelope;

describe("targetNodeId", () => {
  it("selects the node a target names", () => {
    expect(targetNodeId({ kind: "node", role: "instance", id: "sep1" }, envelope)).toBe(
      "instance:sep1",
    );
    expect(
      targetNodeId(
        { kind: "handle", role: "instance", node: "mix1", port: "feed", index: null },
        envelope,
      ),
    ).toBe("instance:mix1");
    expect(
      targetNodeId({ kind: "parameter", node: "p1", name: "outlet_pressure" }, envelope),
    ).toBe("instance:p1");
    expect(targetNodeId({ kind: "node", role: "input", id: "feed_1" }, envelope)).toBe(
      "input:feed_1",
    );
  });

  it("follows an edge to the node that produces the stream", () => {
    const target: Target = { kind: "edge", from: "p1.outlet", to: "hx1.inlet" };
    expect(targetNodeId(target, envelope)).toBe("instance:p1");
  });

  it("selects nothing where the target names no node", () => {
    const nothing: Target[] = [
      // A name nobody declared — the node it implies does not exist to select.
      { kind: "edge", from: "nowhere", to: "mix1.feed" },
      { kind: "endpoint", endpoint: "purge" },
      { kind: "palette", id: "unit_ops.pump", port: null, parameter: null, field: null },
      { kind: "document" },
    ];
    for (const target of nothing) {
      expect(targetNodeId(target, envelope)).toBeNull();
    }
  });
});

describe("selectedNode", () => {
  it("answers with the node, or none", () => {
    expect(selectedNode(envelope, "instance:sep1")?.data.unit).toBe("unit_ops.separator");
    expect(selectedNode(envelope, "instance:nope")).toBeNull();
    expect(selectedNode(envelope, null)).toBeNull();
  });
});
