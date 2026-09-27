/**
 * The boundary's parser, held to the rule that makes it worth having.
 *
 * A decoder that accepts everything is a cast with more lines, and it would pass every test its
 * fixtures pass — `fixtures.test.ts` only ever hands it documents the library *did* write. So the
 * cases here are the other half: a document that is not the shape refuses, and the refusal names
 * the path it was at. Each case is one field moved by one kind.
 *
 * **Why the path matters as much as the refusal.** The failure this replaces is a cast, whose
 * symptom is `undefined` in a widget three files from the field that moved. A sentence naming
 * `the envelope.flowsheet.graph.nodes` is the difference between a boundary that reports and one
 * that hides.
 */

import { describe, expect, it } from "vitest";

import catalogueJson from "./fixtures/catalogue.json";
import envelopeJson from "./fixtures/envelope.json";
import { WireError, asRole, decodeCatalogue, decodeEnvelope, parseEnvelope } from "../src/wire/decode";

/** A copy of the shipped envelope with one field replaced — the shape a field moving in Rust takes. */
function moved(field: string, value: unknown = undefined): unknown {
  const copy = structuredClone(envelopeJson) as Record<string, unknown>;
  copy[field] = value;
  return copy;
}

/** The shipped envelope with one more diagnostic, whose target is the part under test. */
function withDiagnostic(target: Record<string, unknown>): unknown {
  const copy = structuredClone(envelopeJson) as { diagnostics: unknown[] };
  copy.diagnostics.push({
    code: "type_mismatch",
    severity: "error",
    section: "instances",
    path: "p1.ports.inlet",
    target,
    message: "anything at all",
    detail: { anything: [1, 2, { deeper: true }] },
  });
  return copy;
}

/** The refusal a document earns, or a failure saying it was accepted when it should not have been. */
function refused(value: unknown): WireError {
  try {
    decodeEnvelope(value);
  } catch (error) {
    if (error instanceof WireError) {
      return error;
    }
    throw error;
  }
  throw new Error("the decoder accepted a document that is not the shape");
}

describe("the decode", () => {
  it("parses the documents the library actually wrote", () => {
    // The whole point of the file: the same parse a browser runs, over the CLI's own output.
    expect(decodeEnvelope(envelopeJson).flowsheet.graph.nodes).toHaveLength(6);
    expect(decodeCatalogue(catalogueJson).unit_ops).toHaveLength(29);
  });

  it("refuses a document whose array is gone, and names the path", () => {
    const value = structuredClone(envelopeJson) as { flowsheet: { graph: Record<string, unknown> } };
    value.flowsheet.graph["nodes"] = undefined;
    const error = refused(value);
    expect(error.path).toBe("the envelope.flowsheet.graph.nodes");
    expect(error.message).toContain("wanted an array");
    expect(error.message).toContain("nothing at all");
  });

  it("refuses a value that is the wrong kind, and names the field", () => {
    const error = refused(moved("dirty", "no"));
    expect(error.path).toBe("the envelope.dirty");
    expect(error.message).toContain("wanted a boolean");
  });

  it("refuses an enum member the library does not have", () => {
    // A field added in Rust is a string this build does not know, which a cast would pass through
    // and a panel would then switch on and match nothing.
    const error = refused(moved("execution_order", "sideways"));
    expect(error.path).toBe("the envelope.execution_order");
    expect(error.message).toContain("wanted one of insertion, topological");
  });

  it("descends into a nested shape, so a stream's quantity is a quantity", () => {
    const value = structuredClone(envelopeJson) as {
      session: { streams: Record<string, { P: { magnitude_si: unknown } }> };
    };
    value.session.streams["p1.outlet"]!.P.magnitude_si = "500000";
    const error = refused(value);
    expect(error.path).toBe("the envelope.session.streams.p1.outlet.P.magnitude_si");
  });

  it("refuses a diagnostic target it cannot mark", () => {
    // The target is what a panel puts a mark *on*, so a kind nobody can draw is a refusal and not
    // a diagnostic drawn nowhere.
    const value = withDiagnostic({ kind: "somewhere_else" });
    const error = refused(value);
    expect(error.path).toBe("the envelope.diagnostics[0].target");
    expect(error.message).toContain("a target kind a canvas can mark");
  });

  it("carries the fields the library declares open, and decodes none of them", () => {
    // `detail`, a node's `parameters` and a published result are `Record<string, unknown>` by
    // declaration: a shape fixed here would be a second opinion about them rather than a reading.
    const value = withDiagnostic({ kind: "document" }) as {
      session: { results: Record<string, Record<string, unknown>> };
    };
    value.session.results["hx1"] = { a_shape: { this_file_does_not_know: 1 } };
    const parsed = decodeEnvelope(value);
    expect(parsed.diagnostics[0]?.detail).toEqual({ anything: [1, 2, { deeper: true }] });
    expect(parsed.session?.results["hx1"]).toEqual({ a_shape: { this_file_does_not_know: 1 } });
  });

  it("refuses the catalogue's own shapes as well", () => {
    const value = structuredClone(catalogueJson) as {
      unit_ops: { runnable: unknown }[];
    };
    value.unit_ops[0]!.runnable = "yes";
    expect(() => decodeCatalogue(value)).toThrow(WireError);
  });

  it("refuses a string that is not JSON at all", () => {
    expect(() => parseEnvelope("{oh no")).toThrow(WireError);
  });
});

describe("the role", () => {
  it("reads the three the projection makes, and nothing else", () => {
    // The parse `nodes.ts` used to cast: a fourth string is not a fourth role, and the switch
    // below it has no branch for one because this returns `null` instead of inventing it.
    expect(asRole("instance")).toBe("instance");
    expect(asRole("input")).toBe("input");
    expect(asRole("product")).toBe("product");
    expect(asRole("recycle")).toBeNull();
    expect(asRole("")).toBeNull();
  });
});
