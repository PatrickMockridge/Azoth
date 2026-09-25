/**
 * The TypeScript mirrors, held to the documents the library actually emits.
 *
 * **This is the gate `gen_stub.py` is for the Python side.** `src/wire/types.ts` is hand-written
 * from the Rust structs, and a hand-written mirror drifts — so the fixtures under `test/fixtures/`
 * are the CLI's own output, captured with `azoth forms --tools` and `azoth edit --json`, and this
 * file reads them as the types. A field renamed in Rust fails here rather than showing up as
 * `undefined` in a widget, which is what the same mistake looks like in a browser.
 *
 * The fixtures are committed, so the check runs with no Rust toolchain and no wasm build.
 */

import { describe, expect, it } from "vitest";

import type { Catalogue, Envelope } from "../src/wire/types";
import brokenJson from "./fixtures/broken.json";
import catalogueJson from "./fixtures/catalogue.json";
import envelopeJson from "./fixtures/envelope.json";
import { controlFor } from "../src/wire/field";

// **Imported, not read.** `vite` resolves a JSON import and TypeScript infers its shape, so the
// fixtures need no `@types/node` and the test says what it expects of them explicitly anyway.
const catalogue = catalogueJson as unknown as Catalogue;
const solved = envelopeJson as unknown as Envelope;
const refused = brokenJson as unknown as Envelope;

describe("the catalogue", () => {
  it("is the palette, with a model for all but two and a kernel for all but three", () => {
    expect(catalogue.unit_ops).toHaveLength(29);
    expect(catalogue.unit_ops.filter((entry) => entry.model !== null)).toHaveLength(27);
    expect(catalogue.unit_ops.filter((entry) => entry.runnable)).toHaveLength(26);
    expect(catalogue.tools).toHaveLength(15);
  });

  it("carries the kind a field needs, which is nowhere in the palette", () => {
    const pump = catalogue.unit_ops.find((entry) => entry.id === "unit_ops.pump");
    expect(pump).toBeDefined();
    const pressure = pump?.parameters.find((p) => p.name === "outlet_pressure");
    expect(pressure?.kind).toBe("quantity");
    expect(pressure?.dimension).toBe("pressure");
    expect(pressure?.unit).toBe("Pa");
    expect(pressure?.required).toBe(true);
    expect(controlFor(pressure?.kind ?? "unknown")).toBe("number");

    // The two entries no model declares are the only ones a field cannot choose a widget for.
    const unknown = catalogue.unit_ops.flatMap((entry) =>
      entry.parameters.filter((parameter) => parameter.kind === "unknown"),
    );
    expect(new Set(unknown.map((parameter) => parameter.name)).size).toBeGreaterThan(0);
    expect(
      catalogue.unit_ops
        .filter((entry) => entry.model === null)
        .flatMap((entry) => entry.parameters)
        .every((parameter) => parameter.kind === "unknown"),
    ).toBe(true);
  });

  it("carries a model's bound with the sentence that explains it", () => {
    const pump = catalogue.unit_ops.find((entry) => entry.id === "unit_ops.pump");
    const efficiency = pump?.parameters.find((p) => p.name === "isentropic_efficiency");
    const range = efficiency?.range[0];
    expect(range?.min).toBe(0);
    expect(range?.min_inclusive).toBe(false);
    expect(range?.max).toBe(1);
    expect(range?.max_inclusive).toBe(true);
    expect(range?.rationale).not.toBe("");
  });
});

describe("the envelope", () => {
  it("is the graph in the editor's own shape", () => {
    const graph = solved.flowsheet.graph;
    expect(graph.nodes.map((node) => node.id)).toEqual([
      "instance:mix1",
      "instance:p1",
      "instance:hx1",
      "instance:sep1",
      "input:feed_1",
      "product:vapour_product",
    ]);
    expect(graph.edges).toHaveLength(6);

    // A handle's id is the stream's own path, which is what makes a gesture and a value agree.
    const loop = graph.edges.find((edge) => edge.data.kind === "recycle");
    expect(loop?.sourceHandle).toBe("sep1.liquid");
    expect(loop?.data.path).toBe("recycle_1");

    // **A tear carries its own seven settings, and a silence is not a default.** The shipped
    // document states none, so every one is `null` - a panel that showed the class's defaults and
    // wrote them back would turn a silence into a pinned number.
    expect(loop?.data.settings).toEqual({
      flow_tolerance: null,
      composition_tolerance: null,
      temperature_tolerance: null,
      pressure_tolerance: null,
      max_iterations: null,
      minimum_flow: null,
      acceleration_method: null,
    });
    // And a connection has none, which is what makes the field a fact about a tear.
    expect(
      graph.edges.find((edge) => edge.data.kind === "connection")?.data.settings,
    ).toBeUndefined();

    // The feed's record is inline, in the units the schema declares.
    const feed = graph.nodes.find((node) => node.id === "input:feed_1");
    expect(feed?.data.input?.components).toEqual(["methane", "n-butane"]);
    expect(feed?.data.input?.P).toBe(5e5);
  });

  it("carries the session, the paths and the verdict", () => {
    expect(solved.ok).toBe(true);
    expect(solved.dirty).toBe(false);
    expect(solved.run_error).toBeNull();
    expect(solved.diagnostics).toEqual([]);
    expect(solved.session?.converged).toBe(true);
    expect(solved.session?.iterations).toBe(2);
    // The order the next run takes is the class's default, and it is in the envelope so a top bar
    // reads it from the session rather than holding a copy.
    expect(solved.execution_order).toBe("insertion");
    expect(solved.paths).toContain("p1.outlet.P");
    // The position the edit set is in the document *and* in the graph, one value read twice.
    expect(solved.flowsheet.document).toContain("[layout.instances]");
    expect(
      solved.flowsheet.graph.nodes.find((node) => node.id === "instance:sep1")?.position,
    ).toEqual({ x: 1234, y: 56 });
  });

  it("carries a diagnostic with a target a canvas can select", () => {
    expect(refused.ok).toBe(false);
    expect(refused.diagnostics.map((diagnostic) => diagnostic.code)).toEqual([
      "underfed_port",
      "underfed_port",
    ]);
    expect(refused.diagnostics[0]?.target).toEqual({
      kind: "handle",
      role: "instance",
      node: "p1",
      port: "outlet",
      index: null,
    });
    expect(refused.diagnostics[0]?.severity).toBe("error");
    expect(refused.diagnostics[0]?.section).toBe("instances");
    expect(refused.diagnostics[0]?.path).toBe("p1.ports.outlet");
    expect(refused.diagnostics[0]?.message).not.toBe("");
  });
});
