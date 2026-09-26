/**
 * The TypeScript mirrors, held to the documents the library actually emits.
 *
 * **This is the gate `gen_stub.py` is for the Python side.** `src/wire/types.ts` is hand-written
 * from the Rust structs, and a hand-written mirror drifts — so the fixtures under `test/fixtures/`
 * are the CLI's own output, and this file reads them as the types. A field renamed in Rust fails
 * here rather than showing up as `undefined` in a widget, which is what the same mistake looks like
 * in a browser.
 *
 * **How to recapture them, exactly** — the commands are the whole pin, because a fixture captured
 * some other way asserts something about a document nobody can reproduce:
 *
 * ```bash
 * cargo run -p azoth-cli -- forms --tools > ui/test/fixtures/catalogue.json
 * cargo run -p azoth-cli -- edit --flowsheet specs/flowsheets/demo.toml \
 *     --command '{"command":"set_position","node":"instance:sep1","x":1234,"y":56}' \
 *     --run --json > ui/test/fixtures/envelope.json
 * ```
 *
 * The edit is a **position**, which is why the assertions below can say that one gesture is in the
 * document *and* in the graph: it is a command that changes the layout and nothing else, so the run
 * it reports is the shipped demo's own.
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

  it("says which family each entry was filed under", () => {
    // **The grouping a palette panel draws, and the only place it is stated.** An id is
    // `unit_ops.<leaf>`, so a reader that split the id took the leaf for a family and drew one
    // group called `other` - and `source` is NeqSim's taxonomy rather than this palette's, where
    // `cooler` is a `two_port` entry inside NeqSim's `heatexchanger/` directory.
    const families = new Set(catalogue.unit_ops.map((entry) => entry.family));
    for (const entry of catalogue.unit_ops) {
      expect(entry.family, `${entry.id} names no family`).not.toBeNull();
    }
    expect([...families].sort()).toEqual([
      "column",
      "heat_exchanger",
      "mixer",
      "reactor",
      "separator",
      "two_port",
      "utility",
    ]);
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

  it("carries what the run computed beside the streams, under the model's own names", () => {
    // **The heater's duty is the number no outlet stream carries**, and it is why the executor
    // had to start publishing results at all: the kernel computed it and the dispatcher dropped it.
    const results = solved.session?.results ?? {};
    expect(Object.keys(results)).toEqual(["hx1"]);
    expect(results.hx1?.outlet_duty).toMatchObject({ unit: "W" });
    // A mixer's whole answer *is* its stream, so it publishes nothing rather than a restatement.
    expect(results).not.toHaveProperty("mix1");
    expect(results).not.toHaveProperty("sep1");
  });

  it("carries the three fields a stream record could not derive", () => {
    // Mass flow and molar mass are the library's own (a component without a molar mass is `null`,
    // which is why the mirror admits it), and the vapour fraction is what the *flash* established -
    // the one field of the three that no record can work out for itself.
    const vapour = solved.session?.streams["sep1.vapour"];
    expect(vapour?.mass_flow?.unit).toBe("kg/s");
    expect(vapour?.molar_mass?.unit).toBe("kg/mol");
    expect(vapour?.vapour_fraction).toBe(1);
    // **And the endpoints are the phase's, not the flash's extrapolated root.** The demo's feed is
    // a gas at 300 K and 5 bar whose Rachford-Rice root converges to `1.9847`; a fraction outside
    // `[0, 1]` is a split that does not exist, and a separator's two outlets are the two endpoints.
    expect(solved.session?.streams["feed_1"]?.vapour_fraction).toBe(1);
    expect(solved.session?.streams["sep1.liquid"]?.vapour_fraction).toBe(0);
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
