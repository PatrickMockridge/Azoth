// @vitest-environment jsdom
/**
 * What a unit operation rests on, as the window shows it.
 *
 * **The fixture is the library's own document**, so the assertions here are about the values
 * `azoth` actually publishes for the shipped demo rather than about values arranged for a test:
 * the heater's source line is the NeqSim class the palette carries, its verification is the
 * status its model declares, and its two hashes are the spec and the kernel the answer came
 * from. A panel that showed something else would be a second record of the same fact.
 */

import { render, screen } from "@testing-library/react";
import { describe, expect, it } from "vitest";

import envelopeJson from "./fixtures/envelope.json";
import { ProvenanceSheet } from "../src/components/property/ProvenanceSheet";
import { provenanceOf, shortHash } from "../src/state/provenance";
import { decodeEnvelope } from "../src/wire/decode";

const envelope = decodeEnvelope(envelopeJson);

describe("what a document rests on", () => {
  it("reads an instance's entry out of the envelope it came in", () => {
    const heater = provenanceOf(envelope, "hx1");
    expect(heater?.unit).toBe("unit_ops.heater");
    expect(heater?.model).toBe("process.heater");
    // The provenance line the palette carries: the class the port names.
    expect(heater?.source).toContain("NeqSim");
    // **The library's own status word**, spelled as the model declares it.
    expect(heater?.verification).toBe("partially_verified");
    expect(heater?.runnable).toBe(true);
    expect(heater?.refusal).toBeNull();
    // Two digests and not one: a spec can change without the kernel changing, and the entry
    // says which of the two moved.
    expect(heater?.spec_sha256).toMatch(/^[0-9a-f]{64}$/);
    expect(heater?.rust_sha256).toMatch(/^[0-9a-f]{64}$/);
  });

  it("has an entry per instantiated unit operation, and none for a name nobody declared", () => {
    expect(envelope.provenance.units.map((entry) => entry.instance)).toEqual([
      "mix1",
      "p1",
      "hx1",
      "sep1",
    ]);
    // A name the document does not carry is `undefined` rather than a fabricated entry: the
    // record is written per instance, and an instance that is not there has nothing to say.
    expect(provenanceOf(envelope, "nosuch")).toBeUndefined();
    // And the library itself, which is what "azoth" is on this sheet.
    expect(envelope.provenance.library).toBe("azoth");
    expect(envelope.provenance.version).not.toBe("");
  });

  it("shows eight characters of a digest and the whole one in the title", () => {
    // A hash is shown shortened because it is for recognising, and the whole of it is a hover
    // away because comparing is what it is *for*.
    expect(shortHash("7307b12a6e89888287c2bc672d75e2fc5ae18723aa450bc44e68072a1553bff9")).toBe(
      "7307b12a",
    );
    // **`null` stays null.** A palette entry with no model has no hashes, and eight characters
    // of nothing would be a value where the library wrote an absence.
    expect(shortHash(null)).toBeNull();
  });

  it("renders the source and the verification, in the library's words", () => {
    render(<ProvenanceSheet envelope={envelope} instance="hx1" />);
    expect(screen.getByText("NeqSim process/equipment/heatexchanger/Heater.java")).toBeDefined();
    expect(screen.getByText("partially_verified")).toBeDefined();
    expect(screen.getByText("process.heater")).toBeDefined();
  });

  it("says so, rather than showing nothing, where the document has no entry", () => {
    // A feed has no provenance entry — the record is per instantiated unit operation — and a
    // sheet that rendered an empty table would look like a panel that had failed to load.
    render(<ProvenanceSheet envelope={envelope} instance="feed_1" />);
    expect(screen.getByText(/no provenance for/)).toBeDefined();
  });
});
