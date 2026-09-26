// @vitest-environment jsdom
/**
 * The one multiplication a display performs.
 *
 * **What this holds is the arithmetic and the fallbacks**, both of which are decisions rather than
 * readings: a unit whose set has no entry is left alone, an unknown unit is left alone, and a
 * remembered set the catalogue no longer carries is the first set rather than no conversion. The
 * factors themselves are the library's - the fixture is the CLI's own catalogue capture, so the
 * `psi` below is `si_factor("psi")` and not a number written here.
 */

import { describe, expect, it } from "vitest";

import { displayOf, inDisplayUnit, readUnitSet, rememberUnitSet, unitsOf } from "../src/state/units";
import fixture from "./fixtures/catalogue.json";
import type { Catalogue } from "../src/wire/types";

const catalogue = fixture as unknown as Catalogue;

describe("the unit reader", () => {
  it("reads a quantity in the set's unit, and divides by the library's factor", () => {
    const field = unitsOf(catalogue, "field");

    // 1 bar is 100 000 Pa, and a field set reads a pressure in psi: 100000 / 6894.757… = 14.5038.
    expect(displayOf(field, "Pa")).toEqual({ unit: "psi", factor: 6894.757293168361 });
    expect(inDisplayUnit(field, 100000, "Pa")).toBeCloseTo(14.5037738, 6);

    // And the SI set is the identity for the same quantity, which is what makes switching a
    // change of unit rather than a change of value.
    const si = unitsOf(catalogue, "si");
    expect(displayOf(si, "Pa")).toEqual({ unit: "Pa", factor: 1 });
    expect(inDisplayUnit(si, 100000, "Pa")).toBe(100000);
  });

  it("finds the set's unit through the dimension the catalogue gives the quantity's unit", () => {
    const field = unitsOf(catalogue, "field");
    // `kmol/h` is the field spelling of a molar flow, and `kJ/mol` the one it leaves alone: a
    // dimension the sets do not name falls through to the library's own unit.
    expect(displayOf(field, "mol/s").unit).toBe("kmol/h");
    expect(displayOf(field, "J/mol").unit).toBe("kJ/mol");
    // **A dimension no set names is unchanged**, and this is the case that keeps a strange unit
    // readable: the attraction parameter is read in its own units or in none.
    expect(displayOf(field, "Pa*m**6/mol**2")).toEqual({ unit: "Pa*m**6/mol**2", factor: 1 });
  });

  it("leaves a unit it has never heard of alone rather than guessing", () => {
    const field = unitsOf(catalogue, "field");
    expect(displayOf(field, "furlong/fortnight")).toEqual({
      unit: "furlong/fortnight",
      factor: 1,
    });
    // And a catalogue that has not loaded yet, which is the state the first render is in.
    const none = unitsOf(null, "field");
    expect(displayOf(none, "Pa")).toEqual({ unit: "Pa", factor: 1 });
    expect(none.sets).toEqual([]);
  });

  it("falls back to the first set when the remembered one is not the library's", () => {
    // A store written by an older build, or a vocabulary that lost a set: a reader who asked for
    // field units and gets SI is better served than one who gets a number with no unit at all.
    expect(unitsOf(catalogue, "hectares").active).toBe("si");
    expect(unitsOf(catalogue, null).active).toBe("si");
    expect(unitsOf(catalogue, "field").active).toBe("field");
  });

  it("refuses a factor that is not a unit", () => {
    // Neither of these is a conversion the library declares, and dividing by one would be a
    // division by zero rather than a unit.
    const broken = unitsOf(
      {
        unit_ops: [],
        units: { Pa: { dimension: "pressure", factor: 0 }, K: { dimension: null, factor: null } },
        unit_sets: [{ id: "si", name: "SI", units: { pressure: "Pa", thermodynamic_temperature: "K" } }],
      },
      "si",
    );
    expect(displayOf(broken, "Pa")).toEqual({ unit: "Pa", factor: 1 });
    expect(displayOf(broken, "K")).toEqual({ unit: "K", factor: 1 });
  });

  it("remembers the choice under one key", () => {
    window.localStorage.clear();
    expect(readUnitSet()).toBeNull();
    rememberUnitSet("field");
    expect(readUnitSet()).toBe("field");
    expect(window.localStorage.getItem("azoth.units")).toBe("field");
  });
});
