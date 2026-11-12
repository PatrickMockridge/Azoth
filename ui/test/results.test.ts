/**
 * Reading a result without knowing which unit operation produced it.
 *
 * **The shape on the wire is the library's**, and this front end deliberately enumerates none of
 * it: a result is whatever that operation's model declares. So the classifiers are what have to be
 * right, and the case that matters is the one a per-unit-op table would have made impossible - a
 * field this file has never seen, of a value kind it has.
 *
 * The fixture is the CLI's own capture of the shipped demo, so the heater's duty below is the
 * library's number and not one written here.
 */

import { describe, expect, it } from "vitest";

import {
  outputRowsOf,
  resultOf,
  resultTabs,
  scalarsOf,
  seriesOf,
  warningsOf,
} from "../src/state/results";
import { unitsOf } from "../src/state/units";
import catalogue from "./fixtures/catalogue.json";
import fixture from "./fixtures/envelope.json";
import type { Catalogue, Envelope, UnitResult } from "../src/wire/types";

const envelope = fixture as unknown as Envelope;
const forms = catalogue as unknown as Catalogue;
const field = unitsOf(forms, "field");
const si = unitsOf(forms, "si");
const heater = resultOf(envelope, "hx1") as UnitResult;
const heaterForm = forms.unit_ops.find((entry) => entry.id === "unit_ops.heater");
const heaterOutputs = heaterForm?.outputs ?? [];

describe("the result reader", () => {
  it("finds the result an instance published, and none where it published none", () => {
    expect(resultOf(envelope, "hx1")).toBeDefined();
    // A mixer's whole answer is its stream, so it published nothing rather than a restatement.
    expect(resultOf(envelope, "mix1")).toBeUndefined();
    expect(resultOf(envelope, "nonesuch")).toBeUndefined();
  });

  it("reads a quantity as its magnitude and the unit the model named", () => {
    const duty = scalarsOf(heater).find((scalar) => scalar.key === "outlet_duty");
    expect(duty?.unit).toBe("W");
    expect(duty?.numeric).toBe(true);
    // **And a bare number is a row too.** A model's flow is a `float` in its own record, which is
    // the shape its Python dataclass publishes - so this reads it as a number with no unit rather
    // than as a quantity with an empty one.
    const flow = scalarsOf(heater).find((scalar) => scalar.key === "outlet_n");
    expect(flow?.unit).toBeNull();
    expect(flow?.text).not.toBe("");
    // The `warnings` field is a list and has its own reader, so it is not a scalar here.
    expect(scalarsOf(heater).map((scalar) => scalar.key)).not.toContain("warnings");
  });

  it("reads a vector of quantities and a vector of numbers alike", () => {
    expect(seriesOf(heater)).toEqual([]);

    const column: UnitResult = {
      tray_temperature: [
        { magnitude_si: 300, unit: "K" },
        { magnitude_si: 310, unit: "K" },
      ],
      tray_gas_n: [0.5, 0.25],
      converged: true,
    };
    const series = seriesOf(column);
    expect(series.map((entry) => entry.key)).toEqual(["tray_temperature", "tray_gas_n"]);
    expect(series[0]?.unit).toBe("K");
    expect(series[0]?.values).toEqual([300, 310]);
    // A vector of bare numbers carries no unit, and one of one entry is still a series - while a
    // scalar beside it is not one.
    expect(series[1]?.unit).toBeNull();
    expect(series.map((entry) => entry.key)).not.toContain("converged");
  });

  it("reads a warning as the library's own sentence", () => {
    expect(warningsOf(heater)).toEqual([]);
    const warned: UnitResult = {
      warnings: [
        { code: "SOLVER_NOT_CONVERGED", message: "the column hit its cap", field: "reboiler" },
        { code: "TRIVIAL_SOLUTION", message: "x = y = z" },
      ],
    };
    expect(warningsOf(warned)).toEqual([
      { code: "SOLVER_NOT_CONVERGED", message: "the column hit its cap", field: "reboiler" },
      { code: "TRIVIAL_SOLUTION", message: "x = y = z", field: null },
    ]);
    // A result with no warnings key at all, which is what a model without any publishes.
    expect(warningsOf({ outlet_n: 1 })).toEqual([]);
  });

  it("divides a scalar into the set's unit, and leaves the unit it cannot", () => {
    // The heater's duty is a power in watts, and the field set reads a power in horsepower.
    const watts = scalarsOf(heater).find((scalar) => scalar.key === "outlet_duty");
    const horse = scalarsOf(heater, field).find((scalar) => scalar.key === "outlet_duty");
    expect(watts?.unit).toBe("W");
    expect(horse?.unit).toBe("hp");
    // 745.6998715822702 W is one horsepower: the library's own number, and the same one
    // `test_units_cross_library.py` holds to `pint` on the Rust side.
    expect(Number(watts?.text) / Number(horse?.text)).toBeCloseTo(745.6998715822702, 4);
    // And a field with no unit is not converted, because there is no unit to convert it to.
    const flow = scalarsOf(heater, field).find((scalar) => scalar.key === "outlet_n");
    expect(flow?.unit).toBeNull();
  });

  it("offers each sheet only where it has something to draw", () => {
    // The heater published scalars and no vector.
    expect(resultTabs(heater, heaterOutputs)).toEqual(["results"]);
    // A column offers both, profile first.
    expect(
      resultTabs({ condenser_duty: { magnitude_si: 1, unit: "W" }, tray_temperature: [300, 310] }),
    ).toEqual(["stages", "results"]);
    // And a document that has not run has no result, so no sheets.
    expect(resultTabs(undefined)).toEqual([]);
    // A result whose every field is absent has nothing to draw either, so it offers no sheet.
    expect(resultTabs({ outlet_duty: null, instruments: [] }, heaterOutputs)).toEqual([]);
  });

  it("draws only the fields the run published, in the model's declared order", () => {
    const rows = outputRowsOf(heater, si, heaterOutputs);
    // hx1 declares 76 fields and publishes 8: an absent optional field is not a dash row.
    expect(rows).toHaveLength(8);
    expect(rows.map((row) => row.key)).toEqual([
      "outlet_n",
      "outlet_p",
      "outlet_t",
      "outlet_h",
      "outlet_duty",
      "max_utilization",
      "capacity_exceeded",
      "hard_limit_exceeded",
    ]);
    expect(rows.every((row) => row.text !== "—")).toBe(true);
    // A bare float takes the unit its declaration gives it.
    const flow = rows.find((row) => row.key === "outlet_n");
    expect(flow?.unit).toBe("mol/s");
    expect(flow?.text).toBe("1");
    // The human label and the library's sentence.
    expect(rows.find((row) => row.key === "outlet_duty")?.label).toBe("Outlet duty");
    expect(rows.find((row) => row.key === "outlet_duty")?.description).toContain("inlet_n");
    // `false` is a value and stays; an empty record list is not.
    expect(rows.find((row) => row.key === "capacity_exceeded")?.text).toBe("no");
    expect(rows.some((row) => row.key === "instruments")).toBe(false);
  });

  it("appends a field the form does not declare rather than dropping it", () => {
    const rows = outputRowsOf({ outlet_n: 1, mystery: 2 }, si, []);
    expect(rows.map((row) => row.key)).toEqual(["outlet_n", "mystery"]);
    expect(rows[1]?.kind).toBeNull();
    expect(rows[1]?.description).toBeNull();
  });

  it("reports a published record list as a count, not a table", () => {
    const rows = outputRowsOf(
      { instruments: [{ tag: "TIC-1" }, { tag: "TIC-2" }] },
      si,
      [{ name: "instruments", kind: "record_list", unit: null, dimension: null, values: [], description: "one row per instrument.", optional: false }],
    );
    expect(rows[0]?.text).toBe("2 records");
    expect(rows[0]?.records).toBe(2);
  });
});
