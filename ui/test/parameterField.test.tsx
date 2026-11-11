// @vitest-environment jsdom
/**
 * A parameter field commits when the typing stops.
 *
 * The defect this guards is not a wrong number: it is three commands and two intermediate
 * documents for one value, which on a served door is three round trips and on either door is a
 * checker state the reader never meant to create.
 */

import { readFileSync } from "node:fs";
import { resolve } from "node:path";
import { cleanup, fireEvent, render, screen } from "@testing-library/react";
import { afterEach, describe, expect, it, vi } from "vitest";

import { ParameterField } from "../src/components/ParameterField";
import { unitsOf } from "../src/state/units";
import { parseCatalogue } from "../src/wire/decode";

const catalogue = parseCatalogue(readFileSync(resolve("test/fixtures/catalogue.json"), "utf8"));
const units = unitsOf(catalogue, "si");
const parameter = catalogue.unit_ops
  .find((entry) => entry.id === "unit_ops.pump")
  ?.parameters.find((entry) => entry.name === "isentropic_efficiency");

afterEach(cleanup);

function field(onChange: (value: unknown) => void) {
  if (parameter === undefined) {
    throw new Error("the pump's efficiency is not in the catalogue");
  }
  return render(
    <ParameterField parameter={parameter} value={0.75} units={units} onChange={onChange} />,
  );
}

describe("a parameter field", () => {
  it("sends nothing while the reader is typing, and one value on blur", () => {
    const onChange = vi.fn();
    const { container } = field(onChange);
    const input = container.querySelector("input");
    expect(input).not.toBeNull();

    fireEvent.change(input as HTMLInputElement, { target: { value: "0.8" } });
    expect(onChange).not.toHaveBeenCalled();

    fireEvent.blur(input as HTMLInputElement);
    expect(onChange).toHaveBeenCalledOnce();
    expect(onChange).toHaveBeenCalledWith(0.8);
  });

  it("commits on Enter and abandons on Escape", () => {
    const onChange = vi.fn();
    const { container } = field(onChange);
    const input = container.querySelector("input") as HTMLInputElement;

    fireEvent.change(input, { target: { value: "0.9" } });
    fireEvent.keyDown(input, { key: "Escape" });
    fireEvent.blur(input);
    expect(onChange).not.toHaveBeenCalled();
    expect(input.value).toBe("0.75");

    fireEvent.change(input, { target: { value: "0.8" } });
    fireEvent.keyDown(input, { key: "Enter" });
    expect(onChange).toHaveBeenCalledWith(0.8);
  });

  it("puts the explanation in a row as the control's visible label", () => {
    const rho = catalogue.calcs
      ?.find((entry) => entry.id === "hydraulics.pump_power")
      ?.inputs.find((entry) => entry.name === "rho");
    if (rho === undefined) {
      throw new Error("the pump's density input is not in the catalogue");
    }
    const { container } = render(
      <ParameterField layout="row" parameter={rho} value={998} units={units} onChange={() => {}} />,
    );

    const description = screen.getByText(/density of the pumped fluid/i);
    const input = container.querySelector("input");
    expect(description.className).toBe("description");
    expect(input?.getAttribute("aria-describedby")).toContain(description.id);
    expect(screen.getByText(/checked limits/i)).toBeDefined();
    // **The unit is beside the value, not in the text column.** Its display spelling is the
    // engineering one; the raw catalogue string stays on the title for a reader who needs it.
    const unit = container.querySelector(".field-control .field-value .unit");
    expect(unit?.textContent).toBe("kg/m³");
    expect(unit?.getAttribute("title")).toBe("kg/m**3");
    expect(container.querySelector(".field-text .unit")).toBeNull();
    expect(container.querySelector(".status.required")?.textContent).toBe("required");
  });

  it("marks an optional field persistently", () => {
    const separator = catalogue.unit_ops?.find((entry) => entry.id === "unit_ops.separator");
    const optional = separator?.parameters.find((entry) => entry.name === "heat_input");
    if (optional === undefined) {
      throw new Error("the separator's heat input is not in the catalogue");
    }
    const { container } = render(
      <ParameterField
        layout="row"
        parameter={optional}
        value={0}
        units={units}
        onChange={() => {}}
      />,
    );
    expect(container.querySelector(".status.optional")?.textContent).toBe("optional");
  });

  it("shows no unit label for a dimensionless value", () => {
    if (parameter === undefined) {
      throw new Error("the pump's efficiency is not in the catalogue");
    }
    const { container } = render(
      <ParameterField layout="row" parameter={parameter} value={0.75} units={units} onChange={() => {}} />,
    );
    expect(container.querySelector(".field-value .unit")).toBeNull();
  });

  it("keeps the flowsheet's stacked field unit in its head", () => {
    const pump = catalogue.unit_ops?.find((entry) => entry.id === "unit_ops.pump");
    const pressure = pump?.parameters.find((entry) => entry.name === "outlet_pressure");
    if (pressure === undefined) {
      throw new Error("the pump's outlet pressure is not in the catalogue");
    }
    const { container } = render(
      <ParameterField
        parameter={pressure}
        value={1_000_000}
        units={units}
        onChange={() => {}}
      />,
    );
    expect(container.querySelector(".field .head .unit")?.textContent).toBe("Pa");
    expect(container.querySelector(".field-control")).toBeNull();
  });
});
