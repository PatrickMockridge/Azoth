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
import { cleanup, fireEvent, render } from "@testing-library/react";
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
});
