/**
 * The field logic: a kind chooses a control, a bound marks a value, and a command carries the kind.
 *
 * What is worth testing here is the half that would otherwise be discovered in a browser: an
 * efficiency bounded `(0, 1]` must not flag one — the reversible limit the capture itself measures
 * — and a vector must reach the document as a list rather than as the string a text input holds.
 */

import { describe, expect, it } from "vitest";

import { boundsOf, commandValue, controlFor, fieldText, formatQuantity, violates } from "../src/wire/field";
import type { FormParameter } from "../src/wire/types";

function parameter(overrides: Partial<FormParameter> = {}): FormParameter {
  return {
    name: "isentropic_efficiency",
    kind: "quantity",
    required: true,
    unit: "dimensionless",
    dimension: "dimensionless",
    values: [],
    description: "",
    range: [],
    ...overrides,
  };
}

describe("controlFor", () => {
  it("maps the five kinds a palette declares, and admits the two it cannot", () => {
    expect(controlFor("quantity")).toBe("number");
    expect(controlFor("vector")).toBe("list");
    expect(controlFor("boolean")).toBe("switch");
    expect(controlFor("enum")).toBe("select");
    expect(controlFor("string")).toBe("text");
    // A kind no model declares: the field says so rather than guessing at a number box.
    expect(controlFor("unknown")).toBe("unknown");
    expect(controlFor("matrix")).toBe("unknown");
  });
});

describe("boundsOf", () => {
  it("keeps each endpoint's own inclusivity", () => {
    const bounds = boundsOf(
      parameter({
        range: [
          {
            min: 0,
            min_inclusive: false,
            max: 1,
            max_inclusive: true,
            equals: null,
            band: "outside",
            severity: "error",
            code: "OUT_OF_VALID_RANGE",
            rationale: "the reversible step divides by it",
          },
        ],
      }),
    );
    expect(bounds.minExclusive).toBe(true);
    expect(bounds.maxExclusive).toBe(false);
    expect(bounds.rationale).toBe("the reversible step divides by it");

    // **One is the reversible limit and it is legal** — a single flag for both ends would mark
    // the capture's own efficiency-of-one row as a violation.
    expect(violates(bounds, 1)).toBe(false);
    expect(violates(bounds, 0.75)).toBe(false);
    expect(violates(bounds, 0)).toBe(true);
    expect(violates(bounds, 1.01)).toBe(true);
  });

  it("marks nothing where the model states no bound", () => {
    const bounds = boundsOf(parameter());
    expect(bounds.severity).toBeNull();
    expect(violates(bounds, -12)).toBe(false);
  });

  it("ignores a value that is not a number", () => {
    const bounds = boundsOf(
      parameter({
        range: [
          {
            min: 0,
            min_inclusive: false,
            max: null,
            max_inclusive: true,
            equals: null,
            band: "outside",
            severity: "error",
            code: "OUT_OF_VALID_RANGE",
            rationale: "an absolute temperature",
          },
        ],
      }),
    );
    // Text is what the checker refuses as a `parameter_kind`; a range is not the rule for it.
    expect(violates(bounds, "high")).toBe(false);
  });
});

describe("commandValue", () => {
  it("shapes a value the way its declaration describes", () => {
    expect(commandValue("quantity", "2.0e6")).toBe(2_000_000);
    expect(commandValue("vector", "0.5, 0.3, 0.2")).toEqual([0.5, 0.3, 0.2]);
    expect(commandValue("boolean", "true")).toBe(true);
    expect(commandValue("boolean", "false")).toBe(false);
    expect(commandValue("enum", "duty")).toBe("duty");
    expect(commandValue("string", "methane")).toBe("methane");
  });

  it("passes text through where a number was expected, for the checker to refuse", () => {
    // The field does not decide that the value is wrong; `parameter_kind` is the library's rule
    // and this is what makes it reachable.
    expect(commandValue("quantity", "high")).toBe("high");
    expect(commandValue("quantity", "")).toBe("");
  });

  it("reads a value the document already holds back into the field", () => {
    expect(fieldText(2_000_000)).toBe("2000000");
    expect(fieldText([0.5, 0.3])).toBe("0.5, 0.3");
    expect(fieldText(true)).toBe("true");
    expect(fieldText(undefined)).toBe("");
  });
});

describe("formatQuantity", () => {
  it("shortens a magnitude without inventing digits", () => {
    expect(formatQuantity(2_000_000, "Pa", 3)).toBe("2000000 Pa");
    expect(formatQuantity(416.0115678, "K")).toBe("416 K");
    expect(formatQuantity(0.0833, "mol/s")).toBe("0.0833 mol/s");
  });
});
