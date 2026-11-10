// @vitest-environment jsdom
/**
 * The feed form: a picker, an explicit sum, and one Apply.
 *
 * The document is not in this test; the form's contract is the record it hands back, which
 * `app.test.tsx` sends through the real module.
 */

import { readFileSync } from "node:fs";
import { resolve } from "node:path";
import { cleanup, fireEvent, render, screen } from "@testing-library/react";
import { afterEach, describe, expect, it, vi } from "vitest";

import { FeedEditor } from "../src/components/feed/FeedEditor";
import { emptyDraft } from "../src/state/composition";
import { unitsOf } from "../src/state/units";
import { parseCatalogue } from "../src/wire/decode";

const catalogue = parseCatalogue(readFileSync(resolve("test/fixtures/catalogue.json"), "utf8"));
const units = unitsOf(catalogue, "si");

afterEach(cleanup);

function pick(name: string) {
  fireEvent.change(screen.getByLabelText("Search substances"), { target: { value: name } });
  fireEvent.click(screen.getByRole("button", { name }));
}

describe("the feed editor", () => {
  it("picks substances and sends one record", () => {
    const onApply = vi.fn();
    render(
      <FeedEditor
        catalogue={catalogue}
        units={units}
        initial={emptyDraft("feed_1")}
        editing={false}
        onApply={onApply}
        onCancel={() => {}}
      />,
    );

    pick("methane");
    fireEvent.click(screen.getByRole("button", { name: "Add the feed" }));

    expect(onApply).toHaveBeenCalledOnce();
    const record = onApply.mock.calls[0]?.[0] as {
      name: string;
      components: string[];
      z: number[];
    };
    expect(record.name).toBe("feed_1");
    expect(record.components).toEqual(["methane"]);
    expect(record.z).toEqual([1]);
  });

  it("normalises only when the reader asks", () => {
    const onApply = vi.fn();
    render(
      <FeedEditor
        catalogue={catalogue}
        units={units}
        initial={emptyDraft("feed_1")}
        editing={false}
        onApply={onApply}
        onCancel={() => {}}
      />,
    );

    pick("methane");
    pick("ethane");
    fireEvent.change(screen.getByLabelText("ethane fraction"), { target: { value: "1" } });
    expect(screen.getByRole("button", { name: "Add the feed" })).toHaveProperty("disabled", true);

    fireEvent.click(screen.getByRole("button", { name: "Normalise" }));
    fireEvent.click(screen.getByRole("button", { name: "Add the feed" }));

    const record = onApply.mock.calls[0]?.[0] as { z: number[] };
    expect(record.z[0]).toBeCloseTo(0.5, 12);
    expect(record.z[1]).toBeCloseTo(0.5, 12);
  });

  it("converts mass fractions to the mole fractions the model receives", () => {
    const onApply = vi.fn();
    render(
      <FeedEditor
        catalogue={catalogue}
        units={units}
        initial={emptyDraft("feed_1")}
        editing={false}
        onApply={onApply}
        onCancel={() => {}}
      />,
    );

    pick("methane");
    pick("ethane");
    fireEvent.change(screen.getByLabelText("methane fraction"), { target: { value: "0.5" } });
    fireEvent.change(screen.getByLabelText("ethane fraction"), { target: { value: "0.5" } });
    fireEvent.click(screen.getByRole("button", { name: "mass fractions" }));
    fireEvent.click(screen.getByRole("button", { name: "Add the feed" }));

    const record = onApply.mock.calls[0]?.[0] as { z: number[] };
    expect(record.z).toHaveLength(2);
    expect(record.z[0]).toBeGreaterThan(record.z[1] as number);
    expect((record.z[0] as number) + (record.z[1] as number)).toBeCloseTo(1, 12);
  });

  it("removes a row and disables Apply when nothing is left", () => {
    render(
      <FeedEditor
        catalogue={catalogue}
        units={units}
        initial={emptyDraft("feed_1")}
        editing={false}
        onApply={() => {}}
        onCancel={() => {}}
      />,
    );
    pick("methane");
    fireEvent.click(screen.getByRole("button", { name: "remove methane" }));
    expect(screen.getByRole("button", { name: "Add the feed" })).toHaveProperty("disabled", true);
    expect(screen.getByText("add at least one substance")).toBeDefined();
  });
});
