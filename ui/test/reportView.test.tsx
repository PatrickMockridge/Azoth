// @vitest-environment jsdom
/**
 * The printable view and its drawing: every object is on the page, and Print is the browser's.
 */

import { readFileSync } from "node:fs";
import { resolve } from "node:path";
import { cleanup, fireEvent, render, screen } from "@testing-library/react";
import { afterEach, describe, expect, it, vi } from "vitest";

import { Pfd } from "../src/report/pfd";
import { ReportView } from "../src/components/report/ReportView";
import { unitsOf } from "../src/state/units";
import { parseCatalogue, parseEnvelope } from "../src/wire/decode";

const catalogue = parseCatalogue(readFileSync(resolve("test/fixtures/catalogue.json"), "utf8"));
const envelope = parseEnvelope(readFileSync(resolve("test/fixtures/envelope.json"), "utf8"));
const units = unitsOf(catalogue, "si");

afterEach(cleanup);

describe("the report drawing", () => {
  it("draws every object and names every connection", () => {
    const { container } = render(<Pfd envelope={envelope} />);
    for (const name of ["feed_1", "mix1", "p1", "hx1", "sep1", "vapour_product"]) {
      // A feed's name is also an edge's `<title>`, so the question is whether the drawing has
      // it, not whether it has it once.
      expect(screen.getAllByText(name).length).toBeGreaterThan(0);
    }
    const titles = [...container.querySelectorAll("title")].map((node) => node.textContent);
    expect(titles).toContain("mix1.product");
    expect(titles).toContain("recycle_1");
  });
});

describe("the report view", () => {
  it("renders the tables and prints through the browser", () => {
    const print = vi.fn();
    Object.defineProperty(window, "print", { value: print, configurable: true });
    const onClose = vi.fn();
    render(<ReportView envelope={envelope} catalogue={catalogue} units={units} onClose={onClose} />);

    expect(screen.getByRole("heading", { name: "Feed, pump, heat, separate, recycle" })).toBeDefined();
    expect(screen.getByRole("heading", { name: "Streams" })).toBeDefined();
    expect(screen.getAllByText("hx1.outlet").length).toBeGreaterThan(0);
    expect(screen.getByText(/not design deliverables/)).toBeDefined();

    fireEvent.click(screen.getByRole("button", { name: "Print / Save as PDF" }));
    expect(print).toHaveBeenCalledOnce();

    fireEvent.click(screen.getByRole("button", { name: "Close" }));
    expect(onClose).toHaveBeenCalledOnce();
  });
});
