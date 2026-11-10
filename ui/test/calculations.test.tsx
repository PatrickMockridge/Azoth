// @vitest-environment jsdom
/**
 * The Calculations workspace, driven through the real wasm module.
 *
 * **The same discipline as `app.test.tsx`.** `fetch` serves the module this build produced, the
 * demo book is the shipped `specs/calcbooks/demo.toml`, and the assertions are about the DOM. The
 * worked examples are the calculations' own, so the number on the screen is the number the spec
 * pins.
 */

import { readFileSync } from "node:fs";
import { resolve } from "node:path";
import { cleanup, fireEvent, render, screen, waitFor, within } from "@testing-library/react";
import { afterEach, beforeAll, describe, expect, it, vi } from "vitest";

import { only, stubDownload, textOf, unstubDownload } from "./download";
import { nodeIds, stubXyflowEnvironment } from "./xyflow-env";

beforeAll(stubXyflowEnvironment);

function stubFetch(): void {
  const wasm = readFileSync(resolve("src/wasm/pkg/azoth_wasm_bg.wasm"));
  globalThis.fetch = (input: RequestInfo | URL): Promise<Response> => {
    expect(String(input)).toContain("azoth_wasm_bg.wasm");
    return Promise.resolve(
      new Response(wasm, { headers: { "Content-Type": "application/wasm" } }),
    );
  };
}

async function open() {
  window.localStorage.setItem("azoth.welcome", "seen");
  window.localStorage.removeItem("azoth.units");
  const { App } = await import("../src/App");
  const rendered = render(<App />);
  await waitFor(() => {
    expect(rendered.container.querySelectorAll(".react-flow__node")).not.toHaveLength(0);
  });
  fireEvent.click(screen.getByRole("button", { name: "Calculations" }));
  await waitFor(() => {
    expect(screen.getAllByRole("tab").map((tab) => tab.textContent)).toEqual([
      "Pump shaft power",
      "Darcy-Weisbach pressure drop",
    ]);
  });
  return rendered;
}

afterEach(() => {
  cleanup();
  vi.unstubAllGlobals();
  unstubDownload();
});

describe("the calculations workspace", () => {
  it("opens the shipped worked-example book and calculates its first sheet", async () => {
    stubFetch();
    await open();

    const sheet = screen.getByRole("tabpanel");
    // **The sheet has a height to scroll within.** Its panel is the flex child that constrains
    // it; without that the content grew past the clipped main column.
    expect(sheet.className).toBe("calc-sheet-panel");
    expect(sheet.querySelector(".calc-sheet")).not.toBeNull();
    // **The explanation is on the sheet before the run.** A reader meets what the calculation is
    // for, what each field means, what the output will be, and when the answer applies.
    expect(
      screen.getAllByText(/shaft power a pump needs from flow, head and efficiency/i).length,
    ).toBeGreaterThan(0);
    expect(within(sheet).getByText(/density of the pumped fluid/i)).toBeDefined();
    expect(within(sheet).getByText(/shaft power the pump must be supplied with/i)).toBeDefined();
    expect(within(sheet).getByText("When this applies")).toBeDefined();
    const manual = within(sheet).getByRole("link", { name: /read the full manual/i });
    expect(manual.getAttribute("href")).toMatch(/\/manual\/hydraulics\/pump_power\.html$/);

    fireEvent.click(within(sheet).getByRole("button", { name: "Calculate" }));
    await waitFor(() => expect(within(sheet).getByText(/3914\.8/)).toBeDefined());
    // No warning on the worked example, and the model check is present.
    expect(within(sheet).queryByText("Warnings")).toBeNull();
    fireEvent.click(within(sheet).getByText("Model check"));
    expect(within(sheet).getAllByText("hydraulics.pump_power").length).toBeGreaterThan(0);
  });

  it("lists the table-shaped calculations and says why they cannot be added", async () => {
    stubFetch();
    await open();
    const fitting = screen.getByRole("button", {
      name: /fitting resistance coefficients/i,
    }) as HTMLButtonElement;
    expect(fitting.disabled).toBe(true);
    expect(screen.getAllByText(/table of rows/i).length).toBeGreaterThan(0);
  });

  it("adds a sheet, removes it, and undo restores it", async () => {
    stubFetch();
    await open();
    fireEvent.click(screen.getAllByRole("button", { name: /reynolds number/i })[0]!);
    await waitFor(() => {
      const fault = document.querySelector(".calc-body .fault");
      if (fault !== null) {
        throw new Error(fault.textContent ?? "a fault with no sentence");
      }
      expect(screen.getAllByRole("tab")).toHaveLength(3);
    });

    const sheet = screen.getByRole("tabpanel");
    fireEvent.click(within(sheet).getByRole("button", { name: "Remove" }));
    await waitFor(() => expect(screen.getAllByRole("tab")).toHaveLength(2));

    fireEvent.click(screen.getByRole("button", { name: "Undo" }));
    await waitFor(() => expect(screen.getAllByRole("tab")).toHaveLength(3));
  });

  it("saves the book under its name, and the flowsheet survives a workspace switch", async () => {
    stubFetch();
    const downloads = stubDownload();
    const rendered = await open();
    const before = nodeIds(rendered.container).sort();

    fireEvent.click(screen.getByRole("button", { name: "Save" }));
    await waitFor(() => expect(downloads.names).toHaveLength(1));
    expect(downloads.names[0]).toBe("Worked_examples.toml");
    expect(await textOf(only(downloads.blobs))).toContain('calculation = "hydraulics.pump_power"');

    fireEvent.click(screen.getByRole("button", { name: "Flowsheet" }));
    await waitFor(() => expect(nodeIds(rendered.container).sort()).toEqual(before));
    fireEvent.click(screen.getByRole("button", { name: "Calculations" }));
    expect(screen.getAllByRole("tab")).toHaveLength(2);
  });
});
