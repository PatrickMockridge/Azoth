// @vitest-environment jsdom
/**
 * The Macintosh desktop's windows: the regions are titled, can be moved and raised, and the
 * reset command drops the inline geometry.
 *
 * **The geometry is mocked.** jsdom does not lay out boxes, so `getBoundingClientRect` is what
 * this test supplies; the code under test is the drag/resize arithmetic and the variable names,
 * not the browser's measurement. The real screenshot is the visual check.
 */

import { readFileSync } from "node:fs";
import { resolve } from "node:path";

import { cleanup, fireEvent, render, screen, waitFor } from "@testing-library/react";
import { afterEach, beforeAll, expect, it, vi } from "vitest";

import { stubXyflowEnvironment } from "./xyflow-env";

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

/** A DOMRect the window manager can read, because jsdom lays out nothing. */
function rect(x: number, y: number, width: number, height: number): DOMRect {
  return {
    x,
    y,
    width,
    height,
    top: y,
    left: x,
    right: x + width,
    bottom: y + height,
    toJSON: () => ({}),
  } as DOMRect;
}

async function openLayout(layout: "xp" | "mac" | "cyberpunk"): Promise<HTMLElement> {
  window.localStorage.setItem("azoth.welcome", "seen");
  window.localStorage.setItem(
    "azoth.appearance",
    JSON.stringify({ layout, theme: "light" }),
  );
  const { App } = await import("../src/App");
  const rendered = render(<App />);
  await waitFor(() => {
    expect(rendered.container.querySelectorAll(".react-flow__node")).not.toHaveLength(0);
  });
  return rendered.container;
}

const openMac = () => openLayout("mac");

const windowOf = (container: HTMLElement, id: string): HTMLElement => {
  const element = container.querySelector<HTMLElement>(`[data-window="${id}"]`);
  if (element === null) {
    throw new Error(`no ${id} window`);
  }
  return element;
};

afterEach(() => {
  cleanup();
  vi.restoreAllMocks();
});

it("titles the regions as separate windows and hides the dock splitter", async () => {
  stubFetch();
  const container = await openMac();

  expect(document.documentElement.dataset.layout).toBe("mac");
  const titles = ["palette", "document", "inspector", "workbook", "status", "commands"];
  for (const id of titles) {
    const title = windowOf(container, id).querySelector(".window-title");
    expect(title?.textContent?.trim(), id).not.toBe("");
  }
  expect(windowOf(container, "palette").querySelector(".window-title")?.textContent).toBe(
    "Navigator",
  );
  expect(windowOf(container, "workbook").querySelector(".window-title")?.textContent).toBe(
    "Workbook",
  );
  // The dock's splitter is the XP/Cyberpunk resizer; the Mac grip replaces it.
  expect(container.querySelector(".dock .splitter")).toBeNull();
});

it("drags from the title, resizes from the grip, raises, and resets", async () => {
  stubFetch();
  const container = await openMac();
  const app = container.querySelector<HTMLElement>(".app");
  const palette = windowOf(container, "palette");
  const documentWindow = windowOf(container, "document");
  if (app === null) {
    throw new Error("no app");
  }

  vi.spyOn(app, "getBoundingClientRect").mockReturnValue(rect(0, 0, 1200, 800));
  vi.spyOn(palette, "getBoundingClientRect").mockReturnValue(rect(100, 100, 300, 300));
  vi.spyOn(documentWindow, "getBoundingClientRect").mockReturnValue(rect(400, 100, 600, 400));

  const title = palette.querySelector(".window-title");
  if (title === null) {
    throw new Error("no title");
  }
  fireEvent.mouseDown(title, { button: 0, clientX: 120, clientY: 120 });
  fireEvent.mouseMove(window, { clientX: 170, clientY: 170 });
  fireEvent.mouseUp(window);
  expect(palette.style.getPropertyValue("--palette-x")).toBe("150px");
  expect(palette.style.getPropertyValue("--palette-y")).toBe("150px");

  const grip = palette.querySelector(".window-grip");
  if (grip === null) {
    throw new Error("no grip");
  }
  fireEvent.mouseDown(grip, { button: 0, clientX: 400, clientY: 400 });
  fireEvent.mouseMove(window, { clientX: 450, clientY: 450 });
  fireEvent.mouseUp(window);
  expect(palette.style.getPropertyValue("--palette-w")).toBe("350px");
  expect(palette.style.getPropertyValue("--palette-h")).toBe("350px");

  // Clicking another window raises it above the palette.
  fireEvent.mouseDown(documentWindow.querySelector(".window-title") as Element, {
    button: 0,
    clientX: 500,
    clientY: 110,
  });
  expect(documentWindow.style.getPropertyValue("--wz")).not.toBe("");
  expect(Number(documentWindow.style.getPropertyValue("--wz"))).toBeGreaterThan(
    Number(palette.style.getPropertyValue("--wz") || "0"),
  );

  fireEvent.click(screen.getByRole("button", { name: "View" }));
  fireEvent.click(screen.getByRole("menuitem", { name: "Reset window positions" }));
  expect(palette.style.getPropertyValue("--palette-x")).toBe("");
  expect(palette.style.getPropertyValue("--palette-w")).toBe("");
});

it("gives the Cyberpunk workbook a real resizer", async () => {
  stubFetch();
  const container = await openLayout("cyberpunk");
  const dock = container.querySelector<HTMLElement>(".dock");
  const splitter = dock?.querySelector<HTMLElement>(".splitter");
  if (dock === null || dock === undefined || splitter === null || splitter === undefined) {
    throw new Error("the Cyberpunk workbook has no splitter");
  }
  // The one control the HUD used to lose: its height is the reader's, by pointer or keyboard.
  expect(splitter.getAttribute("role")).toBe("separator");
  expect(dock.style.height).toBe("200px");
  fireEvent.keyDown(splitter, { key: "ArrowUp" });
  expect(dock.style.height).toBe("220px");
  // And the palette's pointer entry says what it opens.
  expect(
    screen
      .getAllByRole("button", { name: "Open the command palette" })
      .some((button) => button.textContent === "Commands"),
  ).toBe(true);
});

it("keeps the workspace switch reachable in the Macintosh and Cyberpunk shells", async () => {
  for (const layout of ["mac", "cyberpunk"] as const) {
    stubFetch();
    const container = await openLayout(layout);
    const flowsheet = container.querySelector<HTMLElement>(".app:not(.calculations)");
    if (flowsheet === null) {
      throw new Error("the flowsheet should be mounted");
    }
    const toCalculations = [...flowsheet.querySelectorAll<HTMLButtonElement>(
      ".workspace-switch button",
    )].find((button) => button.textContent === "Calculations");
    if (toCalculations === undefined) {
      throw new Error("the workspace switch has no Calculations button");
    }
    // The switch is the visible path; the View menu is the same commands, not a substitute.
    expect(toCalculations.textContent).toBe("Calculations");
    fireEvent.click(toCalculations);
    const calculations = container.querySelector<HTMLElement>(".app.calculations");
    if (calculations === null) {
      throw new Error("the calculation workspace should open");
    }
    expect(calculations.hidden).toBe(false);
    expect(flowsheet.hidden).toBe(true);

    const back = calculations.querySelectorAll<HTMLButtonElement>(
      ".workspace-switch button",
    )[0];
    if (back === undefined) {
      throw new Error("the switch has no Flowsheet button");
    }
    fireEvent.click(back);
    expect(flowsheet.hidden).toBe(false);
    expect(calculations.hidden).toBe(true);
    cleanup();
  }
});
