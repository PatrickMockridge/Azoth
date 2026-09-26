// @vitest-environment jsdom
/**
 * The editor, driven the way a person drives it.
 *
 * **Nothing here is mocked.** `fetch` is stubbed to the *real* wasm bytes this repo just built, so
 * the real glue instantiates the real module, the real `Editor` opens the real
 * `specs/flowsheets/demo.toml`, and the assertions are about the DOM the real components render.
 * What that buys is the layer neither of the other two suites can reach: `crates/azoth-wasm/test/
 * smoke.mjs` proves the JS boundary and the three `test/*.test.ts` files prove pure functions over
 * JSON, and neither can tell whether a *gesture* emits the command the UI claims.
 *
 * The environment is jsdom, which is missing three things xyflow calls. Each polyfill names the
 * library that needs it, so the next reader does not have to rediscover it — and a fourth would
 * fail as a `ReferenceError` rather than silently.
 */

import { readFileSync } from "node:fs";
import { resolve } from "node:path";
import { cleanup, fireEvent, render, screen, waitFor } from "@testing-library/react";
import { afterEach, beforeAll, describe, expect, it, vi } from "vitest";

import { only, stubDownload, textOf, unstubDownload } from "./download";
import { drag, drawnPosition } from "./drag";
import { idPill, nodeIds, pill, stubXyflowEnvironment } from "./xyflow-env";

beforeAll(stubXyflowEnvironment);

/** The module the build just produced, served the way a static host serves it. */
function stubFetch(): void {
  // From the project root, because `vitest` runs with `ui/` as its cwd — and because
  // `import.meta.url` under the jsdom environment is not a `file:` URL, so resolving beside this
  // file is not available here.
  const wasm = readFileSync(resolve("src/wasm/pkg/azoth_wasm_bg.wasm"));
  globalThis.fetch = (input: RequestInfo | URL): Promise<Response> => {
    // The glue must ask for the module beside itself. A stub that answered anything would hide a
    // broken asset URL, which is the one failure this test is in a position to notice.
    expect(String(input)).toContain("azoth_wasm_bg.wasm");
    return Promise.resolve(
      new Response(wasm, { headers: { "Content-Type": "application/wasm" } }),
    );
  };
}

/** The app, with the module loaded and the shipped document open. */
async function open() {
  const { App } = await import("../src/App");
  const rendered = render(<App />);
  // The demo's first node only exists once the module has loaded and the document was opened.
  await waitFor(() => {
    expect(rendered.container.querySelectorAll(".react-flow__node")).not.toHaveLength(0);
  });
  return rendered;
}

afterEach(() => {
  // **Not automatic here.** React Testing Library unmounts between tests only when `afterEach`
  // is a *global*, and this suite imports it from `vitest` with `globals` off — so without this
  // line every test's tree stays in the document and a query finds the previous test's buttons.
  cleanup();
  vi.unstubAllGlobals();
  unstubDownload();
});

describe("the editor", () => {
  it("renders the palette the module carries, and draws the shipped document", async () => {
    stubFetch();
    const { container } = await open();

    // The palette: 29 entries of which 3 are refused by the executor, read off the DOM rather
    // than off the catalogue - which is what the smoke driver already checked.
    const entries = container.querySelectorAll("aside.side .entry");
    expect(entries.length).toBeGreaterThanOrEqual(29);
    expect(screen.getByRole("button", { name: /pump/i })).toBeDefined();

    // And the document, with the ids the projection gives it.
    expect(nodeIds(container).sort()).toEqual([
      "input:feed_1",
      "instance:hx1",
      "instance:mix1",
      "instance:p1",
      "instance:sep1",
      "product:vapour_product",
    ]);
  });

  it("starts stale, and Solve makes it solved", async () => {
    stubFetch();
    const { container } = await open();
    expect(pill(container)).toBe("stale");

    fireEvent.click(screen.getByRole("button", { name: "Solve" }));
    await waitFor(() => expect(pill(container)).toBe("solved"));
  });

  /**
   * **New opens a document, and what it leaves out is the checker's to say.** The blank flowsheet is
   * `specs/flowsheets/blank.toml`, the two keys a flowsheet cannot omit; everything a flow sheet
   * still needs reaches the editor as a diagnostic rather than as a page that refused to open.
   */
  it("opens a blank document rather than a fault", async () => {
    stubFetch();
    const { container } = await open();
    expect(idPill(container)).toBe("flowsheets.demo");

    fireEvent.click(screen.getByRole("button", { name: "New" }));

    await waitFor(() => expect(nodeIds(container)).toEqual([]));
    // A document, not a refusal: the id is the blank one's and there is no fault to show.
    expect(idPill(container)).toBe("flowsheets.blank");
    expect(container.querySelector(".fault")).toBeNull();
    // It has not run, so its values are older than it.
    expect(pill(container)).toBe("stale");
    // And the palette is the module's, which no document changes.
    expect(container.querySelectorAll("aside.side .entry").length).toBeGreaterThanOrEqual(29);
  });

  /** A gesture, the command it emits, and the value that can only come back from the library. */
  it("edits a parameter and the canvas shows what the run reached", async () => {
    stubFetch();
    const { container } = await open();
    fireEvent.click(screen.getByRole("button", { name: "Solve" }));
    await waitFor(() => expect(pill(container)).toBe("solved"));

    // Select the heater, and read the value it has now.
    const heater = container.querySelector<HTMLElement>('[data-id="instance:hx1"]');
    expect(heater).not.toBeNull();
    fireEvent.click(heater as HTMLElement);
    await waitFor(() => expect(screen.getAllByText("outlet_temperature").length).toBeGreaterThan(0));

    // A field's `onChange` is what the app listens to; the value goes out as `set_parameter`.
    const field = screen.getByDisplayValue("320");
    fireEvent.change(field, { target: { value: "340" } });

    // The edit landed, so the values are older than the document...
    await waitFor(() => expect(pill(container)).toBe("stale"));
    // ...and solving again recomputes them: 340 K can only appear if the command reached the
    // document, the library ran it and the envelope came back.
    fireEvent.click(screen.getByRole("button", { name: "Solve" }));
    await waitFor(() => expect(pill(container)).toBe("solved"));
    // Scoped to the heater's own node: the first `.readout` in the document belongs to the
    // mixer, and asserting on it would pass or fail for a reason the test is not about.
    await waitFor(() => {
      const readouts = [
        ...container.querySelectorAll('[data-id="instance:hx1"] .readout'),
      ].map((node) => node.textContent ?? "");
      expect(readouts.join(" | ")).toContain("340 K");
    });
  });

  /**
   * **The regression this part exists for.** xyflow's Backspace goes through `onNodesChange`,
   * which updates local state and sends no command — and the next envelope re-derives the nodes,
   * so the node came back. A delete that silently undoes itself is the defect; a node that stays
   * gone, with its two endpoints reported unfed, is the fix.
   */
  it("deletes a node with a command rather than silently undoing the gesture", async () => {
    stubFetch();
    const { container } = await open();

    const heater = container.querySelector<HTMLElement>('[data-id="instance:hx1"]');
    fireEvent.click(heater as HTMLElement);
    fireEvent.keyDown(document, { key: "Backspace" });

    await waitFor(() => expect(nodeIds(container)).toHaveLength(5));
    expect(nodeIds(container)).not.toContain("instance:hx1");
    // The two ports the heater was wired between, named by the checker's own record.
    await waitFor(() => {
      const codes = [...container.querySelectorAll(".diagnostic .code")].map(
        (node) => node.textContent,
      );
      expect(codes).toContain("underfed_port");
    });
  });

  /**
   * **The one wire command no other gesture test reaches.** `set_position` goes out of
   * `onNodeDragStop`, and a drag is the one gesture a synthetic event cannot produce cheaply:
   * xyflow drags through `d3-drag`, which listens for *mouse* events and attaches its move handler
   * to `event.view`. `./drag.ts` states what that costs and how it is paid.
   *
   * **Three things have to hold, and each fails on its own.** The gesture has to emit the command
   * (the pill goes stale, which no local state can do); the canvas has to move (the drawn
   * position); and the *document* has to carry what the canvas draws — a position xyflow moved in
   * its own state and never sent would satisfy the first two.
   */
  it("moves a node with a command, and the document carries the position", async () => {
    stubFetch();
    const { container } = await open();
    fireEvent.click(screen.getByRole("button", { name: "Solve" }));
    await waitFor(() => expect(pill(container)).toBe("solved"));

    const before = drawnPosition(container, "instance:hx1");
    const heater = container.querySelector<HTMLElement>('[data-id="instance:hx1"]');
    drag(heater as HTMLElement, [200, 200], [260, 240]);

    await waitFor(() => expect(pill(container)).toBe("stale"));
    const after = drawnPosition(container, "instance:hx1");
    expect(after).not.toEqual(before);

    // **The document, which is only reachable through a save.** Its layout table is what the
    // library holds, and it has to be the numbers the canvas is drawing.
    const download = stubDownload();
    fireEvent.click(screen.getByRole("button", { name: "Save" }));
    await waitFor(() => expect(download.names).toEqual(["flowsheets.demo.toml"]));
    expect(await textOf(only(download.blobs))).toContain(`hx1 = [${after[0]}.0, ${after[1]}.0]`);
  });

  it("sets the order the next run takes, from the session rather than a copy", async () => {
    stubFetch();
    const { container } = await open();
    const order = container.querySelector<HTMLSelectElement>("select.order");
    expect(order?.value).toBe("insertion");

    fireEvent.change(order as HTMLSelectElement, { target: { value: "topological" } });
    await waitFor(() =>
      expect(container.querySelector<HTMLSelectElement>("select.order")?.value).toBe(
        "topological",
      ),
    );
    // An order is not a change to the document, so the values are still current.
    expect(pill(container)).toBe("stale");
  });
});
