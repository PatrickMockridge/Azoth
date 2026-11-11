// @vitest-environment jsdom
/**
 * The editor, driven the way a person drives it.
 *
 * **Nothing here is mocked.** `fetch` is stubbed to the *real* wasm bytes this repo just built, so
 * the real glue instantiates the real module, the real `Editor` opens the real
 * `specs/flowsheets/examples/two_feed_mixer.toml`, and the assertions are about the DOM the real
 * components render.
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
import { cleanup, fireEvent, render, screen, waitFor, within } from "@testing-library/react";
import { afterEach, beforeAll, describe, expect, it, vi } from "vitest";

import { PFD_EDGE, connectable } from "../src/components/Flowsheet";

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

/**
 * The app, with the module loaded and the shipped document open.
 *
 * The appearance is cleared each time - a choice from an earlier test is a reader's, not this
 * case's starting state - and may be preset by key for the tests about remembering one.
 */
async function open(stored: Record<string, string> = {}) {
  // The first-run card is another test's subject; these tests are about the editor behind it.
  window.localStorage.setItem("azoth.welcome", "seen");
  // Each case starts in SI; a unit choice from an earlier test is a reader's, not this case's.
  window.localStorage.removeItem("azoth.units");
  window.localStorage.removeItem("azoth.appearance");
  window.localStorage.removeItem("azoth.theme");
  for (const [key, value] of Object.entries(stored)) {
    window.localStorage.setItem(key, value);
  }
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
      "input:feed_2",
      "instance:hx1",
      "instance:mix1",
      "instance:p1",
      "instance:sep1",
      "product:vapour_product",
    ]);
    // **The node plate is not the unit label.** An unscoped `.unit` rule put a 168px plate
    // around every field's unit string; the plate is `.unit-plate` and no node carries `.unit`.
    expect(container.querySelectorAll(".react-flow__node .unit-plate").length).toBeGreaterThan(0);
    expect(container.querySelector(".react-flow__node .unit")).toBeNull();
  });

  it("labels every connection point with its direction and stream", async () => {
    stubFetch();
    const { container } = await open();
    const labels = [...container.querySelectorAll(".port-label")].map(
      (label) => label.textContent,
    );
    // A single-inlet machine and a boundary stream do not have to be inferred either.
    expect(labels).toEqual(expect.arrayContaining(["IN", "OUT"]));
    expect(labels).toEqual(expect.arrayContaining(["OUT feed", "IN product"]));

    // **The heat exchanger is the four-handle case.** It is not in the demo, so it is added from
    // the palette; the labels appear before anything is wired.
    fireEvent.click(screen.getByRole("button", { name: "Heat exchanger" }));
    await waitFor(() => {
      const added = [...container.querySelectorAll(".port-label")].map(
        (label) => label.textContent,
      );
      expect(added).toEqual(
        expect.arrayContaining(["IN hot", "IN cold", "OUT hot", "OUT cold"]),
      );
    });
  });

  it("refuses a second line into an occupied connection point", () => {
    const occupied = new Set(["mix1.feed[0]"]);
    expect(connectable(occupied, "feed_1", "mix1.feed[0]")).toBe(false);
    expect(connectable(occupied, "feed_1", "mix1.feed[1]")).toBe(true);
    expect(connectable(occupied, null, "mix1.feed[1]")).toBe(false);
    expect(connectable(occupied, "feed_1", undefined)).toBe(false);
  });

  it("uses the HYSYS connector shape", () => {
    // **The canvas cannot draw a path in jsdom** (it measures handle rectangles), so this is the
    // configuration the real canvas consumes: a smoothstep route, not the default bezier, with a
    // 4 px corner. The visual check is the drawing itself.
    expect(PFD_EDGE).toEqual({ type: "smoothstep", borderRadius: 4 });
  });

  it("opens a node overview on double-click and hands it back to the bar", async () => {
    stubFetch();
    const { container } = await open();
    fireEvent.click(screen.getByRole("button", { name: "Solve" }));
    await waitFor(() => expect(pill(container)).toBe("solved"));

    const heater = container.querySelector<HTMLElement>('[data-id="instance:hx1"]');
    fireEvent.doubleClick(heater as HTMLElement);
    const dialog = await screen.findByRole("dialog");
    expect(within(dialog).getByRole("heading", { name: "Heater" })).toBeDefined();
    expect(within(dialog).getByText(/streams/i)).toBeDefined();
    expect(within(dialog).getByText("hx1.outlet")).toBeDefined();

    fireEvent.keyDown(window, { key: "Escape" });
    await waitFor(() => expect(screen.queryByRole("dialog")).toBeNull());

    fireEvent.doubleClick(heater as HTMLElement);
    const again = await screen.findByRole("dialog");
    fireEvent.click(within(again).getByRole("button", { name: /edit in the property bar/i }));
    await waitFor(() => expect(screen.queryByRole("dialog")).toBeNull());
    expect(container.querySelector(".prop-header")?.textContent).toContain("hx1");
  });

  it("uses a distinct inlet point for every mixer feed", async () => {
    stubFetch();
    const { container } = await open();
    const labels = [...container.querySelectorAll('[data-id="instance:mix1"] .port-label')].map(
      (label) => label.textContent ?? "",
    );
    expect(labels).toEqual(expect.arrayContaining(["IN feed 1", "IN feed 2", "IN feed 3"]));
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
    expect(idPill(container)).toBe("flowsheets.two_feed_mixer");

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

    // A field commits when the typing stops; the value goes out as one `set_parameter`.
    const field = screen.getByDisplayValue("300");
    fireEvent.change(field, { target: { value: "340" } });
    fireEvent.blur(field);

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

    await waitFor(() => expect(nodeIds(container)).toHaveLength(6));
    expect(nodeIds(container)).not.toContain("instance:hx1");
    // The two ports the heater was wired between, named by the checker's own record.
    fireEvent.click(screen.getByRole("tab", { name: "Messages" }));
    await waitFor(() => {
      const codes = [...container.querySelectorAll(".message code")].map(
        (node) => node.textContent ?? "",
      );
      expect(codes.join(" ")).toContain("underfed_port");
    });
  });

  it("adds a palette entry selected, and deletes it with Delete or the Remove button", async () => {
    stubFetch();
    const { container } = await open();

    const heaterEntry = [...container.querySelectorAll<HTMLElement>("aside.side .entry")].find(
      (entry) => entry.textContent?.includes("Heater"),
    );
    expect(heaterEntry).toBeDefined();
    fireEvent.click(heaterEntry as HTMLElement);

    await waitFor(() => expect(nodeIds(container)).toContain("instance:heater_1"));
    // The new machine is selected, so its window is open and the reader can see what happened.
    expect(screen.getByRole("button", { name: "Remove" })).toBeDefined();
    // It has no connections, so the document is incomplete and cannot be run.
    expect(screen.getByRole("button", { name: "Solve" })).toHaveProperty("disabled", true);

    fireEvent.keyDown(document, { key: "Delete" });
    await waitFor(() => expect(nodeIds(container)).not.toContain("instance:heater_1"));
    expect(screen.getByRole("button", { name: "Solve" })).toHaveProperty("disabled", false);

    fireEvent.click(heaterEntry as HTMLElement);
    await waitFor(() => expect(nodeIds(container)).toContain("instance:heater_1"));
    fireEvent.click(screen.getByRole("button", { name: "Remove" }));
    await waitFor(() => expect(nodeIds(container)).not.toContain("instance:heater_1"));
  });

  it("turns a connection into one tear, not a tear plus a connection", async () => {
    stubFetch();
    const { container } = await open();

    const pump = container.querySelector<HTMLElement>('[data-id="instance:p1"]');
    expect(pump).not.toBeNull();
    fireEvent.click(pump as HTMLElement);
    fireEvent.click(screen.getByRole("tab", { name: "Connections" }));
    const outlet = screen
      .getAllByRole("button", { name: "p1.outlet" })
      .find((button) => button.closest(".side.right") !== null);
    expect(outlet).toBeDefined();
    fireEvent.click(outlet as HTMLElement);
    await waitFor(() => expect(screen.getByRole("button", { name: "Make it a tear" })).toBeDefined());

    fireEvent.click(screen.getByRole("button", { name: "Make it a tear" }));
    // **The connection is gone rather than left beside the recycle.** A one inlet would be
    // overfed, and the Messages panel would say so; the conversion leaves it clean.
    fireEvent.click(screen.getByRole("tab", { name: "Messages" }));
    await waitFor(() => expect(screen.getByText("No problems found.")).toBeDefined());
    expect(screen.queryByText(/connections; it takes one/)).toBeNull();
  });

  it("undoes and redoes an edit through a document snapshot", async () => {
    stubFetch();
    const { container } = await open();

    const heaterEntry = [...container.querySelectorAll<HTMLElement>("aside.side .entry")].find(
      (entry) => entry.textContent?.includes("Heater"),
    );
    fireEvent.click(heaterEntry as HTMLElement);
    await waitFor(() => expect(nodeIds(container)).toContain("instance:heater_1"));
    await waitFor(() =>
      expect(screen.getByRole("button", { name: "Undo" })).toHaveProperty("disabled", false),
    );

    fireEvent.click(screen.getByRole("button", { name: "Undo" }));
    await waitFor(() => expect(nodeIds(container)).not.toContain("instance:heater_1"));
    await waitFor(() =>
      expect(screen.getByRole("button", { name: "Redo" })).toHaveProperty("disabled", false),
    );

    fireEvent.click(screen.getByRole("button", { name: "Redo" }));
    await waitFor(() => expect(nodeIds(container)).toContain("instance:heater_1"));
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
    await waitFor(() => expect(download.names).toEqual(["Two_feeds_mix_heat_separate_recycle.toml"]));
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

  it("offers the menu bar over the same commands without duplicating their buttons", async () => {
    stubFetch();
    const { container } = await open();

    // The visible command is still one button; the menu item is a menuitem, not a second button.
    expect(screen.getAllByRole("button", { name: "Solve" })).toHaveLength(1);
    fireEvent.click(screen.getByRole("button", { name: "Run" }));
    expect(screen.getByRole("menuitem", { name: "Solve" })).toBeDefined();
    expect(screen.getAllByRole("button", { name: "Solve" })).toHaveLength(1);

    fireEvent.click(screen.getByRole("menuitem", { name: "Solve" }));
    await waitFor(() => expect(pill(container)).toBe("solved"));
    expect(screen.queryByRole("menu")).toBeNull();
  });

  it("closes a menu on Escape and returns focus to its trigger", async () => {
    stubFetch();
    await open();

    const trigger = screen.getByRole("button", { name: "Case" });
    fireEvent.click(trigger);
    expect(screen.getByRole("menu", { name: "Case" })).toBeDefined();

    fireEvent.keyDown(window, { key: "Escape" });
    await waitFor(() => expect(screen.queryByRole("menu")).toBeNull());
    expect(document.activeElement).toBe(trigger);
  });

  it("opens the Cyberpunk command palette, filters and runs a command", async () => {
    stubFetch();
    const { container } = await open({
      "azoth.appearance": JSON.stringify({ layout: "cyberpunk", theme: "dark" }),
    });
    expect(document.documentElement.dataset.layout).toBe("cyberpunk");

    fireEvent.keyDown(window, { key: "k", ctrlKey: true });
    const dialog = await screen.findByRole("dialog", { name: "Command palette" });
    const search = within(dialog).getByRole("combobox", { name: "Search commands" });

    fireEvent.change(search, { target: { value: "solve" } });
    fireEvent.keyDown(search, { key: "Enter" });

    await waitFor(() => expect(pill(container)).toBe("solved"));
    expect(screen.queryByRole("dialog", { name: "Command palette" })).toBeNull();
  });

  it("leaves Ctrl+K to the browser outside Cyberpunk", async () => {
    stubFetch();
    await open();
    fireEvent.keyDown(window, { key: "k", ctrlKey: true });
    expect(screen.queryByRole("dialog", { name: "Command palette" })).toBeNull();
  });

  it("hides and shows the Cyberpunk HUD rails and opens the palette from the strip", async () => {
    stubFetch();
    const { container } = await open({
      "azoth.appearance": JSON.stringify({ layout: "cyberpunk", theme: "dark" }),
    });
    const app = container.querySelector(".app");
    expect(app?.getAttribute("data-hud-palette")).toBe("on");

    fireEvent.click(screen.getByRole("button", { name: "Hide Navigator" }));
    expect(app?.getAttribute("data-hud-palette")).toBe("off");
    fireEvent.click(screen.getByRole("button", { name: "Show Navigator" }));
    expect(app?.getAttribute("data-hud-palette")).toBe("on");

    fireEvent.click(screen.getByRole("button", { name: "Hide Workbook" }));
    expect(app?.getAttribute("data-hud-workbook")).toBe("off");

    fireEvent.click(screen.getByRole("button", { name: "Open the command palette" }));
    expect(screen.getByRole("dialog", { name: "Command palette" })).toBeDefined();
  });

  it("opens in Windows XP light, and remembers the layout and colour chosen", async () => {
    stubFetch();
    await open();

    // **The attributes are overrides and not what makes the default.** `:root` carries the XP
    // light tokens, so a document with no attributes is already the default - which is what makes
    // it a default rather than a value somebody set.
    expect(document.documentElement.dataset.layout).toBe("xp");
    expect(document.documentElement.dataset.theme).toBe("light");

    fireEvent.click(screen.getByRole("button", { name: "Appearance" }));
    fireEvent.click(screen.getByRole("button", { name: "Macintosh" }));
    await waitFor(() => expect(document.documentElement.dataset.layout).toBe("mac"));

    // The panel stays open so the two axes can be compared before it closes.
    fireEvent.click(screen.getByRole("button", { name: "Dark" }));
    await waitFor(() => expect(document.documentElement.dataset.theme).toBe("dark"));
    expect(
      JSON.parse(window.localStorage.getItem("azoth.appearance") ?? "{}"),
    ).toEqual({ layout: "mac", theme: "dark" });

    // Escape closes it and hands focus back to the trigger, which is where a keyboard user was.
    fireEvent.keyDown(window, { key: "Escape" });
    await waitFor(() => expect(screen.queryByRole("dialog", { name: "Appearance" })).toBeNull());
    expect(document.activeElement).toBe(screen.getByRole("button", { name: "Appearance" }));
  });

  it("adopts a remembered appearance, and migrates the old theme key", async () => {
    stubFetch();
    await open({
      "azoth.appearance": JSON.stringify({ layout: "cyberpunk", theme: "neutral" }),
    });
    await waitFor(() =>
      expect(document.documentElement.dataset.layout).toBe("cyberpunk"),
    );
    expect(document.documentElement.dataset.theme).toBe("neutral");

    cleanup();
    // A build before the layout axis remembered only a colour; the layout is the one that build
    // most resembled.
    await open({ "azoth.theme": "dark" });
    await waitFor(() => expect(document.documentElement.dataset.theme).toBe("dark"));
    expect(document.documentElement.dataset.layout).toBe("xp");
  });

  it("shows the duty a heater reached, under the field name its own model gives it", async () => {
    stubFetch();
    const { container } = await open();
    fireEvent.click(screen.getByRole("button", { name: "Solve" }));
    await waitFor(() => expect(pill(container)).toBe("solved"));

    fireEvent.click(container.querySelector<HTMLElement>('[data-id="instance:hx1"]') as HTMLElement);
    await waitFor(() => expect(screen.getAllByText("outlet_temperature").length).toBeGreaterThan(0));

    // **The sheet is offered because the run published a result**, and it says the field by the
    // name the operation's model declares - which is the name a case file and the Python half use.
    // The number is the library's, in the unit the model named: this is the whole point of carrying
    // a result on the wire, since the duty is on no outlet stream.
    fireEvent.click(screen.getByRole("tab", { name: "Results" }));
    const rows = [...container.querySelectorAll("#sheet-results tbody tr")].map(
      (row) => row.textContent ?? "",
    );
    const duty = rows.find((text) => text.startsWith("outlet_duty"));
    expect(duty).toBeDefined();
    expect(duty).toContain("W");
    // And the heater here is *cooling* - the pump's outlet is above the 320 K the document states -
    // so the duty is negative, which is the model's own sign and not this front end's.
    expect(duty).toContain("-");
  });

  it("lists every stream the document wires, and Solve changes the values and not the rows", async () => {
    stubFetch();
    const { container } = await open();

    // The workbook is not the tab the dock opens on - messages are, because a document that does
    // not check is the state a person most needs told about.
    fireEvent.click(screen.getByRole("tab", { name: "Workbook" }));
    const rowIds = () =>
      [...container.querySelectorAll(".dock .grid tbody tr")].map((row) =>
        row.getAttribute("data-id"),
      );

    // **The rows are the *producing* endpoints**, so no stream is listed twice: `mix1.product` and
    // `p1.inlet` are one stream, and only the producer carries its value.
    const before = rowIds();
    expect([...before].sort()).toEqual(
      [
        "feed_1",
        "feed_2",
        "hx1.outlet",
        "mix1.product",
        "p1.outlet",
        "sep1.liquid",
        "sep1.vapour",
        "vapour_product",
      ].sort(),
    );

    // And a run fills the cells rather than adding rows, which is what makes it a worksheet.
    fireEvent.click(screen.getByRole("button", { name: "Solve" }));
    await waitFor(() => expect(pill(container)).toBe("solved"));
    expect(rowIds()).toEqual(before);
  });

  it("groups the palette by the family the library filed each entry under", async () => {
    stubFetch();
    const { container } = await open();

    const headings = [...container.querySelectorAll("aside.side .group > .label")].map(
      (label) => label.textContent,
    );
    // **Eight, and not the one group called `other` this drew until the wire carried the family.**
    // The fault was here and the missing fact was in the wire: every id is `unit_ops.<leaf>`, so a
    // split on `.` took the leaf for a family.
    expect(headings).toEqual([
      // The boundary is first now: a blank case starts with a feed, and it was below the whole
      // palette before.
      "Boundary",
      "column",
      "heat exchanger",
      "machine",
      "mixer",
      "reactor",
      "separator",
      "two port",
      "utility",
    ]);
    // And every entry is in one of them, so a group cannot be a heading with nothing under it.
    const grouped = container.querySelectorAll("aside.side .group .entry");
    expect(grouped.length).toBeGreaterThanOrEqual(29);
    // A refused entry says why in the page, not only in a disabled button's title.
    expect(screen.getByText(/fixed-point loop over MDEA/)).toBeDefined();
  });

  /**
   * **A field reads and writes in the set, and the document keeps its own unit.** This is the half
   * the display-only switcher could not do, and the failure it guards against is invisible in the
   * editor: a field that showed a temperature in °F and put the typed number straight into the
   * document would show the reader exactly what they typed while the run used a state 273 K colder.
   * So the assertion is not on the box - it is on what the document holds once the box is out of
   * the way, which is why the set is switched back to SI to read it.
   */
  it("accepts a value in the set's unit and stores the document's", async () => {
    stubFetch();
    window.localStorage.clear();
    const { container } = await open();
    fireEvent.click(screen.getByRole("button", { name: "Solve" }));
    await waitFor(() => expect(pill(container)).toBe("solved"));

    const chooser = () => container.querySelector<HTMLSelectElement>("select.units");
    fireEvent.change(chooser() as HTMLSelectElement, { target: { value: "field" } });
    expect(container.querySelector(".status-bar")?.textContent).toContain("Field");
    fireEvent.click(container.querySelector<HTMLElement>('[data-id="instance:hx1"]') as HTMLElement);
    await waitFor(() => expect(screen.getAllByText("outlet_temperature").length).toBeGreaterThan(0));

    // The field says which unit it accepts, and it is the set's rather than the spec's: the heater
    // declares `outlet_temperature` in K and a field set reads a temperature in °F.
    const head = screen.getByText("outlet_temperature").parentElement;
    expect(head?.textContent).toContain("°F");
    // The document's 300 K, shown as a reader would say it: 80.33 °F.
    const field = screen.getByDisplayValue("80.33");

    fireEvent.change(field, { target: { value: "200" } });
    fireEvent.blur(field);
    await waitFor(() => expect(pill(container)).toBe("stale"));

    // **Back to SI to read what was stored.** 200 °F is 366.483 K, so a document that received the
    // typed number would show 200 here.
    fireEvent.change(chooser() as HTMLSelectElement, { target: { value: "si" } });
    // The box shows the number the document holds, in its own unit - which is the full float,
    // because SI is the identity conversion and an identity is not something to round.
    await waitFor(() => {
      const stored = screen.getByDisplayValue(/^366\.48/) as HTMLInputElement;
      expect(Number(stored.value)).toBeCloseTo(366.4833, 3);
    });

    // And the run uses it: the outlet reaches the temperature the document holds.
    fireEvent.click(screen.getByRole("button", { name: "Solve" }));
    await waitFor(() => expect(pill(container)).toBe("solved"));
    await waitFor(() => {
      const readouts = [
        ...container.querySelectorAll('[data-id="instance:hx1"] .readout'),
      ].map((node) => node.textContent ?? "");
      expect(readouts.join(" | ")).toContain("366.5");
    });
  });

  /**
   * **The navigator and the canvas are one selection**, because both call the app's own
   * `setSelected`. So the proof is the object's window opening — which nothing but the selection
   * can produce — and then the same thing in the other direction: a click on the drawing moves the
   * row rather than lighting a second one.
   */
  it("selects from the navigator what the canvas selects, in both directions", async () => {
    stubFetch();
    const { container } = await open();

    const rows = () => [...container.querySelectorAll<HTMLElement>(".tree .tree-row")];
    const rowFor = (label: string) =>
      rows().find((row) => row.querySelector(".name")?.textContent === label);
    const opened = () => container.querySelector(".prop-header .note")?.textContent;

    // One row per object: the four unit operations, the six streams the edges carry, and the two
    // boundary names — with the recycle among the streams, said as a tear.
    expect(rows()).toHaveLength(14);
    expect(rowFor("recycle_1")?.querySelector(".tag")?.textContent).toBe("recycle");
    expect(rowFor("mix1.product")?.querySelector(".tag")).toBeNull();

    fireEvent.click(rowFor("hx1") as HTMLElement);
    await waitFor(() => expect(opened()).toBe("hx1"));
    expect(rowFor("hx1")?.getAttribute("aria-current")).toBe("true");

    // The drawing is the same selection, so this click moves the mark instead of adding one.
    fireEvent.click(
      container.querySelector<HTMLElement>('[data-id="instance:sep1"]') as HTMLElement,
    );
    await waitFor(() => expect(opened()).toBe("sep1"));
    expect(rowFor("sep1")?.getAttribute("aria-current")).toBe("true");
    expect(rowFor("hx1")?.getAttribute("aria-current")).toBeNull();
  });

  /**
   * **The unit set is a reading of the run and not a second run.** The sets come from the library's
   * catalogue, the factors from the same conversions a calculation runs, and switching one converts
   * what is drawn and nothing else — so a column's header and its cells move together, and a
   * dimension no set names does not move at all.
   */
  it("reads the workbook in the set the reader chose, and remembers it", async () => {
    stubFetch();
    window.localStorage.clear();
    const { container } = await open();
    fireEvent.click(screen.getByRole("button", { name: "Solve" }));
    await waitFor(() => expect(pill(container)).toBe("solved"));
    fireEvent.click(screen.getByRole("tab", { name: "Workbook" }));

    const chooser = () => container.querySelector<HTMLSelectElement>("select.units");
    const headers = () =>
      [...container.querySelectorAll(".dock .grid thead th")].map((th) => th.textContent);
    // A quantity cell reads as a number *and* its unit - `formatQuantity`'s own shape - so a
    // value is parsed off the front of it rather than the cell being read as a number.
    const cell = (row: string, index: number) =>
      Number.parseFloat(
        container.querySelectorAll(`.dock .grid tbody tr[data-id="${row}"] td`)[index]
          ?.textContent ?? "",
      );

    // The sets are the library's, in its own order, and SI is the one a document opens in.
    expect([...(chooser()?.options ?? [])].map((option) => option.textContent)).toEqual([
      "SI",
      "Field",
    ]);
    expect(chooser()?.value).toBe("si");
    expect(headers()).toEqual([
      "flow mol/s",
      "mass flow kg/s",
      "M kg/mol",
      "P Pa",
      "T K",
      "h J/mol",
      "VF",
    ]);
    const pascals = cell("feed_1", 3);

    fireEvent.change(chooser() as HTMLSelectElement, { target: { value: "field" } });
    await waitFor(() => expect(headers()[3]).toBe("P psi"));

    // Five of the seven move, and the two that do not are the ones no set names: a molar mass in
    // `lb/lbmol` would need a pound-mole, which the vocabulary does not carry, and a vapour
    // fraction has no unit at all.
    expect(headers()).toEqual([
      "flow kmol/h",
      "mass flow lb/h",
      "M kg/mol",
      "P psi",
      "T °F",
      "h kJ/mol",
      "VF",
    ]);
    // And the number is the library's, in the unit the header names: one psi is 6894.757… Pa, so
    // the same pressure read twice agrees to the four significant figures the column shows.
    const pounds = cell("feed_1", 3);
    expect(pounds * 6894.757293168361).toBeGreaterThan(pascals * 0.99);
    expect(pounds * 6894.757293168361).toBeLessThan(pascals * 1.01);

    // **And the temperature is a temperature, not a relabeled kelvin.** The demo's mixer outlet is
    // at 300 K, which is 80.33 °F - so the affine inverse is `si / factor - offset` and not a
    // division, and a front end that dropped the offset would show 540 here. This is the one cell
    // in the document that can tell the two apart.
    const temperature = cell("feed_1", 4);
    expect(temperature).toBeCloseTo(80.33, 1);
    // And the definition inverts it: `si = (value + offset) * factor`, which is uom's own form -
    // `(80.33 + 459.67) * 5/9` is the 300 K the record holds.
    expect((temperature + 459.67) * 0.5555555555555556).toBeCloseTo(300, 0);

    // A reader's choice, remembered the way the theme is - and not on the wire, because it is not
    // a fact about the document.
    expect(window.localStorage.getItem("azoth.units")).toBe("field");
    // The run did not change: only what is drawn did. The workbook's rows are the same rows.
    expect(container.querySelectorAll(".dock .grid tbody tr")).toHaveLength(8);
  });

  it("selects a feed instead of deleting it, and edits it through the feed form", async () => {
    stubFetch();
    const { container } = await open();

    const row = container.querySelector<HTMLElement>(".side .boundary-row .entry");
    expect(row).not.toBeNull();
    fireEvent.click(row as HTMLElement);

    // The click selected the feed, which is what the old row did not do: it removed the feed and
    // every connection that drew from it.
    expect(nodeIds(container)).toContain("input:feed_1");
    await waitFor(() => expect(screen.getByRole("button", { name: "Edit feed" })).toBeDefined());

    fireEvent.click(screen.getByRole("button", { name: "Edit feed" }));
    const conditions = container.querySelectorAll<HTMLInputElement>(
      ".feed-editor .conditions input",
    );
    expect(conditions).toHaveLength(3);
    fireEvent.change(conditions[2] as HTMLInputElement, { target: { value: "340" } });
    fireEvent.click(screen.getByRole("button", { name: "Apply changes" }));

    await waitFor(() => expect(pill(container)).toBe("stale"));
    expect(nodeIds(container)).toContain("input:feed_1");
  });

  it("names the case, saves by that name, and marks it unsaved until it is saved", async () => {
    stubFetch();
    const download = stubDownload();
    const { container } = await open();

    const name = screen.getByLabelText("Case name");
    fireEvent.change(name, { target: { value: "My case" } });
    fireEvent.blur(name);
    await waitFor(() => expect(screen.getByDisplayValue("My case")).toBeDefined());
    expect(screen.getByText("unsaved")).toBeDefined();

    fireEvent.click(screen.getByRole("button", { name: "Save" }));
    await waitFor(() => expect(download.names).toEqual(["My_case.toml"]));
    await waitFor(() => expect(screen.queryByText("unsaved")).toBeNull());

    // An edit is unsaved again; New asks first, and a no leaves the document in place.
    fireEvent.click(screen.getByRole("button", { name: "Solve" }));
    await waitFor(() => expect(pill(container)).toBe("solved"));
    fireEvent.click(container.querySelector<HTMLElement>('[data-id="instance:hx1"]') as HTMLElement);
    const field = await screen.findByDisplayValue("300");
    fireEvent.change(field, { target: { value: "330" } });
    fireEvent.blur(field);
    await waitFor(() => expect(screen.getByText("unsaved")).toBeDefined());

    Object.defineProperty(window, "confirm", { value: () => false, configurable: true });
    fireEvent.click(screen.getByRole("button", { name: "New" }));
    expect(idPill(container)).toBe("flowsheets.two_feed_mixer");
  });

  it("prints a solved case and exports the workbook", async () => {
    stubFetch();
    const download = stubDownload();
    // A previous case may have left the field set remembered; this export is read in SI.
    window.localStorage.removeItem("azoth.units");
    const { container } = await open();
    fireEvent.click(screen.getByRole("button", { name: "Solve" }));
    await waitFor(() => expect(pill(container)).toBe("solved"));

    fireEvent.click(screen.getByRole("button", { name: "Report" }));
    expect(screen.getByRole("dialog", { name: "Case report" })).toBeDefined();
    fireEvent.click(screen.getByRole("button", { name: "Close" }));

    fireEvent.click(screen.getByRole("tab", { name: "Workbook" }));
    fireEvent.click(screen.getByRole("button", { name: "Export CSV" }));
    await waitFor(() => expect(download.names).toHaveLength(1));
    expect(download.names[0]).toContain(".csv");
    const text = await textOf(only(download.blobs));
    expect(text).toContain("feed_1");
    expect(text).toContain("T [K]");
  });
});
