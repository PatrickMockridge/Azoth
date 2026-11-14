// @vitest-environment jsdom
/**
 * The Drawing workspace: import a file, read the coverage, and pass the scaffold on.
 *
 * **The door is a fake**, so the test asserts the workspace's own behaviour rather than the
 * module's: the view is the committed fixture, the export is a refusal with a diagnostic, and
 * "Open as flowsheet" is a callback the test sees called with the scaffold's TOML.
 */

import { cleanup, fireEvent, render, screen, waitFor } from "@testing-library/react";
import { afterEach, describe, expect, it, vi } from "vitest";

import { DexpiWorkspace } from "../src/components/dexpi/DexpiWorkspace";
import { decodeDexpiView } from "../src/wire/decode";
import type { ExchangeDoor } from "../src/wire/session";
import type { DexpiExport, DexpiViewResponse } from "../src/wire/types";
import viewJson from "./fixtures/dexpi-view.json";

const view = decodeDexpiView({ ok: true, view: viewJson }).view;
if (view === undefined) {
  throw new Error("the fixture is a view");
}

const REFUSED: DexpiExport = {
  ok: false,
  xml: null,
  diagnostics: [
    {
      severity: "error",
      code: "refused_mapping",
      subject: "f1",
      message: "`unit_ops.flare` has no reviewed DEXPI Process type",
    },
  ],
  steps: 0,
  streams: 0,
  items: 0,
  segments: 0,
  diagram_groups: 0,
  diagram_primitives: 0,
};

function door(overrides: Partial<ExchangeDoor> = {}): ExchangeDoor {
  return {
    inspect: vi.fn(async (): Promise<DexpiViewResponse> => ({ ok: true, view: view! })),
    exportFlowsheet: vi.fn(async (): Promise<DexpiExport> => REFUSED),
    ...overrides,
  };
}

async function importFile(container: HTMLElement): Promise<void> {
  const input = container.querySelector<HTMLInputElement>('input[type="file"]');
  if (input === null) {
    throw new Error("the workspace has a file input");
  }
  // **A file-shaped object, not a `File`.** The workspace reads `name` and `text()`, and a
  // jsdom `File` need not implement the latter; the test supplies exactly what is read.
  const file = {
    name: "drawing.xml",
    text: async () => "<PlantModel/>",
  } as unknown as File;
  fireEvent.change(input, { target: { files: [file] } });
}

afterEach(cleanup);

describe("the Drawing workspace", () => {
  it("shows the drawing, the coverage and the scaffold", async () => {
    const { container } = render(
      <DexpiWorkspace
        exchange={door()}
        document="id = 'flowsheets.x'"
        openable
        hidden={false}
        onOpenFlowsheet={() => true}
        onWorkspace={() => {}}
      />,
    );
    await importFile(container);
    await screen.findAllByText("Modelled");
    expect(screen.getAllByText("Drawing only").length).toBeGreaterThan(0);
    expect(container.querySelector("svg.dexpi-svg")).not.toBeNull();
    expect(screen.getByText(/This drawing would become/)).toBeTruthy();
  });

  it("hands the scaffold to the editor", async () => {
    const open = vi.fn((_toml: string) => true);
    const workspace = vi.fn();
    const { container } = render(
      <DexpiWorkspace
        exchange={door()}
        document="id = 'flowsheets.x'"
        openable
        hidden={false}
        onOpenFlowsheet={open}
        onWorkspace={workspace}
      />,
    );
    await importFile(container);
    await screen.findAllByText("Modelled");
    fireEvent.click(screen.getByRole("button", { name: "Open as flowsheet" }));
    await waitFor(() => expect(open).toHaveBeenCalled());
    expect(open.mock.calls[0]?.[0]).toBe(view.scaffold?.toml);
    expect(workspace).toHaveBeenCalledWith("flowsheet");
  });

  it("renders an export refusal's diagnostics", async () => {
    const { container } = render(
      <DexpiWorkspace
        exchange={door()}
        document="id = 'flowsheets.x'"
        openable
        hidden={false}
        onOpenFlowsheet={() => true}
        onWorkspace={() => {}}
      />,
    );
    await importFile(container);
    await screen.findAllByText("Modelled");
    fireEvent.click(screen.getByRole("button", { name: "Export PFD" }));
    expect(await screen.findByText("refused_mapping")).toBeTruthy();
    expect(screen.getByText(/The export was refused/)).toBeTruthy();
  });
});
