// @vitest-environment jsdom
/**
 * The messages panel: one row per object, plain sentences, and counts that match what is on screen.
 */

import { readFileSync } from "node:fs";
import { resolve } from "node:path";
import { cleanup, fireEvent, render, screen } from "@testing-library/react";
import { afterEach, describe, expect, it, vi } from "vitest";

import { DiagnosticsPanel } from "../src/components/DiagnosticsPanel";
import { groupDiagnostics } from "../src/state/diagnostics";
import { statusWord } from "../src/state/status";
import { parseEnvelope } from "../src/wire/decode";
import type { Diagnostic, Envelope } from "../src/wire/types";

const solved = parseEnvelope(readFileSync(resolve("test/fixtures/envelope.json"), "utf8"));
const broken = parseEnvelope(readFileSync(resolve("test/fixtures/broken.json"), "utf8"));

afterEach(cleanup);

function withDiagnostic(envelope: Envelope, diagnostic: Diagnostic): Envelope {
  return { ...envelope, ok: false, diagnostics: [diagnostic] };
}

const warning: Diagnostic = {
  code: "unused_feed",
  severity: "warning",
  section: "inputs",
  path: "spare",
  target: { kind: "node", role: "input", id: "spare" },
  message: "feed `spare` is declared and never consumed",
  detail: {},
};

describe("the message groups", () => {
  it("puts a document's facts under the object they are about", () => {
    const groups = groupDiagnostics(broken);
    expect(groups.map((group) => group.title)).toEqual(["Pump p1", "Separator sep1"]);
    expect(groups[0]?.select).toBe("instance:p1");
    expect(groups[0]?.summary).toBe("not connected yet");
    expect(groups[0]?.messages[0]?.text).toBe("The outlet `outlet` is not connected to anything.");
    expect(groups[1]?.messages[0]?.text).toBe("The inlet `feed` is not connected.");
  });

  it("calls a construction-only document incomplete, not refused", () => {
    expect(statusWord(broken)).toBe("incomplete");
    expect(statusWord(solved)).toBe("solved");
    expect(statusWord({ ...solved, dirty: true })).toBe("stale");
    expect(statusWord({ ...solved, run_error: "the flash refused" })).toBe("failed");
    expect(statusWord(withDiagnostic(broken, { ...warning, code: "type_mismatch", severity: "error" }))).toBe(
      "refused",
    );
  });
});

describe("the messages panel", () => {
  it("leads with the object and folds the code into a technical detail", () => {
    const onSelect = vi.fn();
    const { container } = render(
      <DiagnosticsPanel envelope={broken} runError={null} onSelect={onSelect} />,
    );

    expect(screen.getByRole("heading", { name: /Messages · 2 errors/ })).toBeDefined();
    expect(screen.getByText("Pump p1")).toBeDefined();
    expect(screen.getByText("The outlet `outlet` is not connected to anything.")).toBeDefined();
    expect(container.querySelector(".where")).toBeNull();
    expect(container.querySelector(".message code")?.textContent).toContain("underfed_port");

    fireEvent.click(screen.getByRole("button", { name: /Pump p1/ }));
    expect(onSelect).toHaveBeenCalledWith("instance:p1");
  });

  it("counts a warning and a run failure instead of saying zero and no objection", () => {
    const warned: Envelope = { ...solved, ok: true, diagnostics: [warning] };
    render(<DiagnosticsPanel envelope={warned} runError={null} onSelect={() => {}} />);
    expect(screen.getByRole("heading", { name: /0 errors · 1 warning/ })).toBeDefined();
  });

  it("shows the run failure as a message and not under an empty state", () => {
    render(
      <DiagnosticsPanel envelope={solved} runError="the flash is a trivial solution" onSelect={() => {}} />,
    );
    expect(screen.getByRole("heading", { name: /1 run failure/ })).toBeDefined();
    expect(screen.getByText("This came from the run, not the document check.")).toBeDefined();
    expect(screen.queryByText("No problems found.")).toBeNull();
  });
});
