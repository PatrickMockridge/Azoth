/**
 * The calculation book's title strip.
 *
 * **No `.pill`.** The flowsheet's status pill is the only `.pill` on the page because the canvas
 * test reads the last one for the status and the first for the id; a calculation book says its
 * state in a `.chip`, which may live anywhere.
 */

import { DraftField } from "../DraftField";
import { WorkspaceSwitch, type Workspace } from "./WorkspaceSwitch";
import type { Theme } from "../../state/theme";
import type { Units } from "../../state/units";
import type { SheetEnvelope } from "../../wire/types";

export function CalcTitleStrip({
  envelope,
  workspace,
  theme,
  units,
  unsaved,
  version,
  onName,
  onWorkspace,
  onTheme,
  onUnits,
}: {
  envelope: SheetEnvelope | null;
  workspace: Workspace;
  theme: Theme;
  units: Units;
  unsaved: boolean;
  /** The library version, from the flowsheet's provenance record. */
  version: string;
  onName: (name: string) => void;
  onWorkspace: (workspace: Workspace) => void;
  onTheme: () => void;
  onUnits: (id: string) => void;
}) {
  const build = (import.meta.env.VITE_AZOTH_BUILD ?? "").trim();
  const state = status(envelope);
  return (
    <header className="bar">
      <h1>azoth</h1>
      <a
        className="screening"
        href="https://github.com/PatrickMockridge/azoth/blob/main/docs/src/index.md#not-for-design-work-yet"
        target="_blank"
        rel="noreferrer"
        title="the port is incomplete, so these are screening values and not design deliverables"
      >
        not for design work yet
      </a>
      <WorkspaceSwitch workspace={workspace} onWorkspace={onWorkspace} />
      {envelope === null ? null : (
        <div className="case-name">
          <DraftField
            value={envelope.book.name}
            kind="text"
            label="Calculation book name"
            onCommit={(value) => {
              if (typeof value === "string" && value.trim() !== "") {
                onName(value.trim());
              }
            }}
          />
          {unsaved ? (
            <span className="unsaved" title="the book has edits that are not saved">
              unsaved
            </span>
          ) : null}
        </div>
      )}
      <span className="spacer" />
      {envelope === null ? null : (
        <>
          <span className="chip">{envelope.book.id}</span>
          <span className={`chip ${state.className}`}>{state.word}</span>
        </>
      )}
      {units.sets.length === 0 ? null : (
        <select
          className="units"
          value={units.active}
          title="the unit a value is read in; the calculation stores SI and the field shows the set"
          onChange={(event) => onUnits(event.target.value)}
        >
          {units.sets.map((set) => (
            <option key={set.id} value={set.id}>
              {set.name}
            </option>
          ))}
        </select>
      )}
      <span className="build" title="the library version, and the commit this page was built from">
        v{version}
        {build === "" ? "" : ` · ${build.slice(0, 7)}`}
      </span>
      <button
        type="button"
        className="button"
        onClick={onTheme}
        title={`show the editor in its ${theme === "dark" ? "light" : "dark"} theme`}
      >
        {theme === "dark" ? "Light" : "Dark"}
      </button>
    </header>
  );
}

/** The book's state, in the chip vocabulary. */
function status(envelope: SheetEnvelope | null): {
  word: string;
  className: "ok" | "stale" | "warn" | "bad";
} {
  if (envelope === null) {
    return { word: "loading", className: "warn" };
  }
  if (!envelope.ok) {
    return { word: "incomplete", className: "warn" };
  }
  if (envelope.book.sheets.some((sheet) => sheet.error !== null)) {
    return { word: "failed", className: "bad" };
  }
  if (envelope.dirty) {
    return { word: "stale", className: "stale" };
  }
  return { word: "calculated", className: "ok" };
}
