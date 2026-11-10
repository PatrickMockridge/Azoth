import { asExecutionOrder } from "../wire/decode";
import { manualHome } from "../manual";
import { DraftField } from "./DraftField";
import { WorkspaceSwitch, type Workspace } from "./calculations/WorkspaceSwitch";
import { statusClass as classOf, statusWord } from "../state/status";
import type { Envelope, ExecutionOrder } from "../wire/types";
import type { Theme } from "../state/theme";
import type { Units } from "../state/units";

/**
 * The top band: what this is, what it holds, and the two settings that are not commands.
 *
 * **The status pill lives here and not in the status bar**, and the reason is a test rather than a
 * taste: `ui/test/xyflow-env.ts` reads the document's last `.pill` for the status and its first for
 * the flowsheet's id, so the two of them are the only `.pill`s anywhere. The status bar repeats the
 * same word as a `.chip`, which is where a chip may live.
 */
export function TitleStrip({
  envelope,
  session,
  theme,
  units,
  unsaved,
  workspace,
  onWorkspace,
  onName,
  onOrder,
  onTheme,
  onUnits,
}: {
  envelope: Envelope | null;
  /** Whether a session is open, which is what makes the order settable. */
  session: boolean;
  theme: Theme;
  units: Units;
  /** The document differs from the one last opened or saved. */
  unsaved: boolean;
  workspace: Workspace;
  onWorkspace: (workspace: Workspace) => void;
  onName: (name: string) => void;
  onOrder: (order: ExecutionOrder) => void;
  onTheme: () => void;
  onUnits: (id: string) => void;
}) {
  // The commit the page was built from, where the deploy stamps it. A cached page can then be
  // told apart from a fresh one without opening the devtools network panel.
  const build = (import.meta.env.VITE_AZOTH_BUILD ?? "").trim();
  return (
    <header className="bar">
      <h1>azoth</h1>
      {/* **The status, before any number is believable.** The port is incomplete and every
          result ships `unverified` or `partially_verified`; a process engineer must meet that
          fact here rather than in a footnote. The link is the book's own front-door statement. */}
      <a
        className="screening"
        href={`${manualHome()}index.html#not-for-design-work-yet`}
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
            value={envelope.flowsheet.name}
            kind="text"
            label="Case name"
            onCommit={(value) => {
              if (typeof value === "string") {
                onName(value);
              }
            }}
          />
          {unsaved ? (
            <span className="unsaved" title="the document has edits that are not saved">
              unsaved
            </span>
          ) : null}
        </div>
      )}
      <span className="spacer" />
      {envelope === null ? null : (
        <>
          <span className="pill">{envelope.flowsheet.id}</span>
          <span className={`pill ${statusClass(envelope)}`}>{status(envelope)}</span>
        </>
      )}
      {envelope === null ? null : (
        <select
          className="order"
          value={envelope.execution_order}
          disabled={!session}
          title="which order the next run takes; ProcessSystem.useGraphBasedExecution is a flag, so there are two"
          onChange={(event) => {
            const chosen = asExecutionOrder(event.target.value);
            if (chosen !== null) {
              onOrder(chosen);
            }
          }}
        >
          <option value="insertion">insertion</option>
          <option value="topological">topological</option>
        </select>
      )}
      {/* **A set converts what is drawn and nothing else.** It is not a command, not a document
          fact and not on the wire: the numbers a run reached are the same numbers, read in another
          unit, and Solve has no opinion about which. The sets are the library's, so a set added to
          the vocabulary appears here without an edit. */}
      {units.sets.length === 0 ? null : (
        <select
          className="units"
          value={units.active}
          title="the unit a value is read in; the library declares the sets, and this changes the display and not the document"
          onChange={(event) => onUnits(event.target.value)}
        >
          {units.sets.map((set) => (
            <option key={set.id} value={set.id}>
              {set.name}
            </option>
          ))}
        </select>
      )}
      <span
        className="build"
        title="the library version, and the commit this page was built from where the deploy names it"
      >
        v{envelope?.provenance.version ?? "—"}
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

/** The word a status pill shows; the rule itself is `state/status.ts`. */
export function status(envelope: Envelope): string {
  return statusWord(envelope);
}

/** The colour class for the word. */
export function statusClass(envelope: Envelope): string {
  return classOf(envelope);
}
