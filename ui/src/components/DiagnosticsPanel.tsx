/**
 * What the checker and the run said, grouped by the object a reader can act on.
 *
 * **The object leads and the code follows.** The panel used to print `underfed_port` and
 * `instances/p1.ports.outlet` as the first line and a sentence third; a process engineer was
 * reading a type system. A group is one machine, feed or connection; its facts are sentences; the
 * code and dotted path are behind a technical disclosure.
 *
 * **A run failure is a message of the same kind.** It counted zero in the tab and sat under an
 * empty “no objection” state, which was two screens disagreeing about one document.
 */

import { groupDiagnostics } from "../state/diagnostics";
import type { Envelope } from "../wire/types";

export interface DiagnosticsPanelProps {
  envelope: Envelope;
  runError: string | null;
  /** The node or edge id a row selects; the panel has already resolved it. */
  onSelect: (id: string | null) => void;
}

export function DiagnosticsPanel({ envelope, runError, onSelect }: DiagnosticsPanelProps) {
  const groups = groupDiagnostics(envelope);
  const errors = groups.filter((group) => group.severity === "error").length;
  const warnings = groups.filter((group) => group.severity === "warning").length;
  const hasRun = runError !== null;

  return (
    <>
      <h2>
        Messages
        {groups.length === 0 && !hasRun
          ? ""
          : ` · ${String(errors)} error${errors === 1 ? "" : "s"}` +
            (warnings === 0 ? "" : ` · ${String(warnings)} warning${warnings === 1 ? "" : "s"}`) +
            (hasRun ? " · 1 run failure" : "")}
      </h2>

      {groups.length === 0 && !hasRun ? (
        <p className="note">No problems found.</p>
      ) : null}

      {groups.map((group) => (
        <section className={`message ${group.severity}`} key={group.key}>
          {group.select === null ? (
            <div className="message-head">
              <span className="title">{group.title}</span>
              <span className="summary">{group.summary}</span>
            </div>
          ) : (
            <button
              type="button"
              className="message-head"
              title="select this object"
              onClick={() => onSelect(group.select)}
            >
              <span className="title">{group.title}</span>
              <span className="summary">{group.summary}</span>
            </button>
          )}
          <ul className="message-list">
            {group.messages.map((message, index) => (
              <li key={`${message.code}-${String(index)}`}>
                <span>{message.text}</span>
                <details>
                  <summary>Technical detail</summary>
                  <code>{message.detail}</code>
                </details>
              </li>
            ))}
          </ul>
        </section>
      ))}

      {hasRun ? (
        <section className="message error run">
          <div className="message-head">
            <span className="title">The run</span>
            <span className="summary">failed</span>
          </div>
          <p className="fault">{runError}</p>
          <p className="note">This came from the run, not the document check.</p>
        </section>
      ) : null}

      {groups.length === 0 && !hasRun ? (
        <p className="note">Select an object on the canvas or in the navigator to inspect it.</p>
      ) : null}
    </>
  );
}
