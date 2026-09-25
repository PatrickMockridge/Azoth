/**
 * What the checker said, each entry clickable to the thing it is about.
 *
 * **A diagnostic is a record, not a line of text**, and this panel is why that matters: the code
 * is what a front-end switches on, the target is what a click selects, and the message is what a
 * person reads — three readings of one object rather than three ways of parsing a string.
 *
 * A run that failed is shown here too, and flagged as what it is: a fact the checker *cannot*
 * know, because whether a substance exists is the databank's answer and not a fact about the
 * document.
 */

import type { Diagnostic, Target } from "../wire/types";

export interface DiagnosticsPanelProps {
  diagnostics: Diagnostic[];
  runError: string | null;
  /**
   * What a click selects, as the *target* rather than as a node id.
   *
   * The checker derives a target from the variant's own fields; turning one into the node a canvas
   * can select is `state/selection.ts`'s job and happens once, for this panel and for a click on
   * the canvas alike.
   */
  onSelect: (target: Target) => void;
}

export function DiagnosticsPanel({ diagnostics, runError, onSelect }: DiagnosticsPanelProps) {
  const errors = diagnostics.filter((diagnostic) => diagnostic.severity === "error").length;

  return (
    <>
      <h2>
        Diagnostics
        {diagnostics.length === 0 ? "" : ` · ${errors} error${errors === 1 ? "" : "s"}`}
      </h2>
      {diagnostics.length === 0 ? (
        <p className="note">the checker has no objection to this document</p>
      ) : (
        diagnostics.map((diagnostic, index) => (
          <button
            key={`${diagnostic.code}-${diagnostic.path}-${index}`}
            type="button"
            className={`diagnostic ${diagnostic.severity}`}
            onClick={() => onSelect(diagnostic.target)}
          >
            <span className="code">{diagnostic.code}</span>{" "}
            <span className="where">{where(diagnostic)}</span>
            <p>{diagnostic.message}</p>
          </button>
        ))
      )}
      {runError === null ? null : (
        <>
          <h2>The run</h2>
          <div className="fault">{runError}</div>
          <p className="note">
            the checker cannot see this one: it holds a document to what the document decides, and
            whether a substance exists is the databank's answer
          </p>
        </>
      )}
    </>
  );
}

/** Where in the document it is about, in the checker's own spelling. */
function where(diagnostic: Diagnostic): string {
  return diagnostic.path === ""
    ? diagnostic.section
    : `${diagnostic.section}/${diagnostic.path}`;
}
