import { status, statusClass } from "./TitleStrip";
import type { Envelope } from "../wire/types";

/**
 * The bottom band: what the last run reached, and what the document holds.
 *
 * **Every word here is a reading rather than a control.** A status bar that offered an action
 * would be a second place to look for it, and the one reading that is *not* yet a choice is the
 * unit set: a value arrives in its canonical SI unit and the middleware converts nothing on the
 * way out, so "units: SI" is a statement of fact - see `docs/src/architecture/middleware.md` on
 * what a unit set would take.
 */
export function StatusBar({ envelope }: { envelope: Envelope | null }) {
  const instances =
    envelope?.flowsheet.graph.nodes.filter((node) => node.role === "instance").length ?? 0;
  const recycles =
    envelope?.flowsheet.graph.edges.filter((edge) => edge.data.kind === "recycle").length ?? 0;
  const streams = envelope?.session === null || envelope?.session === undefined
    ? 0
    : Object.keys(envelope.session.streams).length;
  const converged = envelope?.session?.converged;

  return (
    <footer className="status-bar" aria-label="status">
      {envelope === null ? (
        <span>no document</span>
      ) : (
        <>
          <span className={`chip ${statusClass(envelope)}`}>{status(envelope)}</span>
          <span>
            passes <span className="value">{envelope.session?.iterations ?? 0}</span>
          </span>
          <span>
            converged <span className="value">{converged === undefined ? "—" : String(converged)}</span>
          </span>
          <span>
            order <span className="value">{envelope.execution_order}</span>
          </span>
          <span>
            unit ops <span className="value">{instances}</span> · streams{" "}
            <span className="value">{streams}</span> · recycles{" "}
            <span className="value">{recycles}</span>
          </span>
        </>
      )}
      <span className="spacer" />
      <span>
        units <span className="value">SI</span>
      </span>
      <span>
        document{" "}
        <span className="value">{envelope === null ? "—" : envelope.flowsheet.id}</span>
      </span>
    </footer>
  );
}
