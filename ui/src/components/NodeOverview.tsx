/**
 * A node's whole story, in one read-only card.
 *
 * **Double-click is the gesture the canvas had free.** Single click selects and opens the property
 * bar; this is the summary a reader wants before deciding which tab to edit. Everything here is
 * read from the one envelope the editor already holds - no second store, no second query, and no
 * second editor: "Edit in the property bar" hands the object back to the surface that owns writes.
 */

import { useEffect, useRef } from "react";

import { groupDiagnostics } from "../state/diagnostics";
import { provenanceOf } from "../state/provenance";
import { outputRowsOf, resultOf, warningsOf } from "../state/results";
import { runNote, statusClass, statusWord } from "../state/status";
import type { Units } from "../state/units";
import { formatQuantity } from "../wire/field";
import type { Catalogue, Envelope, GraphNode } from "../wire/types";

export function NodeOverview({
  envelope,
  catalogue,
  units,
  node,
  onClose,
  onEdit,
}: {
  envelope: Envelope;
  catalogue: Catalogue | null;
  units: Units;
  node: GraphNode;
  onClose: () => void;
  onEdit: () => void;
}) {
  const dialog = useRef<HTMLElement>(null);
  useEffect(() => {
    const previous = document.activeElement;
    dialog.current?.querySelector<HTMLElement>("button")?.focus();
    return () => {
      if (previous instanceof HTMLElement) {
        previous.focus();
      }
    };
  }, []);

  useEffect(() => {
    const onKey = (event: KeyboardEvent) => {
      if (event.key === "Escape") {
        onClose();
      }
    };
    window.addEventListener("keydown", onKey);
    return () => window.removeEventListener("keydown", onKey);
  }, [onClose]);

  const form = catalogue?.unit_ops.find((entry) => entry.id === node.data.unit);
  const values = node.data.parameters ?? {};
  const provenance = node.role === "instance" ? provenanceOf(envelope, node.data.name) : undefined;
  const result = node.role === "instance" ? resultOf(envelope, node.data.name) : undefined;
  const messages = groupDiagnostics(envelope).find((group) => group.select === node.id);
  const title = node.data.unit_name ?? node.data.name;
  const status = statusWord(envelope);
  const streams = envelope.session?.streams ?? {};
  const ports = node.data.ports;
  const rows = [
    ...(ports?.inlets ?? []).flatMap((port) =>
      port.handles.map((path) => ({ path, side: "in" as const, port: port.name })),
    ),
    ...(ports?.outlets ?? []).flatMap((port) =>
      port.handles.map((path) => ({ path, side: "out" as const, port: port.name })),
    ),
  ];
  const outputRows =
    result === undefined ? [] : outputRowsOf(result, units, form?.outputs ?? []);
  const warnings = result === undefined ? [] : warningsOf(result);
  const setParameters =
    form === undefined
      ? []
      : form.parameters.filter((parameter) => {
          const value = values[parameter.name];
          return value !== undefined && value !== null && value !== "";
        });

  return (
    <div
      className="overview-backdrop"
      onMouseDown={(event) => {
        if (event.target === event.currentTarget) {
          onClose();
        }
      }}
    >
      <section
        ref={dialog}
        className="overview"
        role="dialog"
        aria-modal="true"
        aria-labelledby="overview-title"
        onKeyDown={(event) => {
          if (event.key !== "Tab") {
            return;
          }
          const buttons = [...(dialog.current?.querySelectorAll<HTMLElement>("button") ?? [])];
          const first = buttons[0];
          const last = buttons[buttons.length - 1];
          if (first === undefined || last === undefined) {
            return;
          }
          if (event.shiftKey && document.activeElement === first) {
            last.focus();
            event.preventDefault();
          } else if (!event.shiftKey && document.activeElement === last) {
            first.focus();
            event.preventDefault();
          }
        }}
      >
        <header className="overview-head">
          <div>
            <h2 id="overview-title">{title}</h2>
            <p className="note">
              {node.data.name}
              {node.data.unit === undefined ? "" : ` · ${node.data.unit}`}
              {form?.family === undefined || form.family === null ? "" : ` · ${form.family}`}
            </p>
          </div>
          <div className="overview-actions">
            <button type="button" onClick={onClose}>
              Close
            </button>
            <button type="button" className="primary" onClick={onEdit}>
              Edit in the property bar
            </button>
          </div>
        </header>

        <div className="overview-body">
          <section className="overview-section">
            <h3>Run</h3>
            <p>
              <span className={`chip ${statusClass(envelope)}`}>{status}</span>
              {messages === undefined ? null : (
                <span className={`overview-severity ${messages.severity}`}>
                  {messages.summary}
                </span>
              )}
            </p>
            <p className="note">{runNote(status, outputRows.length > 0)}</p>
            {provenance === undefined ? (
              <p className="note">
                {node.role === "instance"
                  ? "the document has no provenance for this instance"
                  : "a boundary stream has no model check"}
              </p>
            ) : (
              <dl className="overview-facts">
                <dt>model</dt>
                <dd>{provenance.model ?? "—"}</dd>
                <dt>source</dt>
                <dd>{provenance.source ?? "—"}</dd>
                <dt>data authority</dt>
                <dd>{envelope.provenance.keycard}</dd>
                <dt>kernel</dt>
                <dd>{provenance.runnable ? "available" : (provenance.refusal ?? "not available")}</dd>
              </dl>
            )}
            <p className="note">the spec and kernel hashes are on the Provenance tab</p>
          </section>

          {node.data.input === undefined ? null : (
            <section className="overview-section">
              <h3>Feed</h3>
              <dl className="overview-facts">
                <dt>flow</dt>
                <dd>{formatQuantity(node.data.input.n, "mol/s", 4, units)}</dd>
                <dt>P</dt>
                <dd>{formatQuantity(node.data.input.P, "Pa", 4, units)}</dd>
                <dt>T</dt>
                <dd>{formatQuantity(node.data.input.T, "K", 4, units)}</dd>
                <dt>components</dt>
                <dd>
                  {node.data.input.components
                    .map((name, index) => `${name} ${node.data.input?.z[index] ?? 0}`)
                    .join(", ")}
                </dd>
              </dl>
            </section>
          )}

          {rows.length === 0 ? null : (
            <section className="overview-section">
              <h3>Streams</h3>
              <table className="grid">
                <thead>
                  <tr>
                    <th scope="col">stream</th>
                    <th scope="col">side</th>
                    <th scope="col">n</th>
                    <th scope="col">T</th>
                    <th scope="col">P</th>
                    <th scope="col">VF</th>
                  </tr>
                </thead>
                <tbody>
                  {rows.map(({ path, side, port }) => {
                    const stream = streams[path];
                    return (
                      <tr key={path}>
                        <th scope="row" title={port}>
                          {path}
                        </th>
                        <td>{side}</td>
                        <td data-num>
                          {stream === undefined
                            ? "—"
                            : formatQuantity(stream.n.magnitude_si, stream.n.unit, 3, units)}
                        </td>
                        <td data-num>
                          {stream === undefined
                            ? "—"
                            : formatQuantity(stream.T.magnitude_si, stream.T.unit, 4, units)}
                        </td>
                        <td data-num>
                          {stream === undefined
                            ? "—"
                            : formatQuantity(stream.P.magnitude_si, stream.P.unit, 3, units)}
                        </td>
                        <td data-num>
                          {stream?.vapour_fraction === null || stream?.vapour_fraction === undefined
                            ? "—"
                            : stream.vapour_fraction.toFixed(4)}
                        </td>
                      </tr>
                    );
                  })}
                </tbody>
              </table>
            </section>
          )}

          {setParameters.length === 0 ? null : (
            <section className="overview-section">
              <h3>Design</h3>
              <dl className="overview-facts">
                {setParameters.map((parameter) => (
                  <div className="overview-fact" key={parameter.name}>
                    <dt title={parameter.description}>{parameter.name}</dt>
                    <dd>{parameterText(values[parameter.name], parameter.unit, units)}</dd>
                  </div>
                ))}
              </dl>
            </section>
          )}

          {outputRows.length === 0 && warnings.length === 0 ? null : (
            <section className="overview-section">
              <h3>Results</h3>
              {outputRows.length === 0 ? null : (
                <table className="grid result-table">
                  <thead>
                    <tr>
                      <th scope="col">output</th>
                      <th scope="col">value</th>
                      <th scope="col">unit</th>
                    </tr>
                  </thead>
                  <tbody>
                    {outputRows.map((row) => (
                      <tr key={row.key}>
                        <th scope="row">
                          <span className="result-label">{row.label}</span>
                          {row.description === null ? null : (
                            <p className="result-description">{row.description}</p>
                          )}
                        </th>
                        <td data-num={row.numeric ? "true" : undefined}>{row.text}</td>
                        <td className="unit">{row.unit ?? ""}</td>
                      </tr>
                    ))}
                  </tbody>
                </table>
              )}
              {warnings.length === 0 ? null : (
                <>
                  <h4>Caveats</h4>
                  {warnings.map((warning, index) => (
                    <p key={`${warning.code}-${String(index)}`}>
                      <strong>{warning.code}</strong> {warning.message}
                    </p>
                  ))}
                </>
              )}
            </section>
          )}

          {messages === undefined ? null : (
            <section className="overview-section">
              <h3>Messages</h3>
              {messages.messages.map((message) => (
                <p key={message.detail} title={message.detail}>
                  {message.text}
                </p>
              ))}
            </section>
          )}
        </div>
      </section>
    </div>
  );
}

/** A declared parameter, in the unit set in force where it has one. */
function parameterText(value: unknown, unit: string | null, units: Units): string {
  if (value === undefined || value === null || value === "") {
    return "—";
  }
  if (typeof value === "number" && unit !== null && unit !== "dimensionless") {
    return formatQuantity(value, unit, 6, units);
  }
  return String(value);
}
