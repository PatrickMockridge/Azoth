/**
 * A connection, or a tear — one panel, because a tear *is* a connection with a name.
 *
 * **This is where the tear's convergence becomes reachable.** The seven settings live on the edge
 * (`graph::EdgeData::settings`), project from the document, and are sent back one at a time as
 * `set_recycle`, so the panel holds no copy: it reads the envelope and writes a command.
 *
 * **`acceleration_method` is a select, and its options come from the wire.** The names live once,
 * in `recycle::ACCELERATION_NAMES`, and the tool schema publishes them; this reads them
 * (`accelerationNames`) rather than typing them again, so the widget cannot offer a name the parser
 * would refuse. One of the names is refused by *this port* — with the measurement its own capture
 * makes — so it is offered, and the checker reports it as `acceleration` on the tear rather than the
 * widget hiding it. A schema that published nothing leaves the field a text box, which is the
 * honest fallback rather than an empty dropdown.
 */

import type { EditorCommand, RecycleField } from "../wire/commands";
import { accelerationNames } from "../wire/field";
import type { Catalogue, Envelope, GraphEdge } from "../wire/types";

export interface EdgePanelProps {
  catalogue: Catalogue | null;
  envelope: Envelope;
  edge: GraphEdge;
  onCommand: (command: EditorCommand) => void;
}

/** The seven settings, in the order `Recycle` declares them, with what each one says. */
const SETTINGS: readonly { field: RecycleField; kind: "number" | "text" | "names"; hint: string }[] = [
  { field: "flow_tolerance", kind: "number", hint: "the flow residual the tear solves to" },
  { field: "composition_tolerance", kind: "number", hint: "the composition residual" },
  { field: "temperature_tolerance", kind: "number", hint: "the temperature residual" },
  { field: "pressure_tolerance", kind: "number", hint: "the pressure residual" },
  { field: "max_iterations", kind: "number", hint: "the pass cap before the tear gives up" },
  { field: "minimum_flow", kind: "number", hint: "below it the tear is switched off outright" },
  {
    field: "acceleration_method",
    kind: "names",
    // The names are not written here: the schema publishes them and this reads them. One of them
    // is refused by this port, and which one - and why - is the checker's `acceleration`
    // diagnostic, carrying the class's own measurement.
    hint: "the checker reports the one this port refuses, with the class's own measurement",
  },
];

export function EdgePanel({ catalogue, envelope, edge, onCommand }: EdgePanelProps) {
  const names = accelerationNames(catalogue);
  const isTear = edge.data.kind === "recycle";
  const settings = edge.data.settings;

  return (
    <>
      <h2>{isTear ? `Tear ${edge.data.path}` : "Connection"}</h2>
      <p className="note">
        {edge.data.from} → {edge.data.to}
      </p>

      <div className="field">
        {isTear ? (
          <button
            type="button"
            className="entry"
            onClick={() => onCommand({ command: "remove_recycle", stream: edge.data.path })}
          >
            Remove the tear
          </button>
        ) : (
          <>
            <button
              type="button"
              className="entry"
              onClick={() =>
                onCommand({
                  command: "disconnect",
                  from: edge.data.from,
                  to: edge.data.to,
                })
              }
            >
              Remove the connection
            </button>
            <button
              type="button"
              className="entry"
              onClick={() =>
                onCommand({
                  command: "add_recycle",
                  stream: nextTear(envelope),
                  from: edge.data.from,
                  to: edge.data.to,
                })
              }
            >
              Make it a tear
            </button>
          </>
        )}
      </div>

      {!isTear || settings === undefined ? null : (
        <>
          <h2>Convergence</h2>
          <p className="note">
            an unstated setting is the class&apos;s own default — the field is empty, and writing a
            number in would pin it
          </p>
          {SETTINGS.map(({ field, kind, hint }) => {
            const value = settings[field];
            return (
              <div className="field" key={field}>
                <div className="head">
                  <span className="name">{field}</span>
                </div>
                {kind === "names" && names.length > 0 ? (
                  // **The options are the schema's**, so this cannot offer a name the parser would
                  // refuse. The empty option is the state the document already has: unstated.
                  <select
                    value={value === null ? "" : String(value)}
                    onChange={(event) => {
                      const raw = event.target.value;
                      if (raw === "") {
                        return;
                      }
                      onCommand({
                        command: "set_recycle",
                        stream: edge.data.path,
                        field,
                        value: raw,
                      });
                    }}
                  >
                    <option value="">the class&apos;s default</option>
                    {names.map((name) => (
                      <option key={name} value={name}>
                        {name}
                      </option>
                    ))}
                  </select>
                ) : (
                  <input
                    value={value === null ? "" : String(value)}
                    inputMode={kind === "number" ? "decimal" : "text"}
                    placeholder="the class's default"
                    onChange={(event) => {
                      const raw = event.target.value;
                      if (raw === "") {
                        // The command model has `set_recycle` and no inverse, so a stated setting
                        // cannot be returned to silence. The field says so by sending nothing.
                        return;
                      }
                      onCommand({
                        command: "set_recycle",
                        stream: edge.data.path,
                        field,
                        value: kind === "number" ? Number(raw) : raw,
                      });
                    }}
                  />
                )}
                <div className="hint">{hint}</div>
              </div>
            );
          })}
        </>
      )}
    </>
  );
}

/** A fresh tear name, over the ones the document already holds. */
function nextTear(envelope: Envelope): string {
  const taken = new Set(envelope.flowsheet.graph.edges.map((edge) => edge.data.path));
  for (let index = 1; index < 100; index += 1) {
    const candidate = `recycle_${index}`;
    if (!taken.has(candidate)) {
      return candidate;
    }
  }
  return "recycle_x";
}
