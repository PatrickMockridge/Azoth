/**
 * The editor: a palette, a canvas, and the panels that read whatever the last call answered with.
 *
 * **The state is one envelope.** There is no second store and no reducer over the graph: an edit
 * goes out as a command, the library answers with the whole document, and this sets it. That is
 * what keeps the canvas, the form and the diagnostics from disagreeing — they are three readings
 * of one object rather than three copies of one state.
 */

import { useCallback, useEffect, useRef, useState } from "react";
import "@xyflow/react/dist/style.css";

import demoDocument from "../../specs/flowsheets/demo.toml?raw";
import { BoundaryPanel } from "./components/BoundaryPanel";
import { DiagnosticsPanel } from "./components/DiagnosticsPanel";
import { EdgePanel } from "./components/EdgePanel";
import { Flowsheet } from "./components/Flowsheet";
import { InputsPanel } from "./components/InputsPanel";
import { Palette } from "./components/Palette";
import { UnitOpPanel } from "./components/UnitOpPanel";
import { selectedEdge, selectedNode, targetNodeId } from "./state/selection";
import type { Middleware, Session } from "./wire/client";
import { wasmMiddleware } from "./wire/client";
import type { Catalogue, Command, Envelope, ExecutionOrder } from "./wire/types";

export function App() {
  const [session, setSession] = useState<Session | null>(null);
  const [catalogue, setCatalogue] = useState<Catalogue | null>(null);
  const [envelope, setEnvelope] = useState<Envelope | null>(null);
  const [selected, setSelected] = useState<string | null>(null);
  const [fault, setFault] = useState<string | null>(null);

  const middleware = useRef<Middleware | null>(null);

  /** Open a document, from anywhere: the demo on load, the Demo button, or a file. */
  const open = useCallback((text: string) => {
    const loaded = middleware.current;
    if (loaded === null) {
      return;
    }
    try {
      const opened = loaded.open(text);
      setSession(opened);
      setEnvelope(opened.envelope());
      setSelected(null);
      setFault(null);
    } catch (error: unknown) {
      // A document the schema cannot read is a fault, which is what the editor already shows for
      // a refusal. The session that was open stays open, because the new one never became one.
      setFault(String(error));
    }
  }, []);

  useEffect(() => {
    let live = true;
    wasmMiddleware()
      .then((loaded) => {
        if (!live) {
          return;
        }
        middleware.current = loaded;
        setCatalogue(loaded.catalogue(true));
        open(demoDocument);
      })
      .catch((error: unknown) => {
        if (live) {
          setFault(String(error));
        }
      });
    return () => {
      live = false;
    };
  }, [open]);

  /** One edit, one answer. A command the library refuses is a fault rather than a silent no-op. */
  const send = useCallback(
    (command: Command) => {
      if (session === null) {
        return;
      }
      try {
        setEnvelope(session.apply(command));
        setFault(null);
      } catch (error: unknown) {
        setFault(String(error));
      }
    },
    [session],
  );

  const solve = useCallback(() => {
    if (session === null) {
      return;
    }
    try {
      setEnvelope(session.run());
      setFault(null);
    } catch (error: unknown) {
      setFault(String(error));
    }
  }, [session]);

  const setOrder = useCallback(
    (order: ExecutionOrder) => {
      if (session === null) {
        return;
      }
      try {
        setEnvelope(session.setOrder(order));
        setFault(null);
      } catch (error: unknown) {
        setFault(String(error));
      }
    },
    [session],
  );

  const save = useCallback(() => {
    if (session === null) {
      return;
    }
    // `text`, not `document`: the global `document` is what the download anchor is created from.
    const text = session.document();
    const url = URL.createObjectURL(new Blob([text], { type: "text/plain" }));
    const link = document.createElement("a");
    link.href = url;
    link.download = `${envelope?.flowsheet.id ?? "flowsheet"}.toml`;
    link.click();
    URL.revokeObjectURL(url);
  }, [session, envelope]);

  const node = envelope === null ? null : selectedNode(envelope, selected);
  const edge = envelope === null ? null : selectedEdge(envelope, selected);

  return (
    <div className="app">
      <header className="bar">
        <h1>azoth</h1>
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
            disabled={session === null}
            title="which order the next run takes; ProcessSystem.useGraphBasedExecution is a flag, so there are two"
            onChange={(event) => setOrder(event.target.value as ExecutionOrder)}
          >
            <option value="insertion">insertion</option>
            <option value="topological">topological</option>
          </select>
        )}
        <button type="button" onClick={solve} disabled={session === null}>
          Solve
        </button>
        <label className="button">
          Open
          <input
            type="file"
            accept=".toml,text/plain"
            onChange={(event) => {
              const file = event.target.files?.[0];
              if (file !== undefined) {
                void file.text().then(open);
              }
              // So that choosing the same file twice fires again.
              event.target.value = "";
            }}
          />
        </label>
        <button type="button" onClick={() => open(demoDocument)} disabled={session === null}>
          Demo
        </button>
        <button type="button" onClick={save} disabled={session === null}>
          Save
        </button>
      </header>

      <div className="body">
        <Palette
          forms={catalogue?.unit_ops ?? []}
          onAdd={(unit) =>
            send({
              command: "add_instance",
              id: nextId(unit, envelope),
              unit,
              parameters: {},
            })
          }
        >
          {envelope === null ? null : (
            <BoundaryPanel envelope={envelope} onCommand={send} />
          )}
        </Palette>

        <div className="canvas">
          {envelope === null ? (
            <p className="note">{fault ?? "loading the middleware…"}</p>
          ) : (
            <Flowsheet
              catalogue={catalogue}
              envelope={envelope}
              selected={selected}
              onSelect={setSelected}
              onCommand={send}
            />
          )}
        </div>

        <aside className="side right">
          {fault !== null ? <div className="fault">{fault}</div> : null}
          {envelope !== null && node !== null && node.role === "instance" ? (
            <UnitOpPanel
              catalogue={catalogue}
              envelope={envelope}
              node={node}
              onCommand={send}
            />
          ) : null}
          {envelope !== null && node !== null && node.role !== "instance" ? (
            <InputsPanel envelope={envelope} node={node} onCommand={send} />
          ) : null}
          {envelope !== null && edge !== null ? (
            <EdgePanel envelope={envelope} edge={edge} onCommand={send} />
          ) : null}
          {envelope !== null ? (
            <DiagnosticsPanel
              diagnostics={envelope.diagnostics}
              runError={envelope.run_error}
              onSelect={(target) => setSelected(targetNodeId(target, envelope))}
            />
          ) : null}
        </aside>
      </div>
    </div>
  );
}

/** The word a status pill shows. */
function status(envelope: Envelope): string {
  if (!envelope.ok) {
    return "refused";
  }
  if (envelope.run_error !== null) {
    return "failed";
  }
  return envelope.dirty ? "stale" : "solved";
}

function statusClass(envelope: Envelope): string {
  if (!envelope.ok || envelope.run_error !== null) {
    return "bad";
  }
  return envelope.dirty ? "stale" : "ok";
}

/**
 * A fresh instance id, from the entry's own name and what the document already holds.
 *
 * A front-end has to choose an id because the document's ids are its own; this is the one place
 * the editor invents a name, and it checks the document first so that a second click on the same
 * palette entry does not replace the first node.
 */
function nextId(unit: string, envelope: Envelope | null): string {
  const stem = unit.replace("unit_ops.", "").replace(/[^a-z0-9]+/g, "_");
  const taken = new Set(envelope?.flowsheet.graph.nodes.map((node) => node.id) ?? []);
  for (let index = 1; index < 100; index += 1) {
    const candidate = `${stem}_${index}`;
    if (!taken.has(`instance:${candidate}`)) {
      return candidate;
    }
  }
  return `${stem}_x`;
}
