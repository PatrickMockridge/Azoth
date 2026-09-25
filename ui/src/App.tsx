/**
 * The editor: a palette, a canvas, and the panels that read whatever the last call answered with.
 *
 * **The state is one envelope.** There is no second store and no reducer over the graph: an edit
 * goes out as a command, the library answers with the whole document, and this sets it. That is
 * what keeps the canvas, the form and the diagnostics from disagreeing — they are three readings
 * of one object rather than three copies of one state.
 */

import { useCallback, useEffect, useState } from "react";
import "@xyflow/react/dist/style.css";

import demoDocument from "../../specs/flowsheets/demo.toml?raw";
import { DiagnosticsPanel } from "./components/DiagnosticsPanel";
import { Flowsheet } from "./components/Flowsheet";
import { InputsPanel } from "./components/InputsPanel";
import { Palette } from "./components/Palette";
import { UnitOpPanel } from "./components/UnitOpPanel";
import { selectedNode, targetNodeId } from "./state/selection";
import type { Session } from "./wire/client";
import { wasmMiddleware } from "./wire/client";
import type { Catalogue, Command, Envelope } from "./wire/types";

export function App() {
  const [session, setSession] = useState<Session | null>(null);
  const [catalogue, setCatalogue] = useState<Catalogue | null>(null);
  const [envelope, setEnvelope] = useState<Envelope | null>(null);
  const [selected, setSelected] = useState<string | null>(null);
  const [fault, setFault] = useState<string | null>(null);

  useEffect(() => {
    let live = true;
    wasmMiddleware()
      .then((loaded) => {
        if (!live) {
          return;
        }
        const opened = loaded.open(demoDocument);
        setCatalogue(loaded.catalogue(true));
        setSession(opened);
        setEnvelope(opened.envelope());
      })
      .catch((error: unknown) => {
        if (live) {
          setFault(String(error));
        }
      });
    return () => {
      live = false;
    };
  }, []);

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
        <button type="button" onClick={solve} disabled={session === null}>
          Solve
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
        />

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
