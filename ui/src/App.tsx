/**
 * The editor: a palette, a canvas, and the panels that read whatever the last call answered with.
 *
 * **The state is one envelope.** There is no second store and no reducer over the graph: an edit
 * goes out as a command, the library answers with the whole document, and this sets it. That is
 * what keeps the canvas, the form and the diagnostics from disagreeing — they are three readings
 * of one object rather than three copies of one state.
 *
 * **The door is chosen once, from the URL** (`./wire/door.ts`), and what happens after that is the
 * same for both: a call, and an envelope back. The two differ in one place only — the controls that
 * *hand a document to the library* exist on a door that can be handed one, and a door that hosts a
 * single document for the life of its process is not one.
 */

import { useCallback, useEffect, useRef, useState } from "react";
import "@xyflow/react/dist/style.css";

import blankDocument from "../../specs/flowsheets/blank.toml?raw";
import demoDocument from "../../specs/flowsheets/demo.toml?raw";
import { BoundaryPanel } from "./components/BoundaryPanel";
import { DiagnosticsPanel } from "./components/DiagnosticsPanel";
import { Dock } from "./components/dock/Dock";
import { Flowsheet } from "./components/Flowsheet";
import { Navigation } from "./components/Navigation";
import { Palette } from "./components/Palette";
import { PropertyView } from "./components/property/PropertyView";
import { RibbonGroup } from "./components/RibbonGroup";
import { StatusBar } from "./components/StatusBar";
import { TitleStrip } from "./components/TitleStrip";
import { Workbook } from "./components/workbook/Workbook";
import { saveDocument } from "./save";
import { selectedEdge, selectedNode, targetNodeId } from "./state/selection";
import { applyTheme, readTheme, type Theme } from "./state/theme";
import { door } from "./wire/door";
import type { Door, Opened, Session } from "./wire/session";
import type { Catalogue, Command, Envelope, ExecutionOrder } from "./wire/types";

/**
 * Why a served door has no New, Open or Demo.
 *
 * In the tool bar as a `title` rather than in a paragraph somewhere: the control is where a person
 * asks the question, so the control is where the answer belongs.
 */
const HOSTED =
  "azoth serve holds one document for the life of its process — restart it with --flowsheet F to serve another";

export function App() {
  const [session, setSession] = useState<Session | null>(null);
  const [catalogue, setCatalogue] = useState<Catalogue | null>(null);
  const [envelope, setEnvelope] = useState<Envelope | null>(null);
  const [selected, setSelected] = useState<string | null>(null);
  const [fault, setFault] = useState<string | null>(null);
  const [openable, setOpenable] = useState(false);
  const [theme, setTheme] = useState<Theme>(readTheme);

  /**
   * The dock's view state, which is the editor's and not the document's.
   *
   * **Messages is the tab it opens on**, because a document that does not check is the state a
   * person most needs told about and the dock is where the telling happens.
   */
  const [docked, setDocked] = useState<"workbook" | "messages">("messages");
  const [dockSize, setDockSize] = useState(200);

  const middleware = useRef<Door | null>(null);

  /**
   * **The newest call's number, so an answer to an older one cannot land on top of it.** The module
   * answers in order because it answers at once; a served door answers when the socket says so, and
   * two edits in flight can come back the other way round — which would leave the editor drawing
   * the older document. Every call takes a ticket, and only the newest one is allowed to set state.
   */
  const issued = useRef(0);

  /** One call's answer, or the sentence that refused it. */
  const answer = useCallback((pending: Promise<Envelope>) => {
    const ticket = (issued.current += 1);
    void pending.then(
      (next) => {
        if (ticket === issued.current) {
          setEnvelope(next);
          setFault(null);
        }
      },
      (error: unknown) => {
        if (ticket === issued.current) {
          setFault(String(error));
        }
      },
    );
  }, []);

  /** A session and the envelope it opened at, as the editor's state. */
  const land = useCallback((opening: Promise<Opened>) => {
    const ticket = (issued.current += 1);
    void opening.then(
      (opened) => {
        if (ticket !== issued.current) {
          return;
        }
        setSession(opened.session);
        setEnvelope(opened.envelope);
        setSelected(null);
        setFault(null);
      },
      (error: unknown) => {
        // A document the schema cannot read is a fault, which is what the editor already shows for
        // a refusal. The session that was open stays open, because the new one never became one.
        if (ticket === issued.current) {
          setFault(String(error));
        }
      },
    );
  }, []);

  /** The theme, in force and remembered. Dark is CSS's default; this only overrides it. */
  useEffect(() => {
    applyTheme(theme);
  }, [theme]);

  /** Hand a document to the library: the Demo button, New, or a file. */
  const open = useCallback(
    (text: string) => {
      const loaded = middleware.current;
      if (loaded === null || loaded.kind !== "openable") {
        return;
      }
      land(loaded.open(text));
    },
    [land],
  );

  useEffect(() => {
    let live = true;
    door(window.location.search)
      .then((loaded) => {
        if (!live) {
          return;
        }
        middleware.current = loaded;
        setOpenable(loaded.kind === "openable");
        void loaded.catalogue(true).then(
          (palette) => {
            if (live) {
              setCatalogue(palette);
            }
          },
          (error: unknown) => {
            if (live) {
              setFault(String(error));
            }
          },
        );
        // **The door says how to open on it.** A served one already holds a document and has no
        // call that would replace it, so asking it to open the demo would be a call it does not
        // have — the kind of thing that reads as a bug in the editor rather than in the request.
        land(loaded.kind === "openable" ? loaded.open(demoDocument) : loaded.attach());
      })
      .catch((error: unknown) => {
        if (live) {
          setFault(String(error));
        }
      });
    return () => {
      live = false;
    };
  }, [land]);

  /** One edit, one answer. A command the library refuses is a fault rather than a silent no-op. */
  const send = useCallback(
    (command: Command) => {
      if (session !== null) {
        answer(session.apply(command));
      }
    },
    [session, answer],
  );

  const solve = useCallback(() => {
    if (session !== null) {
      answer(session.run());
    }
  }, [session, answer]);

  const setOrder = useCallback(
    (order: ExecutionOrder) => {
      if (session !== null) {
        answer(session.setOrder(order));
      }
    },
    [session, answer],
  );

  /**
   * Write the document where a person says.
   *
   * **Off the envelope rather than off the session**, because every envelope carries
   * `flowsheet.document` — which is what makes a save no call at all, on either door.
   */
  const save = useCallback(() => {
    if (envelope === null) {
      return;
    }
    void saveDocument(`${envelope.flowsheet.id}.toml`, envelope.flowsheet.document).catch(
      (error: unknown) => setFault(String(error)),
    );
  }, [envelope]);

  const node = envelope === null ? null : selectedNode(envelope, selected);
  const edge = envelope === null ? null : selectedEdge(envelope, selected);

  return (
    <div className="app">
      <TitleStrip
        envelope={envelope}
        session={session !== null}
        theme={theme}
        onOrder={setOrder}
        onTheme={() => setTheme(theme === "dark" ? "light" : "dark")}
      />

      <div className="ribbon">
        <RibbonGroup label="Run">
          <button
            type="button"
            className="primary"
            onClick={solve}
            disabled={session === null}
            title="run the physics; a check happens on every edit and a run is when you ask"
          >
            Solve
          </button>
        </RibbonGroup>
        <RibbonGroup label="Case">
          <button
            type="button"
            onClick={() => open(blankDocument)}
            disabled={!openable}
            title={openable ? "a document with nothing in it" : HOSTED}
          >
            New
          </button>
          <label className={`button${openable ? "" : " off"}`} title={openable ? undefined : HOSTED}>
            Open
            {openable ? (
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
            ) : null}
          </label>
          <button
            type="button"
            onClick={() => open(demoDocument)}
            disabled={!openable}
            title={openable ? "the shipped demo flowsheet" : HOSTED}
          >
            Demo
          </button>
          <button
            type="button"
            onClick={save}
            disabled={session === null}
            title="write the document where you say; it rides every envelope, so this is no call"
          >
            Save
          </button>
        </RibbonGroup>
      </div>

      <div className="body">
        {/* **Three readings of one document, in one column.** The navigator lists what this
            flowsheet *has*, the palette what can be added, and the boundary what it is between —
            and each selects or commands through the one envelope the app holds. */}
        <aside className="side">
          {envelope === null ? null : (
            <Navigation
              catalogue={catalogue}
              envelope={envelope}
              selected={selected}
              onSelect={setSelected}
            />
          )}
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
          {envelope === null ? null : (
            <BoundaryPanel envelope={envelope} onCommand={send} />
          )}
        </aside>

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
          {envelope === null ? null : (
            <PropertyView
              catalogue={catalogue}
              envelope={envelope}
              node={node}
              edge={edge}
              onCommand={send}
              onSelect={setSelected}
            />
          )}
        </aside>
      </div>

      {/* **The two panels read one envelope from two sides**: the workbook is what the run reached,
          the messages are what the document is. Neither holds a copy of either. */}
      <Dock
        tabs={[
          { id: "workbook", label: "Workbook" },
          {
            id: "messages",
            label: "Messages",
            count: envelope?.diagnostics.length ?? 0,
          },
        ]}
        active={docked}
        onTab={(id) => setDocked(id === "workbook" ? "workbook" : "messages")}
        size={dockSize}
        onSize={setDockSize}
      >
        {envelope === null ? (
          <p className="note">no document</p>
        ) : docked === "workbook" ? (
          <Workbook envelope={envelope} onSelect={setSelected} />
        ) : (
          <DiagnosticsPanel
            diagnostics={envelope.diagnostics}
            runError={envelope.run_error}
            onSelect={(target) => setSelected(targetNodeId(target, envelope))}
          />
        )}
      </Dock>

      <StatusBar envelope={envelope} />
    </div>
  );
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
