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

import { useCallback, useEffect, useMemo, useRef, useState } from "react";
import "@xyflow/react/dist/style.css";

import blankDocument from "../../specs/flowsheets/blank.toml?raw";
import demoDocument from "../../specs/flowsheets/examples/two_feed_mixer.toml?raw";
import { BoundaryPanel } from "./components/BoundaryPanel";
import { CommandIcon } from "./components/CommandIcon";
import { CommandPalette } from "./components/CommandPalette";
import { menusOf, type CommandSpec, type MenuSpec } from "./commands/model";
import { CalculationsApp } from "./components/calculations/CalculationsApp";
import type { Workspace } from "./components/calculations/WorkspaceSwitch";
import { DiagnosticsPanel } from "./components/DiagnosticsPanel";
import { Dock } from "./components/dock/Dock";
import { Flowsheet } from "./components/Flowsheet";
import { NodeOverview } from "./components/NodeOverview";
import { Help } from "./components/Help";
import { Navigation } from "./components/Navigation";
import { Palette } from "./components/Palette";
import { PropertyView } from "./components/property/PropertyView";
import { ReportView } from "./components/report/ReportView";
import { RibbonGroup } from "./components/RibbonGroup";
import { WindowChrome } from "./components/shell/WindowChrome";
import { StarterHints } from "./components/StarterHints";
import { StatusBar } from "./components/StatusBar";
import { TitleStrip } from "./components/TitleStrip";
import { Welcome } from "./components/Welcome";
import { Workbook } from "./components/workbook/Workbook";
import { EXAMPLES } from "./examples";
import { manualHome } from "./manual";
import { saveDocument } from "./save";
import { useWindowManager } from "./hooks/useWindowManager";
import { selectedEdge, selectedNode } from "./state/selection";
import {
  applyAppearance,
  isPreview,
  readAppearance,
  type Appearance,
} from "./state/appearance";
import { readUnitSet, rememberUnitSet, unitsOf } from "./state/units";
import { readWelcome, rememberWelcome } from "./state/welcome";
import { door } from "./wire/door";
import type { Door, Opened, Session, SheetDoor } from "./wire/session";
import type { Catalogue, Command, Envelope, ExecutionOrder } from "./wire/types";

/**
 * Why a served door has no New, Open or Demo.
 *
 * In the tool bar as a `title` rather than in a paragraph somewhere: the control is where a person
 * asks the question, so the control is where the answer belongs.
 */
const HOSTED =
  "azoth serve holds one document for the life of its process — restart it with --flowsheet F to serve another";

/** One document the editor has stood at, for undo and redo. */
interface Snapshot {
  document: string;
  selected: string | null;
}

/** The file stem for a case: its name, or its id when the name has nothing usable. */
function fileStem(name: string, id: string): string {
  const stem = name.trim().replace(/[^A-Za-z0-9._-]+/g, "_");
  return stem === "" ? id : stem;
}

export function App() {
  const [session, setSession] = useState<Session | null>(null);
  const [catalogue, setCatalogue] = useState<Catalogue | null>(null);
  const [envelope, setEnvelope] = useState<Envelope | null>(null);
  const [selected, setSelected] = useState<string | null>(null);
  const [overview, setOverview] = useState<string | null>(null);
  const [fault, setFault] = useState<string | null>(null);
  const [openable, setOpenable] = useState(false);
  const [initialAppearance] = useState(readAppearance);
  const [previewed] = useState(isPreview);
  const [appearance, setAppearance] = useState<Appearance>(initialAppearance);
  /** Which of the two documents is in front; both stay mounted. */
  const [workspace, setWorkspace] = useState<Workspace>("flowsheet");
  /** The calculation book is loaded on first visit and then kept, so switching back keeps it. */
  const [calculationsSeen, setCalculationsSeen] = useState(false);
  const [sheetDoor, setSheetDoor] = useState<SheetDoor | null>(null);

  /**
   * The unit set a reader wants values in, which is remembered like the theme and for the same
   * reason: it is a reader's choice and not a document fact, so reopening a flowsheet finds the
   * units the person looking at it prefers.
   */
  const [unitSet, setUnitSet] = useState<string | null>(readUnitSet);

  /**
   * The dock's view state, which is the editor's and not the document's.
   *
   * **Messages is the tab it opens on**, because a document that does not check is the state a
   * person most needs told about and the dock is where the telling happens.
   */
  const [docked, setDocked] = useState<"workbook" | "messages" | "help">("workbook");
  const dockTouched = useRef(false);
  const [dockSize, setDockSize] = useState(200);

  /**
   * The first-run card, remembered like the theme and the unit set.
   *
   * **Shown over the demo, not instead of it.** A newcomer can press Solve without reading
   * anything, and the card is one click from gone.
   */
  const [welcome, setWelcome] = useState(readWelcome);

  /** The document text as it was last opened or saved, so "unsaved" is a fact and not a guess. */
  const [savedDocument, setSavedDocument] = useState<string | null>(null);

  /**
   * **Undo is a history of documents, not a stack of inverse commands.** Every envelope carries
   * `flowsheet.document`, so a snapshot is the save format; a restore opens it through the same
   * door an example does. The strings are not a second live document and nothing derives from
   * them except the undo itself.
   */
  const [history, setHistory] = useState<Snapshot[]>([]);
  const [future, setFuture] = useState<Snapshot[]>([]);
  const lastDocument = useRef<string | null>(null);
  const restoring = useRef(false);

  /** The printable report, a view over the last run and not a second call. */
  const [reporting, setReporting] = useState(false);
  /** The command palette is the Cyberpunk shell's; it is closed by any other layout. */
  const [paletteOpen, setPaletteOpen] = useState(false);
  /** The Cyberpunk HUD's visible rails; the Mac and XP shells ignore the attributes. */
  const [hudState, setHudState] = useState<Record<string, boolean>>({
    palette: true,
    inspector: true,
    workbook: true,
    commands: true,
  });
  const toggleHud = useCallback((id: string) => {
    setHudState((current) => ({ ...current, [id]: current[id] === false }));
  }, []);
  const canReport =
    envelope !== null &&
    envelope.ok &&
    envelope.run_error === null &&
    !envelope.dirty &&
    envelope.session !== null;

  /**
   * **A refused document is not run.** The button used to attempt the physics and return the
   * executor's own sentence beside the checker's list; the Messages panel is where a fix starts.
   */
  const canSolve = session !== null && envelope !== null && envelope.ok;
  const canUndo = openable && history.length > 0;
  const canRedo = openable && future.length > 0;

  // **A finding opens the dock that says what it is.** A newcomer's first edit can leave the
  // document incomplete; opening on Workbook hides the sentences that explain why. Once the
  // reader has chosen a tab, their choice wins.
  useEffect(() => {
    if (!dockTouched.current && envelope !== null && envelope.diagnostics.length > 0) {
      setDocked("messages");
    }
  }, [envelope]);

  const middleware = useRef<Door | null>(null);
  const openInput = useRef<HTMLInputElement>(null);
  const appRoot = useRef<HTMLDivElement>(null);

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
        setSavedDocument(opened.envelope.flowsheet.document);
        // A newly opened case has no edits behind it.
        setHistory([]);
        setFuture([]);
        lastDocument.current = opened.envelope.flowsheet.document;
        restoring.current = false;
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

  /** Record the document an edit is about to leave, so Undo has somewhere to return. */
  useEffect(() => {
    if (envelope === null) {
      return;
    }
    const document = envelope.flowsheet.document;
    if (lastDocument.current === null) {
      lastDocument.current = document;
      return;
    }
    if (document === lastDocument.current) {
      return;
    }
    if (restoring.current) {
      restoring.current = false;
      lastDocument.current = document;
      return;
    }
    const previous = lastDocument.current;
    lastDocument.current = document;
    setHistory((entries) => [...entries.slice(-49), { document: previous, selected: null }]);
    setFuture([]);
  }, [envelope]);

  /** Restore a snapshot through the door, without resetting the history it came from. */
  const restore = useCallback((snapshot: Snapshot) => {
    const loaded = middleware.current;
    if (loaded === null || loaded.kind !== "openable") {
      return;
    }
    restoring.current = true;
    void loaded.open(snapshot.document).then(
      (opened) => {
        setSession(opened.session);
        setEnvelope(opened.envelope);
        setSelected(snapshot.selected);
        setFault(null);
      },
      (error: unknown) => {
        restoring.current = false;
        setFault(String(error));
      },
    );
  }, []);

  const undo = useCallback(() => {
    if (history.length === 0 || envelope === null) {
      return;
    }
    const previous = history[history.length - 1];
    if (previous === undefined) {
      return;
    }
    setFuture((entries) => [
      { document: envelope.flowsheet.document, selected },
      ...entries.slice(0, 49),
    ]);
    setHistory(history.slice(0, -1));
    restore(previous);
  }, [history, envelope, selected, restore]);

  const redo = useCallback(() => {
    if (future.length === 0 || envelope === null) {
      return;
    }
    const next = future[0];
    if (next === undefined) {
      return;
    }
    setHistory((entries) => [
      ...entries.slice(-49),
      { document: envelope.flowsheet.document, selected },
    ]);
    setFuture(future.slice(1));
    restore(next);
  }, [future, envelope, selected, restore]);

  // **The keyboard path is the conventional one**, and it stays out of a text field so the
  // browser's own undo is not stolen from a draft.
  useEffect(() => {
    const onKey = (event: KeyboardEvent) => {
      if (!(event.metaKey || event.ctrlKey) || event.key.toLowerCase() !== "z") {
        return;
      }
      const target = event.target;
      if (
        target instanceof HTMLElement &&
        target.closest("input, textarea, select, [contenteditable='true']") !== null
      ) {
        return;
      }
      event.preventDefault();
      if (event.shiftKey) {
        redo();
      } else {
        undo();
      }
    };
    window.addEventListener("keydown", onKey);
    return () => window.removeEventListener("keydown", onKey);
  }, [undo, redo]);

  /**
   * The library's units and sets, resolved against the catalogue and the remembered choice.
   *
   * **Derived and never stored**: the catalogue arrives asynchronously on both doors, and a set
   * captured before it loaded would convert nothing for the rest of the session.
   */
  const units = useMemo(() => unitsOf(catalogue, unitSet), [catalogue, unitSet]);

  /** Choose a set: in force at once, and remembered for the next visit. */
  const chooseUnits = useCallback((id: string) => {
    setUnitSet(id);
    rememberUnitSet(id);
  }, []);

  /**
   * The appearance, in force and remembered. XP light is CSS's default; this only overrides it.
   *
   * **A preview link is shown but not saved** until the reader chooses something themselves, so a
   * shared `?layout=mac&theme=dark` does not silently rewrite a visitor's preference.
   */
  useEffect(() => {
    applyAppearance(appearance, {
      persist: !previewed || appearance !== initialAppearance,
    });
  }, [appearance, initialAppearance, previewed]);

  /**
   * The document has edits that were never saved.
   *
   * **Not the same fact as `dirty`,** which says the values are older than the document. A case
   * can be run and still unsaved, and both have to be said separately.
   */
  const unsaved =
    envelope !== null &&
    savedDocument !== null &&
    envelope.flowsheet.document !== savedDocument;

  /** Hand a document to the library: the Demo button, New, or a file. */
  const open = useCallback(
    (text: string) => {
      const loaded = middleware.current;
      if (loaded === null || loaded.kind !== "openable") {
        return;
      }
      if (
        unsaved &&
        !window.confirm("This case has unsaved changes. Open another and lose them?")
      ) {
        return;
      }
      land(loaded.open(text));
    },
    [land, unsaved],
  );

  /** Move between the flowsheet and the calculation book, loading the book on first visit. */
  const chooseWorkspace = useCallback((next: Workspace) => {
    setWorkspace(next);
    if (next === "calculations") {
      setCalculationsSeen(true);
    }
  }, []);

  /** The card has been dealt with: gone for this visit, and not shown again here. */
  const dismissWelcome = useCallback(() => {
    setWelcome(false);
    rememberWelcome();
  }, []);

  const openDemo = useCallback(() => {
    open(demoDocument);
    dismissWelcome();
  }, [open, dismissWelcome]);

  const openBlank = useCallback(() => {
    open(blankDocument);
    dismissWelcome();
  }, [open, dismissWelcome]);

  const openExample = useCallback(
    (document: string) => {
      open(document);
      dismissWelcome();
    },
    [open, dismissWelcome],
  );

  const showHelp = useCallback(() => {
    setDocked("help");
    dismissWelcome();
  }, [dismissWelcome]);

  useEffect(() => {
    let live = true;
    door(window.location.search)
      .then((loaded) => {
        if (!live) {
          return;
        }
        middleware.current = loaded.flowsheet;
        setSheetDoor(loaded.sheets);
        setOpenable(loaded.flowsheet.kind === "openable");
        void loaded.flowsheet.catalogue(true).then(
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
        land(
          loaded.flowsheet.kind === "openable"
            ? loaded.flowsheet.open(demoDocument)
            : loaded.flowsheet.attach(),
        );
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

  /** Run, and let the first-run card get out of the way; it told the reader to press this. */
  const run = useCallback(() => {
    dismissWelcome();
    solve();
  }, [dismissWelcome, solve]);

  /**
   * Add a palette entry, select it and dismiss the card.
   *
   * **A new node has to be visibly the thing that just happened.** It used to arrive unselected,
   * often under the welcome card, so a click on the palette could look like nothing.
   */
  const addUnit = useCallback(
    (unit: string) => {
      const id = nextId(unit, envelope);
      setSelected(`instance:${id}`);
      dismissWelcome();
      send({ command: "add_instance", id, unit, parameters: {} });
    },
    [envelope, dismissWelcome, send],
  );

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
    const document = envelope.flowsheet.document;
    const stem = fileStem(envelope.flowsheet.name, envelope.flowsheet.id);
    void saveDocument(`${stem}.toml`, document).then(
      () => setSavedDocument(document),
      (error: unknown) => setFault(String(error)),
    );
  }, [envelope]);

  /** Set the case's own name, which the report and the saved file take. */
  const setCaseName = useCallback(
    (name: string) => {
      if (name.trim() !== "") {
        send({ command: "set_name", name: name.trim() });
      }
    },
    [send],
  );

  useEffect(() => {
    if (!unsaved) {
      return;
    }
    const warn = (event: BeforeUnloadEvent) => {
      event.preventDefault();
      event.returnValue = "";
    };
    window.addEventListener("beforeunload", warn);
    return () => window.removeEventListener("beforeunload", warn);
  }, [unsaved]);

  const node = envelope === null ? null : selectedNode(envelope, selected);
  const edge = envelope === null ? null : selectedEdge(envelope, selected);
  const overviewNode = envelope === null ? null : selectedNode(envelope, overview);

  /** The Macintosh windows are live only while this workspace is front and the layout is Mac. */
  const { reset: resetWindows } = useWindowManager(appRoot, {
    active: appearance.layout === "mac" && workspace === "flowsheet",
  });

  /** Open the file picker the ribbon's Open button and the Case menu share. */
  const openFile = useCallback(
    (file: File) => {
      void file.text().then(open);
    },
    [open],
  );

  const commands: CommandSpec[] = [
    { id: "new", label: "New", menu: "case", run: openBlank, enabled: openable },
    {
      id: "open",
      label: "Open…",
      menu: "case",
      run: () => openInput.current?.click(),
      enabled: openable,
    },
    { id: "demo", label: "Demo", menu: "case", run: openDemo, enabled: openable },
    ...EXAMPLES.map((example) => ({
      id: `example-${example.file}`,
      label: `Open “${example.title}”`,
      menu: "case" as const,
      run: () => openExample(example.document),
      enabled: openable,
      keywords: "example case",
    })),
    {
      id: "save",
      label: "Save…",
      menu: "case",
      run: save,
      enabled: session !== null,
      shortcut: "Ctrl+S",
    },
    {
      id: "undo",
      label: "Undo",
      menu: "edit",
      run: undo,
      enabled: canUndo,
      shortcut: "Ctrl+Z",
    },
    {
      id: "redo",
      label: "Redo",
      menu: "edit",
      run: redo,
      enabled: canRedo,
      shortcut: "Ctrl+Y",
    },
    { id: "solve", label: "Solve", menu: "run", run, enabled: canSolve },
    {
      id: "report",
      label: "Report…",
      menu: "run",
      run: () => setReporting(true),
      enabled: canReport,
    },
    {
      id: "workspace-flowsheet",
      label: "Flowsheet",
      menu: "view",
      run: () => chooseWorkspace("flowsheet"),
      enabled: true,
    },
    {
      id: "workspace-calculations",
      label: "Calculations",
      menu: "view",
      run: () => chooseWorkspace("calculations"),
      enabled: true,
    },
    {
      id: "reset-windows",
      label: "Reset window positions",
      menu: "view",
      run: resetWindows,
      enabled: appearance.layout === "mac",
    },
    { id: "guide", label: "Open the guide", menu: "help", run: showHelp, enabled: true },
    {
      id: "manual",
      label: "User manual",
      menu: "help",
      run: () => {
        window.open(manualHome(), "_blank", "noopener");
      },
      enabled: true,
    },
  ];
  const menus: MenuSpec[] = menusOf(
    commands,
    { case: "Case", book: "Book", edit: "Edit", run: "Run", view: "View", help: "Help" },
    ["case", "edit", "run", "view", "help"],
  );

  useEffect(() => {
    if (appearance.layout !== "cyberpunk") {
      setPaletteOpen(false);
    }
  }, [appearance.layout]);

  useEffect(() => {
    const onKey = (event: KeyboardEvent) => {
      if (!(event.ctrlKey || event.metaKey) || event.key.toLowerCase() !== "k") {
        return;
      }
      if (appearance.layout !== "cyberpunk" || workspace !== "flowsheet") {
        return;
      }
      if (overview !== null || reporting) {
        return;
      }
      event.preventDefault();
      setPaletteOpen((open) => !open);
    };
    window.addEventListener("keydown", onKey);
    return () => window.removeEventListener("keydown", onKey);
  }, [appearance.layout, overview, reporting, workspace]);

  return (
    <>
    <div
      className={`app${workspace === "calculations" ? " off" : ""}`}
      hidden={workspace === "calculations"}
      ref={appRoot}
      data-hud-palette={hudState.palette ? "on" : "off"}
      data-hud-inspector={hudState.inspector ? "on" : "off"}
      data-hud-workbook={hudState.workbook ? "on" : "off"}
      data-hud-commands={hudState.commands ? "on" : "off"}
    >
      <input
        ref={openInput}
        className="sr-only"
        type="file"
        accept=".toml,text/plain"
        tabIndex={-1}
        aria-hidden="true"
        onChange={(event) => {
          const file = event.target.files?.[0];
          if (file !== undefined) {
            openFile(file);
          }
          event.target.value = "";
        }}
      />
      <TitleStrip
        envelope={envelope}
        session={session !== null}
        menus={menus}
        appearance={appearance}
        units={units}
        unsaved={unsaved}
        workspace={workspace}
        onWorkspace={chooseWorkspace}
        onName={setCaseName}
        onOrder={setOrder}
        onAppearance={setAppearance}
        onUnits={chooseUnits}
        hud={{
          controls: [
            { id: "palette", label: "Navigator", short: "Nav" },
            { id: "inspector", label: "Inspector", short: "Insp" },
            { id: "workbook", label: "Workbook", short: "WkBk" },
            { id: "commands", label: "Toolbar", short: "Deck" },
          ],
          state: hudState,
          onToggle: toggleHud,
          onCommands: () => setPaletteOpen(true),
        }}
      />

      <div className="ribbon" data-window="commands">
        <WindowChrome title="Tools">
        <RibbonGroup label="Run">
          <button
            type="button"
            className="primary"
            onClick={run}
            disabled={!canSolve}
            title={
              canSolve
                ? "run the physics; a check happens on every edit and a run is when you ask"
                : "fix the messages first; a refused document cannot run"
            }
          >
            <CommandIcon kind="solve" />
            Solve
          </button>
          <button
            type="button"
            onClick={undo}
            disabled={!canUndo}
            title={openable ? "undo the last edit" : "undo is not available on a served document"}
          >
            <CommandIcon kind="undo" />
            Undo
          </button>
          <button
            type="button"
            onClick={redo}
            disabled={!canRedo}
            title={openable ? "redo the last undone edit" : "redo is not available on a served document"}
          >
            <CommandIcon kind="redo" />
            Redo
          </button>
          <button
            type="button"
            onClick={() => setReporting(true)}
            disabled={!canReport}
            title={
              canReport
                ? "print or save the case: the drawing, the streams, the results and their source"
                : "solve the document first; a stale case is not a report"
            }
          >
            <CommandIcon kind="report" />
            Report
          </button>
        </RibbonGroup>
        <RibbonGroup label="Case">
          <button
            type="button"
            onClick={openBlank}
            disabled={!openable}
            title={openable ? "a document with nothing in it" : HOSTED}
          >
            <CommandIcon kind="new" />
            New
          </button>
          <button
            type="button"
            onClick={() => openInput.current?.click()}
            disabled={!openable}
            title={openable ? "open a TOML case" : HOSTED}
          >
            <CommandIcon kind="open" />
            Open
          </button>
          <button
            type="button"
            onClick={openDemo}
            disabled={!openable}
            title={openable ? "the shipped demo flowsheet" : HOSTED}
          >
            <CommandIcon kind="demo" />
            Demo
          </button>
          <select
            className="examples"
            value=""
            disabled={!openable}
            title={openable ? "a worked case to open" : HOSTED}
            aria-label="Load a worked example"
            onChange={(event) => {
              const chosen = EXAMPLES.find((example) => example.file === event.target.value);
              if (chosen !== undefined) {
                openExample(chosen.document);
              }
            }}
          >
            <option value="">Examples…</option>
            {EXAMPLES.map((example) => (
              <option key={example.file} value={example.file}>
                {example.title}
              </option>
            ))}
          </select>
          <button
            type="button"
            onClick={save}
            disabled={session === null}
            title="write the document where you say; it rides every envelope, so this is no call"
          >
            <CommandIcon kind="save" />
            Save
          </button>
        </RibbonGroup>
        </WindowChrome>
      </div>

      <div className="body">
        {/* **Three readings of one document, in one column.** The navigator lists what this
            flowsheet *has*, the palette what can be added, and the boundary what it is between —
            and each selects or commands through the one envelope the app holds. */}
        <aside className="side" data-window="palette">
          <WindowChrome title="Navigator">
          {envelope === null ? null : (
            <Navigation
              catalogue={catalogue}
              envelope={envelope}
              selected={selected}
              onSelect={setSelected}
            />
          )}
          {/* **The first step is the first panel.** A blank case starts with a feed, and the
              Boundary used to sit below the whole palette, off the bottom of the column. */}
          {envelope === null ? null : (
            <BoundaryPanel
              catalogue={catalogue}
              envelope={envelope}
              units={units}
              onCommand={send}
              onSelect={setSelected}
            />
          )}
          <Palette
            forms={catalogue?.unit_ops ?? []}
            onAdd={addUnit}
          />
          </WindowChrome>
        </aside>

        <div className="canvas" data-window="document">
          <WindowChrome title={envelope?.flowsheet.name ?? "Untitled"}>
          {envelope === null ? (
            <p className="note">{fault ?? "loading the engine…"}</p>
          ) : envelope.flowsheet.graph.nodes.length === 0 ? (
            <StarterHints />
          ) : (
            <Flowsheet
              catalogue={catalogue}
              envelope={envelope}
              selected={selected}
              units={units}
              onSelect={setSelected}
              onCommand={send}
              onOverview={setOverview}
            />
          )}
          {/* **Over the drawing, not instead of it.** The demo is already loaded behind the
              card, so the first action can be Solve whether or not the card is read. A served
              door holds one document and cannot open another, so it offers no card. */}
          {welcome && openable && envelope !== null ? (
            <Welcome
              examples={EXAMPLES}
              onDemo={openDemo}
              onBlank={openBlank}
              onExample={openExample}
              onHelp={showHelp}
              onDismiss={dismissWelcome}
            />
          ) : null}
          </WindowChrome>
        </div>

        <aside className="side right" data-window="inspector">
          <WindowChrome title={node?.id ?? edge?.id ?? "Inspector"}>
          {fault !== null ? <div className="fault">{fault}</div> : null}
          {envelope === null ? null : (
            <PropertyView
              catalogue={catalogue}
              envelope={envelope}
              node={node}
              edge={edge}
              units={units}
              onCommand={send}
              onSelect={setSelected}
            />
          )}
          </WindowChrome>
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
            count: (envelope?.diagnostics.length ?? 0) + (envelope?.run_error === null || envelope?.run_error === undefined ? 0 : 1),
          },
          { id: "help", label: "Help" },
        ]}
        active={docked}
        onTab={(id) => {
          dockTouched.current = true;
          setDocked(id === "workbook" ? "workbook" : id === "help" ? "help" : "messages");
        }}
        size={dockSize}
        onSize={setDockSize}
        managed={appearance.layout === "mac" || appearance.layout === "cyberpunk"}
      >
        <div role="tabpanel" id={`sheet-${docked}`}>
        {docked === "help" ? (
          <Help />
        ) : envelope === null ? (
          <p className="note">no document</p>
        ) : docked === "workbook" ? (
          <Workbook
            catalogue={catalogue}
            envelope={envelope}
            units={units}
            onSelect={setSelected}
          />
        ) : (
          <DiagnosticsPanel
            envelope={envelope}
            runError={envelope.run_error}
            onSelect={setSelected}
          />
        )}
        </div>
      </Dock>

      <StatusBar envelope={envelope} units={units} />

      {paletteOpen && appearance.layout === "cyberpunk" ? (
        <CommandPalette commands={commands} onClose={() => setPaletteOpen(false)} />
      ) : null}

      {overviewNode !== null && envelope !== null ? (
        <NodeOverview
          envelope={envelope}
          catalogue={catalogue}
          units={units}
          node={overviewNode}
          onClose={() => setOverview(null)}
          onEdit={() => {
            setSelected(overviewNode.id);
            setOverview(null);
          }}
        />
      ) : null}

      {reporting && envelope !== null ? (
        <ReportView
          envelope={envelope}
          catalogue={catalogue}
          units={units}
          onClose={() => setReporting(false)}
        />
      ) : null}
    </div>
    {calculationsSeen ? (
      <CalculationsApp
        door={sheetDoor}
        catalogue={catalogue}
        units={units}
        appearance={appearance}
        version={envelope?.provenance.version ?? "—"}
        workspace={workspace}
        hidden={workspace !== "calculations"}
        onWorkspace={chooseWorkspace}
        onAppearance={setAppearance}
        onUnits={chooseUnits}
      />
    ) : null}
    </>
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
