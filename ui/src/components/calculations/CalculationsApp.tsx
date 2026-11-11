/**
 * The Calculations workspace, as its own document and its own envelope.
 *
 * **Mounted beside the flowsheet and hidden when it is not in use.** Switching workspaces must not
 * throw the book away, so both are rendered and one is `display: none`; each still reads only its
 * own envelope.
 */

import { useCallback, useEffect, useRef, useState } from "react";

import blankBook from "../../../../specs/calcbooks/blank.toml?raw";
import demoBook from "../../../../specs/calcbooks/demo.toml?raw";
import { saveText } from "../../save";
import { declaredToSI, type Units } from "../../state/units";
import type { Appearance } from "../../state/appearance";
import { RibbonGroup } from "../RibbonGroup";
import { WindowChrome } from "../shell/WindowChrome";
import type { SheetDoor, SheetOpened, SheetSession } from "../../wire/session";
import type { SheetCommand } from "../../wire/sheetCommands";
import type { CalcForm, Catalogue, SheetEnvelope } from "../../wire/types";
import { CalcReportView } from "./CalcReportView";
import { blockedReason } from "./CalcSheet";
import { CalcTitleStrip } from "./CalcTitleStrip";
import { CalcsWorkspace } from "./CalcsWorkspace";
import type { Workspace } from "./WorkspaceSwitch";

/** One book the workspace has stood at, for undo and redo. */
interface Snapshot {
  document: string;
  active: string | null;
}

function fileStem(name: string, id: string): string {
  const stem = name.trim().replace(/[^A-Za-z0-9._-]+/g, "_");
  return stem === "" ? id : stem;
}

export function CalculationsApp({
  door,
  catalogue,
  units,
  appearance,
  version,
  workspace,
  hidden,
  onWorkspace,
  onAppearance,
  onUnits,
}: {
  door: SheetDoor | null;
  catalogue: Catalogue | null;
  units: Units;
  appearance: Appearance;
  version: string;
  workspace: Workspace;
  hidden: boolean;
  onWorkspace: (workspace: Workspace) => void;
  onAppearance: (appearance: Appearance) => void;
  onUnits: (id: string) => void;
}) {
  const [session, setSession] = useState<SheetSession | null>(null);
  const [envelope, setEnvelope] = useState<SheetEnvelope | null>(null);
  const [savedDocument, setSavedDocument] = useState<string | null>(null);
  const [active, setActive] = useState<string | null>(null);
  const [fault, setFault] = useState<string | null>(null);
  const [reporting, setReporting] = useState(false);
  const [history, setHistory] = useState<Snapshot[]>([]);
  const [future, setFuture] = useState<Snapshot[]>([]);

  const issued = useRef(0);
  const lastDocument = useRef<string | null>(null);
  const activeRef = useRef<string | null>(null);

  useEffect(() => {
    activeRef.current = active;
  }, [active]);

  const openable = door?.kind === "openable";
  const unsaved =
    envelope !== null &&
    savedDocument !== null &&
    envelope.book.document !== savedDocument;

  /** One call's answer, or the sentence that refused it. */
  const answer = useCallback((pending: Promise<SheetEnvelope>, record: boolean) => {
    const ticket = (issued.current += 1);
    void pending.then(
      (next) => {
        if (ticket !== issued.current) {
          return;
        }
        if (record && lastDocument.current !== null && lastDocument.current !== next.book.document) {
          // **The values are captured before the updater runs.** A functional `setHistory` may run
          // during a later render, after `lastDocument` has already advanced - which is how two
          // queued edits both recorded the same document and Undo restored the wrong one.
          const previous = lastDocument.current;
          const where = activeRef.current;
          setHistory((entries) =>
            [...entries, { document: previous, active: where }].slice(-50),
          );
          setFuture([]);
        }
        lastDocument.current = next.book.document;
        setEnvelope(next);
        setFault(null);
      },
      (error: unknown) => {
        if (ticket === issued.current) {
          setFault(String(error));
        }
      },
    );
  }, []);

  /** A book opened from a file, the shipped demo, or blank. */
  const land = useCallback((opening: Promise<SheetOpened>) => {
    const ticket = (issued.current += 1);
    void opening.then(
      (opened) => {
        if (ticket !== issued.current) {
          return;
        }
        setSession(opened.session);
        setEnvelope(opened.envelope);
        setSavedDocument(opened.envelope.book.document);
        lastDocument.current = opened.envelope.book.document;
        setHistory([]);
        setFuture([]);
        setActive(opened.envelope.book.sheets[0]?.id ?? null);
        setFault(null);
      },
      (error: unknown) => {
        if (ticket === issued.current) {
          setFault(String(error));
        }
      },
    );
  }, []);

  useEffect(() => {
    if (door === null) {
      return;
    }
    land(door.kind === "openable" ? door.open(demoBook) : door.attach());
  }, [door, land]);

  const send = useCallback(
    (command: SheetCommand) => {
      if (session !== null) {
        answer(session.apply(command), true);
      }
    },
    [session, answer],
  );

  const calculate = useCallback(
    (sheet: string | null) => {
      if (session !== null) {
        answer(session.calculate(sheet), false);
      }
    },
    [session, answer],
  );

  const addSheet = useCallback(
    (form: CalcForm) => {
      const id = nextSheetId(form, envelope);
      setActive(id);
      send({ command: "add_sheet", id, calculation: form.id, name: form.name });
    },
    [envelope, send],
  );

  const removeSheet = useCallback(
    (id: string) => {
      if (envelope !== null) {
        const rest = envelope.book.sheets.filter((sheet) => sheet.id !== id);
        setActive(rest[0]?.id ?? null);
      }
      send({ command: "remove_sheet", id });
    },
    [envelope, send],
  );

  const workedExample = useCallback(
    (sheet: string, form: CalcForm) => {
      if (session === null || envelope === null) {
        return;
      }
      const commands: SheetCommand[] = [];
      for (const [name, value] of Object.entries(form.worked_example)) {
        const input = form.inputs.find((entry) => entry.name === name);
        if (input === undefined) {
          continue;
        }
        if (typeof value === "number" && input.unit !== null) {
          commands.push({
            command: "set_input",
            sheet,
            input: name,
            value: declaredToSI(units, value, input.unit),
          });
        } else if (typeof value === "number" || typeof value === "boolean" || typeof value === "string") {
          commands.push({ command: "set_input", sheet, input: name, value });
        }
      }
      if (commands.length === 0) {
        return;
      }
      const pending = commands.reduce<Promise<SheetEnvelope>>(
        (previous, command) => previous.then(() => session.apply(command)),
        Promise.resolve(envelope),
      );
      answer(pending, true);
    },
    [session, envelope, units, answer],
  );

  const save = useCallback(() => {
    if (envelope === null) {
      return;
    }
    const stem = fileStem(envelope.book.name, envelope.book.id);
    void saveText(`${stem}.toml`, envelope.book.document).then(
      () => setSavedDocument(envelope.book.document),
      (error: unknown) => setFault(String(error)),
    );
  }, [envelope]);

  const open = useCallback(
    (text: string) => {
      if (door === null || door.kind !== "openable") {
        return;
      }
      if (
        unsaved &&
        !window.confirm("This calculation book has unsaved changes. Open another and lose them?")
      ) {
        return;
      }
      land(door.open(text));
    },
    [door, land, unsaved],
  );

  /** Restore a snapshot through the same open path a file uses. */
  const restore = useCallback(
    (snapshot: Snapshot) => {
      if (door === null || door.kind !== "openable") {
        return;
      }
      const ticket = (issued.current += 1);
      void door.open(snapshot.document).then(
        (opened) => {
          if (ticket !== issued.current) {
            return;
          }
          setSession(opened.session);
          setEnvelope(opened.envelope);
          lastDocument.current = snapshot.document;
          const exists = opened.envelope.book.sheets.some((sheet) => sheet.id === snapshot.active);
          setActive(exists ? snapshot.active : (opened.envelope.book.sheets[0]?.id ?? null));
          setFault(null);
        },
        (error: unknown) => {
          if (ticket === issued.current) {
            setFault(String(error));
          }
        },
      );
    },
    [door],
  );

  const undo = useCallback(() => {
    if (!openable || history.length === 0) {
      return;
    }
    const snapshot = history[history.length - 1];
    if (snapshot === undefined) {
      return;
    }
    const current = lastDocument.current ?? "";
    const where = activeRef.current;
    setHistory((entries) => entries.slice(0, -1));
    setFuture((entries) => [...entries, { document: current, active: where }].slice(-50));
    restore(snapshot);
  }, [history, openable, restore]);

  const redo = useCallback(() => {
    if (!openable || future.length === 0) {
      return;
    }
    const snapshot = future[future.length - 1];
    if (snapshot === undefined) {
      return;
    }
    const current = lastDocument.current ?? "";
    const where = activeRef.current;
    setFuture((entries) => entries.slice(0, -1));
    setHistory((entries) => [...entries, { document: current, active: where }].slice(-50));
    restore(snapshot);
  }, [future, openable, restore]);

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

  const sheets = envelope?.book.sheets ?? [];
  const current = sheets.find((sheet) => sheet.id === active) ?? sheets[0];
  const blocked = current === undefined ? null : blockedReason(current);
  const canReport = envelope !== null && envelope.ok && !envelope.dirty;
  const canUndo = openable && history.length > 0;
  const canRedo = openable && future.length > 0;

  return (
    <div className={`app calculations${hidden ? " off" : ""}`} hidden={hidden}>
      <CalcTitleStrip
        envelope={envelope}
        workspace={workspace}
        appearance={appearance}
        units={units}
        unsaved={unsaved}
        version={version}
        onName={(name) => send({ command: "set_book_name", name })}
        onWorkspace={onWorkspace}
        onAppearance={onAppearance}
        onUnits={onUnits}
      />

      <div className="ribbon" data-window="commands">
        <WindowChrome title="Tools">
        <RibbonGroup
          label="Run"
          note={
            blocked === null ? undefined : (
              <span id="ribbon-calculate-reason">{blocked}</span>
            )
          }
        >
          <button
            type="button"
            className="primary"
            onClick={() => {
              if (current !== undefined) {
                calculate(current.id);
              }
            }}
            disabled={current === undefined || !current.ready}
            aria-describedby={blocked === null ? undefined : "ribbon-calculate-reason"}
            title={
              current === undefined
                ? "add a calculation sheet first"
                : blocked ?? "run this sheet"
            }
          >
            Calculate
          </button>
          <button
            type="button"
            onClick={() => calculate(null)}
            disabled={!envelope?.dirty}
            title="run every sheet that is ready and older than its inputs"
          >
            Calculate all
          </button>
          <button
            type="button"
            onClick={undo}
            disabled={!canUndo}
            title={openable ? "undo the last edit" : "undo is not available on a served book"}
          >
            Undo
          </button>
          <button
            type="button"
            onClick={redo}
            disabled={!canRedo}
            title={openable ? "redo the last undone edit" : "redo is not available on a served book"}
          >
            Redo
          </button>
          <button
            type="button"
            onClick={() => setReporting(true)}
            disabled={!canReport}
            title={
              canReport
                ? "print the book: inputs, answers, warnings and the model check"
                : "calculate every sheet first; a stale book is not a report"
            }
          >
            Report
          </button>
        </RibbonGroup>
        <RibbonGroup label="Book">
          <button
            type="button"
            onClick={() => open(blankBook)}
            disabled={!openable}
            title={openable ? "a book with no sheets" : "a served book cannot be replaced"}
          >
            New
          </button>
          <label className={`button${openable ? "" : " off"}`} title={openable ? undefined : "a served book cannot be replaced"}>
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
                  event.target.value = "";
                }}
              />
            ) : null}
          </label>
          <button
            type="button"
            onClick={() => open(demoBook)}
            disabled={!openable}
            title={openable ? "the shipped worked-example book" : "a served book cannot be replaced"}
          >
            Demo
          </button>
          <button
            type="button"
            onClick={save}
            disabled={session === null}
            title="write the book where you say; it rides every envelope, so this is no call"
          >
            Save
          </button>
        </RibbonGroup>
        </WindowChrome>
      </div>

      <div className="calc-body">
        {fault === null ? null : <div className="fault">{fault}</div>}
        <CalcsWorkspace
          forms={catalogue?.calcs ?? []}
          envelope={envelope}
          units={units}
          active={current?.id ?? null}
          onAdd={addSheet}
          onSelect={setActive}
          onCommand={send}
          onCalculate={calculate}
          onWorkedExample={workedExample}
          onName={(id, name) => send({ command: "rename_sheet", id, name })}
          onRemove={removeSheet}
        />
      </div>

      {reporting && envelope !== null && catalogue !== null ? (
        <CalcReportView
          envelope={envelope}
          catalogue={catalogue.calcs ?? []}
          units={units}
          onClose={() => setReporting(false)}
        />
      ) : null}
    </div>
  );
}

function nextSheetId(form: CalcForm, envelope: SheetEnvelope | null): string {
  const leaf = form.id.split(".")[1] ?? "sheet";
  const taken = new Set(envelope?.book.sheets.map((sheet) => sheet.id) ?? []);
  for (let index = 1; index < 100; index += 1) {
    const candidate = `${leaf}_${index}`;
    if (!taken.has(candidate)) {
      return candidate;
    }
  }
  return `${leaf}_100`;
}
