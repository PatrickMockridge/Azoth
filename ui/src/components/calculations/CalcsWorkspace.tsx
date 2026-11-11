/**
 * The Calculations workspace.
 *
 * **A workspace, not a panel.** The flowsheet stays loaded behind the switch; this surface owns
 * its own document (the calculation book) and its own envelope, and holds no part of either.
 */

import { DockTabs } from "../property/DockTabs";
import { WindowChrome } from "../shell/WindowChrome";
import { manualHome } from "../../manual";
import { CalcPicker } from "./CalcPicker";
import { CalcSheet } from "./CalcSheet";
import type { Units } from "../../state/units";
import type { SheetCommand } from "../../wire/sheetCommands";
import type { CalcForm, SheetEnvelope } from "../../wire/types";

export function CalcsWorkspace({
  forms,
  envelope,
  units,
  active,
  onAdd,
  onSelect,
  onCommand,
  onCalculate,
  onWorkedExample,
  onName,
  onRemove,
}: {
  forms: CalcForm[];
  envelope: SheetEnvelope | null;
  units: Units;
  active: string | null;
  onAdd: (form: CalcForm) => void;
  onSelect: (id: string) => void;
  onCommand: (command: SheetCommand) => void;
  onCalculate: (sheet: string) => void;
  onWorkedExample: (sheetId: string, form: CalcForm) => void;
  onName: (sheetId: string, name: string) => void;
  onRemove: (sheetId: string) => void;
}) {
  if (envelope === null) {
    return (
      <div className="calc-workspace">
        <CalcPicker forms={forms} onAdd={onAdd} />
        <p className="note">loading the calculation book…</p>
      </div>
    );
  }

  const sheets = envelope.book.sheets;
  const current = sheets.find((sheet) => sheet.id === active) ?? sheets[0];
  const form =
    current === undefined ? undefined : forms.find((entry) => entry.id === current.calculation);

  return (
    <div className="calc-workspace">
      <aside className="calc-side" data-window="picker">
        <WindowChrome title="Calculators">
          <CalcPicker forms={forms} onAdd={onAdd} />
        </WindowChrome>
      </aside>
      <section className="calc-main" data-window="worksheet">
        <WindowChrome title={current?.name ?? envelope.book.name}>
          {current === undefined ? (
            <CalcStarterHints />
          ) : (
            <>
              <DockTabs
                tabs={sheets.map((sheet) => ({ id: sheet.id, label: sheet.name }))}
                active={current.id}
                onTab={onSelect}
                label="Calculation sheets"
              />
              <div
                className="calc-sheet-panel"
                role="tabpanel"
                id={`sheet-${current.id}`}
                aria-labelledby={`tab-${current.id}`}
              >
                <CalcSheet
                  form={form}
                  sheet={current}
                  units={units}
                  onCommand={onCommand}
                  onCalculate={() => onCalculate(current.id)}
                  onWorkedExample={() => {
                    if (form !== undefined) {
                      onWorkedExample(current.id, form);
                    }
                  }}
                  onName={(name) => onName(current.id, name)}
                  onRemove={() => onRemove(current.id)}
                />
              </div>
            </>
          )}
        </WindowChrome>
      </section>
    </div>
  );
}

function CalcStarterHints() {
  return (
    <div className="calc-starter">
      <h2>Calculation sheets</h2>
      <p>
        Pick a calculation on the left to add a sheet. Each sheet holds one calculation's inputs
        and its last answer, and a book can hold as many sheets as a case needs.
      </p>
      <ol>
        <li>Choose a calculation, or a common one from the top of the list.</li>
        <li>Fill the inputs; the unit set at the top changes what the fields show.</li>
        <li>Press Calculate. A red sentence is the calculation refusing an input; a warning is a number that is valid but less certain.</li>
        <li>Save the book in TOML, and print or export a sheet when it needs to be a record.</li>
      </ol>
      <p>
        <a href={`${manualHome()}using/calculations.html`} target="_blank" rel="noreferrer">
          Read the calculation-sheet guide
        </a>{" "}
        for the full explanation of a sheet.
      </p>
    </div>
  );
}
