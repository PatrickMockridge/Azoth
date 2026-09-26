/**
 * The palette: the entries a flowsheet can be built from.
 *
 * **Grouped by the family the entry was filed under**, which is the directory its spec lives in —
 * `two_port`, `separator`, `column`, and so on — and which the *library* names, in `Form.family`.
 *
 * **The grouping is not read off the id, and the difference is a defect this file used to have.**
 * Every id is `unit_ops.<leaf>`, so a split on `.` takes the leaf for a family: all twenty-nine
 * entries came out in one group called `other`. `Form.source` cannot supply it either — that is
 * NeqSim's taxonomy, where `cooler` is a `two_port` entry inside NeqSim's `heatexchanger/`
 * directory.
 *
 * An entry the executor refuses is listed and not offered: it is part of the palette, its refusal
 * is a measurement rather than an absence, and hiding it would make the palette look smaller than
 * the declaration it is generated from.
 *
 * **The column is the app's and this is a section of it.** The navigator sits above, this in the
 * middle and the boundary below, and all three are one `<aside>`: they are the readings of a
 * document that are neither the drawing nor the selected object's own window.
 */

import { familyLabel, familyOfForm } from "../state/navigation";
import type { Form } from "../wire/types";

export interface PaletteProps {
  forms: Form[];
  onAdd: (unit: string) => void;
}

export function Palette({ forms, onAdd }: PaletteProps) {
  const families = new Map<string, Form[]>();
  for (const form of forms) {
    const key = familyOfForm(form);
    families.set(key, [...(families.get(key) ?? []), form]);
  }

  return (
    <>
      <h2>Palette</h2>
      {[...families.entries()].map(([name, entries]) => (
        <div className="group" key={name}>
          <div className="label">{familyLabel(name)}</div>
          {entries.map((form) => (
            <button
              key={form.id}
              type="button"
              className={`entry${form.runnable ? "" : " refused"}`}
              title={form.refusal ?? form.source ?? form.id}
              disabled={!form.runnable}
              onClick={() => onAdd(form.id)}
            >
              <span className="name">{form.name}</span>
              <span className="tag">{form.parameters.length}</span>
            </button>
          ))}
        </div>
      ))}
    </>
  );
}
