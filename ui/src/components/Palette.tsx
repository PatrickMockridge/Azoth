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
 */

import type { Form } from "../wire/types";

export interface PaletteProps {
  forms: Form[];
  onAdd: (unit: string) => void;
  /**
   * What else belongs in this column.
   *
   * **The boundary sits under the palette** rather than in a column of its own: feeds and products
   * are what a flowsheet is *between*, and a fourth panel would be a fourth place to look for
   * something a user adds once. The palette owns the `<aside>`, so this is how it is handed one.
   */
  children?: React.ReactNode;
}

/**
 * The family an entry is filed under, or `unfiled` where the catalogue did not say.
 *
 * **`unfiled` is a visible group rather than a silent merge into `other`**: a palette entry that
 * arrives without its family is a wire that stopped carrying it, and a person should see that as
 * its own heading rather than as one more entry in a group that looks like a family.
 */
function family(form: Form): string {
  return form.family ?? "unfiled";
}

export function Palette({ forms, onAdd, children }: PaletteProps) {
  const families = new Map<string, Form[]>();
  for (const form of forms) {
    const key = family(form);
    families.set(key, [...(families.get(key) ?? []), form]);
  }

  return (
    <aside className="side">
      <h2>Palette</h2>
      {[...families.entries()].map(([name, entries]) => (
        <div className="group" key={name}>
          <div className="label">{name.replace("_", " ")}</div>
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
      {children}
    </aside>
  );
}
