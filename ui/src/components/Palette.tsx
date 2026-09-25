/**
 * The palette: the entries a flowsheet can be built from.
 *
 * **Grouped by the entry's own family**, which is the directory the spec lives in — `two_port`,
 * `separator`, `column`, and so on. The grouping is read off the id rather than declared a second
 * time, so a palette entry cannot arrive in a family the palette does not know about.
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

/** The family a palette id names, e.g. `unit_ops.two_port.pump` → `two_port`. */
function family(id: string): string {
  const parts = id.split(".");
  return parts.length > 2 ? (parts[1] ?? "other") : "other";
}

export function Palette({ forms, onAdd, children }: PaletteProps) {
  const families = new Map<string, Form[]>();
  for (const form of forms) {
    const key = family(form.id);
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
