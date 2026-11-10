/**
 * The calculation picker.
 *
 * **A catalogue, grouped the way the registry is.** The four calculations with table-shaped inputs
 * are shown and disabled with the reason visible: a picker that hid them would leave a reader
 * wondering whether the calculation exists at all.
 */

import { useState } from "react";

import type { CalcForm } from "../../wire/types";

const NAMESPACES: { id: string; label: string }[] = [
  { id: "hydraulics", label: "Hydraulics" },
  { id: "thermal", label: "Heat transfer" },
  { id: "eos", label: "Equations of state and properties" },
  { id: "characterization", label: "Characterization" },
  { id: "reactions", label: "Reactions" },
];

/** The calculations a first-time reader is most likely to want. */
const COMMON = [
  "hydraulics.pump_power",
  "hydraulics.darcy_weisbach",
  "hydraulics.reynolds_number",
  "thermal.lmtd_correction_factor",
  "eos.pr_z_factor",
];

export function CalcPicker({
  forms,
  onAdd,
}: {
  forms: CalcForm[];
  onAdd: (form: CalcForm) => void;
}) {
  const [query, setQuery] = useState("");
  const needle = query.trim().toLowerCase();
  const matches = forms.filter(
    (form) =>
      needle === "" ||
      form.name.toLowerCase().includes(needle) ||
      form.id.toLowerCase().includes(needle),
  );
  const common = COMMON.map((id) => forms.find((form) => form.id === id)).filter(
    (form): form is CalcForm => form !== undefined,
  );

  return (
    <div className="calc-picker">
      <input
        className="calc-search"
        value={query}
        aria-label="Search calculations"
        placeholder="Search calculations"
        onChange={(event) => setQuery(event.target.value)}
      />
      <div className="calc-list">
        {needle === "" ? (
          <>
            <h3>Common</h3>
            {common.map((form) => (
              <Entry key={form.id} form={form} onAdd={onAdd} />
            ))}
          </>
        ) : null}
        {NAMESPACES.map((namespace) => {
          const entries = matches.filter((form) => form.namespace === namespace.id);
          return entries.length === 0 ? null : (
            <section key={namespace.id}>
              <h3>{namespace.label}</h3>
              {entries.map((form) => (
                <Entry key={form.id} form={form} onAdd={onAdd} />
              ))}
            </section>
          );
        })}
        {matches.length === 0 ? <p className="note">nothing matches “{query}”</p> : null}
      </div>
    </div>
  );
}

function Entry({ form, onAdd }: { form: CalcForm; onAdd: (form: CalcForm) => void }) {
  return (
    <div className="calc-entry">
      <button
        type="button"
        disabled={!form.available}
        title={form.available ? `add ${form.name} as a sheet` : form.refusal ?? "not available"}
        onClick={() => onAdd(form)}
      >
        {form.name}
      </button>
      {form.available ? null : <p className="refusal">{form.refusal}</p>}
    </div>
  );
}
