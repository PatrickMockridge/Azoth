import { useState } from "react";

import type { EditorCommand } from "../../wire/commands";
import { displayUnitOf, inDisplayUnit, toDisplayUnit, type Units } from "../../state/units";
import { formatQuantity } from "../../wire/field";
import type { Catalogue, Envelope, GraphNode, StreamRecord } from "../../wire/types";
import { FeedEditor } from "../feed/FeedEditor";
import { draftFromInput } from "../../state/composition";

/**
 * A boundary stream's conditions: the record a feed declares, or what a product reached.
 *
 * **`h` is the field a user does not write.** It is a state function of `T`, `P` and `z`, so the
 * kernel derives it and a written one is refused rather than ignored — which is why the feed's
 * form has four fields and the product's readout has five.
 *
 * **A feed is edited through the same editor the boundary panel uses.** The numbers here are the
 * document's, shown in the unit set in force; changing a composition means replacing the record,
 * which is one `add_input`, so the Edit control opens that form rather than pretending a setter
 * per field exists.
 *
 * **A typed condition commits on blur or Enter.** Writing `340` used to be three commands and
 * two states the checker could see; it is now one.
 */
export function ConditionsSheet({
  catalogue,
  envelope,
  node,
  units,
  onCommand,
}: {
  catalogue: Catalogue | null;
  envelope: Envelope;
  node: GraphNode;
  /** The unit a reader wants the values in; a catalogue that has not loaded converts nothing. */
  units: Units;
  onCommand: (command: EditorCommand) => void;
}) {
  const [editing, setEditing] = useState(false);
  const record = node.data.input;
  const reached: StreamRecord | undefined = envelope.session?.streams[node.data.name];

  if (editing && record !== undefined) {
    return (
      <div className="sheet" id="sheet-conditions" role="tabpanel" aria-labelledby="tab-conditions">
        <FeedEditor
          catalogue={catalogue}
          units={units}
          initial={draftFromInput(node.data.name, record, units)}
          editing
          onApply={(prepared) => {
            onCommand({ command: "add_input", ...prepared });
            setEditing(false);
          }}
          onCancel={() => setEditing(false)}
        />
      </div>
    );
  }

  const write = (field: "n" | "P" | "T", value: number) => {
    if (record === undefined || !Number.isFinite(value)) {
      return;
    }
    const unit = field === "n" ? "mol/s" : field === "P" ? "Pa" : "K";
    const stored = toDisplayUnit(units, value, unit);
    onCommand({
      command: "add_input",
      name: node.data.name,
      components: record.components,
      n: field === "n" ? stored : record.n,
      z: record.z,
      P: field === "P" ? stored : record.P,
      T: field === "T" ? stored : record.T,
    });
  };

  return (
    <div className="sheet" id="sheet-conditions" role="tabpanel" aria-labelledby="tab-conditions">
      {record === undefined ? null : (
        <>
          <div className="feed-head">
            <span className="name">the feed</span>
            <button
              type="button"
              className="button"
              onClick={() => setEditing(true)}
              title="change the fluid, its fractions or its conditions; connections survive"
            >
              Edit feed
            </button>
          </div>
          <div className="field">
            <div className="head">
              <span className="name">components</span>
            </div>
            <input value={record.components.join(", ")} readOnly />
            <div className="hint">mole fractions: {record.z.join(", ")}</div>
          </div>
          {(["n", "P", "T"] as const).map((field) => {
            const unit = field === "n" ? "mol/s" : field === "P" ? "Pa" : "K";
            return (
              <DraftNumber
                key={field}
                label={field}
                unit={displayUnitOf(units, unit) ?? unit}
                value={inDisplayUnit(
                  units,
                  field === "n" ? record.n : field === "P" ? record.P : record.T,
                  unit,
                )}
                onCommit={(value) => write(field, value)}
              />
            );
          })}
        </>
      )}

      {reached === undefined ? (
        <p className="note">no values yet — solve to see what this stream reached</p>
      ) : (
        <>
          <h2>What it reached</h2>
          <table className="grid">
            <thead>
              <tr>
                <th scope="col">field</th>
                <th scope="col">value</th>
                <th scope="col">unit</th>
              </tr>
            </thead>
            <tbody>
              {(
                [
                  ["n", reached.n],
                  ["mass flow", reached.mass_flow],
                  ["molar mass", reached.molar_mass],
                  ["P", reached.P],
                  ["T", reached.T],
                  ["h", reached.h],
                ] as const
              ).map(([label, quantity]) =>
                quantity === undefined || quantity === null ? null : (
                  <tr key={label}>
                    <th scope="row">{label}</th>
                    <td data-num>
                      {formatQuantity(quantity.magnitude_si, quantity.unit, 6, units)}
                    </td>
                    <td className="unit">{shownUnit(units, quantity.unit)}</td>
                  </tr>
                ),
              )}
              {/* **The vapour fraction is the one field a stream's record cannot derive** — it is
                  what the flash said — so a run that reached one shows it and a run that did not
                  shows a dash rather than a zero. */}
              <tr>
                <th scope="row">vapour fraction</th>
                <td data-num>
                  {reached.vapour_fraction === null || reached.vapour_fraction === undefined
                    ? "—"
                    : reached.vapour_fraction.toFixed(6)}
                </td>
                <td className="unit">—</td>
              </tr>
            </tbody>
          </table>
        </>
      )}
    </div>
  );
}

/** One condition, edited as a draft and sent once when the typing stops. */
function DraftNumber({
  label,
  unit,
  value,
  onCommit,
}: {
  label: string;
  unit: string;
  value: number;
  onCommit: (value: number) => void;
}) {
  const [draft, setDraft] = useState<string | null>(null);
  const shown = draft ?? String(value);
  const commit = () => {
    if (draft !== null) {
      onCommit(Number(draft));
      setDraft(null);
    }
  };
  return (
    <div className="field">
      <div className="head">
        <span className="name">{label}</span>
        <span className="unit">{unit}</span>
      </div>
      <input
        value={shown}
        aria-label={label}
        inputMode="decimal"
        onChange={(event) => setDraft(event.target.value)}
        onBlur={commit}
        onKeyDown={(event) => {
          if (event.key === "Enter") {
            commit();
          } else if (event.key === "Escape") {
            setDraft(null);
          }
        }}
      />
    </div>
  );
}

/** The unit a cell is in, which is the set's and not the record's. */
function shownUnit(units: Units, unit: string): string {
  return displayUnitOf(units, unit) ?? unit;
}
