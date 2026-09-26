import type { EditorCommand } from "../../wire/commands";
import { displayOf, type Units } from "../../state/units";
import { formatQuantity } from "../../wire/field";
import type { Envelope, GraphNode, StreamRecord } from "../../wire/types";

/**
 * A boundary stream's conditions: the record a feed declares, or what a product reached.
 *
 * **`h` is the field a user does not write.** It is a state function of `T`, `P` and `z`, so the
 * kernel derives it and a written one is refused rather than ignored — which is why the feed's form
 * has four fields and the product's readout has five.
 *
 * **A feed that is also a product is one object with one record.** A feed's conditions come from
 * the document (`node.data.input`) and are editable through `add_input`, which replaces the whole
 * record; a product's come from the run. Where both exist - a product that is also consumed - the
 * document's record is what is drawn, because that is the one a command can change.
 */
export function ConditionsSheet({
  envelope,
  node,
  units,
  onCommand,
}: {
  envelope: Envelope;
  node: GraphNode;
  /** The unit a reader wants the values in; a catalogue that has not loaded converts nothing. */
  units: Units;
  onCommand: (command: EditorCommand) => void;
}) {
  const record = node.data.input;
  const reached: StreamRecord | undefined = envelope.session?.streams[node.data.name];

  return (
    <div className="sheet" id="sheet-conditions" role="tabpanel" aria-labelledby="tab-conditions">
      {record === undefined ? null : (
        <>
          <div className="field">
            <div className="head">
              <span className="name">components</span>
            </div>
            <input
              value={record.components.join(", ")}
              onChange={() => {
                // A substance list is not a parameter: changing it means a different fluid, and
                // the command model has `add_input` for the whole record rather than a setter per
                // field. The field is drawn read-only so that is visible rather than surprising.
              }}
              readOnly
            />
            <div className="hint">the fluid — change it with `add_input`, which replaces the record</div>
          </div>
          {(["n", "P", "T"] as const).map((field) => (
            <div className="field" key={field}>
              <div className="head">
                <span className="name">{field}</span>
                <span className="unit">{field === "n" ? "mol/s" : field === "P" ? "Pa" : "K"}</span>
              </div>
              <input
                value={String(record[field])}
                inputMode="decimal"
                onChange={(event) =>
                  onCommand({
                    command: "add_input",
                    name: node.data.name,
                    components: record.components,
                    n: field === "n" ? Number(event.target.value) : record.n,
                    z: record.z,
                    P: field === "P" ? Number(event.target.value) : record.P,
                    T: field === "T" ? Number(event.target.value) : record.T,
                  })
                }
              />
            </div>
          ))}
          <div className="field">
            <div className="head">
              <span className="name">z</span>
            </div>
            <input value={record.z.join(", ")} readOnly />
            <div className="hint">mole fractions, one per substance</div>
          </div>
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

/** The unit a cell is in, which is the set's and not the record's. */
function shownUnit(units: Units, unit: string): string {
  return displayOf(units, unit).unit;
}
