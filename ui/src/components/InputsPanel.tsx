/**
 * A boundary stream: the record a feed declares, or what a product reached.
 *
 * **`h` is the field a user does not write.** It is a state function of `T`, `P` and `z`, so the
 * kernel derives it and a written one is refused rather than ignored — which is why the feed's
 * form has four fields and the product's readout has five.
 */

import { formatQuantity } from "../wire/field";
import type { Envelope, GraphNode } from "../wire/types";

export interface InputsPanelProps {
  envelope: Envelope;
  node: GraphNode;
  onCommand: (command: { command: string } & Record<string, unknown>) => void;
}

export function InputsPanel({ envelope, node, onCommand }: InputsPanelProps) {
  const stream = envelope.session?.streams[node.data.name];
  const record = node.data.input;

  return (
    <>
      <h2>{node.role === "input" ? "Feed" : "Product"}</h2>
      <p className="note">{node.data.name}</p>

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
            <div className="hint">
              the fluid — change it with `add_input`, which replaces the record
            </div>
          </div>
          {(["n", "P", "T"] as const).map((field) => (
            <div className="field" key={field}>
              <div className="head">
                <span className="name">{field}</span>
                <span className="unit">
                  {field === "n" ? "mol/s" : field === "P" ? "Pa" : "K"}
                </span>
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

      {stream === undefined ? (
        <p className="note">no values yet — solve to see what this stream reached</p>
      ) : (
        <>
          <h2>What it reached</h2>
          {(
            [
              ["n", stream.n],
              ["P", stream.P],
              ["T", stream.T],
              ["h", stream.h],
            ] as const
          ).map(([field, quantity]) => (
            <div className="field" key={field}>
              <div className="head">
                <span className="name">{field}</span>
              </div>
              <div className="readout">
                {formatQuantity(quantity.magnitude_si, quantity.unit, 6)}
              </div>
            </div>
          ))}
        </>
      )}
    </>
  );
}
