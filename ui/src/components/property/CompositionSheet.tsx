import { componentAxis, fraction } from "../../state/columns";
import { displayOf, type Units } from "../../state/units";
import type { Envelope, GraphNode } from "../../wire/types";

/**
 * A stream's substances, one row each.
 *
 * **The names come from the document and the numbers from the run**, and where the two disagree
 * about the fluid the axis says so rather than picking a feed - see `state/columns.ts`.
 *
 * **One derived column, and it is labelled as one.** `n · z` is a mole flow, and it is two fields
 * of *one record* multiplied: no kernel's arithmetic is repeated here. A mass column would need
 * each substance's own molar mass, which no record carries, so it is absent rather than guessed.
 */
export function CompositionSheet({
  envelope,
  node,
  units,
}: {
  envelope: Envelope;
  node: GraphNode;
  /** The unit a reader wants the values in; a catalogue that has not loaded converts nothing. */
  units: Units;
}) {
  const axis = componentAxis(envelope);
  const record = node.data.input;
  const reached = envelope.session?.streams[node.data.name];
  const z = record?.z ?? reached?.z ?? [];
  // A feed's flow is stated in the document's own unit, which is the vocabulary's; a solved
  // stream's is the run's. Both are SI magnitudes by the time they reach here, so one division
  // covers the two - and `mol/s` is the dimension, not the feed's spelling.
  const n = (record?.n ?? reached?.n.magnitude_si ?? 0) / flowFactor(units);
  const names = axis.names.length >= z.length ? axis.names : z.map((_, index) => `z[${index}]`);

  return (
    <div className="sheet" id="sheet-composition" role="tabpanel" aria-labelledby="tab-composition">
      {axis.note === null ? null : <p className="note">{axis.note}</p>}
      {z.length === 0 ? (
        <p className="note">nothing to show — this stream has no composition yet</p>
      ) : (
        <table className="grid">
          <thead>
            <tr>
              <th scope="col">component</th>
              <th scope="col">mole fraction</th>
              <th scope="col">mole flow {flowUnit(units)}</th>
            </tr>
          </thead>
          <tbody>
            {z.map((value, index) => (
              <tr key={names[index] ?? index}>
                <th scope="row">{names[index] ?? `z[${index}]`}</th>
                <td data-num>{fraction(value)}</td>
                <td data-num>{fraction(value * n)}</td>
              </tr>
            ))}
          </tbody>
        </table>
      )}
      <p className="note">
        mole flow is `n · z`, two fields of one record — a mass column would need each substance&apos;s
        own molar mass, which the record does not carry
      </p>
    </div>
  );
}

/** The unit the derived mole-flow column is read in. */
function flowUnit(units: Units): string {
  return displayOf(units, "mol/s").unit;
}

/** The factor that column's cells are divided by, which is the same unit's. */
function flowFactor(units: Units): number {
  return displayOf(units, "mol/s").factor;
}
