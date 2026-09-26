import { componentAxis, fraction, molarMasses } from "../../state/columns";
import { displayOf, type Units } from "../../state/units";
import type { Catalogue, Envelope, GraphNode } from "../../wire/types";

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
  catalogue,
  envelope,
  node,
  units,
}: {
  catalogue: Catalogue | null;
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
  // **The databank's molar masses, and nothing derived from them that the library computes.**
  // A stream's record carries the mixture's molar mass and not each substance's, so the two mass
  // columns below are the reason the catalogue has a components region at all.
  const masses = molarMasses(catalogue);
  const weights = names.slice(0, z.length).map((name) => masses[name] ?? null);
  const total = weights.every((mass) => mass !== null)
    ? z.reduce((sum, value, index) => sum + value * (weights[index] as number), 0)
    : null;

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
              <th scope="col">mass fraction</th>
              <th scope="col">mass flow {massUnit(units)}</th>
            </tr>
          </thead>
          <tbody>
            {z.map((value, index) => (
              <tr key={names[index] ?? index}>
                <th scope="row">{names[index] ?? `z[${index}]`}</th>
                <td data-num>{fraction(value)}</td>
                <td data-num>{fraction(value * n)}</td>
                <td data-num>
                  {total === null || total === 0 || weights[index] === null
                    ? "—"
                    : fraction((value * (weights[index] as number)) / total)}
                </td>
                <td data-num>
                  {/* A dash where the databank cannot weigh the substance, rather than a zero:
                      the record carries the mixture's molar mass and not each one's. */}
                  {weights[index] === null
                    ? "—"
                    : fraction((value * n * (weights[index] as number)) / massFactor(units))}
                </td>
              </tr>
            ))}
          </tbody>
        </table>
      )}
      <p className="note">
        mole flow is `n · z`, mass flow is `n · z · M` — two fields of one record times the
        substance&apos;s molar mass from the databank, and mass fraction is the same numerator over
        the sum. A substance the databank cannot weigh shows a dash in both rather than a guess.
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

/** The unit the mass columns are read in, and the factor their cells are divided by. */
function massUnit(units: Units): string {
  return displayOf(units, "kg/s").unit;
}

function massFactor(units: Units): number {
  return displayOf(units, "kg/s").factor;
}
