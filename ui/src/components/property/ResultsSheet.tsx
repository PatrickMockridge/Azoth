import { scalarsOf, warningsOf, type Scalar } from "../../state/results";
import type { Units } from "../../state/units";
import type { UnitResult } from "../../wire/types";

/**
 * Everything a unit operation reached that is on no outlet stream.
 *
 * **The labels are the model's own field names**, which is the point: a result is the operation's
 * registered record, so `condenser_duty` and `reboiler_duty` appear under the names the Python half
 * publishes them under and the names its case file declares. Nothing here writes a second label for
 * a number the library already named.
 *
 * The sheet is therefore the same one for a heater and for a distillation column, and a field added
 * to a model appears the next time that operation runs.
 */
export function ResultsSheet({
  result,
  units,
}: {
  result: UnitResult;
  /** The unit a reader wants the values in; a catalogue that has not loaded converts nothing. */
  units: Units;
}) {
  const scalars = scalarsOf(result, units);
  const warnings = warningsOf(result);

  return (
    <div className="sheet" id="sheet-results" role="tabpanel" aria-labelledby="tab-results">
      {scalars.length === 0 ? (
        <p className="note">this unit operation published no scalars</p>
      ) : (
        <table className="grid">
          <thead>
            <tr>
              <th scope="col">field</th>
              <th scope="col">value</th>
              <th scope="col">unit</th>
            </tr>
          </thead>
          <tbody>
            {scalars.map((scalar: Scalar) => (
              <tr key={scalar.key}>
                <th scope="row">{scalar.key}</th>
                <td data-num={scalar.numeric ? "true" : undefined}>{scalar.text}</td>
                <td className="unit">{scalar.unit ?? "—"}</td>
              </tr>
            ))}
          </tbody>
        </table>
      )}

      {/* **The library's own sentences, and not a summary of them.** A column that hit its
          iteration cap says so in its own words, with the code beside it - and under a class of its
          own, because a kernel's caveat is not a checker's diagnostic and the messages dock is
          where those live. */}
      {warnings.length === 0 ? null : (
        <>
          <h2>Caveats</h2>
          {warnings.map((warning, index) => (
            <div className="caveat" key={`${warning.code}-${String(index)}`}>
              <span className="code">{warning.code}</span>
              {warning.field === null ? null : <span className="where"> {warning.field}</span>}
              <p>{warning.message}</p>
            </div>
          ))}
        </>
      )}
    </div>
  );
}
