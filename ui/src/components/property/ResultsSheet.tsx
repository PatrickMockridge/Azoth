import { outputRowsOf, warningsOf } from "../../state/results";
import { runNote, type StatusWord } from "../../state/status";
import type { Units } from "../../state/units";
import type { CalcOutput, UnitResult } from "../../wire/types";

/**
 * What a unit operation reached that is on no outlet stream.
 *
 * **Two layers.** The top of the sheet is for a person: a plain label, the value and unit in
 * the unit set in force, and the library's own sentence about how the number is obtained -
 * often the equation behind it, which is what lets a reader cross-check the answer without
 * reading the kernel. It draws only fields the run published; a field that was never configured
 * is not a dash row. The collapsed block at the bottom is the same values in the model's own
 * field names and the library's units, for a reader or an agent comparing records.
 *
 * **Not a caveat sheet.** The verification status is a fact about the model, not about this
 * answer, and it lives on the Provenance tab. What stays here is the library's own warnings
 * under their own heading, because a warning is about a number on this sheet.
 */
export function ResultsSheet({
  result,
  units,
  outputs,
  status,
}: {
  result: UnitResult;
  /** The unit a reader wants the values in; a catalogue that has not loaded converts nothing. */
  units: Units;
  /** The model's declared outputs, which supply the unit of a bare float and the sentence. */
  outputs: readonly CalcOutput[];
  /** The document's run state, which decides whether these values are from this run. */
  status: StatusWord;
}) {
  const rows = outputRowsOf(result, units, outputs);
  const machineRows = outputRowsOf(result, undefined, outputs);
  const warnings = warningsOf(result);

  return (
    <div className="sheet" id="sheet-results" role="tabpanel" aria-labelledby="tab-results">
      <p className="note">{runNote(status, rows.length > 0)}</p>

      {rows.length === 0 ? (
        <p className="note">this unit operation published no result values</p>
      ) : (
        <table className="grid result-table">
          <thead>
            <tr>
              <th scope="col">output</th>
              <th scope="col">value</th>
              <th scope="col">unit</th>
            </tr>
          </thead>
          <tbody>
            {rows.map((row) => (
              <tr key={row.key}>
                <th scope="row">
                  <span className="result-label">{row.label}</span>
                  {row.description === null ? null : (
                    <p className="result-description">{row.description}</p>
                  )}
                </th>
                <td data-num={row.numeric ? "true" : undefined}>{row.text}</td>
                <td className="unit">{row.unit ?? ""}</td>
              </tr>
            ))}
          </tbody>
        </table>
      )}

      {machineRows.length === 0 ? null : (
        <details className="result-machine">
          <summary>Machine-readable fields</summary>
          <p className="note">
            The same values in the model&apos;s own field names and the library&apos;s units.
          </p>
          <table className="grid">
            <thead>
              <tr>
                <th scope="col">field</th>
                <th scope="col">value</th>
                <th scope="col">unit</th>
              </tr>
            </thead>
            <tbody>
              {machineRows.map((row) => (
                <tr key={row.key}>
                  <th scope="row">
                    <code>{row.key}</code>
                  </th>
                  <td data-num={row.numeric ? "true" : undefined}>{row.text}</td>
                  <td className="unit">{row.unit ?? ""}</td>
                </tr>
              ))}
            </tbody>
          </table>
        </details>
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
