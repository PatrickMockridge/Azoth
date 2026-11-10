/**
 * One calculation sheet: inputs, one Calculate, and the answer.
 *
 * **The inputs are the unit-op form's fields.** `CalcInput` is structurally `FormParameter`, so
 * the same control chooses a number, a select, a switch or a text box and converts the chosen
 * unit set on the way in. A sheet holds SI magnitudes; the field shows the set in force.
 *
 * **The calculation is explicit.** Editing marks the result stale; Calculate is the only call
 * that runs physics, because a flash or a column rating is not free.
 */

import { ParameterField } from "../ParameterField";
import { DraftField } from "../DraftField";
import type { Units } from "../../state/units";
import type { SheetCommand } from "../../wire/sheetCommands";
import { formatQuantity } from "../../wire/field";
import type { CalcForm, CalcInput, CalcOutput, Sheet, SheetValue } from "../../wire/types";

export function CalcSheet({
  form,
  sheet,
  units,
  onCommand,
  onCalculate,
  onWorkedExample,
  onName,
  onRemove,
}: {
  form: CalcForm | undefined;
  sheet: Sheet;
  units: Units;
  onCommand: (command: SheetCommand) => void;
  onCalculate: () => void;
  onWorkedExample: () => void;
  onName: (name: string) => void;
  onRemove: () => void;
}) {
  const write = (input: CalcInput, value: unknown) => {
    if (value === "") {
      onCommand({ command: "unset_input", sheet: sheet.id, input: input.name });
      return;
    }
    if (input.kind === "quantity" && typeof value !== "number") {
      return;
    }
    if ((input.kind === "enum" || input.kind === "string") && typeof value !== "string") {
      return;
    }
    if (input.kind === "boolean" && typeof value !== "boolean") {
      return;
    }
    onCommand({
      command: "set_input",
      sheet: sheet.id,
      input: input.name,
      value: value as number | boolean | string,
    });
  };

  return (
    <article className="calc-sheet">
      <header className="calc-sheet-head">
        <div className="calc-sheet-name">
          <DraftField
            value={sheet.name}
            kind="text"
            label="Sheet name"
            onCommit={(value) => {
              if (typeof value === "string" && value.trim() !== "") {
                onName(value.trim());
              }
            }}
          />
          <span className="calculation">{sheet.calculation}</span>
        </div>
        <div className="calc-sheet-actions">
          {form?.available && form.worked_example !== undefined ? (
            <button type="button" onClick={onWorkedExample}>
              Load worked example
            </button>
          ) : null}
          <button
            type="button"
            className="primary"
            onClick={onCalculate}
            disabled={!sheet.ready}
            title={
              sheet.ready
                ? "run this calculation"
                : `fill the required inputs first: ${sheet.missing.join(", ")}`
            }
          >
            Calculate
          </button>
          <button type="button" onClick={onRemove}>
            Remove
          </button>
        </div>
      </header>

      {form === undefined ? null : (
        <>
          <div className="calc-inputs">
            {form.inputs.map((input) => (
              <div className="calc-input" key={input.name}>
                <ParameterField
                  parameter={input}
                  value={sheet.inputs[input.name] ?? ""}
                  units={units}
                  onChange={(value) => write(input, value)}
                />
                {input.required ? null : (
                  <button
                    type="button"
                    className="clear"
                    title="leave this input unset and let the calculation use its own default"
                    onClick={() =>
                      onCommand({ command: "unset_input", sheet: sheet.id, input: input.name })
                    }
                  >
                    clear
                  </button>
                )}
              </div>
            ))}
          </div>

          {sheet.diagnostics.length > 0 ? (
            <div className="calc-diagnostics">
              {sheet.diagnostics.map((diagnostic) => (
                <p key={`${diagnostic.code}:${diagnostic.field ?? ""}`}>{diagnostic.message}</p>
              ))}
            </div>
          ) : null}

          {sheet.error === null ? null : (
            <div className="fault">
              {sheet.error.message}
              {sheet.error.field === null ? null : (
                <span className="where"> · {sheet.error.field}</span>
              )}
            </div>
          )}

          {sheet.result === null && sheet.error === null ? (
            <p className="note">no answer yet; press Calculate</p>
          ) : null}

          {sheet.result === null ? null : (
            <>
              <table className="grid calc-results">
                <thead>
                  <tr>
                    <th scope="col">output</th>
                    <th scope="col">value</th>
                    <th scope="col">unit</th>
                  </tr>
                </thead>
                <tbody>
                  {form.outputs.map((output) => (
                    <tr key={output.name}>
                      <th scope="row" title={output.description}>
                        {output.name}
                      </th>
                      <td data-num>{shown(sheet.result?.outputs[output.name], output, units)}</td>
                      <td className="unit">{unitOf(sheet.result?.outputs[output.name], output)}</td>
                    </tr>
                  ))}
                </tbody>
              </table>
              {sheet.result.warnings.length === 0 ? null : (
                <div className="calc-warnings">
                  <h3>Warnings</h3>
                  {sheet.result.warnings.map((warning, index) => (
                    <p key={`${warning.code}:${index}`}>
                      <strong>{warning.code}</strong> {warning.message}
                    </p>
                  ))}
                </div>
              )}
            </>
          )}

          {form.equation === null ? null : (
            <details className="calc-details">
              <summary>Equation</summary>
              <pre>{form.equation}</pre>
            </details>
          )}
          {form.assumptions.length === 0 ? null : (
            <details className="calc-details">
              <summary>Assumptions</summary>
              <ul>
                {form.assumptions.map((assumption) => (
                  <li key={assumption}>{assumption}</li>
                ))}
              </ul>
            </details>
          )}
          <details className="calc-details">
            <summary>Model check</summary>
            <dl>
              <dt>model</dt>
              <dd>{form.id}</dd>
              <dt>verification</dt>
              <dd>{form.provenance.verification}</dd>
              <dt>spec hash</dt>
              <dd>{form.provenance.spec_sha256 ?? "—"}</dd>
              <dt>rust hash</dt>
              <dd>{form.provenance.rust_sha256 ?? "—"}</dd>
              <dt>tests</dt>
              <dd>
                {form.provenance.tests_active} active, {form.provenance.tests_skipped} skipped
              </dd>
            </dl>
          </details>
        </>
      )}
    </article>
  );
}

/** One output as a reader sees it, in the unit set in force. */
function shown(value: SheetValue | undefined, output: CalcOutput, units: Units): string {
  if (value === undefined || value === null) {
    return "—";
  }
  if (Array.isArray(value)) {
    return value.join(", ");
  }
  if (typeof value === "object") {
    return formatQuantity(value.magnitude_si, value.unit, 6, units);
  }
  if (typeof value === "number" && output.kind === "quantity") {
    return output.unit === null || output.unit === "dimensionless"
      ? String(value)
      : formatQuantity(value, output.unit, 6, units);
  }
  return String(value);
}

/** The unit a value is in, where it has one. */
function unitOf(value: SheetValue | undefined, output: CalcOutput): string {
  if (value === null || value === undefined) {
    return output.unit ?? "—";
  }
  if (typeof value === "object" && !Array.isArray(value)) {
    return value.unit;
  }
  return output.kind === "quantity" ? (output.unit ?? "—") : "—";
}
