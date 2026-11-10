/**
 * One calculation sheet: what it is for, what it needs, what it returns, and when to believe it.
 *
 * **The inputs are the unit-op form's fields.** `CalcInput` is structurally `FormParameter`, so
 * the same control chooses a number, a select, a switch or a text box and converts the chosen
 * unit set on the way in. A sheet holds SI magnitudes; the field shows the set in force. The row
 * layout is the sheet's, because a reader needs the explanation beside the control.
 *
 * **The calculation is explicit.** Editing marks the result stale; Calculate is the only call
 * that runs physics, because a flash or a column rating is not free.
 */

import { ParameterField } from "../ParameterField";
import { DraftField } from "../DraftField";
import type { Units } from "../../state/units";
import type { SheetCommand } from "../../wire/sheetCommands";
import { quantityParts } from "../../wire/field";
import { calcManualUrl } from "../../manual";
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
          <div className="calculation">
            {form === undefined ? null : <span className="calc-name">{form.name}</span>}
            <span className="calc-id">{sheet.calculation}</span>
            {form === undefined ? null : (
              <span className={`calc-badge ${form.provenance.verification}`}>
                {verificationLabel(form.provenance.verification)}
              </span>
            )}
          </div>
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
          <section className="calc-about">
            <p className="calc-summary">{form.summary}</p>
            <div className="calc-about-links">
              <a href={calcManualUrl(form.id)} target="_blank" rel="noreferrer">
                Read the full manual and validity notes
              </a>
            </div>
            {form.equation === null ? null : (
              <details className="calc-details">
                <summary>Equation</summary>
                <pre>{form.equation}</pre>
              </details>
            )}
            {form.references.length === 0 ? null : (
              <details className="calc-details">
                <summary>Source and references</summary>
                <ul>
                  {form.references.map((reference) => (
                    <li key={reference}>{reference}</li>
                  ))}
                </ul>
              </details>
            )}
          </section>

          <section className="calc-section">
            <h3>Inputs</h3>
            <div className="calc-inputs">
              {form.inputs.map((input) => (
                <ParameterField
                  key={input.name}
                  layout="row"
                  parameter={input}
                  value={sheet.inputs[input.name] ?? ""}
                  units={units}
                  onChange={(value) => write(input, value)}
                  action={
                    input.required ? null : (
                      <button
                        type="button"
                        className="clear"
                        title="leave this input unset and let the calculation use its own default"
                        onClick={() =>
                          onCommand({
                            command: "unset_input",
                            sheet: sheet.id,
                            input: input.name,
                          })
                        }
                      >
                        clear
                      </button>
                    )
                  }
                />
              ))}
            </div>
          </section>

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

          {form.assumptions.length === 0 ? null : (
            <section className="calc-section calc-validity">
              <h3>When this applies</h3>
              <p className="calc-section-note">
                These conditions are not checked at runtime. The result means what it says only
                when they hold; the checked limits beside each input are enforced.
              </p>
              <ul>
                {form.assumptions.map((assumption) => (
                  <li key={assumption}>{assumption}</li>
                ))}
              </ul>
            </section>
          )}

          <section className="calc-section">
            <h3>Results</h3>
            {sheet.result === null ? (
              <p className="note">
                {sheet.error === null ? "not calculated yet; press Calculate" : "not calculated"}
              </p>
            ) : null}
            <div className="calc-outputs">
              {form.outputs.map((output) => {
                const value = sheet.result?.outputs[output.name];
                return (
                  <div className="calc-output" key={output.name}>
                    <div className="calc-output-text">
                      <div className="head">
                        <span className="name">{output.name}</span>
                        {output.optional ? <span className="optional">optional</span> : null}
                      </div>
                      <p className="description">{output.description}</p>
                    </div>
                    <div className="calc-output-value">
                      <span className="value">
                        {sheet.result === null ? "—" : shown(value, output, units)}
                      </span>
                      <span className="unit">
                        {sheet.result === null ? "" : unitOf(value, output, units)}
                      </span>
                    </div>
                  </div>
                );
              })}
            </div>
            {sheet.result !== null && sheet.result.warnings.length > 0 ? (
              <div className="calc-warnings">
                <h4>Warnings</h4>
                {sheet.result.warnings.map((warning, index) => (
                  <p key={`${warning.code}:${index}`}>
                    <strong>{warning.code}</strong> {warning.message}
                  </p>
                ))}
              </div>
            ) : null}
          </section>

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
            <p>
              <a href={`${calcManualUrl(form.id)}#notes`} target="_blank" rel="noreferrer">
                Implementation notes in the manual
              </a>
            </p>
          </details>
        </>
      )}
    </article>
  );
}

/** The verification status as a reader's words rather than an enum's spelling. */
function verificationLabel(status: string): string {
  switch (status) {
    case "verified":
      return "verified";
    case "partially_verified":
      return "partially verified";
    case "source_needed":
      return "source needed";
    default:
      return "unverified";
  }
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
    return quantityParts(value.magnitude_si, value.unit, 6, units).value;
  }
  if (typeof value === "number" && output.kind === "quantity") {
    return output.unit === null || output.unit === "dimensionless"
      ? String(value)
      : quantityParts(value, output.unit, 6, units).value;
  }
  return String(value);
}

/** The unit a value is in, in the unit set in force, or `—` where it has none. */
function unitOf(value: SheetValue | undefined, output: CalcOutput, units: Units): string {
  if (value === null || value === undefined || Array.isArray(value)) {
    return output.unit === null || output.unit === "dimensionless" ? "—" : output.unit;
  }
  if (typeof value === "object") {
    return quantityParts(value.magnitude_si, value.unit, 6, units).unit;
  }
  if (typeof value === "number" && output.kind === "quantity") {
    return output.unit === null || output.unit === "dimensionless"
      ? "—"
      : quantityParts(value, output.unit, 6, units).unit;
  }
  return "—";
}
