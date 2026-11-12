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

import { useId } from "react";

import { ParameterField } from "../ParameterField";
import { DraftField } from "../DraftField";
import { engineeringUnit, type Units } from "../../state/units";
import { definitionOf } from "../../state/verification";
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
  const missingId = useId();
  const blocked = blockedReason(sheet);
  const staleResult = sheet.stale && sheet.result !== null;

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
    <article className="calc-sheet" data-stale={staleResult ? "true" : undefined}>
      <header className="calc-sheet-head">
        <div className="calc-sheet-name">
          <h2 className="calc-name">{form?.name ?? sheet.calculation}</h2>
          <div className="calculation">
            <span className="calc-id">{sheet.calculation}</span>
            {/* **The run, not the model's verification.** A badge beside the name says whether
                this sheet has an answer and whether it is from before the last edit; how far the
                model itself is checked is the provenance record below, where a reader looks for
                it. */}
            <span
              className={`calc-badge ${
                sheet.result === null
                  ? "not-calculated"
                  : sheet.stale
                    ? "out-of-date"
                    : "calculated"
              }`}
            >
              {sheet.result === null
                ? "not calculated"
                : sheet.stale
                  ? "out of date"
                  : "calculated"}
            </span>
          </div>
          <label className="calc-sheet-rename">
            <span className="calc-sheet-rename-label">Sheet name</span>
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
          </label>
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
            aria-describedby={blocked === null ? undefined : missingId}
            title={blocked ?? "run this calculation"}
          >
            Calculate
          </button>
          <button type="button" onClick={onRemove}>
            Remove
          </button>
        </div>
      </header>

      {blocked === null ? null : (
        <p className="calc-missing" id={missingId} role="status">
          {blocked}
        </p>
      )}

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
                        aria-label={`Clear ${input.name}`}
                        title="leave this input unset; the calculation uses its declared default or omits what needs it"
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
            {staleResult ? (
              <p className="calc-stale" role="status">
                These are the last calculated values; the inputs changed after the last Calculate.
              </p>
            ) : null}
            <div className="calc-outputs">
              {form.outputs.map((output) => {
                const value = sheet.result?.outputs[output.name];
                const unit = sheet.result === null ? "" : unitOf(value, output, units);
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
                      <span className="unit" title={unit === "" ? undefined : unit}>
                        {unit === "" ? "" : engineeringUnit(unit)}
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
              <dd>
                <code>{form.provenance.verification}</code>
                <p className="result-description">
                  {definitionOf(form.provenance.verification).definition}
                </p>
              </dd>
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

/**
 * Why Calculate is unavailable, in the words the sheet and the ribbon can both show.
 *
 * **`ready` is not "missing inputs".** A sheet can be unready because a value's kind is refused,
 * with `missing` empty; the fallback keeps that state explainable instead of rendering an empty
 * list after a colon.
 */
export function blockedReason(sheet: Sheet): string | null {
  if (sheet.ready) {
    return null;
  }
  if (sheet.missing.length > 0) {
    return `Required before Calculate: ${sheet.missing.join(", ")}.`;
  }
  return sheet.diagnostics[0]?.message ?? "This sheet cannot run until its inputs are valid.";
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
    const unit = quantityParts(value.magnitude_si, value.unit, 6, units).unit;
    return unit === "dimensionless" ? "—" : unit;
  }
  if (typeof value === "number" && output.kind === "quantity") {
    return output.unit === null || output.unit === "dimensionless"
      ? "—"
      : quantityParts(value, output.unit, 6, units).unit;
  }
  return "—";
}
