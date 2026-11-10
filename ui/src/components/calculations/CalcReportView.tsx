/**
 * A printable calculation book: one sheet per section.
 *
 * **A reading of the envelope, not a second run.** It prints the inputs, the answers, the
 * warnings and the model check; it does not call the library again.
 */

import { useEffect, useRef } from "react";

import { quantityParts } from "../../wire/field";
import { displayUnitOf, type Units } from "../../state/units";
import type {
  CalcForm,
  CalcInput,
  CalcOutput,
  SheetEnvelope,
  SheetValue,
} from "../../wire/types";

export function CalcReportView({
  envelope,
  catalogue,
  units,
  onClose,
}: {
  envelope: SheetEnvelope;
  catalogue: CalcForm[];
  units: Units;
  onClose: () => void;
}) {
  const close = useRef<HTMLButtonElement>(null);
  useEffect(() => {
    close.current?.focus();
    const escape = (event: KeyboardEvent) => {
      if (event.key === "Escape") {
        onClose();
      }
    };
    window.addEventListener("keydown", escape);
    return () => window.removeEventListener("keydown", escape);
  }, [onClose]);

  return (
    <div className="report" role="dialog" aria-modal="true" aria-label="Calculation book">
      <header className="report-head">
        <div>
          <h2>{envelope.book.name}</h2>
          <p className="report-summary">
            {envelope.book.sheets.length} sheet(s) · {envelope.book.id}
          </p>
        </div>
        <div className="report-actions">
          <button type="button" className="primary" onClick={() => window.print()}>
            Print / PDF
          </button>
          <button type="button" ref={close} onClick={onClose}>
            Close
          </button>
        </div>
      </header>

      {envelope.book.sheets.map((sheet) => {
        const form = catalogue.find((entry) => entry.id === sheet.calculation);
        return (
          <section className="calc-report-sheet" key={sheet.id}>
            <h2>{sheet.name}</h2>
            <p className="report-summary">
              {sheet.calculation}
              {form === undefined ? "" : ` · ${form.name}`}
            </p>
            {form === undefined ? null : <p>{form.summary}</p>}

            {form === undefined ? (
              <p className="note">this calculation is not in this build</p>
            ) : (
              <>
                <h3>Inputs</h3>
                <table className="grid">
                  <thead>
                    <tr>
                      <th scope="col">input</th>
                      <th scope="col">meaning</th>
                      <th scope="col">value</th>
                      <th scope="col">unit</th>
                    </tr>
                  </thead>
                  <tbody>
                    {form.inputs.map((input) => (
                      <tr key={input.name}>
                        <th scope="row">{input.name}</th>
                        <td>{input.description}</td>
                        <td data-num>{inputText(sheet.inputs[input.name], input, units)}</td>
                        <td className="unit">{unitText(input.unit, units)}</td>
                      </tr>
                    ))}
                  </tbody>
                </table>

                <h3>Outputs</h3>
                {sheet.result === null ? (
                  <p className="note">
                    {sheet.error === null ? "not calculated" : sheet.error.message}
                  </p>
                ) : (
                  <table className="grid">
                    <thead>
                      <tr>
                        <th scope="col">output</th>
                        <th scope="col">meaning</th>
                        <th scope="col">value</th>
                        <th scope="col">unit</th>
                      </tr>
                    </thead>
                    <tbody>
                      {form.outputs.map((output) => (
                        <tr key={output.name}>
                          <th scope="row">{output.name}</th>
                          <td>{output.description}</td>
                          <td data-num>
                            {outputText(sheet.result?.outputs[output.name], output, units)}
                          </td>
                          <td className="unit">
                            {outputUnit(sheet.result?.outputs[output.name], output, units)}
                          </td>
                        </tr>
                      ))}
                    </tbody>
                  </table>
                )}

                {form.assumptions.length === 0 ? null : (
                  <>
                    <h3>When this applies</h3>
                    <p>
                      These conditions are not checked at runtime; the result means what it says
                      only when they hold.
                    </p>
                    <ul>
                      {form.assumptions.map((assumption) => (
                        <li key={assumption}>{assumption}</li>
                      ))}
                    </ul>
                  </>
                )}

                {form.references.length === 0 ? null : (
                  <>
                    <h3>Source and references</h3>
                    <ul>
                      {form.references.map((reference) => (
                        <li key={reference}>{reference}</li>
                      ))}
                    </ul>
                  </>
                )}

                {sheet.result !== null && sheet.result.warnings.length > 0 ? (
                  <>
                    <h3>Warnings</h3>
                    {sheet.result.warnings.map((warning, index) => (
                      <p key={`${warning.code}:${index}`}>
                        <strong>{warning.code}</strong> {warning.message}
                      </p>
                    ))}
                  </>
                ) : null}

                <h3>Model check</h3>
                <table className="grid">
                  <tbody>
                    <tr>
                      <th scope="row">model</th>
                      <td>{form.id}</td>
                    </tr>
                    <tr>
                      <th scope="row">verification</th>
                      <td>{form.provenance.verification}</td>
                    </tr>
                    <tr>
                      <th scope="row">spec hash</th>
                      <td>{form.provenance.spec_sha256 ?? "—"}</td>
                    </tr>
                    <tr>
                      <th scope="row">rust hash</th>
                      <td>{form.provenance.rust_sha256 ?? "—"}</td>
                    </tr>
                    <tr>
                      <th scope="row">tests</th>
                      <td>
                        {form.provenance.tests_active} active ·{" "}
                        {form.provenance.tests_skipped} skipped
                      </td>
                    </tr>
                  </tbody>
                </table>
              </>
            )}
          </section>
        );
      })}
    </div>
  );
}

function unitText(unit: string | null, units: Units): string {
  if (unit === null || unit === "dimensionless") {
    return "—";
  }
  return displayUnitOf(units, unit) ?? unit;
}

function inputText(
  value: number | boolean | string | undefined,
  input: CalcInput,
  units: Units,
): string {
  if (value === undefined) {
    return "—";
  }
  if (typeof value === "number" && input.unit !== null) {
    return input.unit === "dimensionless"
      ? String(value)
      : quantityParts(value, input.unit, 6, units).value;
  }
  return String(value);
}

function outputText(value: SheetValue | undefined, output: CalcOutput, units: Units): string {
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

function outputUnit(
  value: SheetValue | undefined,
  output: CalcOutput,
  units: Units,
): string {
  if (typeof value === "object" && value !== null && !Array.isArray(value)) {
    return displayUnitOf(units, value.unit) ?? value.unit;
  }
  return unitText(output.unit, units);
}
