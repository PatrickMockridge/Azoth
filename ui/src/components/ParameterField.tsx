/**
 * One parameter, as a control.
 *
 * **The control is chosen by the kind the library declared**, and the value goes back shaped by
 * the same kind — so a field cannot write a document the checker will refuse for its type. What a
 * field *does* do is mark a value the model's own bound excludes, with the model's own sentence
 * beside it: the field reports, the run decides.
 *
 * **A unit-bearing field reads and writes in the unit set in force**, which is what makes the set
 * a set of units rather than a way of reading: the document keeps the unit its spec declares, the
 * field shows and accepts the one a reader chose, and the two conversions are `wire/field.ts`'s.
 * The bounds stay in the document's unit, because the model's range is a statement about the
 * physics and not about the display.
 *
 * **A row layout is the calculation sheet's.** The flowsheet stacks a field because its window is
 * narrow; a sheet has the width to put the explanation beside the control, where it can be read
 * before the value is typed. The control, the draft and the conversion are the same in both.
 */

import { useId, useState, type ReactNode } from "react";

import { displayUnitOf, inDisplayUnit, unitLabelOf, type Units } from "../state/units";
import {
  boundsOf,
  commandValueIn,
  controlFor,
  fieldTextIn,
  rangeTextIn,
  violationOf,
} from "../wire/field";
import type { FormParameter } from "../wire/types";

export interface ParameterFieldProps {
  parameter: FormParameter;
  value: unknown;
  /** The unit a reader wants values in, or `null` where there is no catalogue yet. */
  units: Units | null;
  onChange: (value: unknown) => void;
  /** Stacked is the flowsheet's field; row is a calculation sheet's labelled line. */
  layout?: "stacked" | "row";
  /** A control-adjacent action, such as a calculation sheet's clear button. */
  action?: ReactNode;
}

export function ParameterField({
  parameter,
  value,
  units,
  onChange,
  layout = "stacked",
  action,
}: ParameterFieldProps) {
  const control = controlFor(parameter.kind);
  const bounds = boundsOf(parameter);
  const text = fieldTextIn(value, parameter.unit, units);
  const inputId = useId();

  /**
   * **The field is a draft until the typing stops.**
   *
   * An edit used to be one command per keystroke: `340` was three commands, two of them states
   * the checker could see, and a served door made three round trips. The draft is form state and
   * not a document copy; Apply is blur or Enter, Escape throws it away, and a select still
   * commits at once because choosing an option is a finished gesture.
   */
  const [draft, setDraft] = useState<string | null>(null);
  const shownText = draft ?? text;
  const pending =
    draft === null ? value : commandValueIn(parameter.kind, draft, parameter.unit, units);
  // **The violation is checked in the document's unit**, which is the value this holds: the
  // model's range is stated against the physics, so a temperature typed in `°F` has already
  // become a kelvin by the time it is compared with the bound.
  const violation = violationOf(parameter.range, pending);
  const outside = violation !== null;
  const empty = (draft ?? text) === "";
  // The unit the field shows, and the bounds in it - so an input's own `min`/`max` are the
  // numbers a reader sees rather than the ones the document holds. The raw spelling stays on the
  // title, because that is the string the catalogue and a CSV carry.
  const rawUnit = displayUnitOf(units, parameter.unit);
  const shown = unitLabelOf(units, parameter.unit);
  const at = (bound: number | null) =>
    bound === null || parameter.unit === null || units === null
      ? bound
      : inDisplayUnit(units, bound, parameter.unit);
  const limits = rangeTextIn(parameter.range, parameter.unit, units);
  const write = (raw: string) =>
    onChange(commandValueIn(parameter.kind, raw, parameter.unit, units));
  const commit = () => {
    if (draft !== null) {
      write(draft);
      setDraft(null);
    }
  };

  const descriptionId = `${inputId}-description`;
  const limitsId = `${inputId}-limits`;
  const violationId = `${inputId}-violation`;
  const describedBy = [
    parameter.kind === "unknown" ? null : descriptionId,
    limits === null ? null : limitsId,
    outside ? violationId : null,
  ]
    .filter((entry): entry is string => entry !== null)
    .join(" ");

  const controlElement =
    control === "select" ? (
      <select
        id={inputId}
        value={text}
        aria-label={layout === "row" ? undefined : parameter.name}
        aria-required={parameter.required}
        aria-invalid={outside ? true : undefined}
        aria-describedby={describedBy === "" ? undefined : describedBy}
        onChange={(event) => write(event.target.value)}
      >
        <option value="">—</option>
        {parameter.values.map((option) => (
          <option key={option} value={option}>
            {option}
          </option>
        ))}
      </select>
    ) : control === "switch" ? (
      <select
        id={inputId}
        value={text}
        aria-label={layout === "row" ? undefined : parameter.name}
        aria-required={parameter.required}
        aria-invalid={outside ? true : undefined}
        aria-describedby={describedBy === "" ? undefined : describedBy}
        onChange={(event) => write(event.target.value)}
      >
        <option value="">—</option>
        <option value="true">true</option>
        <option value="false">false</option>
      </select>
    ) : control === "unknown" ? (
      // A kind the library cannot name. The field says so rather than guessing: this is the
      // entry `UNRUNNABLE` already names, and the palette's own refusal explains it.
      <div className="hint bad">no declaration says what kind of value this is</div>
    ) : (
      <input
        id={inputId}
        className={outside ? "violating" : ""}
        aria-label={layout === "row" ? undefined : parameter.name}
        aria-required={parameter.required}
        aria-invalid={outside ? true : undefined}
        aria-describedby={describedBy === "" ? undefined : describedBy}
        value={shownText}
        inputMode={control === "number" ? "decimal" : "text"}
        placeholder={control === "list" ? "0.5, 0.3, 0.2" : ""}
        onChange={(event) => setDraft(event.target.value)}
        onBlur={commit}
        onKeyDown={(event) => {
          if (event.key === "Enter") {
            commit();
          } else if (event.key === "Escape") {
            setDraft(null);
          }
        }}
        {...(at(bounds.min) === null ? {} : { min: at(bounds.min) as number })}
        {...(at(bounds.max) === null ? {} : { max: at(bounds.max) as number })}
      />
    );

  const status = (
    <span
      className={`status ${parameter.required ? "required" : "optional"}`}
      title={
        parameter.required
          ? "a kernel cannot run without it"
          : "left unset, the calculation uses its declared default or omits what needs it"
      }
    >
      {parameter.required ? "required" : "optional"}
    </span>
  );

  const head = (withStatus: boolean, withUnit: boolean) => (
    <div className="head">
      <span className="name">{parameter.name}</span>
      {withStatus ? (
        status
      ) : parameter.required && empty ? (
        <span className="required" title="a kernel cannot run without it">
          *
        </span>
      ) : null}
      {withUnit && shown !== null ? (
        <span className="unit" title={rawUnit ?? undefined}>
          {shown}
        </span>
      ) : null}
    </div>
  );

  const explanation =
    parameter.kind === "unknown" ? null : (
      <label className="description" htmlFor={inputId} id={descriptionId}>
        {parameter.description}
      </label>
    );

  const violationHint =
    violation === null ? null : (
      <div
        className={`hint ${violation.severity === "error" ? "bad" : "warn"}`}
        id={violationId}
      >
        {violation.rationale || "outside the range this model is validated in"}
      </div>
    );

  if (layout === "row") {
    return (
      <div className="field field-row">
        <div className="field-text">
          {head(true, false)}
          {explanation}
          {limits === null ? null : (
            <div className="limits" id={limitsId}>
              Checked limits: {limits}
            </div>
          )}
          {violationHint}
        </div>
        <div className="field-control">
          <div className="field-value">
            {controlElement}
            {shown === null ? null : (
              <span className="unit" title={rawUnit ?? undefined}>
                {shown}
              </span>
            )}
          </div>
          {action}
        </div>
      </div>
    );
  }

  return (
    <div className="field">
      {head(false, true)}
      {controlElement}
      {violationHint}
      {explanation}
    </div>
  );
}
