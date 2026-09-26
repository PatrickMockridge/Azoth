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
 */

import { displayUnitOf, inDisplayUnit, type Units } from "../state/units";
import {
  boundsOf,
  commandValueIn,
  controlFor,
  fieldTextIn,
  violates,
} from "../wire/field";
import type { FormParameter } from "../wire/types";

export interface ParameterFieldProps {
  parameter: FormParameter;
  value: unknown;
  /** The unit a reader wants values in, or `null` where there is no catalogue yet. */
  units: Units | null;
  onChange: (value: unknown) => void;
}

export function ParameterField({ parameter, value, units, onChange }: ParameterFieldProps) {
  const control = controlFor(parameter.kind);
  const bounds = boundsOf(parameter);
  const text = fieldTextIn(value, parameter.unit, units);
  // **The violation is checked in the document's unit**, which is the value this holds: the
  // model's range is stated against the physics, so a temperature typed in `°F` has already
  // become a kelvin by the time it is compared with the bound.
  const outside = violates(bounds, value);
  const empty = value === undefined || value === null;
  // The unit the field shows, and the bounds in it - so an input's own `min`/`max` are the
  // numbers a reader sees rather than the ones the document holds.
  const shown = displayUnitOf(units, parameter.unit);
  const at = (bound: number | null) =>
    bound === null || parameter.unit === null || units === null
      ? bound
      : inDisplayUnit(units, bound, parameter.unit);
  const write = (raw: string) =>
    onChange(commandValueIn(parameter.kind, raw, parameter.unit, units));

  return (
    <div className="field">
      <div className="head">
        <span className="name">{parameter.name}</span>
        {parameter.required && empty ? (
          <span className="required" title="a kernel cannot run without it">
            *
          </span>
        ) : null}
        {shown === null ? null : <span className="unit">{shown}</span>}
      </div>

      {control === "select" ? (
        <select
          value={text}
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
          value={text}
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
          className={outside ? "violating" : ""}
          value={text}
          inputMode={control === "number" ? "decimal" : "text"}
          placeholder={control === "list" ? "0.5, 0.3, 0.2" : ""}
          onChange={(event) => write(event.target.value)}
          {...(at(bounds.min) === null ? {} : { min: at(bounds.min) as number })}
          {...(at(bounds.max) === null ? {} : { max: at(bounds.max) as number })}
        />
      )}

      {outside ? (
        <div className={`hint ${bounds.severity === "error" ? "bad" : "warn"}`}>
          {bounds.rationale ?? "outside the range this model is validated in"}
        </div>
      ) : null}
      {parameter.kind === "unknown" ? null : (
        <div className="hint">{parameter.description}</div>
      )}
    </div>
  );
}
