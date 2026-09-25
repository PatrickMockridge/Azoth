/**
 * One parameter, as a control.
 *
 * **The control is chosen by the kind the library declared**, and the value goes back shaped by
 * the same kind — so a field cannot write a document the checker will refuse for its type. What a
 * field *does* do is mark a value the model's own bound excludes, with the model's own sentence
 * beside it: the field reports, the run decides.
 */

import { boundsOf, commandValue, controlFor, fieldText, violates } from "../wire/field";
import type { FormParameter } from "../wire/types";

export interface ParameterFieldProps {
  parameter: FormParameter;
  value: unknown;
  onChange: (value: unknown) => void;
}

export function ParameterField({ parameter, value, onChange }: ParameterFieldProps) {
  const control = controlFor(parameter.kind);
  const bounds = boundsOf(parameter);
  const text = fieldText(value);
  const outside = violates(bounds, value);
  const empty = value === undefined || value === null;

  return (
    <div className="field">
      <div className="head">
        <span className="name">{parameter.name}</span>
        {parameter.required && empty ? (
          <span className="required" title="a kernel cannot run without it">
            *
          </span>
        ) : null}
        {parameter.unit === null ? null : <span className="unit">{parameter.unit}</span>}
      </div>

      {control === "select" ? (
        <select
          value={text}
          onChange={(event) => onChange(commandValue(parameter.kind, event.target.value))}
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
          onChange={(event) => onChange(commandValue(parameter.kind, event.target.value))}
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
          onChange={(event) => onChange(commandValue(parameter.kind, event.target.value))}
          {...(bounds.min === null ? {} : { min: bounds.min })}
          {...(bounds.max === null ? {} : { max: bounds.max })}
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
