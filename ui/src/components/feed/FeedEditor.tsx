/**
 * A feed, edited as the record it becomes.
 *
 * **One Apply, one `add_input`.** The draft is local text; the moment it becomes the document is
 * the moment the reader presses Apply, and the command carries the whole record so the feed's
 * connections survive an edit. Normalisation is a visible action rather than a silent correction:
 * the sum is shown, and a feed that does not sum to one cannot be applied until it does.
 *
 * **Mass fractions are a view.** A reader may type them, but the model is sent mole fractions;
 * the converted values are shown before Apply, so the conversion is never a surprise.
 */

import { useState } from "react";

import { displayOf, displayUnitOf, inDisplayUnit, type Units } from "../../state/units";
import {
  canUseMass,
  emptyDraft,
  moleFractions,
  molarMassesOf,
  normalise,
  numberText,
  readyToSend,
  recordOf,
  sumOf,
  type FeedDraft,
  type PreparedFeed,
} from "../../state/composition";
import type { Catalogue } from "../../wire/types";
import { ComponentPicker } from "./ComponentPicker";

export interface FeedEditorProps {
  catalogue: Catalogue | null;
  units: Units;
  /** The feed to edit, or a new one. */
  initial: FeedDraft;
  /** Whether the feed already exists in the document. Its name is then fixed. */
  editing: boolean;
  onApply: (record: PreparedFeed) => void;
  onCancel: () => void;
}

export function FeedEditor({
  catalogue,
  units,
  initial,
  editing,
  onApply,
  onCancel,
}: FeedEditorProps) {
  const [draft, setDraft] = useState(initial);
  const sum = sumOf(draft);
  const masses = molarMassesOf(draft, catalogue);
  const massAllowed = canUseMass(draft, catalogue);
  const computed = moleFractions(draft, catalogue);
  const record = recordOf(draft, catalogue, units);
  const ready = readyToSend(draft) && record !== null;

  const change = (next: Partial<FeedDraft>) => setDraft({ ...draft, ...next });

  return (
    <div className="feed-editor">
      <div className="field">
        <div className="head">
          <span className="name">feed name</span>
        </div>
        <input
          value={draft.name}
          readOnly={editing}
          onChange={(event) => (editing ? undefined : change({ name: event.target.value }))}
        />
      </div>

      <h3>Composition</h3>
      <div className="mode" role="group" aria-label="fraction basis">
        <button
          type="button"
          className={draft.mode === "mole" ? "on" : ""}
          onClick={() => change({ mode: "mole" })}
        >
          mole fractions
        </button>
        <button
          type="button"
          className={draft.mode === "mass" ? "on" : ""}
          disabled={!massAllowed}
          title={
            massAllowed
              ? "type mass fractions; the model receives the mole fractions they stand for"
              : "every substance must have a molar mass in the databank to convert"
          }
          onClick={() => change({ mode: "mass" })}
        >
          mass fractions
        </button>
      </div>

      {draft.components.length === 0 ? (
        <p className="hint">no substances yet — search for one below</p>
      ) : (
        <div className="component-rows">
          {draft.components.map((component, index) => (
            <div className="component-row" key={`${component.name}-${String(index)}`}>
              <span className="name">{component.name}</span>
              <span className="mass">
                {masses[index] === null || masses[index] === undefined
                  ? "—"
                  : `${numberText(
                      inDisplayUnit(units, masses[index] as number, "kg/mol"),
                    )} ${displayOf(units, "kg/mol").unit}`}
              </span>
              <input
                value={component.fraction}
                inputMode="decimal"
                aria-label={`${component.name} fraction`}
                onChange={(event) =>
                  change({
                    components: draft.components.map((row, at) =>
                      at === index ? { ...row, fraction: event.target.value } : row,
                    ),
                  })
                }
              />
              <button
                type="button"
                aria-label={`remove ${component.name}`}
                onClick={() =>
                  change({ components: draft.components.filter((_, at) => at !== index) })
                }
              >
                ×
              </button>
            </div>
          ))}
        </div>
      )}

      <ComponentPicker
        catalogue={catalogue}
        chosen={draft.components.map((component) => component.name)}
        onPick={(name) =>
          change({
            components: [...draft.components, { name, fraction: draft.components.length === 0 ? "1" : "0" }],
          })
        }
      />

      <div className="sum">
        <span>
          sum <strong className={Math.abs(sum - 1) < 1e-6 ? "" : "bad"}>{numberText(sum)}</strong>
        </span>
        <button
          type="button"
          disabled={!(sum > 0) || Math.abs(sum - 1) < 1e-6}
          onClick={() => {
            const next = normalise(draft);
            if (next !== null) {
              setDraft(next);
            }
          }}
        >
          Normalise
        </button>
      </div>

      {draft.mode === "mass" && computed !== null ? (
        <p className="hint">
          the model will receive mole fractions{" "}
          {computed.map((value) => numberText(value)).join(", ")}
        </p>
      ) : null}

      <div className="conditions">
        {(
          [
            ["n", "mol/s", draft.n],
            ["P", "Pa", draft.P],
            ["T", "K", draft.T],
          ] as const
        ).map(([field, unit, value]) => (
          <label className="field" key={field}>
            <span className="head">
              <span className="name">{field}</span>
              <span className="unit">{displayUnitOf(units, unit) ?? unit}</span>
            </span>
            <input
              value={value}
              inputMode="decimal"
              onChange={(event) => {
                const typed = event.target.value;
                change(field === "n" ? { n: typed } : field === "P" ? { P: typed } : { T: typed });
              }}
            />
          </label>
        ))}
      </div>

      {ready ? null : (
        <p className="hint bad">
          {draft.components.length === 0
            ? "add at least one substance"
            : draft.mode === "mass" && computed === null
              ? "a substance here cannot be weighed; use mole fractions"
              : sum <= 0
                ? "the fractions must have a positive sum"
                : Math.abs(sum - 1) >= 1e-6
                  ? "the fractions must sum to one; press Normalise or correct them"
                  : "the flow, pressure and temperature must be numbers"}
        </p>
      )}

      <div className="feed-actions">
        <button
          type="button"
          className="primary"
          disabled={!ready}
          onClick={() => {
            if (record !== null) {
              onApply(record);
            }
          }}
        >
          {editing ? "Apply changes" : "Add the feed"}
        </button>
        <button type="button" onClick={onCancel}>
          Cancel
        </button>
      </div>
    </div>
  );
}

/** A new feed with no row and a name the document does not already hold. */
export function newFeedDraft(taken: string[]): FeedDraft {
  for (let index = 1; index < 100; index += 1) {
    const candidate = `feed_${String(index)}`;
    if (!taken.includes(candidate)) {
      return emptyDraft(candidate);
    }
  }
  return emptyDraft("feed_x");
}
