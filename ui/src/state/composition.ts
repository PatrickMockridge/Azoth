/**
 * A feed being edited, before it is a document record.
 *
 * **The draft is text and the record is numbers.** A field a person is halfway through typing is
 * not a composition, and writing each keystroke into the document would make the model check read
 * a state nobody meant. So the editor holds strings, validates them locally, and sends exactly
 * one `add_input` when the reader applies the feed.
 *
 * **Mass fractions are converted here, with the catalogue's molar masses.** There is no
 * conversion table in the front end, and no guessed molar mass: a component the databank cannot
 * weigh disables the mass mode rather than inventing a number.
 */

import { inDisplayUnit, toDisplayUnit, type Units } from "./units";
import type { Catalogue, InputRecord } from "../wire/types";

/** One row: a substance and the fraction the reader typed. */
export interface FeedComponentDraft {
  name: string;
  fraction: string;
}

/** Which measure the typed fractions are. Mole is the first-class one. */
export type FractionMode = "mole" | "mass";

/** A feed in the form, before it is sent. */
export interface FeedDraft {
  name: string;
  components: FeedComponentDraft[];
  mode: FractionMode;
  /** Molar flow in the unit set's unit. */
  n: string;
  /** Pressure in the unit set's unit. */
  P: string;
  /** Temperature in the unit set's unit. */
  T: string;
}

/** A draft from a document record, shown in the unit set in force. */
export function draftFromInput(name: string, record: InputRecord, units: Units): FeedDraft {
  return {
    name,
    components: record.components.map((component, index) => ({
      name: component,
      fraction: numberText(record.z[index] ?? 0),
    })),
    mode: "mole",
    n: numberText(inDisplayUnit(units, record.n, "mol/s")),
    P: numberText(inDisplayUnit(units, record.P, "Pa")),
    T: numberText(inDisplayUnit(units, record.T, "K")),
  };
}

/** A draft with one row, which is the smallest feed the checker accepts. */
export function emptyDraft(name: string, units: Units): FeedDraft {
  return {
    name,
    components: [],
    mode: "mole",
    // **The defaults are physical, and the boxes are in the set in force.** These used to be SI
    // literals labelled as the display unit, so a Field-set feed silently opened at 500000 psi.
    n: numberText(inDisplayUnit(units, 1, "mol/s")),
    P: numberText(inDisplayUnit(units, 5e5, "Pa")),
    T: numberText(inDisplayUnit(units, 300, "K")),
  };
}

/** The fractions as numbers, with a non-numeric entry as zero so the sum still means something. */
export function typedFractions(draft: FeedDraft): number[] {
  return draft.components.map((component) => {
    const value = Number(component.fraction);
    return Number.isFinite(value) ? value : 0;
  });
}

/** The sum the form shows; a feed that does not sum to one is not normalised behind the reader. */
export function sumOf(draft: FeedDraft): number {
  return typedFractions(draft).reduce((sum, value) => sum + value, 0);
}

/** Every row has a name, and every fraction is a finite non-negative number. */
export function readyToSend(draft: FeedDraft): boolean {
  return (
    draft.name.trim() !== "" &&
    draft.components.length > 0 &&
    draft.components.every((component) => component.name.trim() !== "") &&
    draft.components.every((component) => {
      const value = Number(component.fraction);
      return Number.isFinite(value) && value >= 0;
    }) &&
    Math.abs(sumOf(draft) - 1) < 1e-6 &&
    [draft.n, draft.P, draft.T].every((value) => Number.isFinite(Number(value)))
  );
}

/** The same draft with its fractions scaled to one, or null when the sum is not positive. */
export function normalise(draft: FeedDraft): FeedDraft | null {
  const sum = sumOf(draft);
  if (!(sum > 0)) {
    return null;
  }
  return {
    ...draft,
    components: draft.components.map((component) => ({
      ...component,
      fraction: numberText(typedFraction(component, sum)),
    })),
  };
}

/** One row's text after normalisation. */
function typedFraction(component: FeedComponentDraft, sum: number): number {
  const value = Number(component.fraction);
  return Number.isFinite(value) ? value / sum : 0;
}

/** The molar masses of a draft's rows, where the catalogue has them. */
export function molarMassesOf(draft: FeedDraft, catalogue: Catalogue | null): (number | null)[] {
  return draft.components.map(
    (component) =>
      catalogue?.components?.find((entry) => entry.name === component.name)?.molar_mass ?? null,
  );
}

/** Whether every row can be weighed, which is what mass-fraction entry needs. */
export function canUseMass(draft: FeedDraft, catalogue: Catalogue | null): boolean {
  return molarMassesOf(draft, catalogue).every((mass) => mass !== null && mass > 0);
}

/**
 * The mole fractions a draft stands for.
 *
 * A mass-fraction row `w` becomes `x = (w/M) / Σ(w/M)`; the mole mode returns the typed values.
 * Null where a mass row cannot be weighed, because the alternative is a guessed molar mass.
 */
export function moleFractions(draft: FeedDraft, catalogue: Catalogue | null): number[] | null {
  const typed = typedFractions(draft);
  if (draft.mode === "mole") {
    return typed;
  }
  const masses = molarMassesOf(draft, catalogue);
  if (!masses.every((mass) => mass !== null && mass > 0)) {
    return null;
  }
  const moles = typed.map((value, index) => value / (masses[index] as number));
  const total = moles.reduce((sum, value) => sum + value, 0);
  if (!(total > 0)) {
    return null;
  }
  return moles.map((value) => value / total);
}

/** A number as a field shows it, without trailing zeroes a reader did not type. */
export function numberText(value: number): string {
  if (!Number.isFinite(value)) {
    return "";
  }
  return String(Number(value.toPrecision(6)));
}

/** The record an `add_input` command carries. */
export interface PreparedFeed {
  name: string;
  components: string[];
  n: number;
  z: number[];
  P: number;
  T: number;
}

/** The document record a draft would send, or null where it is not ready. */
export function recordOf(
  draft: FeedDraft,
  catalogue: Catalogue | null,
  units: Units,
): PreparedFeed | null {
  if (!readyToSend(draft)) {
    return null;
  }
  const z = moleFractions(draft, catalogue);
  if (z === null) {
    return null;
  }
  return {
    name: draft.name.trim(),
    components: draft.components.map((component) => component.name.trim()),
    z,
    n: toDisplayUnit(units, Number(draft.n), "mol/s"),
    P: toDisplayUnit(units, Number(draft.P), "Pa"),
    T: toDisplayUnit(units, Number(draft.T), "K"),
  };
}
