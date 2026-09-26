/**
 * What a table's columns are, and where the substance axis comes from.
 *
 * **A solved stream's `z` is a bare array**, positional and unnamed: the library's record is
 * `n, z, P, T, h` and the names live on the *feed* that declared the fluid. So the axis has to be
 * recovered from the document, and the rule for that is this module's first content - with the one
 * case it must refuse to guess about.
 */

import type { Envelope } from "../wire/types";

/** The substance axis a composition table is laid out over. */
export interface Axis {
  /** One entry per column. */
  names: readonly string[];
  /**
   * Whether the names are **positions rather than names**.
   *
   * True when the feeds disagree about the fluid - two lists of the same length with different
   * names, or lists of different lengths - because then no feed's own list describes every stream
   * and a header drawn from one of them would be a claim the document does not make.
   */
  positional: boolean;
  /** Why, when it is positional; `null` otherwise. */
  note: string | null;
}

/**
 * The substances every stream of this document is made of, read from the feeds that declare them.
 *
 * **A feed is the only place the names are stated.** A product has no list of its own and a solved
 * stream's `z` is a bare array, so there is nothing else to read - and where the feeds disagree,
 * `positional` says so rather than silently preferring the first one.
 */
export function componentAxis(envelope: Envelope): Axis {
  const declared = envelope.flowsheet.graph.nodes
    .filter((node) => node.role === "input")
    .map((node) => node.data.input?.components ?? [])
    .filter((components) => components.length > 0);

  const first = declared[0];
  if (first === undefined) {
    return {
      names: [],
      positional: true,
      note: "no feed this document declares names its substances, so a composition has no names to show",
    };
  }
  if (declared.length === 1) {
    return { names: first, positional: false, note: null };
  }

  const agrees = declared.every(
    (list) => list.length === first.length && list.every((name, i) => name === first[i]),
  );
  if (agrees) {
    return { names: first, positional: false, note: null };
  }

  const longest = declared.reduce((a, b) => (b.length > a.length ? b : a));
  const lists = declared.map((list) => list.join(" / ")).join("; and ");
  return {
    names: longest,
    positional: true,
    note: `the feeds disagree about the fluid (${lists}), so these columns are positions and not names`,
  };
}

/** How a fractional column's cell reads, to four decimals - enough for a mole fraction. */
export function fraction(value: number): string {
  return value.toFixed(4);
}
