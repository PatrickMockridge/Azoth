/**
 * What one instantiated unit operation rests on.
 *
 * **A reading of the envelope and not a second record.** The library publishes a provenance
 * entry per instance — the palette entry, the registered model, the source line the palette
 * carries, whether the executor runs it, and the spec and kernel hashes behind it — because
 * "what does this answer rest on" is a question a front end has to be able to show. This module
 * is the lookup, and it invents nothing: every value it returns is one the library wrote.
 *
 * **The lookup is by the document's own name for the instance**, which is the same string
 * `session.paths()` spells and the same one the graph's node carries, so a panel, a path and a
 * hash are three readings of one name rather than three that have to be reconciled.
 */

import type { Envelope, UnitProvenance } from "../wire/types";

/** What one instance rests on, or `undefined` where the document has no entry for it. */
export function provenanceOf(
  envelope: Envelope,
  instance: string,
): UnitProvenance | undefined {
  return envelope.provenance.units.find((entry) => entry.instance === instance);
}

/**
 * A hash as a panel shows it: the first eight characters, because the rest is for comparing.
 *
 * **`null` stays `null`**, and that is the whole reason this is a function rather than a
 * `slice` at the call site: a palette entry with no model has no hashes, and a panel that
 * showed eight characters of nothing would be showing a value where the library wrote an
 * absence.
 */
export function shortHash(hash: string | null): string | null {
  return hash === null ? null : hash.slice(0, 8);
}
