/**
 * What the library's verification status means, in one sentence.
 *
 * **A record field, not a verdict on the answer.** The status says how far the *model* is
 * checked against evidence outside its own examples - an independent validation case, a source
 * that resolved, a test that ran. It does not say whether this run calculated, which the run
 * status and the presence of a result already say. It is shown once, on the Provenance tab and
 * in the report's machine appendix, and never beside a value.
 *
 * **An unknown spelling stays visible.** The library's status is derived rather than declared,
 * so a future vocabulary is possible; mapping an unrecognised word to `unverified` would hide a
 * new fact behind an old one. The raw string is printed and the sentence says it is unrecognised.
 */

export interface VerificationInfo {
  /** The library's own spelling, or `null` where the entry has no model. */
  raw: string | null;
  /** The one sentence a reader needs, in the same words on every surface. */
  definition: string;
  /** Whether this editor knows the spelling. Untrue is a finding, not a fallback. */
  known: boolean;
}

const KNOWN: ReadonlyMap<string, string> = new Map([
  [
    "verified",
    "An independent validation case is recorded and agrees with the model.",
  ],
  [
    "partially_verified",
    "The kernels reproduce the examples pinned in the model's own specification; no independent validation case is recorded. This is the normal state for the port.",
  ],
  [
    "source_needed",
    "A cited source was sought and could not be resolved; the model's answer is not attributed to it.",
  ],
  ["unverified", "No test exercises this model, so its accuracy cannot be claimed."],
]);

/** The raw status and its definition, or an explicit statement that the spelling is unknown. */
export function definitionOf(raw: string | null): VerificationInfo {
  if (raw === null) {
    return {
      raw: null,
      definition: "This entry has no registered model.",
      known: true,
    };
  }
  const definition = KNOWN.get(raw);
  if (definition === undefined) {
    return {
      raw,
      definition:
        "The library published a verification status this editor does not recognise; do not read it as verified.",
      known: false,
    };
  }
  return { raw, definition, known: true };
}
