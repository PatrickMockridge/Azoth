/**
 * The one word a document's state gets, and the one class that colours it.
 *
 * **A half-built case is not a broken one.** `underfed_port` and `missing_parameter` are facts
 * about a document still being built; the checker rightly refuses to run it, but calling that
 * "refused" makes adding a palette entry look like a model error. Everything else that fails the
 * check is refused, a run that stopped is failed, and an edited solved document is stale.
 */

import type { Envelope } from "../wire/types";

/** The construction facts: a document with only these is incomplete, not refused. */
const CONSTRUCTION = new Set(["underfed_port", "missing_parameter"]);

export type StatusWord = "solved" | "stale" | "incomplete" | "refused" | "failed";

/** What the document is, from the envelope's own flags and diagnostics. */
export function statusWord(envelope: Envelope): StatusWord {
  // A failed run is cleared by an edit in the library, so a failure that is still here describes
  // the document in hand.
  if (envelope.run_error !== null) {
    return "failed";
  }
  if (!envelope.ok) {
    const errors = envelope.diagnostics.filter((diagnostic) => diagnostic.severity === "error");
    const construction =
      errors.length > 0 && errors.every((diagnostic) => CONSTRUCTION.has(diagnostic.code));
    return construction ? "incomplete" : "refused";
  }
  return envelope.dirty ? "stale" : "solved";
}

/**
 * Whether a number is calculated, in a sentence rather than a status word.
 *
 * **The question a reader asks first.** The status chip answers "what is the document", and a
 * result sheet has to answer "is this number from this run, or from before the last edit" in
 * the same breath. Kept here so the sheet, the overview and any report sentence agree.
 */
export function runNote(status: StatusWord, published: boolean): string {
  if (!published) {
    switch (status) {
      case "solved":
        return "No separate result: this operation's answer is its outlet streams.";
      case "stale":
        return "Not calculated: the document has changed since the last solve.";
      case "incomplete":
        return "Not calculated: the document is incomplete.";
      case "refused":
        return "Not calculated: the model check refused the document.";
      case "failed":
        return "The run failed; there is no result.";
    }
  }
  switch (status) {
    case "solved":
      return "Calculated, from this run.";
    case "stale":
      return "Calculated, from the last run; the document has changed since.";
    case "incomplete":
      return "Calculated before the document was left incomplete.";
    case "refused":
      return "Calculated before the document was refused.";
    case "failed":
      return "Calculated before the last run failed.";
  }
}

/** The colour class for the word: green, amber, or red. */
export function statusClass(envelope: Envelope): "ok" | "stale" | "warn" | "bad" {
  switch (statusWord(envelope)) {
    case "solved":
      return "ok";
    case "stale":
      return "stale";
    case "incomplete":
      return "warn";
    case "refused":
    case "failed":
      return "bad";
  }
}
