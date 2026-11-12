/**
 * The machine-readable run record: what the report's appendix keeps in a PDF, as JSON.
 *
 * **The envelope, not a summary of it.** The document, the session's streams and results exactly
 * as the library serialised them, the provenance and the diagnostics are copied rather than
 * reinterpreted, so an agent holding this file can reproduce the run and compare its own answer
 * field by field. The one derived value is the status word, and the flags it comes from
 * (`ok`, `dirty`, `run_error`) travel beside it rather than being replaced by it.
 *
 * **The front end adds no number.** There is no converted magnitude, no rounded display value
 * and no label here; the units are the library's own and the field names are the record's.
 */

import { statusWord } from "../state/status";
import type { Units } from "../state/units";
import type { Envelope } from "../wire/types";

export interface RunRecord {
  /** The record's own name and version, so a reader can refuse one it does not know. */
  format: "azoth.run-record/1";
  case: {
    id: string;
    name: string;
    /** The case document, so the run can be repeated rather than only read. */
    document: string;
  };
  generated: string;
  /** The reading of the flags below, not a replacement for them. */
  status: string;
  ok: boolean;
  dirty: boolean;
  execution_order: string;
  units: {
    active: string;
    name: string;
  };
  library: {
    name: string;
    version: string;
  };
  keycard: string;
  session: Envelope["session"];
  provenance: Envelope["provenance"];
  diagnostics: Envelope["diagnostics"];
  run_error: string | null;
}

/** Assemble the record. `generated` is passed in so a test can hold it to one document. */
export function buildRunRecord(
  envelope: Envelope,
  units: Units,
  generated: string,
): RunRecord {
  const set = units.sets.find((candidate) => candidate.id === units.active);
  return {
    format: "azoth.run-record/1",
    case: {
      id: envelope.flowsheet.id,
      name: envelope.flowsheet.name,
      document: envelope.flowsheet.document,
    },
    generated,
    status: statusWord(envelope),
    ok: envelope.ok,
    dirty: envelope.dirty,
    execution_order: envelope.execution_order,
    units: {
      active: units.active,
      name: set?.name ?? units.active,
    },
    library: {
      name: envelope.provenance.library,
      version: envelope.provenance.version,
    },
    keycard: envelope.provenance.keycard,
    session: envelope.session,
    provenance: envelope.provenance,
    diagnostics: envelope.diagnostics,
    run_error: envelope.run_error,
  };
}
