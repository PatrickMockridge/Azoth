import { provenanceOf, shortHash } from "../../state/provenance";
import { resultOf } from "../../state/results";
import { runNote, statusWord } from "../../state/status";
import { definitionOf } from "../../state/verification";
import type { Envelope } from "../../wire/types";

/**
 * What this unit operation rests on: where it came from, and how far it is checked.
 *
 * **The machine record, in the one place a reader looks for it.** `source` is the provenance
 * line the palette carries - the NeqSim class a unit operation was ported from - `verification`
 * is the model's derived status, the hashes are of the spec and the Rust kernel the answer came
 * from, and `keycard` is the data authority the run used. Nothing here is a second label for
 * something the library already named; the status is printed as the library spells it, with one
 * sentence saying what that spelling means.
 *
 * **A document with no entry for the object says so.** The record is written per *instantiated*
 * unit operation, so a feed, a product and a connection have none, and an instance the
 * projection drew but the palette cannot resolve has none either - which is a fact about the
 * document rather than an error, and it is rendered as one.
 */
export function ProvenanceSheet({
  envelope,
  instance,
}: {
  envelope: Envelope;
  instance: string;
}) {
  const entry = provenanceOf(envelope, instance);
  if (entry === undefined) {
    return (
      <div className="sheet" id="sheet-provenance" role="tabpanel" aria-labelledby="tab-provenance">
        <p className="note">
          the document has no provenance for <code>{instance}</code>, which is what a unit
          operation the palette cannot resolve reads as
        </p>
      </div>
    );
  }

  const status = statusWord(envelope);
  const verification = definitionOf(entry.verification);
  const calculated = resultOf(envelope, instance) !== undefined;

  return (
    <div className="sheet" id="sheet-provenance" role="tabpanel" aria-labelledby="tab-provenance">
      <table className="grid">
        <tbody>
          <tr>
            <th scope="row">palette entry</th>
            <td>{entry.unit}</td>
          </tr>
          <tr>
            <th scope="row">model</th>
            <td>{entry.model ?? "—"}</td>
          </tr>
          <tr>
            <th scope="row">source</th>
            <td>{entry.source ?? "—"}</td>
          </tr>
          <tr>
            <th scope="row">data authority</th>
            <td>{envelope.provenance.keycard}</td>
          </tr>
          <tr>
            <th scope="row">run</th>
            <td>
              {status} · {runNote(status, calculated)}
            </td>
          </tr>
          <tr>
            <th scope="row">verification</th>
            <td>
              <code>{verification.raw ?? "—"}</code>
              <p className="result-description">{verification.definition}</p>
            </td>
          </tr>
          <tr>
            <th scope="row">spec</th>
            {/* The whole digest in the title, because the eight characters are for reading. */}
            <td title={entry.spec_sha256 ?? undefined}>{shortHash(entry.spec_sha256) ?? "—"}</td>
          </tr>
          <tr>
            <th scope="row">kernel</th>
            <td title={entry.rust_sha256 ?? undefined}>
              {shortHash(entry.rust_sha256) ?? "—"}
            </td>
          </tr>
        </tbody>
      </table>
    </div>
  );
}
