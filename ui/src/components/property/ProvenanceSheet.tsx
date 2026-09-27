import { provenanceOf, shortHash } from "../../state/provenance";
import type { Envelope } from "../../wire/types";

/**
 * What this unit operation rests on: where it came from, and how far it is checked.
 *
 * **Every value here is the library's own word.** `source` is the provenance line the palette
 * carries — the NeqSim class a unit operation was ported from — `verification` is the model's
 * declared status, and the two hashes are of the spec and the Rust kernel the answer came from.
 * Nothing on this sheet is a second label for something the library already named, which is the
 * rule the result sheets follow and the reason this sheet exists at all: a reader asking "how
 * far is this checked" should not have to open a TOML file to find out.
 *
 * **A document with no entry for the object says so.** The record is written per *instantiated*
 * unit operation, so a feed, a product and a connection have none, and an instance the
 * projection drew but the palette cannot resolve has none either — which is a fact about the
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
            <th scope="row">verification</th>
            {/* The library's own status word, and not a sentence about it. */}
            <td>{entry.verification ?? "—"}</td>
          </tr>
          <tr>
            <th scope="row">runs</th>
            <td>{entry.runnable ? "yes" : (entry.refusal ?? "no")}</td>
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
