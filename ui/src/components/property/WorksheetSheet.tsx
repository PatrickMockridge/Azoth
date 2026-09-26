import { formatQuantity } from "../../wire/field";
import type { Envelope, GraphNode } from "../../wire/types";

/**
 * One unit operation's streams, as the row of a workbook its window shows.
 *
 * **The rows come from the *document* and the values from the run**, which is what makes this a
 * worksheet rather than a run log: a port path is a fact about the declaration, so the table has
 * the same shape before and after Solve and a stream the run did not reach is a row of blanks
 * rather than a missing row.
 */
export function WorksheetSheet({ envelope, node }: { envelope: Envelope; node: GraphNode }) {
  const ports = node.data.ports;
  const inlets = ports?.inlets.flatMap((port) => port.handles) ?? [];
  const outlets = ports?.outlets.flatMap((port) => port.handles) ?? [];

  const rows = [
    ...inlets.map((path) => ({ path, side: "in" as const })),
    ...outlets.map((path) => ({ path, side: "out" as const })),
  ];

  return (
    <div className="sheet" id="sheet-worksheet" role="tabpanel" aria-labelledby="tab-worksheet">
      {rows.length === 0 ? (
        <p className="note">this object has no ports, so it carries no streams</p>
      ) : (
        <table className="grid">
          <thead>
            <tr>
              <th scope="col">stream</th>
              <th scope="col">side</th>
              <th scope="col">n mol/s</th>
              <th scope="col">T K</th>
              <th scope="col">P Pa</th>
              <th scope="col">h J/mol</th>
              <th scope="col">VF</th>
            </tr>
          </thead>
          <tbody>
            {rows.map(({ path, side }) => {
              const stream = envelope.session?.streams[path];
              return (
                <tr key={path}>
                  <th scope="row">{path}</th>
                  <td>{side}</td>
                  <td data-num>
                    {stream === undefined ? "—" : formatQuantity(stream.n.magnitude_si, stream.n.unit, 3)}
                  </td>
                  <td data-num>
                    {stream === undefined ? "—" : formatQuantity(stream.T.magnitude_si, stream.T.unit, 4)}
                  </td>
                  <td data-num>
                    {stream === undefined ? "—" : formatQuantity(stream.P.magnitude_si, stream.P.unit, 3)}
                  </td>
                  <td data-num>
                    {stream === undefined ? "—" : formatQuantity(stream.h.magnitude_si, stream.h.unit, 4)}
                  </td>
                  <td data-num>
                    {stream?.vapour_fraction === null || stream?.vapour_fraction === undefined
                      ? "—"
                      : stream.vapour_fraction.toFixed(4)}
                  </td>
                </tr>
              );
            })}
          </tbody>
        </table>
      )}
    </div>
  );
}
