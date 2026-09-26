import type { EditorCommand } from "../../wire/commands";
import type { Envelope, GraphEdge, GraphNode } from "../../wire/types";

/**
 * What is wired to this object, one edge per row.
 *
 * **The rows select, which is what the old list could not do.** It printed the paths; an edge is
 * also an object with its own window - a tear has seven settings - so a row that selects the edge
 * is the way from a unit op to the loop it sits in.
 */
export function ConnectionsSheet({
  envelope,
  node,
  edge,
  onCommand,
  onSelect,
}: {
  envelope: Envelope;
  /** The object whose wiring this is, or `null` when the selection is an edge. */
  node: GraphNode | null;
  /** The selected edge, or `null` when the selection is a node. */
  edge: GraphEdge | null;
  onCommand: (command: EditorCommand) => void;
  onSelect: (id: string | null) => void;
}) {
  const edges =
    node === null
      ? []
      : envelope.flowsheet.graph.edges.filter(
          (candidate) => candidate.source === node.id || candidate.target === node.id,
        );

  return (
    <div className="sheet" id="sheet-connections" role="tabpanel" aria-labelledby="tab-connections">
      {edge === null ? null : (
        <div className="field">
          <div className="head">
            <span className="name">
              {edge.data.from} → {edge.data.to}
            </span>
            <span className="unit">{edge.data.kind}</span>
          </div>
        </div>
      )}

      {node === null ? null : edges.length === 0 ? (
        <p className="note">nothing is wired to it yet</p>
      ) : (
        <table className="grid">
          <thead>
            <tr>
              <th scope="col">stream</th>
              <th scope="col">kind</th>
              <th scope="col">from</th>
              <th scope="col">to</th>
            </tr>
          </thead>
          <tbody>
            {edges.map((candidate) => {
              const leaves = candidate.source === node.id;
              const other = leaves ? candidate.target : candidate.source;
              return (
                <tr key={candidate.id} data-id={candidate.id}>
                  <th scope="row">
                    <button type="button" className="link" onClick={() => onSelect(candidate.id)}>
                      {candidate.data.path}
                    </button>
                  </th>
                  <td>{candidate.data.kind}</td>
                  <td>{leaves ? "this" : other.replace(/^instance:/, "")}</td>
                  <td>{leaves ? other.replace(/^instance:/, "") : "this"}</td>
                </tr>
              );
            })}
          </tbody>
        </table>
      )}

      {edge === null ? null : (
        <div className="field">
          {edge.data.kind === "recycle" ? (
            <button
              type="button"
              className="entry"
              onClick={() => onCommand({ command: "remove_recycle", stream: edge.data.path })}
            >
              Remove the tear
            </button>
          ) : (
            <>
              <button
                type="button"
                className="entry"
                onClick={() =>
                  onCommand({
                    command: "disconnect",
                    from: edge.data.from,
                    to: edge.data.to,
                  })
                }
              >
                Remove the connection
              </button>
              <button
                type="button"
                className="entry"
                onClick={() =>
                  onCommand({
                    command: "add_recycle",
                    stream: nextTear(envelope),
                    from: edge.data.from,
                    to: edge.data.to,
                  })
                }
              >
                Make it a tear
              </button>
            </>
          )}
        </div>
      )}
    </div>
  );
}

/** A fresh tear name, over the ones the document already holds. */
export function nextTear(envelope: Envelope): string {
  const taken = new Set(envelope.flowsheet.graph.edges.map((candidate) => candidate.data.path));
  for (let index = 1; index < 100; index += 1) {
    const candidate = `recycle_${index}`;
    if (!taken.has(candidate)) {
      return candidate;
    }
  }
  return "recycle_x";
}
