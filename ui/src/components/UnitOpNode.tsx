/**
 * A node, drawn from the declaration rather than from a picture of it.
 *
 * **The handles are the ports the projection emitted, and their ids are the stream paths.** That
 * is the payoff of the projection's second rule: an edge's `sourceHandle` is `sep1.liquid`, which
 * is the session's key for the same stream, so a connection needs no translation and a value
 * needs no lookup.
 *
 * A `many` outlet draws one handle per stream it returns — three for a splitter with three
 * factors — and a `many` inlet draws one however many edges reach it, which is what makes a
 * recycle wireable to a mixer.
 */

import { Handle, Position, type NodeProps } from "@xyflow/react";

import type { Form, FormPort, GraphNode } from "../wire/types";

/** One line of a node's readout: what a stream on one of its ports reached. */
export type Readout = { port: string; text: string };

/**
 * What the projection put on a node, plus the form that explains it.
 *
 * The graph document carries the ports and the values; the form carries the names, the units and
 * the kinds. A node with no form is an instance the palette does not carry, which is drawn with
 * no ports because there is no declaration to draw them from.
 */
export type NodePayload = {
  graph: GraphNode;
  form?: Form;
  readout: Readout[];
  bad: boolean;
};

export const UNIT_NODE = "unit_op";
export const STREAM_NODE = "stream";

/** The side a port's handles sit at. */
function side(direction: "in" | "out"): Position {
  return direction === "in" ? Position.Left : Position.Right;
}

/**
 * The handle ids one port draws.
 *
 * **Read from the graph and not computed here.** How many streams a `many` outlet returns is a
 * value — the splitter's own factors, or the positions the wiring names — and the projection is
 * where that was decided, so this takes its answer rather than making a second one.
 */
function handlesOf(port: FormPort, graph: GraphNode): string[] {
  const declared = graph.data.ports;
  const ports = declared?.[port.direction === "in" ? "inlets" : "outlets"];
  return ports?.find((candidate) => candidate.name === port.name)?.handles ?? [];
}

export function UnitOpNode({ data, selected }: NodeProps) {
  // The canvas is handed this payload by `Flowsheet`, which builds it from the graph document.
  const payload = data as unknown as NodePayload;
  const { graph, form, readout, bad } = payload;

  if (graph.role !== "instance") {
    // A boundary stream: one handle on the side the stream leaves or arrives at, named by the
    // stream itself, so a feed's handle id is the feed's own path.
    const leaving = graph.role === "input";
    return (
      <div className={`stream${selected === true ? " selected" : ""}${bad ? " bad" : ""}`}>
        <Handle
          type={leaving ? "source" : "target"}
          position={side(leaving ? "out" : "in")}
          id={graph.data.name}
        />
        {graph.data.name}
        <div className="hint">{graph.role === "input" ? "feed" : "product"}</div>
        {readout.map((line) => (
          <div key={line.port} className="readout">
            {line.text}
          </div>
        ))}
      </div>
    );
  }

  return (
    <div className={`unit${selected === true ? " selected" : ""}${bad ? " bad" : ""}`}>
      {(form?.ports ?? []).flatMap((port) =>
        handlesOf(port, graph).map((handle) => (
          <Handle
            key={handle}
            type={port.direction === "in" ? "target" : "source"}
            position={side(port.direction)}
            id={handle}
            title={port.name}
          />
        )),
      )}
      <header>{graph.data.unit_name ?? graph.data.unit ?? "unknown unit"}</header>
      <div className="id">
        {graph.data.name}
        {form !== undefined && !form.runnable ? " · not runnable" : ""}
      </div>
      {readout.map((line) => (
        <div key={line.port} className="readout">
          {line.text}
        </div>
      ))}
    </div>
  );
}

/** Both node kinds are one component: a boundary is the same box with no ports to draw. */
export const nodeTypes = { [UNIT_NODE]: UnitOpNode, [STREAM_NODE]: UnitOpNode };
