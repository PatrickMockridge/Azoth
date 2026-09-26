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
 * recycle wireable to a mixer. **They are spread along the side rather than stacked**, so a mixer
 * with three feeds shows three places to connect rather than one that three edges leave from.
 *
 * **Which side is `state/glyphs.ts`'s answer and not the direction's.** An inlet is on the left and
 * an outlet on the right for most machinery, but a column's distillate leaves the top and its
 * bottoms the bottom — and a drawing that put those on the sides would be drawing a different
 * machine from the one the declaration describes.
 */

import { Handle, Position, type NodeProps } from "@xyflow/react";

import { glyphFor, isVertical, senseFor, sideOf, type Side } from "../state/glyphs";
import { Glyph } from "./Glyph";
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

/** A side, as xyflow names it. */
const POSITION: Record<Side, Position> = {
  top: Position.Top,
  bottom: Position.Bottom,
  left: Position.Left,
  right: Position.Right,
};

/** One handle: which port it belongs to, which side it sits on, and how far along that side. */
type Slot = {
  handle: string;
  port: string;
  direction: "in" | "out";
  side: Side;
  at: number;
};

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

/** Every handle this node draws, with the place its port's side puts it. */
function slotsOf(form: Form | undefined, graph: GraphNode): Slot[] {
  const unit = graph.data.unit ?? "";
  return (form?.ports ?? []).flatMap((port) =>
    handlesOf(port, graph).map((handle, index, all) => ({
      handle,
      port: port.name,
      direction: port.direction,
      side: sideOf(unit, port.name, port.direction),
      // Evenly along the side, never at an end: a handle at 0% sits on the corner an edge would
      // leave from anyway, and one at 100% reads as belonging to the next side.
      at: (index + 1) / (all.length + 1),
    })),
  );
}

/** The offset along a side a handle sits at, which is the axis the side does not fix. */
function offset(slot: Slot): { left?: string; top?: string } {
  const percent = `${String(slot.at * 100)}%`;
  return isVertical(slot.side) ? { left: percent } : { top: percent };
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
          position={leaving ? Position.Right : Position.Left}
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
      {slotsOf(form, graph).map((slot) => (
        <Handle
          key={slot.handle}
          type={slot.direction === "in" ? "target" : "source"}
          position={POSITION[slot.side]}
          id={slot.handle}
          title={slot.port}
          style={offset(slot)}
        />
      ))}
      <Glyph
        kind={glyphFor(graph.data.unit)}
        sense={senseFor(graph.data.unit)}
        selected={selected === true}
        flagged={bad}
      />
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
