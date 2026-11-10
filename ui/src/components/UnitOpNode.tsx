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

import { Handle, Position, type Node, type NodeProps } from "@xyflow/react";

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
  /** The worst diagnostic about this object, or none. */
  severity: "error" | "warning" | null;
  /** The handles an edge already uses, so a port can show whether it is wired. */
  connected: ReadonlySet<string>;
};

/**
 * The node xyflow holds, typed by its payload.
 *
 * **`NodeProps` is generic over the node**, so naming the payload here is what makes `data` the
 * projection's shape in the component below instead of `Record<string, unknown>`. The alternative
 * was a cast or a re-read per render, and neither is a statement about the type — this is.
 */
export type AppNode = Node<NodePayload, typeof UNIT_NODE | typeof STREAM_NODE>;

export const UNIT_NODE = "unit_op";
export const STREAM_NODE = "stream";

/** A side, as xyflow names it. */
const POSITION: Record<Side, Position> = {
  top: Position.Top,
  bottom: Position.Bottom,
  left: Position.Left,
  right: Position.Right,
};

/** Port names that are only a direction, so a label does not read `IN inlet`. */
const GENERIC = new Set(["in", "inlet", "out", "outlet"]);

/**
 * The visible text beside one handle.
 *
 * **The direction is a word and the port is its name**, so a reader does not have to know that a
 * left-hand dot is an inlet. A port called `hot_in` says `IN hot`; a bare `inlet` says `IN`; a
 * three-way splitter's `products` become `OUT 1`, `OUT 2`, `OUT 3`. The number is only there when
 * a port draws more than one handle, because a number on a single connection is noise.
 */
export function portLabel(slot: Slot): string {
  const prefix = slot.direction === "in" ? "IN" : "OUT";
  let name = slot.port
    .replace(/_(in|out|inlet|outlet)$/i, "")
    .replace(/_/g, " ")
    .trim();
  if (GENERIC.has(name.toLowerCase())) {
    name = "";
  }
  if (name.toLowerCase() === "products") {
    name = slot.count > 1 ? "" : "product";
  }
  const number = slot.count > 1 ? ` ${slot.index + 1}` : "";
  return name === "" ? `${prefix}${number}` : `${prefix} ${name}${number}`;
}

/** One handle: which port it belongs to, which side it sits on, and how far along that side. */
type Slot = {
  handle: string;
  port: string;
  direction: "in" | "out";
  side: Side;
  at: number;
  /** Its place among the handles of that port, and how many there are. */
  index: number;
  count: number;
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

/**
 * Every handle this node draws, with the place its port's side puts it.
 *
 * **The distribution is per side, not per port.** A port's handles used to be spread within that
 * port alone, so two single-handle ports on the same side both landed at 50% and one was drawn on
 * top of the other - the heat exchanger's two inlets looked like one. A side is the space a reader
 * sees, so every handle on one is distributed across it in declaration order.
 */
export function portSlots(form: Form | undefined, graph: GraphNode): Slot[] {
  const unit = graph.data.unit ?? "";
  const slots: Slot[] = (form?.ports ?? []).flatMap((port) =>
    handlesOf(port, graph).map((handle, index, handles) => ({
      handle,
      port: port.name,
      direction: port.direction,
      side: sideOf(unit, port.name, port.direction),
      at: 0,
      index,
      count: handles.length,
    })),
  );

  const sides = new Map<Side, Slot[]>();
  for (const slot of slots) {
    const group = sides.get(slot.side) ?? [];
    group.push(slot);
    sides.set(slot.side, group);
  }
  for (const group of sides.values()) {
    // Evenly along the side, never at an end: a handle at 0% sits on the corner an edge would
    // leave from anyway, and one at 100% reads as belonging to the next side.
    group.forEach((slot, index) => {
      slot.at = (index + 1) / (group.length + 1);
    });
  }
  return slots;
}

/** The offset along a side a handle sits at, which is the axis the side does not fix. */
function offset(slot: Slot): { left?: string; top?: string } {
  const percent = `${String(slot.at * 100)}%`;
  return isVertical(slot.side) ? { left: percent } : { top: percent };
}

export function UnitOpNode({ data, selected }: NodeProps<AppNode>) {
  // **No cast.** The canvas is handed this payload by `Flowsheet`, which builds it from the graph
  // document, and `NodeProps<AppNode>` is what makes that a statement about the type rather than
  // an assertion in the component that reads it.
  const { graph, form, readout, severity, connected } = data;

  if (graph.role !== "instance") {
    // A boundary stream: one handle on the side the stream leaves or arrives at, named by the
    // stream itself, so a feed's handle id is the feed's own path.
    const leaving = graph.role === "input";
    return (
      <div className={`stream${selected === true ? " selected" : ""}${severity === "error" ? " bad" : severity === "warning" ? " warn" : ""}`}>
        <Handle
          type={leaving ? "source" : "target"}
          position={leaving ? Position.Right : Position.Left}
          id={graph.data.name}
          aria-label={leaving ? "outlet feed" : "inlet product"}
          className={`handle${connected.has(graph.data.name) ? " connected" : " unconnected"}`}
        />
        <span
          className={`port-label ${leaving ? "right out" : "left in"}`}
          title={leaving ? "outlet feed" : "inlet product"}
        >
          {leaving ? "OUT feed" : "IN product"}
        </span>
        {graph.data.name}
        {readout.map((line) => (
          <div key={line.port} className="readout">
            {line.text}
          </div>
        ))}
      </div>
    );
  }

  const slots = portSlots(form, graph);
  return (
    <div className={`unit${selected === true ? " selected" : ""}${severity === "error" ? " bad" : severity === "warning" ? " warn" : ""}`}>
      {slots.map((slot) => (
        <Handle
          key={slot.handle}
          type={slot.direction === "in" ? "target" : "source"}
          position={POSITION[slot.side]}
          id={slot.handle}
          title={slot.port}
          aria-label={`${slot.direction === "in" ? "inlet" : "outlet"} ${slot.port}`}
          className={`handle${connected.has(slot.handle) ? " connected" : " unconnected"}`}
          style={offset(slot)}
        />
      ))}
      {slots.map((slot) => (
        <span
          key={`label:${slot.handle}`}
          className={`port-label ${slot.side} ${slot.direction}`}
          style={offset(slot)}
          title={`${slot.direction === "in" ? "inlet" : "outlet"} ${slot.port}`}
        >
          {portLabel(slot)}
        </span>
      ))}
      <Glyph
        kind={glyphFor(graph.data.unit)}
        sense={senseFor(graph.data.unit)}
        selected={selected === true}
        flagged={severity === "error"}
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
