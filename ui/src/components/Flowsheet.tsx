/**
 * The canvas.
 *
 * **The editor's model is the schema's, which is what makes this file short.** The graph document
 * the library emits is already xyflow's node/edge shape — `id`, `type`, `position`, `data`,
 * `sourceHandle` — so drawing it is a pass-through, and every gesture goes back out as one command
 * that returns the whole document again. Nothing here holds a copy of the flowsheet.
 *
 * The one piece of local state is the drag: xyflow applies a position change every frame, and the
 * library is told once, on release. A `set_position` on every frame would be a command per mouse
 * move and a re-projection per frame for a figure that has not moved yet.
 */

import {
  Background,
  BackgroundVariant,
  ConnectionLineType,
  Controls,
  MarkerType,
  ReactFlow,
  useEdgesState,
  useNodesState,
  type Connection,
  type Edge,
  type Node,
} from "@xyflow/react";
import { useCallback, useEffect, useMemo } from "react";

import type { EditorCommand } from "../wire/commands";
import type { Units } from "../state/units";
import { formatQuantity } from "../wire/field";
import { removeCommandFor } from "../wire/nodes";
import { asEdgeData } from "../wire/decode";
import type { Catalogue, Envelope, GraphNode } from "../wire/types";
import { PfdEdge } from "./PfdEdge";
import { STREAM_NODE, UNIT_NODE, nodeTypes, type NodePayload } from "./UnitOpNode";

/** The connector: a route computed around the node plates, not xyflow's handle-to-handle curve. */
export const PFD_EDGE = { type: "pfd" as const };

const edgeTypes = { pfd: PfdEdge };

/** The plate sizes the router treats as obstacles; they are the layout tokens' fixed values. */
const UNIT_OBSTACLE = { width: 190, height: 168 };
const STREAM_OBSTACLE = { width: 160, height: 42 };

/**
 * **One line per point.** A handle already carrying an edge is full; a connection needs a free
 * source and a free target. The projection draws a spare point on every `many` inlet, so a mixer
 * can always take another feed.
 */
export function connectable(
  occupied: ReadonlySet<string>,
  source: string | null | undefined,
  target: string | null | undefined,
): boolean {
  return (
    typeof source === "string" &&
    typeof target === "string" &&
    !occupied.has(source) &&
    !occupied.has(target)
  );
}

export interface FlowsheetProps {
  catalogue: Catalogue | null;
  envelope: Envelope;
  selected: string | null;
  /**
   * The unit a reader wants a number in; a catalogue that has not loaded converts nothing.
   *
   * A readout is the one place a value is drawn *on* the drawing, so it converts like a table
   * does and by the same factor - `formatQuantity`'s.
   */
  units: Units;
  onSelect: (id: string | null) => void;
  onCommand: (command: EditorCommand) => void;
  /** A double-click opens the node's read-only overview. */
  onOverview: (id: string) => void;
}

export function Flowsheet({
  catalogue,
  envelope,
  selected,
  units,
  onSelect,
  onCommand,
  onOverview,
}: FlowsheetProps) {
  const [nodes, setNodes, onNodesChange] = useNodesState<Node>([]);
  const [edges, setEdges, onEdgesChange] = useEdgesState<Edge>([]);

  // **The envelope is the token, and the document is not enough of one.** A node's readout is
  // what the *run* reached, and a run changes no document — so a derivation keyed on the document
  // would keep the values it computed before Solve was pressed, and the canvas would never show a
  // number. The envelope's identity changes on every call, which is exactly the right cadence.
  const derived = useMemo(
    () => derive(envelope, catalogue, selected, units),
    [envelope, catalogue, selected, units],
  );

  useEffect(() => {
    setNodes(derived.nodes);
    setEdges(derived.edges);
  }, [derived, setNodes, setEdges]);

  /**
   * **One line per connection point.** A handle already carrying an edge is full: the next line
   * must find its own point, and the projection draws a spare on every `many` inlet so it can.
   * This is also the front-end half of the boundary rule - a feed with two consumers or a product
   * with two producers would be a hidden split or merge.
   */
  const occupied = useMemo(() => {
    const handles = new Set<string>();
    for (const edge of edges) {
      if (typeof edge.sourceHandle === "string") {
        handles.add(edge.sourceHandle);
      }
      if (typeof edge.targetHandle === "string") {
        handles.add(edge.targetHandle);
      }
    }
    return handles;
  }, [edges]);

  const connect = useCallback(
    (connection: Connection) => {
      // A handle's id *is* the endpoint the document writes, so a gesture needs no translation.
      if (connection.sourceHandle === null || connection.targetHandle === null) {
        return;
      }
      onCommand({
        command: "connect",
        from: connection.sourceHandle,
        to: connection.targetHandle,
      });
    },
    [onCommand],
  );

  return (
    <ReactFlow
      key={envelope.flowsheet.id}
      nodes={nodes}
      edges={edges}
      nodeTypes={nodeTypes}
      edgeTypes={edgeTypes}
      onNodesChange={onNodesChange}
      onEdgesChange={onEdgesChange}
      onConnect={connect}
      defaultEdgeOptions={{ type: PFD_EDGE.type }}
      connectionLineType={ConnectionLineType.SmoothStep}
      isValidConnection={(connection) =>
        connectable(occupied, connection.sourceHandle, connection.targetHandle)
      }
      snapToGrid
      snapGrid={[20, 20]}
      onNodeDoubleClick={(_event, node) => onOverview(node.id)}
      onNodeDragStop={(_event, node) =>
        onCommand({
          command: "set_position",
          node: node.id,
          x: Math.round(node.position.x),
          y: Math.round(node.position.y),
        })
      }
      // **Both keys, not xyflow's Backspace alone.** On a Windows/Linux keyboard the key a
      // reader reaches for is Delete; the default binding silently ignored it.
      deleteKeyCode={["Backspace", "Delete"]}
      onNodeClick={(_event, node) => onSelect(node.id)}
      onEdgeClick={(_event, edge) => onSelect(edge.id)}
      onPaneClick={() => onSelect(null)}
      onNodesDelete={(deleted) => {
        // **The defect this is here for.** xyflow's default Backspace goes through
        // `onNodesChange`, which updates this component's local state and sends no command - and
        // the next envelope re-derives the nodes, so the node reappears. A gesture that looks like
        // a delete and is silently undone is the "silently ignored" class, in the front end.
        for (const node of deleted) {
          const command = removeCommandFor(node.id);
          if (command !== null) {
            onCommand(command);
          }
        }
      }}
      onEdgesDelete={(deleted) => {
        for (const edge of deleted) {
          // A tear goes by its **name**: `disconnect` removes the first connection or recycle that
          // joins a pair, and a tear that shares its endpoints with a connection would take the
          // wrong one. The projection puts a tear's name in `data.path`.
          if (edge.data === undefined) {
            continue;
          }
          // **Read back, not asserted.** xyflow types `data` as `Record<string, unknown>` — it
          // carries the projection's object and does not know its shape — so this reads it again
          // rather than asserting it. An edge whose data is not the projection's is one a canvas
          // cannot make a command from, and dropping it is what `data === undefined` above does.
          const data = asEdgeData(edge.data);
          if (data === null) {
            continue;
          }
          onCommand(
            data.kind === "recycle"
              ? { command: "remove_recycle", stream: data.path }
              : { command: "disconnect", from: data.from, to: data.to },
          );
        }
      }}
      fitView
      proOptions={{ hideAttribution: true }}
    >
      {/* **xyflow's own grid, and not a gradient on the container.** A background painted on the
          pane is fixed to the viewport while the flowsheet moves under it, which is the one thing a
          PFD's grid must not do; this one pans and zooms with the drawing and takes its colour from
          the theme through `--xy-background-pattern-color`. */}
      <Background variant={BackgroundVariant.Dots} gap={20} size={1} />
      <Controls showInteractive={false} />
    </ReactFlow>
  );
}

/** The canvas's nodes and edges, from the graph document and the run. */
function derive(
  envelope: Envelope,
  catalogue: Catalogue | null,
  selected: string | null,
  units: Units,
): { nodes: Node[]; edges: Edge[] } {
  const streams = envelope.session?.streams ?? {};
  const forms = new Map((catalogue?.unit_ops ?? []).map((form) => [form.id, form]));

  // **A node's worst diagnostic and an edge's diagnostic are two readings of the same list.**
  // The canvas used one boolean for a node, so a warning looked like an error; and it did not mark
  // an edge at all. Both are derived here, not held as a second state.
  const severityById = new Map<string, "error" | "warning">();
  const severityByEdge = new Map<string, "error" | "warning">();
  for (const diagnostic of envelope.diagnostics) {
    const id = diagnosticNode(diagnostic.target);
    if (id !== "") {
      const current = severityById.get(id);
      if (current !== "error") {
        severityById.set(id, diagnostic.severity);
      }
    }
    if (diagnostic.target.kind === "edge") {
      severityByEdge.set(
        `${diagnostic.target.from}\u0000${diagnostic.target.to}`,
        diagnostic.severity,
      );
    }
  }

  const connected = new Set<string>();
  for (const edge of envelope.flowsheet.graph.edges) {
    if (edge.sourceHandle !== null) {
      connected.add(edge.sourceHandle);
    }
    if (edge.targetHandle !== null) {
      connected.add(edge.targetHandle);
    }
  }

  const nodes: Node[] = envelope.flowsheet.graph.nodes.map((node: GraphNode) => {
    const form = node.data.unit === undefined ? undefined : forms.get(node.data.unit);
    const payload: NodePayload = {
      graph: node,
      ...(form === undefined ? {} : { form }),
      readout: readoutOf(node, envelope, units),
      severity: severityById.get(node.id) ?? null,
      connected,
    };
    return {
      id: node.id,
      type: node.type === "unit_op" ? UNIT_NODE : STREAM_NODE,
      position: node.position,
      data: payload,
      selected: node.id === selected,
    };
  });

  const obstacles = nodes.map((node) => ({
    x: node.position.x,
    y: node.position.y,
    width: node.type === UNIT_NODE ? UNIT_OBSTACLE.width : STREAM_OBSTACLE.width,
    height: node.type === UNIT_NODE ? UNIT_OBSTACLE.height : STREAM_OBSTACLE.height,
  }));
  let returnLane = 0;

  const edges: Edge[] = envelope.flowsheet.graph.edges.map((edge) => {
    const stream = streams[edge.data.path];
    const recycle = edge.data.kind === "recycle";
    return {
      id: edge.id,
      type: PFD_EDGE.type,
      source: edge.source,
      sourceHandle: edge.sourceHandle,
      target: edge.target,
      targetHandle: edge.targetHandle,
      data: {
        ...edge.data,
        route: {
          obstacles,
          recycle,
          lane: recycle ? returnLane++ : 0,
        },
      },
      markerEnd: {
        type: MarkerType.ArrowClosed,
        width: 14,
        height: 14,
        color: edge.data.kind === "recycle" ? "var(--edge-recycle)" : "var(--edge)",
      },
      // A recycle is drawn moving, which is what a loop is: the tear's value is the previous
      // pass's, and the animation is the honest picture of that.
      animated: edge.data.kind === "recycle",
      // **A tear is stroked in its own hue**, and the class is where the stylesheet learns which
      // edge this is: `animated` is xyflow's own and says nothing about *why* the edge moves, and
      // `domAttributes` is typed as SVG attributes, which `data-*` is not. A diagnostic on the
      // pair adds the same reading in the error/warning colour.
      className: [
        edge.data.kind === "recycle" ? "tear" : "",
        severityByEdge.get(`${edge.data.from}\u0000${edge.data.to}`) === "error"
          ? "bad"
          : severityByEdge.get(`${edge.data.from}\u0000${edge.data.to}`) === "warning"
            ? "warn"
            : "",
      ]
        .filter((part) => part !== "")
        .join(" "),
      ...(stream === undefined
        ? {}
        : {
            label: `${formatQuantity(stream.T.magnitude_si, stream.T.unit, 4, units)} · ${formatQuantity(
              stream.P.magnitude_si,
              stream.P.unit,
              3,
              units,
            )}`,
          }),
    };
  });

  return { nodes, edges };
}

/** The lines a node shows: one per outlet the run has a value for. */
function readoutOf(node: GraphNode, envelope: Envelope, units: Units): NodePayload["readout"] {
  const streams = envelope.session?.streams ?? {};
  const ports = node.data.ports?.outlets ?? [];
  return ports
    .flatMap((port) => port.handles)
    .flatMap((path) => {
      const stream = streams[path];
      if (stream === undefined) {
        return [];
      }
      return [
        {
          port: path,
          text: `${path}  ${formatQuantity(stream.T.magnitude_si, stream.T.unit, 4, units)}  ${formatQuantity(
            stream.n.magnitude_si,
            stream.n.unit,
            3,
            units,
          )}`,
        },
      ];
    });
}

/** The node id a diagnostic is about, where it is about one. */
function diagnosticNode(target: { kind: string; [key: string]: unknown }): string {
  switch (target.kind) {
    case "node":
      return `${String(target.role)}:${String(target.id)}`;
    case "handle":
      return `${String(target.role)}:${String(target.node)}`;
    case "parameter":
      return `instance:${String(target.node)}`;
    default:
      return "";
  }
}
