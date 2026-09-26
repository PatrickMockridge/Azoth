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
  Controls,
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
import type { Catalogue, Envelope, GraphEdge, GraphNode } from "../wire/types";
import { STREAM_NODE, UNIT_NODE, nodeTypes, type NodePayload } from "./UnitOpNode";

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
}

export function Flowsheet({
  catalogue,
  envelope,
  selected,
  units,
  onSelect,
  onCommand,
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
      nodes={nodes}
      edges={edges}
      nodeTypes={nodeTypes}
      onNodesChange={onNodesChange}
      onEdgesChange={onEdgesChange}
      onConnect={connect}
      onNodeDragStop={(_event, node) =>
        onCommand({
          command: "set_position",
          node: node.id,
          x: Math.round(node.position.x),
          y: Math.round(node.position.y),
        })
      }
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
          // xyflow types `data` as `Record<string, unknown>`; this is the projection's own shape,
          // which `GraphEdge` mirrors and `fixtures.test.ts` holds to the emitted document.
          const data = edge.data as unknown as GraphEdge["data"];
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
  const flagged = new Set(
    envelope.diagnostics.map((diagnostic) => diagnosticNode(diagnostic.target)),
  );

  const nodes: Node[] = envelope.flowsheet.graph.nodes.map((node: GraphNode) => {
    const form = node.data.unit === undefined ? undefined : forms.get(node.data.unit);
    const payload: NodePayload = {
      graph: node,
      ...(form === undefined ? {} : { form }),
      readout: readoutOf(node, envelope, units),
      bad: flagged.has(node.id),
    };
    return {
      id: node.id,
      type: node.type === "unit_op" ? UNIT_NODE : STREAM_NODE,
      position: node.position,
      data: payload,
      selected: node.id === selected,
    };
  });

  const edges: Edge[] = envelope.flowsheet.graph.edges.map((edge) => {
    const stream = streams[edge.data.path];
    return {
      id: edge.id,
      source: edge.source,
      sourceHandle: edge.sourceHandle,
      target: edge.target,
      targetHandle: edge.targetHandle,
      // A recycle is drawn moving, which is what a loop is: the tear's value is the previous
      // pass's, and the animation is the honest picture of that.
      animated: edge.data.kind === "recycle",
      // **A tear is stroked in its own hue**, and the class is where the stylesheet learns which
      // edge this is: `animated` is xyflow's own and says nothing about *why* the edge moves, and
      // `domAttributes` is typed as SVG attributes, which `data-*` is not.
      className: edge.data.kind === "recycle" ? "tear" : "",
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
