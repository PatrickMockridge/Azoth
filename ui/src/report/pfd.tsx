/**
 * The flowsheet as a printable drawing.
 *
 * **A schematic from the graph, not a capture of the canvas.** The editor's canvas is a live
 * xyflow surface with handles, hit areas and a viewport; a report needs the same machines and the
 * same connections on a page. This reads `position`, `name` and the edge endpoints and draws
 * boxes and lines, so a print or an SVG export carries no browser interaction and no dependency.
 */

import type { RefObject } from "react";

import type { Envelope } from "../wire/types";

const UNIT_WIDTH = 128;
const UNIT_HEIGHT = 54;
const PAD = 24;

/** The SVG namespace, which is a URL and not a magnitude. */
const SVG_NAMESPACE = "http://www.w3.org/2000/svg"; // ui-ok: an XML namespace is a URL, not a quantity

/** The report's drawing of the document. */
export function Pfd({
  envelope,
  svgRef,
}: {
  envelope: Envelope;
  svgRef?: RefObject<SVGSVGElement>;
}) {
  const nodes = envelope.flowsheet.graph.nodes;
  const positions = new Map(nodes.map((node) => [node.id, node.position]));
  const maxX = nodes.reduce((max, node) => Math.max(max, node.position.x), 0);
  const maxY = nodes.reduce((max, node) => Math.max(max, node.position.y), 0);
  const width = Math.max(maxX + UNIT_WIDTH + PAD, UNIT_WIDTH + PAD * 2);
  const height = Math.max(maxY + UNIT_HEIGHT + PAD, UNIT_HEIGHT + PAD * 2);

  return (
    <svg
      ref={svgRef}
      className="pfd"
      viewBox={`0 0 ${String(width)} ${String(height)}`}
      role="img"
      aria-label="Process flow diagram"
      xmlns={SVG_NAMESPACE}
    >
      <defs>
        <marker id="pfd-arrow" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" markerHeight="6" orient="auto">
          <path d="M 0 0 L 10 5 L 0 10 z" className="arrow" />
        </marker>
      </defs>

      {envelope.flowsheet.graph.edges.map((edge) => {
        const from = positions.get(edge.source);
        const to = positions.get(edge.target);
        if (from === undefined || to === undefined) {
          return null;
        }
        const x1 = from.x + UNIT_WIDTH;
        const y1 = from.y + UNIT_HEIGHT / 2;
        const x2 = to.x;
        const y2 = to.y + UNIT_HEIGHT / 2;
        const mid = (x1 + x2) / 2;
        const path =
          x2 >= x1
            ? `M ${String(x1)} ${String(y1)} L ${String(x2)} ${String(y2)}`
            : `M ${String(x1)} ${String(y1)} L ${String(mid)} ${String(y1)} L ${String(mid)} ${String(y2)} L ${String(x2)} ${String(y2)}`;
        return (
          <g key={edge.id}>
            <title>{edge.data.path}</title>
            <path
              d={path}
              className={edge.data.kind === "recycle" ? "connection recycle" : "connection"}
              markerEnd="url(#pfd-arrow)"
            />
          </g>
        );
      })}

      {nodes.map((node) => (
        <g key={node.id} transform={`translate(${String(node.position.x)} ${String(node.position.y)})`}>
          <rect
            width={UNIT_WIDTH}
            height={UNIT_HEIGHT}
            rx="4"
            className={`machine ${node.role}`}
          />
          <text x={UNIT_WIDTH / 2} y={UNIT_HEIGHT / 2 - 2} textAnchor="middle" className="name">
            {node.data.name}
          </text>
          <text x={UNIT_WIDTH / 2} y={UNIT_HEIGHT / 2 + 14} textAnchor="middle" className="kind">
            {node.data.unit_name ?? node.role}
          </text>
        </g>
      ))}
    </svg>
  );
}
