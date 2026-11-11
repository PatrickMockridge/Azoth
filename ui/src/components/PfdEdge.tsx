/**
 * The PFD edge: a custom orthogonal route rather than xyflow's smooth step.
 *
 * The route is computed by `pfd/route` from the node rectangles the projection hands it; the
 * built-in smooth step knows only the two handle points and will happily run a recycle back over a
 * plate or along its own outbound line. The label is rendered by xyflow's label portal so the
 * stream's temperature and pressure stay where the old edge put them.
 */

import {
  BaseEdge,
  EdgeLabelRenderer,
  getSmoothStepPath,
  type EdgeProps,
  type Position,
} from "@xyflow/react";

import { route, type Rect, type Side } from "../pfd/route";

interface RouteData {
  obstacles: Rect[];
  recycle: boolean;
  lane: number;
}

function isRecord(value: unknown): value is Record<string, unknown> {
  return typeof value === "object" && value !== null;
}

function rectOf(value: unknown): Rect | null {
  if (!isRecord(value)) {
    return null;
  }
  const { x, y, width, height } = value;
  return typeof x === "number" &&
    typeof y === "number" &&
    typeof width === "number" &&
    typeof height === "number"
    ? { x, y, width, height }
    : null;
}

/** The projection's routing metadata, read back rather than asserted. */
function routeDataOf(value: unknown): RouteData | null {
  if (!isRecord(value) || !isRecord(value.route)) {
    return null;
  }
  const raw = value.route;
  const obstacles = Array.isArray(raw.obstacles)
    ? raw.obstacles.map(rectOf).filter((rect): rect is Rect => rect !== null)
    : [];
  if (
    obstacles.length === 0 ||
    typeof raw.recycle !== "boolean" ||
    typeof raw.lane !== "number"
  ) {
    return null;
  }
  return { obstacles, recycle: raw.recycle, lane: raw.lane };
}

function sideOf(position: Position): Side {
  switch (position) {
    case "left":
      return "left";
    case "top":
      return "top";
    case "bottom":
      return "bottom";
    default:
      return "right";
  }
}

export function PfdEdge({
  id,
  sourceX,
  sourceY,
  targetX,
  targetY,
  sourcePosition,
  targetPosition,
  markerEnd,
  style,
  label,
  data,
}: EdgeProps) {
  const routing = routeDataOf(data);
  const routed =
    routing === null
      ? null
      : route({
          sx: sourceX,
          sy: sourceY,
          tx: targetX,
          ty: targetY,
          sourceSide: sideOf(sourcePosition),
          targetSide: sideOf(targetPosition),
          obstacles: routing.obstacles,
          recycle: routing.recycle,
          lane: routing.lane,
        });

  const [fallback, fallbackX, fallbackY] = getSmoothStepPath({
    sourceX,
    sourceY,
    targetX,
    targetY,
    sourcePosition,
    targetPosition,
    borderRadius: 4,
  });
  const path = routed?.path ?? fallback;
  const mid = routed?.mid ?? { x: fallbackX, y: fallbackY };

  return (
    <>
      <BaseEdge
        id={id}
        path={path}
        {...(markerEnd === undefined ? {} : { markerEnd })}
        style={style}
        interactionWidth={22}
      />
      {label === undefined || label === null ? null : (
        <EdgeLabelRenderer>
          <div
            className="pfd-edge-label"
            style={{ transform: `translate(-50%, -50%) translate(${mid.x}px, ${mid.y}px)` }}
          >
            {label}
          </div>
        </EdgeLabelRenderer>
      )}
    </>
  );
}
