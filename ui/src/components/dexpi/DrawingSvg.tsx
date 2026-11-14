/**
 * A read-only drawing: the file's own primitives, at the file's own positions.
 *
 * **Not a canvas and not an editor.** The five primitives the decoder resolves are drawn as
 * SVG, the view can pan and zoom, and a click reports which represented object a group draws.
 * Nothing here edits a symbol or a coordinate; the drawing is a reading of the file.
 */

import { useRef, useState } from "react";

import { coverageIndex, groupKey, groupStatus, viewBoxOf } from "../../state/dexpi";
import type {
  DexpiView,
  Drawing,
  DrawingColor,
  DrawingPoint,
  DrawingPrimitive,
  DrawingStroke,
} from "../../wire/types";

interface Zoom {
  scale: number;
  tx: number;
  ty: number;
}

const FIT: Zoom = { scale: 1, tx: 0, ty: 0 };

export function DrawingSvg({
  drawing,
  view,
  selected,
  onSelect,
}: {
  drawing: Drawing;
  view: DexpiView;
  selected: string | null;
  onSelect: (key: string | null) => void;
}) {
  const [zoom, setZoom] = useState<Zoom>(FIT);
  const dragging = useRef<{ x: number; y: number } | null>(null);
  const index = coverageIndex(view);

  const ordered = drawing.groups
    .map((group, position) => ({ group, position }))
    .sort((left, right) => layer(left.group.kind) - layer(right.group.kind));

  return (
    <div className="dexpi-canvas">
      <div className="dexpi-zoom">
        <button type="button" onClick={() => setZoom(FIT)}>
          Fit
        </button>
        <button
          type="button"
          onClick={() => setZoom((current) => ({ ...current, scale: current.scale * 1.1 }))}
        >
          +
        </button>
        <button
          type="button"
          onClick={() => setZoom((current) => ({ ...current, scale: current.scale * 0.9 }))}
        >
          −
        </button>
      </div>
      <svg
        className="dexpi-svg"
        viewBox={viewBoxOf(drawing)}
        role="img"
        aria-label={drawing.name}
        onWheel={(event) => {
          event.preventDefault();
          setZoom((current) => ({
            ...current,
            scale: current.scale * (event.deltaY < 0 ? 1.1 : 0.9),
          }));
        }}
        onPointerDown={(event) => {
          dragging.current = { x: event.clientX, y: event.clientY };
          event.currentTarget.setPointerCapture(event.pointerId);
        }}
        onPointerMove={(event) => {
          if (dragging.current === null) {
            return;
          }
          const dx = event.clientX - dragging.current.x;
          const dy = event.clientY - dragging.current.y;
          dragging.current = { x: event.clientX, y: event.clientY };
          setZoom((current) => ({ ...current, tx: current.tx + dx, ty: current.ty + dy }));
        }}
        onPointerUp={() => {
          dragging.current = null;
        }}
        onClick={() => onSelect(null)}
      >
        <g transform={`translate(${zoom.tx} ${zoom.ty}) scale(${zoom.scale})`}>
          {ordered.map(({ group, position }) => {
            const key = groupKey(group, position);
            const status = groupStatus(index, group);
            const selection = selected === null ? undefined : index.get(selected);
            const chosen =
              selected === key ||
              (selection !== undefined &&
                (selection.id === key || selection.azoth === key));
            return (
              <g
                key={key}
                className={`dexpi-group status-${status}${chosen ? " selected" : ""}`}
                onClick={(event) => {
                  event.stopPropagation();
                  onSelect(key);
                }}
              >
                {group.primitives.map((primitive, primitiveIndex) => (
                  <PrimitiveShape key={primitiveIndex} primitive={primitive} />
                ))}
              </g>
            );
          })}
        </g>
      </svg>
    </div>
  );
}

function layer(kind: string): number {
  if (kind === "edge") {
    return 0;
  }
  if (kind === "node") {
    return 1;
  }
  return 2;
}

function PrimitiveShape({ primitive }: { primitive: DrawingPrimitive }) {
  switch (primitive.kind) {
    case "polyline":
      return (
        <polyline
          points={points(primitive.points)}
          fill="none"
          {...strokeAttributes(primitive.stroke)}
        />
      );
    case "polygon":
      return (
        <polygon
          points={points(primitive.points)}
          fill={fillColor(primitive.fill)}
          {...strokeAttributes(primitive.stroke)}
        />
      );
    case "ellipse":
      return (
        <ellipse
          cx={primitive.center.x}
          cy={primitive.center.y}
          rx={primitive.radius_x}
          ry={primitive.radius_y}
          fill={fillColor(primitive.fill)}
          transform={rotate(primitive.rotation, primitive.center)}
          {...strokeAttributes(primitive.stroke)}
        />
      );
    case "ellipse_arc":
      return (
        <path
          d={arcPath(primitive)}
          fill="none"
          transform={rotate(primitive.rotation, primitive.center)}
          {...strokeAttributes(primitive.stroke)}
        />
      );
    case "text":
      return (
        <text
          x={primitive.position.x}
          y={primitive.position.y}
          fontSize={primitive.size}
          fill={rgb(primitive.color)}
          textAnchor={textAnchor(primitive.alignment)}
          dominantBaseline={textBaseline(primitive.alignment)}
          transform={rotate(primitive.rotation, primitive.position)}
        >
          {primitive.text}
        </text>
      );
  }
}

function points(values: DrawingPoint[]): string {
  return values.map((point) => `${point.x},${point.y}`).join(" ");
}

function strokeAttributes(stroke: DrawingStroke): {
  stroke: string;
  strokeWidth: number;
  strokeDasharray?: string;
} {
  const dash =
    stroke.dash === "Solid"
      ? undefined
      : stroke.dash === "Dot"
        ? "1 3"
        : stroke.dash === "DashShortDash"
          ? "6 3 2 3"
          : "6 3";
  return {
    stroke: rgb(stroke.color),
    strokeWidth: stroke.width,
    ...(dash === undefined ? {} : { strokeDasharray: dash }),
  };
}

function rgb(color: DrawingColor): string {
  return `rgb(${color.r} ${color.g} ${color.b})`;
}

function fillColor(fill: string): string {
  return fill === "Transparent" ? "none" : "#e5e7eb";
}

function rotate(degrees: number, centre: DrawingPoint): string | undefined {
  return degrees === 0 ? undefined : `rotate(${degrees} ${centre.x} ${centre.y})`;
}

function textAnchor(alignment: string): "start" | "middle" | "end" {
  if (alignment.startsWith("Left")) {
    return "start";
  }
  return alignment.startsWith("Right") ? "end" : "middle";
}

function textBaseline(alignment: string): "hanging" | "middle" | "text-after-edge" {
  if (alignment.endsWith("Top")) {
    return "hanging";
  }
  return alignment.endsWith("Bottom") ? "text-after-edge" : "middle";
}

function arcPath(primitive: {
  center: DrawingPoint;
  radius_x: number;
  radius_y: number;
  start_angle: number;
  end_angle: number;
}): string {
  const start = onEllipse(primitive, primitive.start_angle);
  const end = onEllipse(primitive, primitive.end_angle);
  const span = Math.abs(primitive.end_angle - primitive.start_angle);
  const large = span > 180 ? 1 : 0;
  const sweep = primitive.end_angle >= primitive.start_angle ? 1 : 0;
  return [
    `M ${start.x} ${start.y}`,
    `A ${primitive.radius_x} ${primitive.radius_y} 0 ${large} ${sweep} ${end.x} ${end.y}`,
  ].join(" ");
}

function onEllipse(
  primitive: { center: DrawingPoint; radius_x: number; radius_y: number },
  degrees: number,
): DrawingPoint {
  const radians = (degrees * Math.PI) / 180;
  return {
    x: primitive.center.x + primitive.radius_x * Math.cos(radians),
    y: primitive.center.y + primitive.radius_y * Math.sin(radians),
  };
}
