/**
 * The DEXPI view's readings: what a status means, and how a drawing is framed.
 *
 * **Derived, never held.** Every function takes the view the door answered with; nothing here
 * keeps a second copy of a drawing or a coverage table. The status order is the legend's, and
 * the index is the only join between a drawn group and the object it represents.
 */

import type {
  CoverageStatus,
  DexpiView,
  Drawing,
  DrawingGroup,
  ObjectCoverage,
} from "../wire/types";

/** The order the legend and the filters read in, from modelled to furniture. */
export const STATUS_ORDER: readonly CoverageStatus[] = [
  "modelled",
  "lossy",
  "ambiguous",
  "unmapped",
  "topology",
  "control",
  "drawing_only",
];

const LABELS: Record<CoverageStatus, string> = {
  modelled: "Modelled",
  lossy: "Lossy",
  ambiguous: "Ambiguous",
  unmapped: "Unmapped",
  topology: "Topology",
  control: "Control",
  drawing_only: "Drawing only",
};

const EXPLANATIONS: Record<CoverageStatus, string> = {
  modelled: "becomes the azoth unit named beside it",
  lossy: "maps to the nearest unit, with something lost",
  ambiguous: "several units share this DEXPI type; the importer chooses one",
  unmapped: "no azoth unit, and it stays dexpi.unmapped",
  topology: "carried as a connection, line or boundary, not a unit",
  control: "carried as instrumentation, not simulated",
  drawing_only: "graphics with no conceptual object behind them",
};

export function statusLabel(status: CoverageStatus): string {
  return LABELS[status];
}

export function statusExplanation(status: CoverageStatus): string {
  return EXPLANATIONS[status];
}

/** The `ObjectCoverage` a group draws, by object id or by the azoth key it became. */
export function coverageIndex(view: DexpiView): Map<string, ObjectCoverage> {
  const index = new Map<string, ObjectCoverage>();
  for (const object of view.objects) {
    index.set(object.id, object);
    if (object.azoth !== null) {
      index.set(object.azoth, object);
    }
  }
  return index;
}

export function groupStatus(
  index: Map<string, ObjectCoverage>,
  group: DrawingGroup,
): CoverageStatus {
  return index.get(group.represents ?? "")?.status ?? "drawing_only";
}

/** The key a group is selected by: what it represents, else its own id. */
export function groupKey(group: DrawingGroup, position: number): string {
  return group.represents ?? group.id ?? `drawing-${position}`;
}

/** The padding around a drawing in its viewBox, in drawing units. */
const PAD = 24;

export function viewBoxOf(drawing: Drawing): string {
  const { min_x, min_y, max_x, max_y } = drawing.bounds;
  return `${min_x - PAD} ${min_y - PAD} ${max_x - min_x + PAD * 2} ${max_y - min_y + PAD * 2}`;
}

/** What the drawing would become, in one sentence a reader can act on. */
export function summaryOf(view: DexpiView): string {
  const scaffold = view.scaffold;
  if (scaffold === null) {
    return "This file was refused as a simulation; the drawing and its diagnostics remain.";
  }
  const parts = [
    `${scaffold.instances} instance${scaffold.instances === 1 ? "" : "s"}`,
    `${scaffold.connections} connection${scaffold.connections === 1 ? "" : "s"}`,
    `${scaffold.instruments} instrument${scaffold.instruments === 1 ? "" : "s"}`,
  ];
  if (view.counts.unmapped > 0) {
    parts.push(`${view.counts.unmapped} unmapped`);
  }
  if (scaffold.check_errors > 0) {
    parts.push(`${scaffold.check_errors} check error${scaffold.check_errors === 1 ? "" : "s"}`);
  }
  return `This drawing would become ${parts.join(", ")}.`;
}

/** Whether an object can be exported back to the same profile. */
export function exportableText(object: ObjectCoverage): string {
  if (!object.imported) {
    return "not imported";
  }
  return object.exportable ? "writable" : "read-only";
}
