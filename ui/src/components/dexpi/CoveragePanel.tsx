/**
 * What the drawing would become, as a list a reader can act on.
 *
 * **The statuses are the answer, not a decoration.** A modelled unit, a lossy mapping, an
 * ambiguous type, an unmapped class, a connection, an instrument and a label are seven
 * different facts about a file, and a viewer that showed one colour would be telling a reader
 * "this is a drawing" and nothing else.
 */

import {
  STATUS_ORDER,
  statusExplanation,
  statusLabel,
  summaryOf,
} from "../../state/dexpi";
import type { CoverageCounts, CoverageStatus, DexpiView, ObjectCoverage } from "../../wire/types";

export function CoveragePanel({
  view,
  selected,
  onSelect,
  filter,
  onFilter,
  showNotDrawn,
  onShowNotDrawn,
}: {
  view: DexpiView;
  selected: string | null;
  onSelect: (key: string | null) => void;
  filter: CoverageStatus | "all";
  onFilter: (filter: CoverageStatus | "all") => void;
  showNotDrawn: boolean;
  onShowNotDrawn: (show: boolean) => void;
}) {
  const objects = view.objects.filter((object) => {
    if (!showNotDrawn && !object.drawn) {
      return false;
    }
    return filter === "all" || object.status === filter;
  });

  return (
    <aside className="dexpi-coverage" aria-label="Model coverage">
      <p className="dexpi-summary">{summaryOf(view)}</p>

      <div className="dexpi-counts" role="group" aria-label="Coverage counts">
        <button
          type="button"
          className={filter === "all" ? "active" : ""}
          onClick={() => onFilter("all")}
        >
          <strong>{view.counts.semantic}</strong>
          <span>objects</span>
        </button>
        {STATUS_ORDER.map((status) => (
          <button
            key={status}
            type="button"
            className={`status-${status}${filter === status ? " active" : ""}`}
            title={statusExplanation(status)}
            onClick={() => onFilter(status)}
          >
            <strong>{countOf(view.counts, status)}</strong>
            <span>{statusLabel(status)}</span>
          </button>
        ))}
      </div>

      <label className="dexpi-toggle">
        <input
          type="checkbox"
          checked={showNotDrawn}
          onChange={(event) => onShowNotDrawn(event.target.checked)}
        />
        Show objects the drawing does not draw
      </label>

      {objects.length === 0 ? (
        <p className="dexpi-empty">Nothing in this category.</p>
      ) : (
        <ul className="dexpi-objects">
          {objects.map((object) => (
            <li key={`${object.id}-${object.status}`}>
              <ObjectRow
                object={object}
                selected={
                  selected !== null && (selected === object.id || selected === object.azoth)
                }
                onSelect={() => onSelect(object.id)}
              />
            </li>
          ))}
        </ul>
      )}

      {view.diagnostics.length > 0 ? (
        <details className="dexpi-diagnostics">
          <summary>
            {view.diagnostics.length} diagnostic
            {view.diagnostics.length === 1 ? "" : "s"}
          </summary>
          <ul>
            {view.diagnostics.map((diagnostic, index) => (
              <li key={`${diagnostic.code}-${index}`} className={`severity-${diagnostic.severity}`}>
                <strong>{diagnostic.code}</strong> {diagnostic.message}
                {diagnostic.subject === "" ? null : (
                  <span className="dexpi-subject"> ({diagnostic.subject})</span>
                )}
              </li>
            ))}
          </ul>
        </details>
      ) : null}
    </aside>
  );
}

function ObjectRow({
  object,
  selected,
  onSelect,
}: {
  object: ObjectCoverage;
  selected: boolean;
  onSelect: () => void;
}) {
  return (
    <button
      type="button"
      className={`dexpi-object status-${object.status}${selected ? " selected" : ""}`}
      onClick={onSelect}
    >
      <span className="dexpi-object-label">{object.label}</span>
      <span className={`dexpi-badge status-${object.status}`}>{statusLabel(object.status)}</span>
      {object.unit === null ? null : <span className="dexpi-unit">{object.unit}</span>}
      <span className="dexpi-reason">{object.reason}</span>
    </button>
  );
}

function countOf(counts: CoverageCounts, status: CoverageStatus): number {
  switch (status) {
    case "modelled":
      return counts.modelled;
    case "lossy":
      return counts.lossy;
    case "ambiguous":
      return counts.ambiguous;
    case "unmapped":
      return counts.unmapped;
    case "topology":
      return counts.topology;
    case "control":
      return counts.control;
    case "drawing_only":
      return counts.drawing_only;
  }
}
