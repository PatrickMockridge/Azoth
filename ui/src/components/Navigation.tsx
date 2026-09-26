/**
 * The navigator: every object in the flowsheet, and a way to select it.
 *
 * **It selects through the same `onSelect` the canvas does**, so the two are one selection rather
 * than two that agree — a row here lights the node there, and a click there lights the row here.
 *
 * **This is also the keyboard and screen-reader path to every object, which the canvas is not.**
 * xyflow makes a node focusable and nothing more: a person on the keyboard can select one and
 * cannot wire a port, place one or read a value off it. The rows below are buttons, so every object
 * is reachable by Tab and activated by Enter, and the one that is selected says so with
 * `aria-current` rather than only with a colour.
 */

import { navigationOf } from "../state/navigation";
import type { Catalogue, Envelope } from "../wire/types";

export interface NavigationProps {
  envelope: Envelope;
  catalogue: Catalogue | null;
  /** The selected object's id, which the canvas and this pane agree on. */
  selected: string | null;
  onSelect: (id: string) => void;
}

export function Navigation({ envelope, catalogue, selected, onSelect }: NavigationProps) {
  const groups = navigationOf(envelope, catalogue);

  return (
    <nav className="tree" aria-label="Objects in this flowsheet">
      <h2>Objects</h2>
      {groups.length === 0 ? (
        <p className="note">nothing placed yet</p>
      ) : (
        groups.map((group) => (
          <div className="tree-group" key={group.label}>
            <div className="tree-label">{group.label}</div>
            {group.rows.map((row) => (
              <button
                key={row.id}
                type="button"
                className="tree-row"
                // The id, because the row shows the name and a person debugging wants the address.
                title={row.id}
                aria-current={row.id === selected ? "true" : undefined}
                onClick={() => onSelect(row.id)}
              >
                <span className="name">{row.label}</span>
                {/* The tag is the data and the attribute is the selector, so a new word beside a
                    row needs no class name and cannot collide with one. */}
                {row.tag === null ? null : <span className="tag" data-tag={row.tag}>{row.tag}</span>}
              </button>
            ))}
          </div>
        ))
      )}
    </nav>
  );
}
