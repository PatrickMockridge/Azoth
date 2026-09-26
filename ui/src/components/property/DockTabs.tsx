import type { KeyboardEvent } from "react";

/**
 * A row of tabs, with the keyboard contract `role="tablist"` promises a screen reader.
 *
 * **Left/Right/Home/End move the selection**, and only the selected tab is in the tab order - which
 * is what makes the strip one stop rather than six. A row of plain buttons with `aria-selected` on
 * them would *look* like tabs and be a set of six independent controls to anything reading the DOM.
 */
export function DockTabs({
  tabs,
  active,
  onTab,
  label,
}: {
  tabs: readonly { id: string; label: string; count?: number }[];
  active: string;
  onTab: (id: string) => void;
  /** What the strip is over, for the group's own accessible name. */
  label: string;
}) {
  const move = (event: KeyboardEvent<HTMLDivElement>): void => {
    const here = tabs.findIndex((tab) => tab.id === active);
    if (here === -1 || tabs.length === 0) {
      return;
    }
    const step =
      event.key === "ArrowRight" ? 1 : event.key === "ArrowLeft" ? -1 : 0;
    let next = here;
    if (step !== 0) {
      next = (here + step + tabs.length) % tabs.length;
    } else if (event.key === "Home") {
      next = 0;
    } else if (event.key === "End") {
      next = tabs.length - 1;
    } else {
      return;
    }
    event.preventDefault();
    const chosen = tabs[next];
    if (chosen !== undefined) {
      onTab(chosen.id);
    }
  };

  return (
    <div className="tabs" role="tablist" aria-label={label} onKeyDown={move}>
      {tabs.map((tab) => (
        <button
          key={tab.id}
          type="button"
          role="tab"
          id={`tab-${tab.id}`}
          aria-selected={tab.id === active}
          aria-controls={`sheet-${tab.id}`}
          tabIndex={tab.id === active ? 0 : -1}
          className={`tab${tab.id === active ? " active" : ""}`}
          onClick={() => onTab(tab.id)}
        >
          {tab.label}
          {/* **Out of the accessible name.** A tab is found by its label - `getByRole("tab", {
              name: "Messages" })` - and a count inside the button would make the name "Messages 3". */}
          {tab.count === undefined ? null : (
            <span className="count" aria-hidden="true">
              {tab.count}
            </span>
          )}
        </button>
      ))}
    </div>
  );
}
