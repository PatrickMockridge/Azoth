import type { ReactNode } from "react";

/**
 * One labelled cluster of controls, with a hairline before it.
 *
 * A group and not a tab, and the difference is a test: six of the eleven rendered cases address
 * New, Open, Demo, Save and Solve by their exact accessible name, and a control hidden behind a
 * ribbon tab is a control `getByRole` does not find. The five document commands stay in one
 * always-visible row, so a later tabbed ribbon can hold everything else without moving them.
 */
export function RibbonGroup({ label, children }: { label: string; children: ReactNode }) {
  return (
    <div className="ribbon-group" role="group" aria-label={label}>
      <div className="ribbon-controls">{children}</div>
      <div className="ribbon-label">{label}</div>
    </div>
  );
}
