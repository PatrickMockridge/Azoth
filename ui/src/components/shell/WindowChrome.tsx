/**
 * The title/body split every shell region exposes.
 *
 * **The title is outside the scrolling body.** A Macintosh window and an XP pane both need a
 * caption that stays put while the content under it scrolls; a pseudo-element on the scrolling
 * element would scroll away. The base stylesheet hides the title and makes the body
 * `display:contents`, so the wrapper is inert until a layout opts in, and the existing grid and
 * flex rules see exactly the children they saw before.
 *
 * The grip is the resize handle for the Macintosh window manager; it is hidden everywhere else.
 */

import type { ReactNode } from "react";

export function WindowChrome({ title, children }: { title: string; children: ReactNode }) {
  return (
    <>
      <div className="window-title" aria-hidden="true">
        {title}
      </div>
      <div className="window-body">{children}</div>
      <span className="window-grip" aria-hidden="true" />
    </>
  );
}
