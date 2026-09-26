import { useEffect, useRef, useState, type ReactNode } from "react";

import { DockTabs } from "../property/DockTabs";

/** The shortest a dock may be drawn, and the tallest - both in pixels. */
const MIN = 80;
const MAX = 460;

/** Inside the range, whatever the drag or the key asked for. */
function clamp(size: number): number {
  return Math.min(MAX, Math.max(MIN, size));
}

/**
 * The bottom dock: the workbook and the messages, under the drawing rather than beside it.
 *
 * **It is resizable by pointer *and* by keyboard**, and the second is the one that is easy to
 * forget: a splitter that only drags is a size a keyboard cannot reach, so it carries
 * `role="separator"` with `aria-valuenow` and the arrow keys, which is what that role promises.
 *
 * **The drag listens on the window**, because the pointer leaves the six-pixel strip immediately -
 * the same reason `d3-drag` does it that way for the canvas.
 *
 * **Nothing here measures anything.** The size is `App`'s state and not a `ResizeObserver`'s: the
 * stub for that in jsdom is a no-op, so an observer-driven dock would be nought pixels tall in
 * every test and would prove nothing about the one place it is used.
 */
export function Dock({
  tabs,
  active,
  onTab,
  size,
  onSize,
  children,
}: {
  tabs: readonly { id: string; label: string; count?: number }[];
  active: string;
  onTab: (id: string) => void;
  size: number;
  onSize: (size: number) => void;
  children: ReactNode;
}) {
  const drag = useRef<{ y: number; size: number } | null>(null);
  const [dragging, setDragging] = useState(false);

  useEffect(() => {
    if (!dragging) {
      return;
    }
    const move = (event: MouseEvent): void => {
      const from = drag.current;
      if (from === null) {
        return;
      }
      // Up is taller: the dock grows away from the splitter, which sits above it.
      onSize(clamp(from.size - (event.clientY - from.y)));
    };
    const release = (): void => {
      drag.current = null;
      setDragging(false);
    };
    window.addEventListener("mousemove", move);
    window.addEventListener("mouseup", release);
    return () => {
      window.removeEventListener("mousemove", move);
      window.removeEventListener("mouseup", release);
    };
  }, [dragging, onSize]);

  return (
    <div className={`dock${dragging ? " dragging" : ""}`} style={{ height: size }}>
      <div
        className="splitter"
        role="separator"
        aria-orientation="horizontal"
        aria-label="resize the dock"
        aria-valuenow={size}
        aria-valuemin={MIN}
        aria-valuemax={MAX}
        tabIndex={0}
        onMouseDown={(event) => {
          drag.current = { y: event.clientY, size };
          setDragging(true);
        }}
        onKeyDown={(event) => {
          const step = event.key === "ArrowUp" ? 20 : event.key === "ArrowDown" ? -20 : 0;
          if (step === 0) {
            return;
          }
          event.preventDefault();
          onSize(clamp(size + step));
        }}
      />
      <DockTabs tabs={tabs} active={active} onTab={onTab} label="docked panels" />
      <div className="dock-body">{children}</div>
    </div>
  );
}
