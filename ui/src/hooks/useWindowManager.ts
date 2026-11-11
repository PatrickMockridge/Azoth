/**
 * The Macintosh desktop's window manager.
 *
 * **The stylesheet owns the resting geometry; this only writes a drag.** Each region has a
 * percentage default in `mac.css`, so windows appear without JavaScript, a newly loaded
 * calculation sheet gets a position, and `reset()` is just dropping the inline values.
 *
 * **A drag writes CSS variables, not React state.** Sixty renders a second would re-render the
 * whole editor for a moving shadow; the variables are read by the same stylesheet rule and the
 * ref is committed only on release. This is the Dock splitter's pattern (`Dock.tsx`) applied to
 * four windows.
 *
 * The title is the drag handle and the grip is the resize handle; both are real elements, so a
 * keyboard can reach them (Arrow moves, Shift+Arrow resizes) and a title bar does not scroll away.
 */

import { useCallback, useEffect, useRef, type RefObject } from "react";

/** The desktop top inset: the menu bar's height, read from the stylesheet where possible. */
const MENU_FALLBACK = 24;

/** A window never shrinks past this, so a title and a grip remain reachable. */
const MIN_W = 160;
const MIN_H = 80;

interface Drag {
  region: HTMLElement;
  mode: "move" | "resize";
  startX: number;
  startY: number;
  rect: DOMRect;
}

function clamp(value: number, low: number, high: number): number {
  return Math.min(Math.max(value, low), Math.max(low, high));
}

/** The desktop size, with a jsdom-friendly fallback when no element has a measured box. */
function desktop(root: HTMLElement | null): { width: number; height: number; top: number } {
  const rect = root?.getBoundingClientRect();
  const width = rect !== undefined && rect.width > 0 ? rect.width : window.innerWidth;
  const height = rect !== undefined && rect.height > 0 ? rect.height : window.innerHeight;
  const style = root === null ? "" : getComputedStyle(root).getPropertyValue("--menubar-h");
  const top = Number.parseFloat(style) || MENU_FALLBACK;
  return { width, height, top };
}

/** The variable names one region's geometry lives under; each window has its own defaults. */
function geometry(region: HTMLElement): { x: string; y: string; w: string; h: string } {
  const id = region.dataset.window ?? "window";
  return { x: `--${id}-x`, y: `--${id}-y`, w: `--${id}-w`, h: `--${id}-h` };
}

/** Write a region's geometry into the inline variables the Mac rules read. */
function place(region: HTMLElement, left: number, top: number, width: number, height: number): void {
  const names = geometry(region);
  region.style.setProperty(names.x, `${Math.round(left)}px`);
  region.style.setProperty(names.y, `${Math.round(top)}px`);
  region.style.setProperty(names.w, `${Math.round(width)}px`);
  region.style.setProperty(names.h, `${Math.round(height)}px`);
}

/** Keep a window inside the desktop after a drag, a resize or a viewport change. */
function contain(region: HTMLElement, root: HTMLElement | null): void {
  const rect = region.getBoundingClientRect();
  if (rect.width === 0 || rect.height === 0) {
    return;
  }
  const bounds = desktop(root);
  const width = clamp(rect.width, MIN_W, bounds.width);
  const height = clamp(rect.height, MIN_H, bounds.height - bounds.top);
  place(
    region,
    clamp(rect.left, 0, bounds.width - width),
    clamp(rect.top, bounds.top, bounds.height - height),
    width,
    height,
  );
}

export function useWindowManager(
  root: RefObject<HTMLElement | null>,
  { active }: { active: boolean },
): { reset: () => void } {
  const drag = useRef<Drag | null>(null);
  const z = useRef(10);

  const reset = useCallback(() => {
    const element = root.current;
    if (element === null) {
      return;
    }
    for (const region of element.querySelectorAll<HTMLElement>("[data-window]")) {
      const names = geometry(region);
      for (const name of [names.x, names.y, names.w, names.h, "--wz"]) {
        region.style.removeProperty(name);
      }
    }
  }, [root]);

  useEffect(() => {
    const element = root.current;
    if (!active || element === null) {
      return;
    }

    const onMove = (event: MouseEvent): void => {
      const current = drag.current;
      if (current === null) {
        return;
      }
      const bounds = desktop(element);
      const dx = event.clientX - current.startX;
      const dy = event.clientY - current.startY;
      if (current.mode === "move") {
        const width = current.rect.width;
        const height = current.rect.height;
        place(
          current.region,
          clamp(current.rect.left + dx, 0, bounds.width - width),
          clamp(current.rect.top + dy, bounds.top, bounds.height - height),
          width,
          height,
        );
      } else {
        const width = clamp(current.rect.width + dx, MIN_W, bounds.width - current.rect.left);
        const height = clamp(current.rect.height + dy, MIN_H, bounds.height - current.rect.top);
        place(current.region, current.rect.left, current.rect.top, width, height);
      }
    };

    const onUp = (): void => {
      drag.current = null;
      window.removeEventListener("mousemove", onMove);
      window.removeEventListener("mouseup", onUp);
    };

    const onDown = (event: MouseEvent): void => {
      if (event.button !== 0) {
        return;
      }
      const target = event.target;
      if (!(target instanceof Element)) {
        return;
      }
      const region = target.closest<HTMLElement>("[data-window]");
      if (region === null) {
        return;
      }
      z.current += 1;
      if (z.current > 18) {
        // Stay below the overview/report layer (20): cycle the other windows back down rather
        // than letting a much-clicked window climb over the modal layer.
        z.current = 11;
        for (const other of element.querySelectorAll<HTMLElement>("[data-window]")) {
          if (other !== region) {
            other.style.setProperty("--wz", "10");
          }
        }
      }
      region.style.setProperty("--wz", String(z.current));

      const interactive = target.closest("input, select, button, textarea, a, [contenteditable]");
      if (interactive !== null) {
        return;
      }
      const mode =
        target.closest(".window-grip") !== null
          ? "resize"
          : target.closest(".window-title") !== null
            ? "move"
            : null;
      if (mode === null) {
        return;
      }
      event.preventDefault();
      drag.current = {
        region,
        mode,
        startX: event.clientX,
        startY: event.clientY,
        rect: region.getBoundingClientRect(),
      };
      window.addEventListener("mousemove", onMove);
      window.addEventListener("mouseup", onUp);
    };

    const onKey = (event: KeyboardEvent): void => {
      const target = event.target;
      if (!(target instanceof HTMLElement)) {
        return;
      }
      const region = target.closest<HTMLElement>("[data-window]");
      if (region === null || !event.key.startsWith("Arrow")) {
        return;
      }
      const rect = region.getBoundingClientRect();
      if (rect.width === 0 || rect.height === 0) {
        return;
      }
      event.preventDefault();
      const step = 16;
      const dx = event.key === "ArrowLeft" ? -step : event.key === "ArrowRight" ? step : 0;
      const dy = event.key === "ArrowUp" ? -step : event.key === "ArrowDown" ? step : 0;
      if (event.shiftKey) {
        const bounds = desktop(element);
        const width = clamp(rect.width + dx, MIN_W, bounds.width - rect.left);
        const height = clamp(rect.height + dy, MIN_H, bounds.height - rect.top);
        place(region, rect.left, rect.top, width, height);
      } else {
        const bounds = desktop(element);
        place(
          region,
          clamp(rect.left + dx, 0, bounds.width - rect.width),
          clamp(rect.top + dy, bounds.top, bounds.height - rect.height),
          rect.width,
          rect.height,
        );
      }
    };

    const onResize = (): void => {
      for (const region of element.querySelectorAll<HTMLElement>("[data-window]")) {
        contain(region, element);
      }
    };

    element.addEventListener("mousedown", onDown);
    element.addEventListener("keydown", onKey);
    window.addEventListener("resize", onResize);
    return () => {
      element.removeEventListener("mousedown", onDown);
      element.removeEventListener("keydown", onKey);
      window.removeEventListener("resize", onResize);
      window.removeEventListener("mousemove", onMove);
      window.removeEventListener("mouseup", onUp);
    };
  }, [active, root]);

  return { reset };
}
