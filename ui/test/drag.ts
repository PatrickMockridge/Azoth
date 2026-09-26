/**
 * A node drag, as the canvas receives one.
 *
 * **d3-drag listens for *mouse* events, and it attaches its move handler to `event.view`.** From
 * `node_modules/d3-drag/src/drag.js`: `mousedown.drag` is bound on the node, and the handler then
 * does `select(event.view).on("mousemove.drag", …).on("mouseup.drag", …)`. `@xyflow/system` uses
 * that instance directly (`dist/esm/index.js`, `drag()` and `d3Selection.call(d3DragInstance)`).
 * Both halves are what jsdom under vitest makes awkward, and each was measured:
 *
 * * a `MouseEvent` built with no `view` gives d3 `select(null)`, which throws inside the handler;
 * * `view: globalThis` is refused by jsdom's brand check — *"member view is not of type Window"* —
 *   because vitest redefines `document.defaultView` to return the global proxy rather than a node.
 *
 * **The Window is still reachable through the descriptor vitest shadowed.** `vitest` installs
 * `document.defaultView` with `Object.defineProperty` on the *document instance*, so it is an own
 * property that shadows the prototype accessor rather than replacing it. Calling the prototype's
 * own getter against the document therefore returns the live jsdom `Window` — the object events
 * belong to, and the one d3 can attach a listener to.
 *
 * **A pointer sequence cannot substitute, whatever library produces it.** `d3-drag` never listens
 * for `pointerdown`, and `@testing-library/user-event` builds every event it dispatches with a
 * `view` own-property of `null`, which is the same `select(null)` by a longer road.
 */

import { act } from "@testing-library/react";

/** jsdom's own `Window`, past the property the test environment shadows. */
function realWindow(): Window {
  const descriptor = Object.getOwnPropertyDescriptor(Document.prototype, "defaultView");
  const window: unknown = descriptor?.get?.call(document);
  if (window === undefined || window === null) {
    throw new Error("Document.prototype's own `defaultView` no longer answers a Window");
  }
  return window as Window;
}

/**
 * One event of the sequence.
 *
 * The global `MouseEvent` *is* jsdom's under this environment, so building it with the global and
 * handing it the window is the same object the window would have built — and `view` is the one
 * member of the init that d3-drag reads.
 */
function mouse(type: string, view: Window, x: number, y: number): MouseEvent {
  return new MouseEvent(type, {
    bubbles: true,
    cancelable: true,
    composed: true,
    view,
    clientX: x,
    clientY: y,
    button: 0,
  });
}

/**
 * Drag `node` from one client point to another.
 *
 * **The move has to cross a pixel**, because that is xyflow's node-drag threshold: a gesture that
 * did not would end without a position ever being written, and the test would pass on a canvas that
 * never moved. The moves are dispatched **on the window** because that is where d3-drag put its
 * listener; the press is on the node, because that is where it put the other one.
 */
export function drag(node: HTMLElement, from: [number, number], to: [number, number]): void {
  const view = realWindow();
  const middle: [number, number] = [(from[0] + to[0]) / 2, (from[1] + to[1]) / 2];

  act(() => {
    node.dispatchEvent(mouse("mousedown", view, ...from));
  });
  // Two moves rather than one, because a drag is a gesture: d3's own `mousemoving` flag is set by
  // the first one, and the second is the one that moves the node.
  act(() => {
    view.dispatchEvent(mouse("mousemove", view, ...middle));
  });
  act(() => {
    view.dispatchEvent(mouse("mousemove", view, ...to));
  });
  act(() => {
    view.dispatchEvent(mouse("mouseup", view, ...to));
  });
}

/** Where the canvas is drawing a node, in the document's own coordinates. */
export function drawnPosition(container: HTMLElement, id: string): [number, number] {
  const node = container.querySelector<HTMLElement>(`[data-id="${id}"]`);
  const transform = node?.style.transform ?? "";
  // xyflow positions each node with `translate(xpx, ypx)` in flow coordinates; the viewport's zoom
  // is a separate transform on the pane, which is what makes this readable as a position rather
  // than as a pixel.
  const match = /translate\((-?[\d.]+)px,\s*(-?[\d.]+)px\)/.exec(transform);
  if (match === null) {
    throw new Error(`${id} is drawn at no position: ${transform || "(none)"}`);
  }
  return [Number(match[1]), Number(match[2])];
}
