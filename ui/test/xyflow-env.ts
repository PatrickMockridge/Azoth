/**
 * What jsdom is missing for xyflow to mount.
 *
 * **Shared by the suites that render the editor**, because a polyfill only one of them installed
 * would fail as a `ReferenceError` in the other. Each class names the library that needs it, so the
 * next reader does not have to rediscover it — and a third would fail loudly rather than silently,
 * which is the property worth keeping.
 */

/** xyflow constructs this for every node it mounts, unguarded. */
class ResizeObserverStub {
  observe(): void {}
  unobserve(): void {}
  disconnect(): void {}
}

/** xyflow reads `m22` off the viewport's computed `transform`; jsdom has no such class. */
class DOMMatrixReadOnlyStub {
  m22 = 1;
}

/** Install both, which a suite does once. */
export function stubXyflowEnvironment(): void {
  (globalThis as unknown as { ResizeObserver: unknown }).ResizeObserver = ResizeObserverStub;
  (globalThis as { DOMMatrixReadOnly?: unknown }).DOMMatrixReadOnly = DOMMatrixReadOnlyStub;
}

/** The nodes a rendered editor is drawing, by their ids. */
export function nodeIds(container: HTMLElement): string[] {
  return [...container.querySelectorAll<HTMLElement>(".react-flow__node")].map(
    (node) => node.dataset.id ?? "",
  );
}

/** The word the status pill is showing. */
export function pill(container: HTMLElement): string {
  return [...container.querySelectorAll(".pill")].at(-1)?.textContent ?? "";
}

/** The id of the document the editor has open, off its own pill. */
export function idPill(container: HTMLElement): string {
  return container.querySelector(".pill")?.textContent ?? "";
}
