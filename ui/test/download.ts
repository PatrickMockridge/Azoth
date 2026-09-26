/**
 * What a save wrote, and the stand-in for the parts of a browser that write it.
 *
 * **jsdom has neither `showSaveFilePicker` nor `URL.createObjectURL`**, so a save in a test either
 * goes through a picker the test installed or through a download the test installed — and that is
 * what makes the two branches tell each other apart: a case that installed one and got the other
 * has a recorded list where it expected none.
 *
 * This is shared because two suites need it — `save.test.ts` for the fallback itself, and
 * `app.test.tsx` for the one place the document's *text* is reachable from the DOM, which is what
 * says the canvas and the document carry the same position rather than two that agree by accident.
 */

import { vi } from "vitest";

/**
 * What a jsdom `Blob` carries.
 *
 * jsdom 25's `Blob` has no `text()`, and its blobs are not the ones Node's `Response` understands —
 * so the read goes through jsdom's own `FileReader`, the one thing in the environment that can open
 * one.
 */
export function textOf(blob: Blob): Promise<string> {
  return new Promise((resolve, reject) => {
    const reader = new FileReader();
    reader.onload = () => resolve(String(reader.result));
    reader.onerror = () => reject(reader.error ?? new Error("unreadable"));
    reader.readAsText(blob);
  });
}

/** The one blob a download wrote, or a failure that says there was none. */
export function only(blobs: Blob[]): Blob {
  const [blob] = blobs;
  if (blob === undefined) {
    throw new Error("the download wrote no blob");
  }
  return blob;
}

let installed: (() => void) | null = null;

/**
 * A browser without a picker, where a save is a Blob and an anchor.
 *
 * The blob is kept rather than read here, because reading it is the case's business.
 */
export function stubDownload(): { blobs: Blob[]; names: string[] } {
  const blobs: Blob[] = [];
  const names: string[] = [];
  const url = URL as unknown as {
    createObjectURL(blob: Blob): string;
    revokeObjectURL(url: string): void;
  };
  url.createObjectURL = (blob: Blob): string => {
    blobs.push(blob);
    return "blob:test";
  };
  url.revokeObjectURL = (): void => {};
  // The anchor is what names the file, and it is the only place that name is written down.
  const anchor = vi
    .spyOn(HTMLAnchorElement.prototype, "click")
    .mockImplementation(function (this: HTMLAnchorElement) {
      names.push(this.download);
    });

  installed = (): void => {
    // Installed on `URL` itself, so `vi.unstubAllGlobals` does not know about them.
    Reflect.deleteProperty(url, "createObjectURL");
    Reflect.deleteProperty(url, "revokeObjectURL");
    anchor.mockRestore();
    installed = null;
  };
  return { blobs, names };
}

/** Put the environment back, which a suite does once per case. */
export function unstubDownload(): void {
  installed?.();
}
