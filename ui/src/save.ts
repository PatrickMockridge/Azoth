/**
 * Writing the document out.
 *
 * **The browser's own picker where there is one, a download where there is not.** `showSaveFilePicker`
 * lets a person choose the file and overwrite it; the anchor is the fallback every browser has. What
 * the editor sees is the same either way, which is why this is a module rather than a branch in the
 * tool bar.
 *
 * **A cancelled picker is not a fault.** `AbortError` is somebody deciding not to save, which is a
 * thing they are allowed to do — and it is the one state this module exists to keep off the
 * editor's error line.
 */

/**
 * The one call this uses of the File System Access API.
 *
 * **Declared here because `lib.dom.d.ts` does not have it.** The API is Chromium's and the standard
 * is still a draft, so the DOM types TypeScript ships describe it nowhere; a local declaration is
 * how a call to it is typed without a dependency for one method.
 */
interface SaveFilePicker {
  (options?: { suggestedName?: string }): Promise<{
    createWritable(): Promise<{ write(data: string): Promise<void>; close(): Promise<void> }>;
  }>;
}

/** Whether an error is the picker being dismissed rather than anything having gone wrong. */
function cancelled(error: unknown): boolean {
  return (
    typeof error === "object" &&
    error !== null &&
    "name" in error &&
    (error as { name: unknown }).name === "AbortError"
  );
}

/** The fallback: a Blob and an anchor, which is the download every browser can do. */
function download(name: string, text: string): void {
  const url = URL.createObjectURL(new Blob([text], { type: "text/plain" }));
  const link = document.createElement("a");
  link.href = url;
  link.download = name;
  link.click();
  URL.revokeObjectURL(url);
}

/**
 * Write `text` as `name`.
 *
 * # Errors
 * Whatever the picker or the write refused, so that a real failure reaches the editor. A dismissal
 * is not one, and is answered with nothing.
 */
export async function saveDocument(name: string, text: string): Promise<void> {
  const picker = (window as { showSaveFilePicker?: SaveFilePicker }).showSaveFilePicker;
  if (picker === undefined) {
    download(name, text);
    return;
  }
  // **The dismissal is caught here rather than around the write.** They are different questions: a
  // person saying no is answered with silence, and a write that fails for any other reason is a
  // failure the caller has to hear about.
  const handle = await picker({ suggestedName: name }).catch((error: unknown) => {
    if (cancelled(error)) {
      return undefined;
    }
    throw error;
  });
  if (handle === undefined) {
    return;
  }
  const writable = await handle.createWritable();
  await writable.write(text);
  await writable.close();
}
