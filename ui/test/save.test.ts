// @vitest-environment jsdom
/**
 * Writing the document out.
 *
 * **Two branches, and each case fails if the other one is taken.** jsdom has neither
 * `showSaveFilePicker` nor `URL.createObjectURL`, so both sides are installed by the test; a
 * download case that quietly went through a picker, or the reverse, would pass a weaker assertion
 * than these.
 *
 * The third case is the one the module exists for: a dismissed picker is somebody deciding not to
 * save, and it must not reach the editor as a fault.
 */

import { afterEach, describe, expect, it, vi } from "vitest";

import { saveDocument } from "../src/save";
import { only, stubDownload, textOf, unstubDownload } from "./download";

/** A browser with the picker, recording what a save asked for. */
function stubPicker(refusal?: Error): { written: string[]; suggested: string[] } {
  const written: string[] = [];
  const suggested: string[] = [];
  vi.stubGlobal(
    "showSaveFilePicker",
    (options?: { suggestedName?: string }) => {
      suggested.push(options?.suggestedName ?? "");
      if (refusal !== undefined) {
        return Promise.reject(refusal);
      }
      return Promise.resolve({
        createWritable: () =>
          Promise.resolve({
            write: (data: string) => {
              written.push(data);
              return Promise.resolve();
            },
            close: () => Promise.resolve(),
          }),
      });
    },
  );
  return { written, suggested };
}

afterEach(() => {
  vi.unstubAllGlobals();
  vi.restoreAllMocks();
  unstubDownload();
});

describe("saving the document", () => {
  it("writes through the picker where the browser has one, and downloads nothing", async () => {
    const picker = stubPicker();
    const download = stubDownload();

    await saveDocument("flowsheets.demo.toml", "id = \"flowsheets.demo\"");

    expect(picker.suggested).toEqual(["flowsheets.demo.toml"]);
    expect(picker.written).toEqual(["id = \"flowsheets.demo\""]);
    expect(download.names).toEqual([]);
    expect(download.blobs).toEqual([]);
  });

  it("downloads where the browser has not got one", async () => {
    const download = stubDownload();

    await saveDocument("flowsheets.demo.toml", "id = \"flowsheets.demo\"");

    expect(download.names).toEqual(["flowsheets.demo.toml"]);
    expect(download.blobs).toHaveLength(1);
    expect(await textOf(only(download.blobs))).toBe("id = \"flowsheets.demo\"");
  });

  it("says nothing when a person dismisses the picker", async () => {
    const refusal = new Error("the user aborted a request");
    refusal.name = "AbortError";
    const picker = stubPicker(refusal);
    const download = stubDownload();

    // Not a rejection: a dismissal is an answer, and the answer is that nothing was written.
    await expect(saveDocument("flowsheets.demo.toml", "id = \"x\"")).resolves.toBeUndefined();
    expect(picker.written).toEqual([]);
    expect(download.names).toEqual([]);
  });

  it("reports anything else the picker refused", async () => {
    stubPicker(new Error("no permission to write there"));

    await expect(saveDocument("flowsheets.demo.toml", "id = \"x\"")).rejects.toThrow(
      "no permission to write there",
    );
  });
});
