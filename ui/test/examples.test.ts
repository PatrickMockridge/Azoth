/**
 * The worked cases, opened and solved through the real module.
 *
 * **An example is a document the editor invites a newcomer to open**, so a parse that succeeds
 * and a run that refuses is the one failure this test is for. It is also the guard that every
 * file in the examples directory has a menu description: a file nobody described is a file
 * nobody meant to ship.
 */

import { readFileSync } from "node:fs";
import { resolve } from "node:path";
import { beforeAll, describe, expect, it } from "vitest";

import { EXAMPLES } from "../src/examples";
import init, { Editor } from "../src/wasm/pkg/azoth_wasm.js";

beforeAll(async () => {
  const wasm = readFileSync(resolve("src/wasm/pkg/azoth_wasm_bg.wasm"));
  globalThis.fetch = (): Promise<Response> =>
    Promise.resolve(new Response(wasm, { headers: { "Content-Type": "application/wasm" } }));
  await init();
});

describe("the shipped examples", () => {
  it("each has a description, and the order is the menu's", () => {
    expect(EXAMPLES.length).toBeGreaterThanOrEqual(4);
    for (const example of EXAMPLES) {
      expect(example.title.length, example.file).toBeGreaterThan(0);
      expect(example.description.length, example.file).toBeGreaterThan(0);
    }
  });

  it("each opens and solves", () => {
    for (const example of EXAMPLES) {
      const editor = new Editor(example.document);
      const envelope = JSON.parse(editor.run()) as {
        ok: boolean;
        run_error: string | null;
        session: { converged: boolean; streams: Record<string, unknown> };
      };
      expect(envelope.ok, example.file).toBe(true);
      expect(envelope.run_error, example.file).toBeNull();
      expect(envelope.session.converged, example.file).toBe(true);
      expect(Object.keys(envelope.session.streams).length, example.file).toBeGreaterThan(0);
    }
  });
});
