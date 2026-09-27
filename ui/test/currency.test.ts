/**
 * The fixtures are the library's output **today**, and this is what keeps them so.
 *
 * `fixtures.test.ts` reads three committed documents the CLI once wrote, and everything the
 * front-end believes about the wire is held by them. What that did not hold is the fixtures
 * themselves: they were recaptured by hand, by whoever remembered, with the commands in that
 * file's header — so a field that moved in Rust would fail no build at all. The mirrors would
 * agree with a copy of the library that no longer exists, which is the failure `gen_stub.py`
 * prevents on the Python side and nothing prevented here.
 *
 * **It was already stale, and this case is what said so.** The catalogue fixture named `mol/s`
 * for a column's product-flow target where the shipped spec had been corrected to `mol/hr`, and
 * no test in the tree could see it: the mirror was being held to a document the library had
 * stopped emitting.
 *
 * **So the fixtures are recomputed rather than trusted.** The module built from this tree is asked
 * for both documents: `palette_json` is the same `middleware::catalogue` the CLI's `forms --tools`
 * calls, and `Editor`'s apply-then-run-then-envelope is the same three calls `azoth edit
 * --command … --run --json` makes.
 *
 * **Numbers are compared by tolerance and every other value exactly**, which is the rule the two
 * implementations are already held to: the fixture was captured from a native build and this
 * compares against a wasm one, and the two differ in the last bits of a converged flash
 * (`0.07701182836646382` against `0.07701182836646368`). A byte comparison would fail on
 * arithmetic neither half controls, and a tolerance loose enough to hide a changed string would
 * fail to notice the `mol/hr` above — so the tolerance is on the numbers and nothing else.
 *
 * This runs against the module `npm run wasm` just built in the `ui` job, so it costs no build of
 * its own. The wasm bytes are stubbed over `fetch` for the reason `app.test.tsx` stubs them: the
 * glue asks for the module beside itself, and this is the URL a static host would answer.
 */

import { readFileSync } from "node:fs";
import { resolve } from "node:path";
import { beforeAll, describe, expect, it } from "vitest";

import demo from "../../specs/flowsheets/demo.toml?raw";
import init, { Editor, palette_json } from "../src/wasm/pkg/azoth_wasm.js";

beforeAll(async () => {
  const wasm = readFileSync(resolve("src/wasm/pkg/azoth_wasm_bg.wasm"));
  globalThis.fetch = (): Promise<Response> =>
    Promise.resolve(new Response(wasm, { headers: { "Content-Type": "application/wasm" } }));
  await init();
});

/** One fixture, as the document it is committed as. */
const fixture = (name: string): unknown =>
  JSON.parse(readFileSync(resolve(`test/fixtures/${name}`), "utf8")) as unknown;

/**
 * Where the two documents differ, or `null`.
 *
 * Numbers are compared by a relative tolerance and every other value exactly — a string that
 * changed is the shape the stale `mol/s` had, and it is the one thing a tolerance must not hide.
 * Keys are compared in both directions, so a field the library added is a difference: a mirror
 * that silently dropped it is what this file is for.
 */
function difference(now: unknown, held: unknown, path: string): string | null {
  if (typeof now === "number" && typeof held === "number") {
    const scale = Math.max(1, Math.abs(now), Math.abs(held));
    return Math.abs(now - held) <= 1e-12 * scale
      ? null
      : `${path}: the library says ${now} and the fixture holds ${held}`;
  }
  if (Array.isArray(now) || Array.isArray(held)) {
    if (!Array.isArray(now) || !Array.isArray(held)) {
      return `${path}: an array on one side only`;
    }
    if (now.length !== held.length) {
      return `${path}: ${now.length} entries where the fixture holds ${held.length}`;
    }
    for (let index = 0; index < now.length; index += 1) {
      const found = difference(now[index], held[index], `${path}[${index}]`);
      if (found !== null) {
        return found;
      }
    }
    return null;
  }
  if (typeof now === "object" && now !== null && typeof held === "object" && held !== null) {
    const here = now as Record<string, unknown>;
    const there = held as Record<string, unknown>;
    for (const key of [...new Set([...Object.keys(here), ...Object.keys(there)])].sort()) {
      if (!(key in here)) {
        return `${path}.${key}: the fixture holds it and the library no longer emits it`;
      }
      if (!(key in there)) {
        return `${path}.${key}: the library emits it and the fixture does not hold it`;
      }
      const found = difference(here[key], there[key], `${path}.${key}`);
      if (found !== null) {
        return found;
      }
    }
    return null;
  }
  return Object.is(now, held)
    ? null
    : `${path}: the library says ${JSON.stringify(now)} and the fixture holds ${JSON.stringify(held)}`;
}

/** The library's document and the fixture, held to each other. */
function agree(now: string, fixtureName: string): void {
  expect(difference(JSON.parse(now), fixture(fixtureName), "the document")).toBeNull();
}

describe("the fixtures", () => {
  it("are the catalogue this tree's module carries", () => {
    // `azoth forms --tools`, which is what the fixtures file's header says captured it.
    agree(palette_json(true), "catalogue.json");
  });

  it("are the envelope this tree's editor answers with", () => {
    // `azoth edit --flowsheet specs/flowsheets/demo.toml --command '{…set_position…}' --run --json`,
    // in the three calls the CLI makes: open, apply, run.
    const editor = new Editor(demo);
    editor.apply('{"command":"set_position","node":"instance:sep1","x":1234,"y":56}');
    editor.run();
    agree(editor.envelope(), "envelope.json");
  });

  it("and the refused document is the same edit without the run", () => {
    // The third fixture is the demo with the heater removed and no run: the checker refuses two
    // ports the removal left underfed, which is the state the editor draws rather than raises.
    const editor = new Editor(demo);
    editor.apply('{"command":"remove_instance","id":"hx1"}');
    agree(editor.envelope(), "broken.json");
  });
});
