/**
 * The middleware as a browser reaches it: the wasm module, in the page.
 *
 * **One door of two, and nothing here decides anything.** `open` hands the document to the library,
 * `apply` hands it one command, and both answer with the envelope the library built. The other door
 * is `./http.ts`, and the interface both implement is `./session.ts`.
 */

import init, { Editor, palette_json } from "../wasm/pkg/azoth_wasm.js";
import type { Door, Opened, OpenableDoor, Session } from "./session";
import type { Catalogue, Command, Envelope, ExecutionOrder } from "./types";

const envelope = (json: string): Envelope => JSON.parse(json) as Envelope;

/**
 * The module's session.
 *
 * **`async` and not `Promise.resolve` around the call, because a refusal is thrown synchronously.**
 * A command the library cannot read answers with an exception rather than an envelope, and a method
 * that wrapped the call would throw *before* a promise existed — which is the one shape a caller
 * writing `.catch` cannot see.
 */
class WasmSession implements Session {
  readonly #editor: Editor;

  constructor(editor: Editor) {
    this.#editor = editor;
  }

  async apply(command: Command): Promise<Envelope> {
    return envelope(this.#editor.apply(JSON.stringify(command)));
  }

  async run(): Promise<Envelope> {
    return envelope(this.#editor.run());
  }

  async setOrder(order: ExecutionOrder): Promise<Envelope> {
    return envelope(this.#editor.set_order(order));
  }
}

class WasmDoor implements OpenableDoor {
  readonly kind = "openable" as const;

  async catalogue(withTools = false): Promise<Catalogue> {
    return JSON.parse(palette_json(withTools)) as Catalogue;
  }

  async open(document: string): Promise<Opened> {
    // **Both, from one construction.** The envelope the editor opens on is the one this session
    // already stands at, so asking the module for it afterwards would spend a call on a fact that
    // was in hand — and a door that had to be polled to boot would be a different shape from this
    // one for no reason.
    const editor = new Editor(document);
    return { session: new WasmSession(editor), envelope: envelope(editor.envelope()) };
  }
}

/**
 * Load the module, and the door with it.
 *
 * The wasm is fetched beside the glue `wasm-pack` wrote (`new URL(..., import.meta.url)`), and
 * `vite` emits it as an asset — which is why the module is built with `--target web` rather than
 * `bundler`: no plugin, no second loader, and the same code path a static host serves.
 */
export async function moduleDoor(): Promise<Door> {
  await init();
  return new WasmDoor();
}
