/**
 * The middleware, as the editor reaches it.
 *
 * **One interface, and the wasm module is its only implementation today.** The functions are
 * string-in and string-out, which is what makes a second implementation - a notebook kernel, a
 * hosted session - a matter of another module here rather than a change anywhere else.
 *
 * Nothing in this file decides anything. `open` hands the document to the library, `apply` hands
 * it one command, and both answer with the envelope the library built: the front-end draws what
 * came back and re-derives nothing, which is the whole reason the wire layer exists.
 */

import init, { Editor, palette_json } from "../wasm/pkg/azoth_wasm.js";
import type { Catalogue, Command, Envelope, Graph } from "./types";

/** A live document, as the editor holds it. */
export interface Session {
  /** Apply one command, and answer with the envelope it left. */
  apply(command: Command): Envelope;
  /** Run the document, and answer with the envelope it left. */
  run(): Envelope;
  /** Everything as it stands, without running anything. */
  envelope(): Envelope;
  /** The document as TOML, which is what a save writes. */
  document(): string;
  /** The connection graph alone. */
  graph(): Graph;
  /** One value by path, e.g. `p1.outlet.P`. */
  value(path: string): number;
  readonly ok: boolean;
  readonly dirty: boolean;
}

export interface Middleware {
  /** The palette as one form per entry, and the agent's tools on request. */
  catalogue(withTools?: boolean): Catalogue;
  /** Open a document from its TOML. */
  open(document: string): Session;
}

const envelope = (json: string): Envelope => JSON.parse(json) as Envelope;

/** The wasm module's implementation, once the module has loaded. */
class WasmSession implements Session {
  readonly #editor: Editor;

  constructor(editor: Editor) {
    this.#editor = editor;
  }

  apply(command: Command): Envelope {
    return envelope(this.#editor.apply(JSON.stringify(command)));
  }

  run(): Envelope {
    return envelope(this.#editor.run());
  }

  envelope(): Envelope {
    return envelope(this.#editor.envelope());
  }

  document(): string {
    return this.#editor.document();
  }

  graph(): Graph {
    return JSON.parse(this.#editor.graph()) as Graph;
  }

  value(path: string): number {
    return this.#editor.value(path);
  }

  get ok(): boolean {
    return this.#editor.ok;
  }

  get dirty(): boolean {
    return this.#editor.dirty;
  }
}

class WasmMiddleware implements Middleware {
  catalogue(withTools = false): Catalogue {
    return JSON.parse(palette_json(withTools)) as Catalogue;
  }

  open(document: string): Session {
    return new WasmSession(new Editor(document));
  }
}

/**
 * Load the module, and the middleware with it.
 *
 * The wasm is fetched beside the glue `wasm-pack` wrote (`new URL(..., import.meta.url)`), and
 * `vite` emits it as an asset — which is why the module is built with `--target web` rather than
 * `bundler`: no plugin, no second loader, and the same code path a static host serves.
 */
export async function wasmMiddleware(): Promise<Middleware> {
  await init();
  return new WasmMiddleware();
}
