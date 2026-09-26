/**
 * The middleware as a page reaches it over a network: `azoth serve`'s `POST /rpc`.
 *
 * **The same envelope, and not a second surface.** A call posts a command and the route answers
 * what the library's own codec wrote, so a panel behind this door draws the object the module door
 * draws — the two differ in where the kernels run and in nothing a widget can see.
 *
 * **This door is `hosted` and not `openable`, and that is the server's shape rather than a
 * limitation of this file.** `azoth serve` is started with `--flowsheet` and holds that one
 * document for the life of its process; editing it is what a client can do to it, and the route has
 * no call that hands it another. `attach` takes the document the process already has.
 *
 * **One call is not the module's.** The served session runs the document after every call, so an
 * edit posted without a `run` would pay the physics on every keystroke — the opposite of the
 * editor's cadence, where a check is cheap and happens on each edit and a run is asked for. The
 * edit says `false`; Solve says `true`.
 */

import type { HostedDoor, Opened, Session } from "./session";
import type { Catalogue, Command, Envelope, ExecutionOrder } from "./types";

/**
 * One request body, and the object it answered with.
 *
 * **A refusal is a thrown sentence, and a refusal is about the *call*.** The route says `400` with
 * `{"error": …}` for a request it cannot read or a command that cannot take effect; a call that
 * *landed* and left a document the checker refuses is a `200` carrying `ok: false`, which is an
 * envelope and not an error — the state the editor draws its "refused" pill and diagnostics for.
 * Keeping those two apart is the whole of this function.
 *
 * `Content-Type: application/json` is what makes this a *preflighted* request, which is why the
 * server has to be told `--allow-origin` before a browser will let the answer be read.
 */
async function post(base: string, body: Record<string, unknown>): Promise<unknown> {
  const response = await fetch(`${base}/rpc`, {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify(body),
  });
  const answer: unknown = await response.json();
  if (!response.ok) {
    throw new Error(sentence(answer, response.status));
  }
  return answer;
}

/** The server's own sentence where it wrote one, so a refusal reaches the editor as it was made. */
function sentence(answer: unknown, status: number): string {
  if (typeof answer === "object" && answer !== null && "error" in answer) {
    const { error } = answer as { error: unknown };
    if (typeof error === "string") {
      return error;
    }
  }
  return `the server answered ${status}`;
}

class HttpSession implements Session {
  readonly #base: string;

  constructor(base: string) {
    this.#base = base;
  }

  async apply(command: Command): Promise<Envelope> {
    return (await post(this.#base, { command, run: false })) as Envelope;
  }

  async run(): Promise<Envelope> {
    // `command: null` is the route's "the envelope as it stands", and `run` is what makes it run.
    return (await post(this.#base, { command: null, run: true })) as Envelope;
  }

  async setOrder(order: ExecutionOrder): Promise<Envelope> {
    return (await post(this.#base, { order, command: null })) as Envelope;
  }
}

/** The `azoth serve` door. */
export class HttpDoor implements HostedDoor {
  readonly kind = "hosted" as const;

  readonly #base: string;

  constructor(base: string) {
    this.#base = base;
  }

  async catalogue(withTools = false): Promise<Catalogue> {
    return (await post(this.#base, { catalogue: withTools })) as Catalogue;
  }

  async attach(): Promise<Opened> {
    const envelope = (await post(this.#base, { command: null })) as Envelope;
    return { session: new HttpSession(this.#base), envelope };
  }
}
