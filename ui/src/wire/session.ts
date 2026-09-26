/**
 * The seam between the editor and the middleware, as a type.
 *
 * **One interface, and the door is not part of it.** The calls are string-in and string-out: a
 * command goes out, the envelope the library built comes back, and the front-end draws that and
 * re-derives nothing. A panel therefore never learns which door it is behind, which is what makes a
 * second door another module rather than a change to a component.
 *
 * **Every call is a promise, and that is a fact rather than a style.** The wasm module answers
 * synchronously; a door that is a `fetch` cannot; and one interface both implement is the only
 * version where a panel does not have to know.
 *
 * **Three calls, and the envelope carries the rest.** Every answer *is* the live envelope, so a
 * getter for the document, for the graph or for a value by path would be a second place the same
 * fact is held — and the editor holds no copy of a flowsheet. `flowsheet.document` is in every
 * envelope, which is why a save needs no call of its own.
 */

import type { Catalogue, Command, Envelope, ExecutionOrder } from "./types";

/** A live document, as the editor holds it. */
export interface Session {
  /** Apply one command, and answer with the envelope it left. */
  apply(command: Command): Promise<Envelope>;
  /** Run the document, and answer with the envelope it left. */
  run(): Promise<Envelope>;
  /**
   * Set which of the class's two orders the next run takes, and answer with the envelope.
   *
   * In the *envelope* rather than behind a getter: a control should read the order from the same
   * object as the diagnostics and the values, or the two can disagree after a call that was
   * refused.
   */
  setOrder(order: ExecutionOrder): Promise<Envelope>;
}

/** A session, and the envelope it was opened with. */
export interface Opened {
  session: Session;
  envelope: Envelope;
}

/** What every door offers, whichever kind it is. */
interface AnyDoor {
  /** The palette as one form per entry, and the agent's tools on request. */
  catalogue(withTools?: boolean): Promise<Catalogue>;
}

/** A door that opens documents: handed the TOML, it builds a session of it. The module is one. */
export type OpenableDoor = AnyDoor & {
  readonly kind: "openable";
  open(document: string): Promise<Opened>;
};

/** A door that hosts a document: the one it was started with, for the life of its process. */
export type HostedDoor = AnyDoor & {
  readonly kind: "hosted";
  attach(): Promise<Opened>;
};

/**
 * A door to the middleware: a palette, and a session.
 *
 * **Two kinds, because a door either opens documents or hosts one.** The wasm module is handed the
 * TOML and builds a session from it. `azoth serve` was started with `--flowsheet` and holds that
 * document for the life of its process: editing it is what a client can do to it, and no call hands
 * it another. That is a difference in the *type* rather than a flag to remember, so an editor
 * drawing a control has to say which door it is drawing it for — the Open button, for one, exists
 * only on the openable arm.
 */
export type Door = OpenableDoor | HostedDoor;
