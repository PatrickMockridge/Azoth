// @vitest-environment jsdom
/**
 * The hosted door: `azoth serve`'s `POST /rpc`, as the editor reaches it.
 *
 * **The bodies this posts are the assertion**, and they are the whole of what the door decides. The
 * route reads `{"command": …, "run": …, "order": …}` and every key of it is optional, so which key
 * a call sends *is* which call it is — and one of them is not obvious: the served session runs the
 * document after every call, so an edit that sent no `run` would pay the physics on every
 * keystroke. A stub that answered anything would hide that, so these assert the request rather than
 * the answer.
 *
 * **The answers are the real ones.** `test/fixtures/` holds what `azoth forms --tools` and
 * `azoth edit --json` wrote, and `/rpc` writes the same object from the same codec — so the App
 * case below renders the document the library actually emitted rather than an emulation of it.
 */

import { cleanup, render, screen, waitFor } from "@testing-library/react";
import { afterEach, beforeAll, describe, expect, it, vi } from "vitest";

import catalogueJson from "./fixtures/catalogue.json";
import envelopeJson from "./fixtures/envelope.json";
import { HttpDoor } from "../src/wire/http";
import { serveUrl } from "../src/wire/door";
import { nodeIds, pill, stubXyflowEnvironment } from "./xyflow-env";

beforeAll(stubXyflowEnvironment);

/**
 * A server that answers with one thing, and records what it was asked.
 *
 * Nothing is emulated but the socket: the object it answers with is the fixture, which is the
 * library's own output.
 */
function stubRoute(answer: unknown, status = 200): { bodies: Record<string, unknown>[] } {
  const bodies: Record<string, unknown>[] = [];
  vi.stubGlobal(
    "fetch",
    (_input: RequestInfo | URL, init?: RequestInit): Promise<Response> => {
      bodies.push(JSON.parse(String(init?.body)) as Record<string, unknown>);
      return Promise.resolve(
        new Response(JSON.stringify(answer), {
          status,
          headers: { "Content-Type": "application/json" },
        }),
      );
    },
  );
  return { bodies };
}

afterEach(() => {
  // The suite imports `react` through `App`, so the tree has to go the way `app.test.tsx` takes it.
  cleanup();
  vi.unstubAllGlobals();
  // And the door is chosen from the URL, so a test that set one has to put it back.
  history.replaceState({}, "", "/");
});

describe("the door the URL names", () => {
  it("is the module unless the query string names a server", () => {
    expect(serveUrl("")).toBeNull();
    expect(serveUrl("?other=1")).toBeNull();
    // An empty value names no server rather than a server at "".
    expect(serveUrl("?serve=")).toBeNull();
    expect(serveUrl("?serve=%20")).toBeNull();
    expect(serveUrl("?serve=http://127.0.0.1:4000")).toBe("http://127.0.0.1:4000");
    // A trailing slash is trimmed, because every call is `${base}/rpc`.
    expect(serveUrl("?serve=http://127.0.0.1:4000/")).toBe("http://127.0.0.1:4000");
  });
});

describe("what a call posts", () => {
  it("sends the body the route reads, for each of the four calls", async () => {
    const door = new HttpDoor("http://127.0.0.1:4000");
    const { bodies } = stubRoute(envelopeJson);

    const { session } = await door.attach();
    await door.catalogue(true);
    await session.apply({ command: "remove_instance", id: "hx1" });
    await session.run();
    await session.setOrder("topological");

    expect(bodies).toEqual([
      // The look that opens on the document the process already holds.
      { command: null },
      { catalogue: true },
      // **`run: false`, and it is the point of the case.** The served session runs after every
      // call, so an edit that said nothing would run the flowsheet on every keystroke.
      { command: { command: "remove_instance", id: "hx1" }, run: false },
      { command: null, run: true },
      { order: "topological", command: null },
    ]);
  });

  it("tells a refused call apart from a document the checker refuses", async () => {
    const door = new HttpDoor("http://127.0.0.1:4000");
    let answer: { body: unknown; status: number } = { body: envelopeJson, status: 200 };
    vi.stubGlobal(
      "fetch",
      (): Promise<Response> =>
        Promise.resolve(
          new Response(JSON.stringify(answer.body), {
            status: answer.status,
            headers: { "Content-Type": "application/json" },
          }),
        ),
    );

    const { session } = await door.attach();

    // A `400` is about the *call*: the server's own sentence reaches the editor as it was made.
    answer = { body: { error: "Unknown tool: delete_everything" }, status: 400 };
    await expect(session.apply({ command: "remove_instance", id: "hx1" })).rejects.toThrow(
      "Unknown tool: delete_everything",
    );

    // A call that *landed* and left a document the checker refuses is an envelope, not an error —
    // which is the "refused" pill and the diagnostics panel, not the fault line.
    answer = { body: { ...envelopeJson, ok: false }, status: 200 };
    const refused = await session.apply({ command: "remove_instance", id: "hx1" });
    expect(refused.ok).toBe(false);
  });
});

/**
 * **The whole editor, behind the served door.** Nothing here is the module: the door is chosen from
 * the URL, every answer arrives over `fetch`, and the assertions are about the DOM — so a door that
 * was selected but not *used*, or one whose envelope was not the editor's state, fails.
 */
describe("the editor on a served door", () => {
  it("draws the document the server holds, and offers no control that would replace it", async () => {
    const asked: string[] = [];
    vi.stubGlobal(
      "fetch",
      (input: RequestInfo | URL, init?: RequestInit): Promise<Response> => {
        asked.push(String(input));
        const body = JSON.parse(String(init?.body)) as Record<string, unknown>;
        const answer = "catalogue" in body ? catalogueJson : envelopeJson;
        return Promise.resolve(
          new Response(JSON.stringify(answer), {
            status: 200,
            headers: { "Content-Type": "application/json" },
          }),
        );
      },
    );

    history.replaceState({}, "", "/?serve=http://127.0.0.1:4000");
    const { App } = await import("../src/App");
    const rendered = render(<App />);

    await waitFor(() => {
      expect(rendered.container.querySelectorAll(".react-flow__node")).not.toHaveLength(0);
    });
    // The six nodes of the document the server holds, by the ids the projection gives them.
    expect(nodeIds(rendered.container)).toHaveLength(6);
    // And the pill is the server's own answer, which is the envelope having come over the socket
    // rather than out of the module.
    expect(pill(rendered.container)).toBe("solved");
    // Nothing was fetched but the route: the wasm module was never loaded, which is what hosted
    // means. A stub that had to serve wasm bytes would fail this.
    expect(asked.length).toBeGreaterThan(0);
    expect(asked.every((url) => url.endsWith("/rpc"))).toBe(true);

    // The palette came over the same socket, which is what a unit-op window is drawn from.
    expect(screen.getByRole("button", { name: /pump/i })).toBeDefined();

    // **And the controls that hand a document to the library are refused**, because this door
    // holds one document for the life of its process and the route has no call to replace it.
    const disabled = (name: string): boolean =>
      (screen.getByRole("button", { name }) as HTMLButtonElement).disabled;
    expect(disabled("New")).toBe(true);
    expect(disabled("Demo")).toBe(true);
    // Save is not one of them: every envelope carries the document, so it works here too.
    expect(disabled("Save")).toBe(false);
    // Open is a label rather than a button, because a file input cannot be styled — so what says
    // it is refused is that there is no input inside it to open.
    const open = screen.getByText("Open").closest("label");
    expect(open?.classList.contains("off")).toBe(true);
    expect(open?.querySelector("input")).toBeNull();
  });
});
