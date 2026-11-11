// @vitest-environment jsdom
/**
 * The calculation book over the served door.
 *
 * **`/sheets` is its own route.** The flowsheet's `/rpc` holds one document and the book is a
 * second; this test is that the browser asks the second route and reads the envelope it answers.
 */

import { cleanup, fireEvent, render, screen, waitFor } from "@testing-library/react";
import { afterEach, beforeAll, describe, expect, it, vi } from "vitest";

import { stubXyflowEnvironment } from "./xyflow-env";
import catalogueJson from "./fixtures/catalogue.json";
import envelopeJson from "./fixtures/envelope.json";

beforeAll(stubXyflowEnvironment);

const sheetEnvelope = {
  ok: true,
  dirty: false,
  book: {
    id: "calcbook_1",
    name: "Served book",
    document: 'id = "calcbook_1"\nname = "Served book"\nsheets = []\n',
    sheets: [],
  },
};

afterEach(() => {
  cleanup();
  vi.unstubAllGlobals();
});

describe("the calculation book on a served door", () => {
  it("asks /sheets for the book and draws its empty state", async () => {
    vi.stubGlobal(
      "fetch",
      (input: RequestInfo | URL, init?: RequestInit): Promise<Response> => {
        const url = String(input);
        const body = JSON.parse(String(init?.body)) as Record<string, unknown>;
        const answer = url.endsWith("/sheets")
          ? sheetEnvelope
          : "catalogue" in body
            ? catalogueJson
            : envelopeJson;
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

    fireEvent.click(screen.getByRole("button", { name: "Calculations" }));
    await waitFor(() => {
      expect(screen.getByText("Calculation sheets")).toBeDefined();
    });
    expect(screen.getByDisplayValue("Served book")).toBeDefined();
  });
});
