// @vitest-environment jsdom
/**
 * The first-run card: it offers the three ways in, and it remembers being dealt with.
 *
 * The card is presentation; the document it opens is the editor's own `open`, which
 * `app.test.tsx` exercises. What is tested here is that a newcomer's choices are wired to the
 * callbacks and that the dismissal survives a reload.
 */

import { cleanup, fireEvent, render, screen } from "@testing-library/react";
import { afterEach, beforeEach, describe, expect, it, vi } from "vitest";

import { Welcome } from "../src/components/Welcome";
import { readWelcome, rememberWelcome } from "../src/state/welcome";

beforeEach(() => {
  window.localStorage.clear();
});

afterEach(() => {
  cleanup();
});

const examples = [
  {
    file: "heat_a_gas.toml",
    title: "Heat a gas",
    description: "One feed, one heater and one product.",
    document: "id = \"flowsheets.example\"\nname = \"Heat a gas\"",
  },
];

describe("the welcome card", () => {
  it("offers the demo, a blank page, the guide and the examples", () => {
    const onDemo = vi.fn();
    const onBlank = vi.fn();
    const onHelp = vi.fn();
    const onExample = vi.fn();
    render(
      <Welcome
        examples={examples}
        onDemo={onDemo}
        onBlank={onBlank}
        onExample={onExample}
        onHelp={onHelp}
        onDismiss={() => {}}
      />,
    );

    fireEvent.click(screen.getByRole("button", { name: "Open the demo" }));
    fireEvent.click(screen.getByRole("button", { name: "Start blank" }));
    fireEvent.click(screen.getByRole("button", { name: "Read the 2-minute guide" }));
    fireEvent.click(screen.getByRole("button", { name: "Heat a gas" }));

    expect(onDemo).toHaveBeenCalledOnce();
    expect(onBlank).toHaveBeenCalledOnce();
    expect(onHelp).toHaveBeenCalledOnce();
    expect(onExample).toHaveBeenCalledWith(examples[0]?.document);
  });

  it("remembers that it was dealt with", () => {
    expect(readWelcome()).toBe(true);
    rememberWelcome();
    expect(readWelcome()).toBe(false);
  });
});
