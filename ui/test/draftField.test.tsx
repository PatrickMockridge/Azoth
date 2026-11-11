// @vitest-environment jsdom
/** A draft field sends one value when the typing stops, and never a value it cannot read. */

import { cleanup, fireEvent, render, screen } from "@testing-library/react";
import { afterEach, describe, expect, it, vi } from "vitest";

import { DraftField } from "../src/components/DraftField";

afterEach(cleanup);

describe("a draft field", () => {
  it("sends one number on blur, including scientific notation", () => {
    const onCommit = vi.fn();
    render(<DraftField value="0.01" kind="number" onCommit={onCommit} />);
    const input = screen.getByDisplayValue("0.01");

    fireEvent.change(input, { target: { value: "1e-5" } });
    expect(onCommit).not.toHaveBeenCalled();
    fireEvent.blur(input);

    expect(onCommit).toHaveBeenCalledOnce();
    expect(onCommit).toHaveBeenCalledWith(1e-5);
  });

  it("does not send a value it cannot read", () => {
    const onCommit = vi.fn();
    render(<DraftField value="0.01" kind="number" onCommit={onCommit} />);
    const input = screen.getByDisplayValue("0.01");

    fireEvent.change(input, { target: { value: "not a number" } });
    fireEvent.blur(input);

    expect(onCommit).not.toHaveBeenCalled();
    expect((input as HTMLInputElement).value).toBe("0.01");
  });

  it("sends an empty string as an unset, and Escape abandons the draft", () => {
    const onCommit = vi.fn();
    render(<DraftField value="0.01" kind="number" onCommit={onCommit} />);
    const input = screen.getByDisplayValue("0.01");

    fireEvent.change(input, { target: { value: "0.2" } });
    fireEvent.keyDown(input, { key: "Escape" });
    fireEvent.blur(input);
    expect(onCommit).not.toHaveBeenCalled();

    fireEvent.change(input, { target: { value: "" } });
    fireEvent.blur(input);
    expect(onCommit).toHaveBeenCalledWith("");
  });
});
