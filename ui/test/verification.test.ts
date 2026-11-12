/**
 * The verification vocabulary, as this editor explains it.
 *
 * **A record field is not a verdict on the answer.** The four statuses the library derives say
 * how far the model itself is checked; the tests here hold the editor to printing the library's
 * spelling and one sentence, and to refusing to turn a status it does not know into one it does.
 */

import { describe, expect, it } from "vitest";

import { definitionOf } from "../src/state/verification";

describe("the verification vocabulary", () => {
  it("defines each status the library derives", () => {
    for (const status of [
      "verified",
      "partially_verified",
      "source_needed",
      "unverified",
    ]) {
      const info = definitionOf(status);
      expect(info.raw).toBe(status);
      expect(info.known).toBe(true);
      expect(info.definition.length).toBeGreaterThan(20);
    }
    // The normal case is named as such rather than left to sound like a defect.
    expect(definitionOf("partially_verified").definition).toContain(
      "normal state for the port",
    );
  });

  it("says an entry with no model has no check", () => {
    const info = definitionOf(null);
    expect(info.raw).toBeNull();
    expect(info.known).toBe(true);
    expect(info.definition).toContain("no registered model");
  });

  it("keeps an unknown spelling visible rather than calling it unverified", () => {
    const info = definitionOf("checked_by_a_future_tool");
    expect(info.raw).toBe("checked_by_a_future_tool");
    expect(info.known).toBe(false);
    expect(info.definition).toContain("does not recognise");
    expect(info.definition).not.toBe(definitionOf("unverified").definition);
  });
});
