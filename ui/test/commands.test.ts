/**
 * The command model, held to the schema the library publishes.
 *
 * **This is the "14 of 14" claim.** The editor's `COMMAND_NAMES` and the tool list in the
 * catalogue fixture must be the same fourteen, and that fixture is `azoth forms --tools` output —
 * so the editor's reachable set, the wire's command enum and the agent's tool list meet in the
 * middle. A fifteenth command in Rust would fail here until the editor could send it.
 */

import { describe, expect, it } from "vitest";

import catalogueJson from "./fixtures/catalogue.json";
import { COMMAND_NAMES, type EditorCommand } from "../src/wire/commands";
import { removeCommandFor } from "../src/wire/nodes";
import type { Catalogue } from "../src/wire/types";

const catalogue = catalogueJson as unknown as Catalogue;

describe("the command model", () => {
  it("is the same fourteen the agent's tool schema publishes", () => {
    const tools = (catalogue.tools ?? []).map((tool) => tool.name);
    expect(tools).toHaveLength(14);
    expect([...COMMAND_NAMES].sort()).toEqual([...tools].sort());
  });

  it("is in the order the Rust enum declares, so a diff is readable", () => {
    expect(COMMAND_NAMES[0]).toBe("add_instance");
    expect(COMMAND_NAMES.at(-1)).toBe("set_position");
    expect(new Set(COMMAND_NAMES).size).toBe(COMMAND_NAMES.length);
  });
});

describe("removeCommandFor", () => {
  it("removes each node kind by its own command", () => {
    expect(removeCommandFor("instance:sep1")).toEqual({ command: "remove_instance", id: "sep1" });
    expect(removeCommandFor("input:feed_1")).toEqual({ command: "remove_input", name: "feed_1" });
    expect(removeCommandFor("product:out")).toEqual({ command: "remove_product", name: "out" });
  });

  it("refuses an id whose role is none of the three rather than guessing one", () => {
    // A node the projection did not make is not a node to invent a command for.
    expect(removeCommandFor("widget:sep1")).toBeNull();
    expect(removeCommandFor("sep1")).toBeNull();
    expect(removeCommandFor("")).toBeNull();
  });

  it("keeps a name that contains a colon, which only the first one separates", () => {
    const command: EditorCommand | null = removeCommandFor("instance:a:b");
    expect(command).toEqual({ command: "remove_instance", id: "a:b" });
  });
});
