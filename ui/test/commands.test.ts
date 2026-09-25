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
import { accelerationNames } from "../src/wire/field";
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

describe("the acceleration names", () => {
  it("are read from the schema rather than spelled in a widget", () => {
    // **One list, read.** The names live in `recycle::ACCELERATION_NAMES` and the schema carries
    // them; the panel offers what this returns, so a name the parser would refuse cannot be
    // offered - and the Rust test `the_acceleration_names_a_tool_offers_are_the_recycle_class_s_own`
    // is what holds the schema to the class.
    expect(accelerationNames(catalogue)).toEqual([
      "direct_substitution",
      "wegstein",
      "broyden",
    ]);
    // A catalogue without the tool is an empty list, and the panel falls back to a text field
    // rather than rendering an empty dropdown.
    expect(accelerationNames(null)).toEqual([]);
    expect(accelerationNames({ unit_ops: [] })).toEqual([]);
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
