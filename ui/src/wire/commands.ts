/**
 * Every edit the middleware has, as a type.
 *
 * **This is the check that the editor can reach all fourteen, and it is the TypeScript
 * `deny_unknown_fields`.** The Rust command enum refuses a key it does not declare; a union here
 * refuses one at compile time, so a typo'd field name is a type error rather than a command the
 * wasm boundary rejects at run time with a sentence about a key nobody meant to write.
 *
 * `test/commands.test.ts` holds `COMMAND_NAMES` to the tool schema the library publishes, and that
 * schema is itself held to `Command`'s variants by
 * `every_command_has_a_tool_and_no_tool_is_a_second_surface` in the Rust crate — so the two lists
 * meet in the middle rather than being maintained side by side.
 */

/** One of the seven `[[recycles]]` keys a tear carries. */
export type RecycleField =
  | "flow_tolerance"
  | "composition_tolerance"
  | "temperature_tolerance"
  | "pressure_tolerance"
  | "max_iterations"
  | "minimum_flow"
  | "acceleration_method";

/** One edit, as the command model reads it — the fourteen, and nothing else. */
export type EditorCommand =
  | { command: "add_instance"; id: string; unit: string; parameters: Record<string, unknown> }
  | { command: "remove_instance"; id: string }
  | { command: "connect"; from: string; to: string }
  | { command: "disconnect"; from: string; to: string }
  | { command: "set_parameter"; instance: string; name: string; value: unknown }
  | { command: "unset_parameter"; instance: string; name: string }
  | {
      command: "add_input";
      name: string;
      components: string[];
      n: number;
      z: number[];
      P: number;
      T: number;
    }
  | { command: "remove_input"; name: string }
  | { command: "add_product"; name: string }
  | { command: "remove_product"; name: string }
  | { command: "add_recycle"; stream: string; from: string; to: string }
  | { command: "remove_recycle"; stream: string }
  | { command: "set_recycle"; stream: string; field: RecycleField; value: number | string }
  | { command: "set_position"; node: string; x: number; y: number };

/** The names, which is what a test can compare against a document. */
export type CommandName = EditorCommand["command"];

/** Every command this editor can send, in the order the Rust enum declares them. */
export const COMMAND_NAMES: readonly CommandName[] = [
  "add_instance",
  "remove_instance",
  "connect",
  "disconnect",
  "set_parameter",
  "unset_parameter",
  "add_input",
  "remove_input",
  "add_product",
  "remove_product",
  "add_recycle",
  "remove_recycle",
  "set_recycle",
  "set_position",
];
