/**
 * One edit to a calculation book.
 *
 * **A closed union, as the flowsheet's commands are.** The inputs are SI magnitudes for
 * quantities; an absent optional input is an `unset_input`, never a null. The tag is the
 * library's `Command` enum's own `snake_case` spelling.
 */

export type SheetCommand =
  | { command: "add_sheet"; id: string; calculation: string; name?: string }
  | { command: "remove_sheet"; id: string }
  | { command: "rename_sheet"; id: string; name: string }
  | { command: "set_input"; sheet: string; input: string; value: number | boolean | string }
  | { command: "unset_input"; sheet: string; input: string }
  | { command: "set_book_name"; name: string };

/** Every command's tag, in the order the library declares them. */
export const SHEET_COMMAND_NAMES: readonly SheetCommand["command"][] = [
  "add_sheet",
  "remove_sheet",
  "rename_sheet",
  "set_input",
  "unset_input",
  "set_book_name",
];
