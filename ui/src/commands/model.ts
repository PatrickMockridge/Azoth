/**
 * One command, in the menu, the palette and the toolbar.
 *
 * **The list is built from the handlers, not from the buttons.** A menu and a palette that
 * re-implemented `Solve` would be a second path to the same call with its own enabled state; this
 * model carries the existing callback and the existing `disabled` expression, so the three
 * surfaces cannot disagree.
 */

export type MenuId = "case" | "book" | "edit" | "run" | "view" | "help";

export interface CommandSpec {
  id: string;
  label: string;
  menu: MenuId;
  run: () => void;
  enabled: boolean;
  shortcut?: string;
  hint?: string;
  /** Words that should find this command in the palette even when the label would not. */
  keywords?: string;
  /** A submenu: the parent is a heading, and these are its entries. */
  children?: readonly CommandSpec[];
}

export interface MenuSpec {
  id: MenuId;
  label: string;
  commands: readonly CommandSpec[];
}

/** The menu a list of commands belongs to, in the order the menu bar shows them. */
export function menusOf(
  commands: readonly CommandSpec[],
  labels: Record<MenuId, string>,
  order: readonly MenuId[],
): MenuSpec[] {
  return order
    .map((id) => ({
      id,
      label: labels[id],
      commands: commands.filter((command) => command.menu === id),
    }))
    .filter((menu) => menu.commands.length > 0);
}
