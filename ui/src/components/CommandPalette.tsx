/**
 * The command palette: the Cyberpunk shell's keyboard-first way to the same commands.
 *
 * **The list is options, not buttons.** The rendered cases look for one visible `Solve`, `New`,
 * `Save` and so on; a palette that rendered them as buttons would multiply every query. The
 * overlay exists only while it is open, its rows are `role="option"`, and it calls the same
 * `CommandSpec.run` as the toolbar and the menus.
 *
 * The shortcut is Ctrl/Cmd+K and is wired by the workspace that owns the palette; it does nothing
 * outside Cyberpunk and does not open over a modal.
 */

import { useEffect, useMemo, useRef, useState, type KeyboardEvent } from "react";

import type { CommandSpec } from "../commands/model";

const MENU_LABEL: Record<string, string> = {
  case: "Case",
  book: "Book",
  edit: "Edit",
  run: "Run",
  view: "View",
  help: "Help",
};

export function CommandPalette({
  commands,
  onClose,
}: {
  commands: readonly CommandSpec[];
  onClose: () => void;
}) {
  const [query, setQuery] = useState("");
  const [active, setActive] = useState(0);
  const input = useRef<HTMLInputElement>(null);
  const restore = useRef<Element | null>(null);

  const filtered = useMemo(() => {
    const needle = query.trim().toLowerCase();
    if (needle === "") {
      return commands;
    }
    return commands.filter((command) =>
      `${command.label} ${command.keywords ?? ""} ${command.menu}`
        .toLowerCase()
        .includes(needle),
    );
  }, [commands, query]);

  useEffect(() => {
    restore.current = document.activeElement;
    input.current?.focus();
    return () => {
      if (restore.current instanceof HTMLElement) {
        restore.current.focus();
      }
    };
  }, []);

  useEffect(() => {
    if (active >= filtered.length) {
      setActive(Math.max(0, filtered.length - 1));
    }
  }, [active, filtered.length]);

  const run = (command: CommandSpec | undefined): void => {
    if (command === undefined || !command.enabled) {
      return;
    }
    onClose();
    command.run();
  };

  const onKey = (event: KeyboardEvent<HTMLInputElement>): void => {
    if (event.key === "ArrowDown") {
      event.preventDefault();
      setActive((current) => Math.min(filtered.length - 1, current + 1));
    } else if (event.key === "ArrowUp") {
      event.preventDefault();
      setActive((current) => Math.max(0, current - 1));
    } else if (event.key === "Home") {
      event.preventDefault();
      setActive(0);
    } else if (event.key === "End") {
      event.preventDefault();
      setActive(Math.max(0, filtered.length - 1));
    } else if (event.key === "Enter") {
      event.preventDefault();
      run(filtered[active]);
    } else if (event.key === "Escape") {
      event.stopPropagation();
      onClose();
    } else if (event.key === "Tab") {
      // The combobox pattern keeps focus in the field while the list is open.
      event.preventDefault();
    }
  };

  return (
    <div
      className="command-palette-backdrop"
      role="presentation"
      onMouseDown={onClose}
    >
      <div
        className="command-palette"
        role="dialog"
        aria-modal="true"
        aria-label="Command palette"
        onMouseDown={(event) => event.stopPropagation()}
      >
        <input
          ref={input}
          type="text"
          role="combobox"
          aria-label="Search commands"
          aria-autocomplete="list"
          aria-expanded="true"
          aria-controls="command-palette-list"
          aria-activedescendant={
            filtered[active] === undefined ? undefined : `command-${filtered[active].id}`
          }
          placeholder="Type a command…"
          value={query}
          onChange={(event) => {
            setQuery(event.target.value);
            setActive(0);
          }}
          onKeyDown={onKey}
        />
        <ul id="command-palette-list" className="command-palette-list" role="listbox" aria-label="Commands">
          {filtered.length === 0 ? (
            <li className="command-empty" role="presentation">
              no command matches
            </li>
          ) : (
            filtered.map((command, index) => (
              <li
                key={command.id}
                id={`command-${command.id}`}
                role="option"
                aria-selected={index === active}
                aria-disabled={!command.enabled}
                className={`command-option${command.enabled ? "" : " off"}${index === active ? " active" : ""}`}
                onMouseEnter={() => setActive(index)}
                onMouseDown={() => run(command)}
              >
                <span>{command.label}</span>
                <span className="command-meta">
                  {MENU_LABEL[command.menu] ?? command.menu}
                  {command.shortcut === undefined ? "" : ` · ${command.shortcut}`}
                </span>
              </li>
            ))
          )}
        </ul>
      </div>
    </div>
  );
}
