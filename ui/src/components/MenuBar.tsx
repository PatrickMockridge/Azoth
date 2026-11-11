/**
 * The application menu bar: period chrome over the same commands as the toolbar.
 *
 * **A menu item is not a button.** The rendered cases address `Solve`, `New`, `Open`, `Demo`,
 * `Save`, `Undo`, `Redo`, `Report` and `Appearance` by exact accessible name and expect one
 * visible button; a menu that mirrored them as buttons would multiply every query. The items are
 * `role="menuitem"`, live only while the menu is open, and call the same handlers.
 *
 * **The bar is inside the title strip**, after the wordmark. Macintosh keeps it on the title row;
 * XP gives it a full-width row under the title; Cyberpunk hides it and uses the palette.
 */

import { useEffect, useRef, useState } from "react";

import type { MenuSpec } from "../commands/model";

export function MenuBar({ menus }: { menus: readonly MenuSpec[] }) {
  const [open, setOpen] = useState<string | null>(null);
  const [submenu, setSubmenu] = useState<string | null>(null);
  const triggers = useRef(new Map<string, HTMLButtonElement | null>());
  const panel = useRef<HTMLDivElement | null>(null);

  useEffect(() => {
    if (open === null) {
      return;
    }
    const onKey = (event: KeyboardEvent) => {
      if (event.key === "Escape") {
        event.stopImmediatePropagation();
        if (submenu !== null) {
          setSubmenu(null);
          return;
        }
        setOpen(null);
        triggers.current.get(open)?.focus();
      }
    };
    const onDown = (event: MouseEvent) => {
      const target = event.target;
      if (target instanceof Node && !panel.current?.contains(target)) {
        const trigger = triggers.current.get(open);
        if (trigger === null || !trigger?.contains(target)) {
          setOpen(null);
          setSubmenu(null);
        }
      }
    };
    window.addEventListener("keydown", onKey);
    window.addEventListener("mousedown", onDown);
    return () => {
      window.removeEventListener("keydown", onKey);
      window.removeEventListener("mousedown", onDown);
    };
  }, [open, submenu]);

  if (menus.length === 0) {
    return null;
  }

  return (
    <nav className="menubar" aria-label="Application menu">
      {menus.map((menu) => (
        <div className="menu-wrap" key={menu.id}>
          <button
            type="button"
            className={`menu-trigger${open === menu.id ? " open" : ""}`}
            aria-haspopup="menu"
            aria-expanded={open === menu.id}
            ref={(element) => {
              triggers.current.set(menu.id, element);
            }}
            onClick={() => {
              setOpen(open === menu.id ? null : menu.id);
              setSubmenu(null);
            }}
          >
            {menu.label}
          </button>
          {open === menu.id ? (
            <div
              className="menu"
              role="menu"
              aria-label={menu.label}
              ref={panel}
              onMouseDown={(event) => event.stopPropagation()}
            >
              {menu.commands.map((command) => (
                <div className="menu-subwrap" key={command.id}>
                  <button
                    type="button"
                    role="menuitem"
                    aria-disabled={!command.enabled}
                    aria-haspopup={command.children === undefined ? undefined : "menu"}
                    aria-expanded={command.children === undefined ? undefined : submenu === command.id}
                    className={`menu-item${command.enabled ? "" : " off"}${command.children === undefined ? "" : " has-submenu"}`}
                    title={command.hint}
                    onMouseEnter={() => {
                      if (command.children !== undefined) {
                        setSubmenu(command.id);
                      }
                    }}
                    onClick={() => {
                      if (!command.enabled) {
                        return;
                      }
                      if (command.children !== undefined) {
                        setSubmenu(submenu === command.id ? null : command.id);
                        return;
                      }
                      setOpen(null);
                      setSubmenu(null);
                      command.run();
                    }}
                  >
                    <span>{command.label}</span>
                    {command.shortcut === undefined ? (
                      command.children === undefined ? null : (
                        <span className="menu-shortcut" aria-hidden="true">
                          ▸
                        </span>
                      )
                    ) : (
                      <span className="menu-shortcut" aria-hidden="true">
                        {command.shortcut}
                      </span>
                    )}
                  </button>
                  {command.children !== undefined && submenu === command.id ? (
                    <div className="menu submenu" role="menu" aria-label={command.label}>
                      {command.children.map((child) => (
                        <button
                          key={child.id}
                          type="button"
                          role="menuitem"
                          aria-disabled={!child.enabled}
                          className={`menu-item${child.enabled ? "" : " off"}`}
                          title={child.hint}
                          onClick={() => {
                            if (!child.enabled) {
                              return;
                            }
                            setOpen(null);
                            setSubmenu(null);
                            child.run();
                          }}
                        >
                          <span>{child.label}</span>
                        </button>
                      ))}
                    </div>
                  ) : null}
                </div>
              ))}
            </div>
          ) : null}
        </div>
      ))}
    </nav>
  );
}
