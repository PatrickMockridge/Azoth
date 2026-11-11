/**
 * The appearance control: one layout, one colour, chosen from the title strip.
 *
 * **A dialog rather than a menu.** The choices are two orthogonal groups of toggle buttons, and
 * `role="menu"` would promise arrow-key roving that grouped buttons do not have. Tab is the
 * keyboard the buttons already answer, Escape closes, and focus returns to the trigger.
 *
 * **A swatch reads the real token blocks.** Setting `data-layout` and `data-theme` on a span makes
 * the same CSS that paints the editor paint the preview, so a layout cannot look one way in the
 * chooser and another in the window.
 */

import { useEffect, useRef, useState } from "react";

import type { Appearance, Layout, Theme } from "../state/appearance";

const LAYOUTS: readonly { id: Layout; name: string }[] = [
  { id: "mac", name: "Macintosh" },
  { id: "xp", name: "Windows XP" },
  { id: "cyberpunk", name: "Cyberpunk" },
];

const THEMES: readonly { id: Theme; name: string }[] = [
  { id: "light", name: "Light" },
  { id: "neutral", name: "Neutral" },
  { id: "dark", name: "Dark" },
];

export function AppearanceMenu({
  appearance,
  onChange,
}: {
  appearance: Appearance;
  onChange: (appearance: Appearance) => void;
}) {
  const [open, setOpen] = useState(false);
  const trigger = useRef<HTMLButtonElement>(null);
  const panel = useRef<HTMLDivElement>(null);
  const first = useRef<HTMLButtonElement>(null);

  useEffect(() => {
    if (!open) {
      return;
    }
    first.current?.focus();
    const onKey = (event: KeyboardEvent) => {
      if (event.key === "Escape") {
        // Another Escape listener (a report, the overview) must not also close its layer.
        event.stopImmediatePropagation();
        setOpen(false);
        trigger.current?.focus();
      }
    };
    const onPointer = (event: PointerEvent) => {
      const target = event.target;
      if (
        target instanceof Node &&
        !panel.current?.contains(target) &&
        !trigger.current?.contains(target)
      ) {
        setOpen(false);
      }
    };
    window.addEventListener("keydown", onKey);
    window.addEventListener("pointerdown", onPointer);
    return () => {
      window.removeEventListener("keydown", onKey);
      window.removeEventListener("pointerdown", onPointer);
    };
  }, [open]);

  return (
    <div className="appearance">
      <button
        ref={trigger}
        type="button"
        className="button"
        aria-haspopup="dialog"
        aria-expanded={open}
        title="choose the layout and the colour"
        onClick={() => setOpen((current) => !current)}
      >
        Appearance
      </button>
      {open ? (
        <div ref={panel} className="appearance-panel" role="dialog" aria-label="Appearance">
          <div className="appearance-group" role="group" aria-labelledby="appearance-layout">
            <span className="appearance-label" id="appearance-layout">
              Layout
            </span>
            <div className="appearance-choices">
              {LAYOUTS.map((choice, index) => (
                <button
                  key={choice.id}
                  ref={index === 0 ? first : undefined}
                  type="button"
                  className="appearance-choice"
                  aria-pressed={appearance.layout === choice.id}
                  onClick={() => onChange({ ...appearance, layout: choice.id })}
                >
                  <span
                    className="appearance-swatch"
                    data-layout={choice.id}
                    data-theme={appearance.theme}
                    aria-hidden="true"
                  />
                  <span>{choice.name}</span>
                </button>
              ))}
            </div>
          </div>
          <div className="appearance-group" role="group" aria-labelledby="appearance-colour">
            <span className="appearance-label" id="appearance-colour">
              Colour
            </span>
            <div className="appearance-choices">
              {THEMES.map((choice) => (
                <button
                  key={choice.id}
                  type="button"
                  className="appearance-choice colour"
                  aria-pressed={appearance.theme === choice.id}
                  onClick={() => onChange({ ...appearance, theme: choice.id })}
                >
                  <span
                    className="appearance-swatch colour"
                    data-layout={appearance.layout}
                    data-theme={choice.id}
                    aria-hidden="true"
                  />
                  <span>{choice.name}</span>
                </button>
              ))}
            </div>
          </div>
        </div>
      ) : null}
    </div>
  );
}
