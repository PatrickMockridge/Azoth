/**
 * A layout may hide decoration; it may not hide a control.
 *
 * **The regression this exists for.** The Macintosh workspace switch was hidden to make the menu
 * bar fit; jsdom applies no stylesheet, so every existing test still found the button and the
 * loss was invisible. This reads the layout CSS as text and refuses a `display: none` on anything
 * a reader acts through. The allowlist is deliberately small and each entry is decoration, a
 * caption, or a user-controlled panel toggle - not a control the reader needs to reach a command.
 */

import { describe, expect, it } from "vitest";

import { blocks, read } from "./css";

const LAYOUTS = ["xp", "mac", "cyberpunk"] as const;

/** The parts of a control's selector that a hide must never name. */
const CONTROLS = [
  ".workspace-switch",
  ".splitter",
  ".hud-toggle",
  ".hud-controls",
  ".menu-trigger",
  ".tab",
  ".entry",
  ".tree-row",
  "button",
  "select",
  "input",
];

/** The only things a layout file may remove: information, captions, and toggled panels. */
const DECORATION = [
  ".build",
  ".app.off",
  ".app[hidden]",
  ".window-title",
  ".window-grip",
  ".spacer",
  ".menubar",
  '.ribbon-label',
  '[data-hud-',
];

const hidden = (css: string): string[] =>
  blocks(css)
    .filter((entry) => /display:\s*none/.test(entry.body))
    .map((entry) => entry.selector);

describe("the layouts' hidden selectors", () => {
  for (const layout of LAYOUTS) {
    it(`${layout} hides no control`, () => {
      const selectors = hidden(read(`layouts/${layout}.css`));
      for (const selector of selectors) {
        for (const control of CONTROLS) {
          expect(
            selector.includes(control),
            `${layout}.css hides \`${selector}\`, which names the control \`${control}\``,
          ).toBe(false);
        }
      }
      // The other half: every hide is one of the allowed kinds, so a new hide has to be argued
      // for here rather than appearing silently.
      for (const selector of selectors) {
        expect(
          DECORATION.some((allowed) => selector.includes(allowed)),
          `${layout}.css hides \`${selector}\`, which is not decoration or a panel toggle`,
        ).toBe(true);
      }
    });
  }

  it("keeps every command surface reachable", () => {
    // The menus stay in the Macintosh bar and the XP frame; Cyberpunk replaces them with the
    // palette and the ribbon, and the palette's button is visible in the strip.
    for (const layout of ["xp", "mac"] as const) {
      expect(hidden(read(`layouts/${layout}.css`)).some((s) => s.includes(".menubar"))).toBe(
        false,
      );
    }
    // The HUD's rail switches and palette entry are shown by the shared command stylesheet in
    // Cyberpunk; no layout file may hide them.
    expect(read("commands.css")).toContain('[data-layout="cyberpunk"] .hud-controls');
    expect(read("layouts/cyberpunk.css")).not.toMatch(/\.splitter[^{]*\{[^}]*display:\s*none/);
    expect(read("layouts/cyberpunk.css")).not.toMatch(
      /\.workspace-switch[^{]*\{[^}]*display:\s*none/,
    );
  });

  it("draws no close or zoom box that does nothing", () => {
    // The Macintosh title bar is a drag handle and the grip is the resize handle. A box drawn
    // like a control but wired to nothing is a false affordance, not decoration.
    const mac = read("layouts/mac.css");
    expect(mac).not.toContain(".window-title::before");
    expect(mac).not.toContain(".window-title::after");
  });
});
