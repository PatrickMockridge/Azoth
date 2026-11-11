/**
 * The nine palettes, measured.
 *
 * **Contrast is computed, not asserted in a comment.** A theme is a set of flat colours, and two
 * of the period palettes (XP's grey-on-grey and Mac's Platinum) are exactly the sort a modern eye
 * finds faint; this holds every cell to the same pairs rather than trusting the one that was
 * checked when it was written. Gradients and `color-mix` paints are deliberately not paired - a
 * ratio is defined between two colours, and the title-bar stops are the anchors a gradient is
 * measured at.
 *
 * Thresholds are WCAG 2.1 AA: 4.5:1 for body text and 3:1 for large text, UI boundaries and a
 * focus indicator.
 */

import { expect, it } from "vitest";

import { COMBINATIONS, literal } from "./css";

const BODY = 4.5;
const LARGE = 3;

interface Pair {
  foreground: string;
  background: string;
  minimum: number;
  why: string;
}

const PAIRS: readonly Pair[] = [
  { foreground: "text", background: "bg-panel", minimum: BODY, why: "body text" },
  { foreground: "text-strong", background: "bg-panel", minimum: BODY, why: "strong text" },
  { foreground: "text-dim", background: "bg-panel", minimum: BODY, why: "dim text" },
  { foreground: "text-faint", background: "bg-panel", minimum: BODY, why: "faint text" },
  { foreground: "text", background: "bg-app", minimum: BODY, why: "text on the desktop" },
  { foreground: "text-dim", background: "bg-panel-2", minimum: BODY, why: "dim text on a header" },
  { foreground: "text-dim", background: "bg-raised", minimum: BODY, why: "dim text on a raised face" },
  { foreground: "text-dim", background: "bg-hover", minimum: BODY, why: "dim text under the pointer" },
  { foreground: "selected-text", background: "bg-selected", minimum: BODY, why: "selected row" },
  { foreground: "accent", background: "bg-panel", minimum: BODY, why: "link text" },
  { foreground: "accent-ink", background: "accent", minimum: BODY, why: "primary button" },
  { foreground: "accent-ink", background: "accent-hi", minimum: BODY, why: "primary button hover" },
  { foreground: "ok", background: "bg-panel", minimum: BODY, why: "success text" },
  { foreground: "warn", background: "bg-panel", minimum: BODY, why: "warning text" },
  { foreground: "error", background: "bg-panel", minimum: BODY, why: "error text" },
  { foreground: "ok", background: "pill-face", minimum: BODY, why: "success pill" },
  { foreground: "warn", background: "pill-face", minimum: BODY, why: "warning pill" },
  { foreground: "error", background: "pill-face", minimum: BODY, why: "error pill" },
  { foreground: "titlebar-text", background: "titlebar-a", minimum: BODY, why: "title bar top" },
  { foreground: "titlebar-text", background: "titlebar-b", minimum: BODY, why: "title bar upper" },
  { foreground: "titlebar-text", background: "titlebar-c", minimum: BODY, why: "title bar lower" },
  { foreground: "titlebar-text", background: "titlebar-d", minimum: BODY, why: "title bar bottom" },
  { foreground: "titlebar-muted", background: "titlebar-a", minimum: BODY, why: "muted title text" },
  { foreground: "titlebar-muted", background: "titlebar-b", minimum: BODY, why: "muted title text" },
  { foreground: "titlebar-muted", background: "titlebar-c", minimum: BODY, why: "muted title text" },
  { foreground: "titlebar-muted", background: "titlebar-d", minimum: BODY, why: "muted title text" },
  { foreground: "ctl-text", background: "ctl-face", minimum: BODY, why: "control label" },
  { foreground: "ctl-text", background: "ctl-face-hi", minimum: BODY, why: "control highlight" },
  { foreground: "ctl-text", background: "ctl-face-low", minimum: BODY, why: "control shadow" },
  { foreground: "tab-text", background: "tab-face", minimum: BODY, why: "inactive tab" },
  { foreground: "tab-text-active", background: "tab-face-active", minimum: BODY, why: "active tab" },
  { foreground: "focus", background: "bg-panel", minimum: LARGE, why: "focus ring" },
  { foreground: "focus", background: "bg-raised", minimum: LARGE, why: "focus ring" },
  { foreground: "edge", background: "bg-canvas", minimum: LARGE, why: "connection line" },
  { foreground: "edge-recycle", background: "bg-canvas", minimum: LARGE, why: "recycle line" },
];

/** The WCAG relative luminance of one sRGB channel. */
function channel(value: number): number {
  const v = value / 255;
  return v <= 0.04045 ? v / 12.92 : ((v + 0.055) / 1.055) ** 2.4;
}

function luminance(hex: string): number {
  const value = hex.replace("#", "");
  return (
    0.2126 * channel(parseInt(value.slice(0, 2), 16)) +
    0.7152 * channel(parseInt(value.slice(2, 4), 16)) +
    0.0722 * channel(parseInt(value.slice(4, 6), 16))
  );
}

function ratio(a: string, b: string): number {
  const [high, low] = [luminance(a), luminance(b)].sort((x, y) => y - x);
  return ((high ?? 0) + 0.05) / ((low ?? 0) + 0.05);
}

it("measures every text pair in every layout and colour", () => {
  for (const combination of COMBINATIONS) {
    for (const pair of PAIRS) {
      const foreground = literal(combination, pair.foreground);
      const background = literal(combination, pair.background);
      expect(foreground, `${combination.layout}/${combination.theme}: --${pair.foreground}`).toMatch(
        /^#[0-9a-fA-F]{6}$/,
      );
      expect(background, `${combination.layout}/${combination.theme}: --${pair.background}`).toMatch(
        /^#[0-9a-fA-F]{6}$/,
      );
      const measured = ratio(foreground, background);
      expect(
        measured,
        `${combination.layout}/${combination.theme}: ${pair.why} - ` +
          `--${pair.foreground} ${foreground} on --${pair.background} ${background} ` +
          `is ${measured.toFixed(2)}:1, needs ${pair.minimum}:1`,
      ).toBeGreaterThanOrEqual(pair.minimum);
    }
  }
});
