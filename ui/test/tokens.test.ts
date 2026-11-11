/**
 * The stylesheet's own gate, read as text.
 *
 * **Nothing else tests the CSS at all.** It is imported by `src/main.tsx` and by nothing the tests
 * render, and jsdom applies no stylesheet and reports no colour - so a mistyped token, a colour
 * written into a rule, or a duration that escapes the reduced-motion block would be invisible to
 * every other test in this directory. These readings are what is mechanically checkable about a
 * token layer; the *look* of it is not, and no assertion here claims to have seen it.
 */

import { expect, it } from "vitest";

import {
  COMBINATIONS,
  FILES,
  block,
  blocks,
  cell,
  cellSelector,
  names,
  read,
  referenced,
  resolve,
} from "./css";

/** The colour and chrome values every tone cell must supply, so none leaks from another layout. */
const TONE_TOKENS: readonly string[] = [
  "bg-app",
  "bg-canvas",
  "bg-panel",
  "bg-panel-2",
  "bg-raised",
  "bg-hover",
  "bg-input",
  "bg-selected",
  "line",
  "line-strong",
  "grid-minor",
  "text",
  "text-strong",
  "text-dim",
  "text-faint",
  "accent",
  "accent-hi",
  "accent-ink",
  "focus",
  "ok",
  "warn",
  "error",
  "edge",
  "edge-recycle",
  "selected-text",
  "titlebar-a",
  "titlebar-b",
  "titlebar-c",
  "titlebar-d",
  "titlebar-text",
  "titlebar-muted",
  "ctl-face",
  "ctl-face-hi",
  "ctl-face-low",
  "ctl-text",
  "ctl-border",
  "well",
  "well-border",
  "pill-face",
  "pill-border",
  "tab-face",
  "tab-face-active",
  "tab-text",
  "tab-text-active",
  "tab-border",
  "table-face",
  "table-head-face",
  "table-head-text",
  "table-stripe",
  "node-face",
  "node-head",
  "node-head-text",
  "node-border",
  "pin-a",
  "pin-b",
  "desktop-a",
  "desktop-b",
];

/** The shape and paint values each layout must supply; XP is `:root`. */
const LAYOUT_TOKENS: readonly string[] = [
  "font-ui",
  "font-display",
  "font-mono",
  "fs-2xs",
  "fs-xs",
  "fs-sm",
  "fs-md",
  "fs-lg",
  "fs-xl",
  "lh-tight",
  "lh-body",
  "ls-label",
  "label-transform",
  "row-h",
  "head-h",
  "ctl-h",
  "tab-h",
  "r1",
  "r2",
  "r3",
  "r4",
  "r-pill",
  "col-nav",
  "col-dock",
  "case-name-w",
  "desktop-gap",
  "node-radius",
  "handle-dot",
  "handle-radius",
  "scrollbar-size",
  "scrollbar-radius",
  "splitter-h",
  "bar-pad-y",
  "bar-pad-x",
  "focus-w",
  "focus-style",
  "focus-offset",
  "surface-app",
  "surface-workspace",
  "surface-canvas",
  "surface-panel",
  "surface-panel-2",
  "surface-raised",
  "surface-bar",
  "surface-tools",
  "surface-ribbon",
  "surface-dock",
  "surface-status",
  "surface-tabs",
  "surface-tab-active",
  "surface-menu",
  "surface-report",
  "surface-node",
  "surface-node-head",
  "surface-cell",
  "surface-cell-alt",
  "ctl-face-paint",
  "ctl-face-hover",
  "ctl-bevel",
  "well-shadow",
  "node-shadow",
  "menu-shadow",
  "e2",
  "e3",
  "e-node",
  "e-highlight",
  "frame",
  "frame-w",
  "bevel-light",
  "bevel-mid",
  "bevel-dark",
  "bevel-w",
  "scrollbar-thumb",
  "scrollbar-track",
  "glow-accent",
  "texture",
  "backdrop",
  "command-icon-size",
  "menubar-h",
  "menubar-ink",
  "surface-menubar",
  "window-title-h",
  "window-title-ink",
  "surface-window-title",
  "window-shadow",
  "window-grip",
  "t-fast",
  "t-mid",
  "t-slow",
  "t-blink",
  "t-pulse",
];

/** A stylesheet without its top-level `@media print` block, which is the fixed paper palette. */
function withoutPrint(css: string): string {
  const clean = css.replace(/\/\*[\s\S]*?\*\//g, "");
  const at = clean.indexOf("@media print");
  if (at === -1) {
    return clean;
  }
  const open = clean.indexOf("{", at);
  if (open === -1) {
    return clean.slice(0, at);
  }
  let depth = 1;
  let end = open + 1;
  while (end < clean.length && depth > 0) {
    if (clean[end] === "{") {
      depth += 1;
    } else if (clean[end] === "}") {
      depth -= 1;
    }
    end += 1;
  }
  return clean.slice(0, at) + clean.slice(end);
}

/**
 * **A `var()` whose name is never declared resolves to nothing**, and the property it was standing
 * for is then absent rather than wrong - a 1px border that is not there, a colour inherited from
 * its parent. Nothing renders a warning.
 */
it("declares every token the stylesheets use", () => {
  const declared = names(block(read("tokens.css"), ":root"));
  const missing = [...referenced()].filter((name) => !declared.includes(name));
  expect(missing).toEqual([]);
});

/**
 * **Every tone cell declares the whole tone set.** A missed name would inherit the tone of the
 * layout whose block happened to appear earlier in the cascade - the muddy compromise the nine
 * cells exist to end - so the check is completeness rather than a subset.
 */
it("gives every layout and tone the same tone contract", () => {
  const tokens = read("tokens.css");
  for (const combination of COMBINATIONS) {
    const selector = cellSelector(combination);
    const declared = names(cell(tokens, selector));
    const missing = TONE_TOKENS.filter((name) => !declared.includes(name));
    expect(
      missing,
      `${combination.layout}/${combination.theme} (${selector}) misses ${missing.join(", ")}`,
    ).toEqual([]);
  }
});

/** Every layout owns the same shape and paint names, so a typo cannot pass as a shape. */
it("gives every layout the same shape contract", () => {
  const tokens = read("tokens.css");
  for (const selector of [":root", '[data-layout="mac"]', '[data-layout="cyberpunk"]']) {
    const declared = names(cell(tokens, selector));
    const missing = LAYOUT_TOKENS.filter((name) => !declared.includes(name));
    expect(missing, `${selector} misses ${missing.join(", ")}`).toEqual([]);
  }
});

/** A block may override a name the contract declares; it may not invent one. */
it("overrides only names the contract declares", () => {
  const contract = names(block(read("tokens.css"), ":root"));
  for (const entry of blocks(read("tokens.css"))) {
    if (
      !entry.selector.startsWith(":root") &&
      !entry.selector.startsWith("[data-layout") &&
      !entry.selector.startsWith("[data-theme")
    ) {
      continue;
    }
    const invented = names(entry.body).filter((name) => !contract.includes(name));
    expect(invented, `${entry.selector} invents ${invented.join(", ")}`).toEqual([]);
  }
});

/** Every referenced token resolves in every one of the nine combinations. */
it("resolves every referenced token in all nine combinations", () => {
  const used = [...referenced()];
  for (const combination of COMBINATIONS) {
    const effective = resolve(combination);
    const missing = used.filter((name) => !effective.has(name));
    expect(
      missing,
      `${combination.layout}/${combination.theme} resolves no ${missing.join(", ")}`,
    ).toEqual([]);
  }
});

/**
 * **A colour written into a rule is a value no other cell can reach**, which is the whole point of
 * the layer - and the same goes for a duration, which the reduced-motion block cannot zero: it
 * zeroes the tokens, so a literal `140ms` keeps animating for a user who asked it not to.
 */
it("writes no literal colour or duration outside the token file", () => {
  const colour = /#[0-9a-fA-F]{3,8}\b|\brgba?\(|\bhsla?\(/;
  const duration = /[0-9.]+m?s\b/;
  for (const file of FILES.filter((name) => name !== "tokens.css")) {
    const css = withoutPrint(read(file));
    for (const line of css.split("\n")) {
      // A declaration is where a value lives; a selector line cannot be one, because this
      // stylesheet has no `.` or `#` in front of a number.
      if (!line.includes(":")) {
        continue;
      }
      expect(colour.test(line), `${file}: ${line.trim()}`).toBe(false);
      expect(duration.test(line), `${file}: ${line.trim()}`).toBe(false);
    }
  }
});

it("covers the stylesheets the editor actually imports", () => {
  expect(FILES).toEqual(
    expect.arrayContaining([
      "appearance.css",
      "layouts/mac.css",
      "layouts/xp.css",
      "layouts/cyberpunk.css",
    ]),
  );
  expect(FILES.at(-1)).toBe("report.css");
});

/**
 * **The PFD plate is not `.unit`.** An unscoped `.unit` rule put a 168px node plate around every
 * unit string in a field, a table and a result; the plate is now `.unit-plate` and this keeps the
 * two vocabularies apart.
 */
it("gives the PFD plate its own class", () => {
  const canvas = read("canvas.css");
  expect(canvas).toMatch(/^\.unit-plate\s*\{/m);
  expect(canvas).not.toMatch(/^\.unit(?!-)/m);
});

/** Geometry the layouts vary is a token, not a number in a component rule. */
it("reads the layout's geometry from tokens", () => {
  expect(block(read("chrome.css"), ".body")).toMatch(
    /grid-template-columns:\s*var\(--col-nav\)\s+1fr\s+var\(--col-dock\)/,
  );
  expect(block(read("canvas.css"), ".unit-plate")).toMatch(/min-width:\s*var\(--node-w\)/);
  expect(block(read("canvas.css"), ".unit-plate")).toMatch(/height:\s*var\(--node-h\)/);
  expect(block(read("canvas.css"), ".react-flow__handle")).toMatch(
    /width:\s*var\(--handle-size\)/,
  );
  expect(read("base.css")).toMatch(/scrollbar-color:\s*var\(--scrollbar-thumb\)/);
});

it("gives the sheet name and picker search a themed well", () => {
  expect(block(read("calculations.css"), ".calc-sheet-name input")).toMatch(
    /background:\s*var\(--well\)/,
  );
  expect(block(read("calculations.css"), ".calc-sheet-name input")).toMatch(
    /border:\s*var\(--bevel-w\) solid var\(--well-border\)/,
  );
  expect(block(read("calculations.css"), ".calc-search")).toMatch(/background:\s*var\(--well\)/);
});

it("gives the report an opaque screen background", () => {
  expect(block(read("report.css"), ".report")).toMatch(/background:\s*var\(--surface-report\)/);
});

it("tells the browser which colour scheme is in force", () => {
  const tokens = read("tokens.css");
  expect(cell(tokens, ":root")).toContain("color-scheme: light");
  expect(cell(tokens, '[data-theme="neutral"]')).toContain("color-scheme: light");
  expect(cell(tokens, '[data-theme="dark"]')).toContain("color-scheme: dark");
  expect(cell(tokens, '[data-layout="mac"]')).toContain("color-scheme: light");
  expect(cell(tokens, '[data-layout="mac"][data-theme="dark"]')).toContain("color-scheme: dark");
  expect(cell(tokens, '[data-layout="cyberpunk"]')).toContain("color-scheme: light");
  expect(cell(tokens, '[data-layout="cyberpunk"][data-theme="dark"]')).toContain(
    "color-scheme: dark",
  );
});

/** Print is the one fixed palette, and the layout decoration must not reach it. */
it("keeps paper out of the layouts and above every cell", () => {
  const report = read("report.css");
  expect(report).toMatch(/@media print/);
  expect(report).toContain("[data-layout][data-theme]");
  expect(report).toContain("--surface-report: #fff");
  // A tone's table ink and a layout's uppercase labels must not reach paper.
  expect(report).toContain("--table-head-text: #333");
  expect(report).toContain("--label-transform: none");
  for (const file of FILES.filter((name) => name.startsWith("layouts/"))) {
    expect(read(file), `${file} has a print block`).not.toContain("@media print");
  }
});

/** Reduced motion zeroes every duration, including the decorative ones. */
it("zeroes every motion token under reduced motion", () => {
  const tokens = read("tokens.css");
  const at = tokens.indexOf("@media (prefers-reduced-motion: reduce)");
  expect(at).toBeGreaterThan(-1);
  const reduced = tokens.slice(at);
  for (const name of ["t-fast", "t-mid", "t-slow", "t-blink", "t-pulse"]) {
    expect(reduced, `--${name} is not zeroed`).toContain(`--${name}: 0ms`);
  }
  expect(reduced).toContain("animation: none");
});
