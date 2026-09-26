/**
 * The stylesheet's own gate, read as text.
 *
 * **Nothing else tests the CSS at all.** It is imported by `src/main.tsx` and by nothing the tests
 * render, and jsdom applies no stylesheet and reports no colour - so a mistyped token, a colour
 * written into a rule, or a duration that escapes the reduced-motion block would be invisible to
 * every other test in this directory. These four readings are what is mechanically checkable about
 * a token layer; the *look* of it is not, and no assertion here claims to have seen it.
 */

import { readFileSync } from "node:fs";
import { resolve } from "node:path";

import { expect, it } from "vitest";

/** The stylesheet files, by name, relative to the package root the runner is started in. */
const FILES = ["tokens.css", "base.css", "chrome.css", "canvas.css", "forms.css"] as const;

const read = (name: string): string =>
  readFileSync(resolve(process.cwd(), "src/styles", name), "utf8");

/**
 * The declarations inside one top-level block, by selector.
 *
 * A hand-rolled scanner rather than PostCSS, which would be a dependency for four assertions - and
 * enough for this: it strips comments, then walks top-level blocks by brace depth, so a `:root`
 * inside a `@media` is *not* one of these, which is what the checks below want.
 */
function block(css: string, selector: string): string {
  const clean = css.replace(/\/\*[\s\S]*?\*\//g, "");
  let cursor = 0;
  while (cursor < clean.length) {
    const open = clean.indexOf("{", cursor);
    if (open === -1) {
      break;
    }
    const found = clean.slice(cursor, open).trim();
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
    if (found.endsWith(selector)) {
      return clean.slice(open + 1, end - 1);
    }
    cursor = end;
  }
  throw new Error(`no \`${selector}\` block in the stylesheet`);
}

/** The one capture group of a match, which every pattern here always has. */
const group = (match: RegExpMatchArray): string => match[1] ?? "";

const names = (body: string): string[] =>
  [...body.matchAll(/--([a-z0-9-]+)\s*:/g)].map(group);

/** Every custom property the stylesheets *use*, across all five files. */
function referenced(): Set<string> {
  const used = new Set<string>();
  for (const file of FILES) {
    for (const match of read(file).matchAll(/var\(--([a-z0-9-]+)/g)) {
      used.add(group(match));
    }
  }
  return used;
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
 * **The light theme overrides a subset, and every name in it has to exist.** The non-colour tokens
 * - the spacing scale, the type scale, the durations - are inherited deliberately, so this is a
 * check that a *typo* cannot pass rather than a demand that both blocks be the same length.
 */
it("overrides only names the dark theme declares", () => {
  const dark = names(block(read("tokens.css"), ":root"));
  const light = names(block(read("tokens.css"), '[data-theme="light"]'));
  const invented = light.filter((name) => !dark.includes(name));
  expect(invented).toEqual([]);
  // And a light theme that declared nothing would pass the check above, so it has to be a palette.
  expect(light.length).toBeGreaterThan(20);
});

/**
 * **A colour written into a rule is a value the other theme cannot reach**, which is the whole
 * point of the layer - and the same goes for a duration, which the reduced-motion block cannot
 * zero: it zeroes the tokens, so a literal `140ms` keeps animating for a user who asked it not to.
 */
it("writes no literal colour or duration outside the token file", () => {
  const colour = /#[0-9a-fA-F]{3,8}\b|\brgba?\(|\bhsla?\(/;
  const duration = /[0-9.]+m?s\b/;
  for (const file of FILES.filter((name) => name !== "tokens.css")) {
    const css = read(file).replace(/\/\*[\s\S]*?\*\//g, "");
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
