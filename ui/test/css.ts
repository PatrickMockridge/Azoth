/**
 * The stylesheet, read as text.
 *
 * jsdom applies no stylesheet and reports no colour, so the two tests that read CSS - the token
 * contract and the contrast measurement - share this small hand-rolled scanner. It is not
 * PostCSS: it strips comments and walks top-level blocks by brace depth, which is enough for a
 * token layer and keeps the tests dependency-free.
 */

import { readFileSync } from "node:fs";
import { resolve as resolvePath } from "node:path";

export function read(name: string): string {
  return readFileSync(resolvePath(process.cwd(), "src/styles", name), "utf8");
}

/** Every stylesheet the editor imports, taken from the one list that cannot drift. */
export const FILES: readonly string[] = [
  ...read("index.css").matchAll(/@import\s+"\.\/([^"]+)";/g),
].map((match) => match[1] ?? "");

export interface Block {
  selector: string;
  body: string;
}

/** The top-level blocks of a stylesheet, in source order, without comments. */
export function blocks(css: string): Block[] {
  const clean = css.replace(/\/\*[\s\S]*?\*\//g, "");
  const found: Block[] = [];
  let cursor = 0;
  while (cursor < clean.length) {
    const open = clean.indexOf("{", cursor);
    if (open === -1) {
      break;
    }
    const selector = clean.slice(cursor, open).trim();
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
    found.push({ selector, body: clean.slice(open + 1, end - 1) });
    cursor = end;
  }
  return found;
}

/** The first top-level block with exactly this selector. */
export function block(css: string, selector: string): string {
  const found = blocks(css).find((entry) => entry.selector === selector);
  if (found === undefined) {
    throw new Error(`no \`${selector}\` block in the stylesheet`);
  }
  return found.body;
}

/** Every block with this selector, joined: a cell may be split into shape and tone blocks. */
export function cell(css: string, selector: string): string {
  return blocks(css)
    .filter((entry) => entry.selector === selector)
    .map((entry) => entry.body)
    .join("\n");
}

/** The custom-property declarations of a body, in source order. */
export function declarations(body: string): [string, string][] {
  return [...body.matchAll(/--([a-z0-9-]+)\s*:\s*([^;{}]+);/g)].map((match) => [
    match[1] ?? "",
    (match[2] ?? "").trim(),
  ]);
}

export function names(body: string): string[] {
  return declarations(body).map(([name]) => name);
}

/** Every custom property the stylesheets use. */
export function referenced(): Set<string> {
  const used = new Set<string>();
  for (const file of FILES) {
    for (const match of read(file).matchAll(/var\(--([a-z0-9-]+)/g)) {
      used.add(match[1] ?? "");
    }
  }
  return used;
}

export type Layout = "xp" | "mac" | "cyberpunk";
export type Tone = "light" | "neutral" | "dark";

export interface Combination {
  layout: Layout;
  theme: Tone;
}

export const COMBINATIONS: readonly Combination[] = (
  ["xp", "mac", "cyberpunk"] as const
).flatMap((layout) =>
  (["light", "neutral", "dark"] as const).map((theme) => ({ layout, theme })),
);

/** The selector whose declarations are a cell's own row in the contract table. */
export function cellSelector({ layout, theme }: Combination): string {
  if (layout === "xp") {
    return theme === "light" ? ":root" : `[data-theme="${theme}"]`;
  }
  return theme === "light"
    ? `[data-layout="${layout}"]`
    : `[data-layout="${layout}"][data-theme="${theme}"]`;
}

/** Whether a top-level selector's declarations apply to one combination. */
function applies(selector: string, combination: Combination): boolean {
  if (selector === ":root") {
    return true;
  }
  const layout = /\[data-layout="([a-z]+)"\]/.exec(selector)?.[1];
  const theme = /\[data-theme="([a-z]+)"\]/.exec(selector)?.[1];
  if (layout === undefined && theme === undefined) {
    return false;
  }
  if (selector.replace(/\[data-(?:layout|theme)="[a-z]+"\]/g, "").trim() !== "") {
    return false;
  }
  return (
    (layout === undefined || layout === combination.layout) &&
    (theme === undefined || theme === combination.theme)
  );
}

/** The effective token map of one combination, in CSS source order. */
export function resolve(combination: Combination): Map<string, string> {
  const map = new Map<string, string>();
  for (const entry of blocks(read("tokens.css"))) {
    if (!applies(entry.selector, combination)) {
      continue;
    }
    for (const [name, value] of declarations(entry.body)) {
      map.set(name, value);
    }
  }
  return map;
}

/** Follow `var()` aliases to a literal, refusing anything a contrast measurement cannot read. */
export function literal(combination: Combination, name: string): string {
  const map = resolve(combination);
  const seen = new Set<string>();
  let current = name;
  while (true) {
    if (seen.has(current)) {
      throw new Error(`${combination.layout}/${combination.theme}: var(--${current}) is circular`);
    }
    seen.add(current);
    const value = map.get(current);
    if (value === undefined) {
      throw new Error(`${combination.layout}/${combination.theme}: --${current} is not declared`);
    }
    const alias = /^var\(--([a-z0-9-]+)\)$/.exec(value);
    if (alias !== null) {
      current = alias[1] ?? "";
      continue;
    }
    return value;
  }
}
