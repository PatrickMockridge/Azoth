/**
 * The bundled faces, read as text.
 *
 * **A token that names a family nobody declares is a fallback nobody chose.** The three layouts
 * are supposed to be typographically distinct; if a `@font-face` URL is stale or a token names a
 * family the CSS never loads, the browser silently falls back and two layouts converge again.
 * This reads `fonts.css` and the effective tokens, and checks the files exist.
 */

import { existsSync } from "node:fs";
import { resolve } from "node:path";

import { expect, it } from "vitest";

import { read, resolve as effective, type Combination } from "./css";

const FONT_STYLES = read("fonts.css");

interface Face {
  family: string;
  file: string;
}

const FACES: readonly Face[] = [...FONT_STYLES.matchAll(/@font-face\s*\{([^}]*)\}/g)].map(
  (match) => {
    const body = match[1] ?? "";
    return {
      family: /font-family:\s*"([^"]+)"/.exec(body)?.[1] ?? "",
      file: /url\("\.\.\/fonts\/([^"]+\.(?:woff2|ttf))"\)/.exec(body)?.[1] ?? "",
    };
  },
);

const LIGHT: readonly Combination[] = [
  { layout: "xp", theme: "light" },
  { layout: "mac", theme: "light" },
  { layout: "cyberpunk", theme: "light" },
];

/** The first family in a font stack, with its quotes removed. */
function firstFamily(stack: string): string {
  return (stack.split(",")[0] ?? "").trim().replace(/^["']|["']$/g, "");
}

it("loads a file for every @font-face", () => {
  expect(FACES.length).toBeGreaterThan(0);
  for (const face of FACES) {
    expect(face.family, "a face has no family").not.toBe("");
    expect(face.file, `${face.family} has no url`).not.toBe("");
    expect(
      existsSync(resolve(process.cwd(), "src/fonts", face.file)),
      `${face.family}: src/fonts/${face.file} is missing`,
    ).toBe(true);
  }
});

it("names a declared family first in every layout token", () => {
  const declared = new Set(FACES.map((face) => face.family));
  for (const combination of LIGHT) {
    const tokens = effective(combination);
    for (const name of ["font-ui", "font-display", "font-mono"]) {
      const stack = tokens.get(name) ?? "";
      const first = firstFamily(stack);
      expect(
        declared.has(first),
        `${combination.layout}: --${name} starts with "${first}", which no @font-face declares`,
      ).toBe(true);
    }
  }
});

it("keeps the chrome distinct and the body readable", () => {
  // The body face is the readable one; Macintosh and Cyberpunk deliberately share it. The
  // identity lives in the title/menu face, which must differ in every layout.
  const uis = LIGHT.map((combination) => firstFamily(effective(combination).get("font-ui") ?? ""));
  expect(new Set(uis)).toEqual(new Set(["Wine Tahoma", "IBM Plex Sans"]));
  const displays = LIGHT.map((combination) =>
    firstFamily(effective(combination).get("font-display") ?? ""),
  );
  expect(new Set(displays)).toEqual(new Set(["Wine Tahoma", "Sysfont", "Rajdhani"]));
});
