// @vitest-environment jsdom
/**
 * The appearance axis: what is read, what is written, and the pre-paint copy.
 *
 * **The module and `index.html` are two writers of the same four values**, which is deliberate:
 * the inline script exists so the first paint is the reader's choice rather than the default. A
 * drift between them would show as a flash on every reload, so the last test reads the HTML and
 * holds the copy to the module's key and defaults.
 */

import { existsSync, readFileSync } from "node:fs";
import { resolve } from "node:path";

import { beforeEach, expect, it } from "vitest";

import {
  DEFAULT,
  applyAppearance,
  isPreview,
  readAppearance,
} from "../src/state/appearance";

beforeEach(() => {
  window.localStorage.clear();
  window.history.replaceState({}, "", "/");
  delete document.documentElement.dataset.layout;
  delete document.documentElement.dataset.theme;
});

it("opens in Windows XP light when nothing was chosen", () => {
  expect(readAppearance()).toEqual({ layout: "xp", theme: "light" });
  expect(DEFAULT).toEqual({ layout: "xp", theme: "light" });
});

it("reads a stored choice and falls back field by field", () => {
  window.localStorage.setItem(
    "azoth.appearance",
    JSON.stringify({ layout: "mac", theme: "dark" }),
  );
  expect(readAppearance()).toEqual({ layout: "mac", theme: "dark" });

  window.localStorage.setItem("azoth.appearance", JSON.stringify({ layout: "mac" }));
  expect(readAppearance()).toEqual({ layout: "mac", theme: "light" });

  window.localStorage.setItem(
    "azoth.appearance",
    JSON.stringify({ layout: "windows", theme: "sepia" }),
  );
  expect(readAppearance()).toEqual({ layout: "xp", theme: "light" });
});

it("survives a corrupt store and a store that throws", () => {
  window.localStorage.setItem("azoth.appearance", "{not json");
  expect(readAppearance()).toEqual({ layout: "xp", theme: "light" });

  const original = Storage.prototype.getItem;
  Storage.prototype.getItem = () => {
    throw new Error("storage is refused");
  };
  try {
    expect(readAppearance()).toEqual({ layout: "xp", theme: "light" });
    expect(() => applyAppearance({ layout: "mac", theme: "dark" })).not.toThrow();
  } finally {
    Storage.prototype.getItem = original;
  }
});

it("migrates the old colour-only key onto Windows XP", () => {
  window.localStorage.setItem("azoth.theme", "dark");
  expect(readAppearance()).toEqual({ layout: "xp", theme: "dark" });
});

it("lets the URL preview a combination without rewriting the saved one", () => {
  window.localStorage.setItem(
    "azoth.appearance",
    JSON.stringify({ layout: "xp", theme: "neutral" }),
  );
  window.history.replaceState({}, "", "/?layout=cyberpunk&theme=dark");
  expect(isPreview()).toBe(true);
  expect(readAppearance()).toEqual({ layout: "cyberpunk", theme: "dark" });
  expect(window.localStorage.getItem("azoth.appearance")).toBe(
    JSON.stringify({ layout: "xp", theme: "neutral" }),
  );

  window.history.replaceState({}, "", "/");
  expect(isPreview()).toBe(false);
});

it("puts a choice in force and remembers it, unless asked not to", () => {
  applyAppearance({ layout: "cyberpunk", theme: "dark" });
  expect(document.documentElement.dataset.layout).toBe("cyberpunk");
  expect(document.documentElement.dataset.theme).toBe("dark");
  expect(window.localStorage.getItem("azoth.appearance")).toBe(
    JSON.stringify({ layout: "cyberpunk", theme: "dark" }),
  );

  applyAppearance({ layout: "mac", theme: "light" }, { persist: false });
  expect(document.documentElement.dataset.layout).toBe("mac");
  expect(window.localStorage.getItem("azoth.appearance")).toBe(
    JSON.stringify({ layout: "cyberpunk", theme: "dark" }),
  );
});

it("links a browser icon that the build carries", () => {
  const html = readFileSync(resolve(process.cwd(), "index.html"), "utf8");
  expect(html).toContain('rel="icon"');
  expect(html).toContain("favicon.svg");
  expect(html).toContain("favicon.ico");
  expect(html).toContain("apple-touch-icon.png");
  for (const name of ["favicon.svg", "favicon.ico", "apple-touch-icon.png"]) {
    expect(existsSync(resolve(process.cwd(), "public", name)), name).toBe(true);
  }
});

it("names the same key and defaults in the pre-paint script", () => {
  const html = readFileSync(resolve(process.cwd(), "index.html"), "utf8");
  expect(html).toContain("azoth.appearance");
  expect(html).toContain("azoth.theme");
  expect(html).toContain('layout = "xp"');
  expect(html).toContain('theme = "light"');
  expect(html).toContain('content="#ece9d8"');

  // The script must run before the module, or the first paint is the default.
  const script = html.indexOf("azoth.appearance");
  const module = html.indexOf('src="/src/main.tsx"');
  expect(script).toBeGreaterThan(-1);
  expect(module).toBeGreaterThan(script);
});
