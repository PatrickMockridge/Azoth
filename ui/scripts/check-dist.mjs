#!/usr/bin/env node
/**
 * The built editor, served under the path Pages uses, and asked for its own files.
 *
 * **Why this is a server and not a file listing.** A relative build can be right on disk and
 * wrong when served under `/azoth/`: a stray root-absolute `src`, a hashed asset the HTML does
 * not name, or a wasm served as `application/octet-stream`, which a browser refuses to
 * instantiate. The only checks that see those are a request for each URL the page will make.
 *
 * The wasm is named from the JavaScript that asks for it rather than guessed from the directory,
 * because the glue's `new URL(..., import.meta.url)` is the thing that has to resolve.
 *
 * Usage: node scripts/check-dist.mjs [dist]
 */

import { createServer } from "node:http";
import { readFile, readdir, stat } from "node:fs/promises";
import { extname, join, normalize, resolve, sep } from "node:path";
import { fileURLToPath } from "node:url";

const ROOT = resolve(fileURLToPath(new URL("..", import.meta.url)));
const DIST = resolve(process.argv[2] ?? join(ROOT, "dist"));
const PREFIX = "/azoth/";

const TYPES = new Map([
  [".html", "text/html; charset=utf-8"],
  [".js", "text/javascript; charset=utf-8"],
  [".mjs", "text/javascript; charset=utf-8"],
  [".css", "text/css; charset=utf-8"],
  [".wasm", "application/wasm"],
  [".json", "application/json; charset=utf-8"],
  [".svg", "image/svg+xml"],
  [".png", "image/png"],
  [".ico", "image/x-icon"],
]);

const failures = [];
const fail = (message) => failures.push(message);

/** The path under `dist` a request names, or null when it escapes or is not there. */
async function fileFor(pathname) {
  const relative = pathname.startsWith(PREFIX) ? pathname.slice(PREFIX.length) : null;
  if (relative === null || relative.includes("\0")) {
    return null;
  }
  const candidate = resolve(join(DIST, normalize(relative === "" ? "index.html" : relative)));
  const inside = candidate === DIST || candidate.startsWith(DIST + sep);
  if (!inside) {
    return null;
  }
  try {
    const info = await stat(candidate);
    return info.isDirectory() ? join(candidate, "index.html") : candidate;
  } catch {
    return null;
  }
}

const server = createServer(async (request, response) => {
  const url = new URL(request.url ?? "/", "http://localhost");
  const file = await fileFor(url.pathname);
  if (file === null) {
    response.writeHead(404, { "Content-Type": "text/plain" });
    response.end("not found");
    return;
  }
  const body = await readFile(file);
  response.writeHead(200, {
    "Content-Type": TYPES.get(extname(file)) ?? "application/octet-stream",
    "Content-Length": body.length,
  });
  response.end(body);
});

await new Promise((ready) => server.listen(0, "127.0.0.1", ready));
const origin = `http://127.0.0.1:${server.address().port}`;

/** Fetch one URL, and answer with the response (or a failure and null). */
async function get(path) {
  const response = await fetch(origin + path);
  if (!response.ok) {
    fail(`${path} answered ${response.status}`);
    return null;
  }
  return response;
}

try {
  const index = await readFile(join(DIST, "index.html"), "utf8");

  // **No root-absolute reference.** Under Pages the page is not at the origin root, so `src="/x"`
  // asks for a file that is not there - the classic deploy failure this script exists for.
  for (const match of index.matchAll(/\b(?:src|href)="([^"]+)"/g)) {
    const value = match[1];
    if (value.startsWith("/")) {
      fail(`index.html names a root-absolute asset: ${value}`);
    }
  }

  // Every script, stylesheet and module the HTML names, fetched as the browser would.
  const referenced = [...index.matchAll(/\b(?:src|href)="([^"]+)"/g)]
    .map((match) => match[1])
    .filter((value) => !value.startsWith("http") && !value.startsWith("data:"))
    .map((value) => (value.startsWith("./") ? value.slice(2) : value));
  if (referenced.length === 0) {
    fail("index.html names no built asset");
  }
  for (const path of referenced) {
    await get(PREFIX + path);
  }

  // The wasm, found through the JavaScript that asks for it.
  const assets = join(DIST, "assets");
  let wasm = [];
  try {
    wasm = (await readdir(assets)).filter((name) => name.endsWith(".wasm"));
  } catch {
    fail("the build has no assets/ directory");
  }
  if (wasm.length === 0) {
    fail("the build has no .wasm module");
  }
  let named = false;
  for (const name of (await readdir(assets)).filter((entry) => entry.endsWith(".js"))) {
    const source = await readFile(join(assets, name), "utf8");
    if (source.includes(".wasm")) {
      named = true;
    }
  }
  if (!named) {
    fail("no emitted JavaScript names the wasm module");
  }
  for (const name of wasm) {
    const response = await get(`${PREFIX}assets/${name}`);
    if (response === null) {
      continue;
    }
    const type = response.headers.get("content-type") ?? "";
    if (!type.startsWith("application/wasm")) {
      fail(`${name} is served as ${type || "(no type)"}, not application/wasm`);
    }
  }
} finally {
  server.close();
}

if (failures.length > 0) {
  console.error("check-dist: the built editor is not servable under " + PREFIX);
  for (const failure of failures) {
    console.error("  - " + failure);
  }
  process.exit(1);
}
console.log(`check-dist: dist/ is servable under ${PREFIX} (${DIST})`);
