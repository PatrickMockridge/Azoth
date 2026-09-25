import { fileURLToPath } from "node:url";
import { resolve } from "node:path";

import react from "@vitejs/plugin-react";
// From `vitest/config`, not `vite`: the `test` block is vitest's, and `vite`'s own
// `defineConfig` types reject it. Vitest 4 removed the lenient overload that accepted it.
import { defineConfig } from "vitest/config";

// The editor opens the same `specs/flowsheets/demo.toml` the CLI runs and the Rust tests
// hold to, which means `src/App.tsx` imports a file from *outside* this package. Vite 6
// tightened `server.fs` to the project root and now refuses that read, so the repository
// root is allowed explicitly. It is the one grant this front-end needs, and it is stated
// rather than widened to `/`.
const packageRoot = fileURLToPath(new URL(".", import.meta.url));
const repositoryRoot = resolve(packageRoot, "..");

export default defineConfig({
  plugins: [react()],
  // The wasm module is built by wasm-pack into `src/wasm/pkg`, and `vite` bundles the glue's
  // `new URL("azoth_wasm_bg.wasm", import.meta.url)` as an asset - which is why the target is
  // `web` and not `bundler`: one file to serve, no plugin, and no second loader.
  optimizeDeps: { exclude: ["./src/wasm/pkg/azoth_wasm.js"] },
  server: { fs: { allow: [repositoryRoot] } },
  test: {
    // Node by default, because most of these tests are pure functions over JSON. The one that
    // renders the app carries `// @vitest-environment jsdom` at the top of its own file, so the
    // environment a file needs is visible in the file rather than inferred from a glob.
    environment: "node",
    include: ["test/**/*.test.{ts,tsx}"],
  },
});
