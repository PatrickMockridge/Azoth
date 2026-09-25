import { defineConfig } from "vite";
import react from "@vitejs/plugin-react";

export default defineConfig({
  plugins: [react()],
  // The wasm module is built by wasm-pack into `src/wasm/pkg`, and `vite` bundles the glue's
  // `new URL("azoth_wasm_bg.wasm", import.meta.url)` as an asset - which is why the target is
  // `web` and not `bundler`: one file to serve, no plugin, and no second loader.
  optimizeDeps: { exclude: ["./src/wasm/pkg/azoth_wasm.js"] },
  test: {
    // Node by default, because most of these tests are pure functions over JSON. The one that
    // renders the app carries `// @vitest-environment jsdom` at the top of its own file, so the
    // environment a file needs is visible in the file rather than inferred from a glob.
    environment: "node",
    include: ["test/**/*.test.{ts,tsx}"],
  },
});
