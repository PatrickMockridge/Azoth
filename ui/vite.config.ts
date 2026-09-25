import { defineConfig } from "vite";
import react from "@vitejs/plugin-react";

export default defineConfig({
  plugins: [react()],
  // The wasm module is built by wasm-pack into `src/wasm/pkg`, and `vite` bundles the glue's
  // `new URL("azoth_wasm_bg.wasm", import.meta.url)` as an asset - which is why the target is
  // `web` and not `bundler`: one file to serve, no plugin, and no second loader.
  optimizeDeps: { exclude: ["./src/wasm/pkg/azoth_wasm.js"] },
  test: { environment: "node", include: ["test/**/*.test.ts"] },
});
