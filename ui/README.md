# The azoth editor

The flowsheet editor: palette, PFD, unit-operation windows, workbook and messages, over one
document. The engine is `crates/azoth-wasm`, a browser build of the same middleware every other
door binds; the page is a reading of the envelope each call answers, not a second model.

## Run it

Prerequisites: Rust with the `wasm32-unknown-unknown` target, `wasm-pack` 0.13.1, and Node 22.

```bash
npm ci
npm run dev        # builds the wasm, then starts Vite
```

The development page is the module door: the wasm runs in the tab. The other door is
`azoth serve`:

```bash
azoth serve --flowsheet specs/flowsheets/demo.toml --allow-origin http://localhost:5173
# then open http://localhost:5173/?serve=http://127.0.0.1:7373
```

The served door holds one document for the life of the process, so New, Open and Examples are
disabled there; Save still works because every envelope carries the document.

## Check it

```bash
npm run typecheck
npm test
npm run build
npm run check:dist   # serves dist/ under /Azoth/ and asks for the page's own assets
python tools/check_ui.py
```

`check:dist` is the deploy gate: it catches a root-absolute asset and a wasm served without
`application/wasm`, which a browser refuses to instantiate. `tools/check_ui.py` holds the three
front-end rules — only `src/wire/decode.ts` casts, no measurement literals, and no component
holding part of the document.

The committed fixtures in `test/fixtures/` are the library's output. Regenerate them with
`python tools/gen_ui_fixtures.py`, and commit them with any change that moves the catalogue,
the palette or a shipped flowsheet; `test/currency.test.ts` recomputes them from the built wasm
and fails when they are stale.

## Publish it

`.github/workflows/pages.yml` builds the wasm and `dist/` on every push to `main` that reaches
the editor or the engine, runs the tests and `check:dist`, and deploys to GitHub Pages. The build
uses a relative `base`, so `dist/` is servable under `/Azoth/` or at a site root. Pages must be
enabled for the repository with **Source: GitHub Actions**.
