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
npm run check:dist   # serves dist/ under /azoth/ and asks for the page's own assets
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

## Appearance

`data-layout` on `<html>` is `xp`, `mac` or `cyberpunk`; `data-theme` is `light`, `neutral` or
`dark`. `:root` is Windows XP light, so a document with no attributes already has the default and
the pre-paint script only overrides it. `?layout=mac&theme=dark` previews a combination without
saving it.

**The layouts are shells, not skins.** Each one arranges the same semantic regions — the
Navigator, document, Inspector, Workbook and status — differently:

- **Macintosh** is a desktop. The menu bar is a screen strip and the palette, document,
  inspector, workbook and status are separate windows. `useWindowManager` drags from a title,
  resizes from the grip, raises on click, clamps to the desktop and exposes a **Reset window
  positions** command. The resting geometry is CSS.
- **Windows XP** is one MDI frame: a Luna title bar, a menu row, a tool row, a large-button
  toolbar, Explorer panes and a white document child under a blue caption.
- **Cyberpunk** is a full-bleed drawing with translucent HUD rails and a command deck. The
  telemetry strip's rail switches hide the Navigator, Inspector, Workbook or Toolbar; the canvas
  never hides. `Ctrl/Cmd+K` opens the command palette.

**Typography is bundled**, not borrowed from the OS: Sysfont for Macintosh chrome, the unmodified
Wine Tahoma for XP, Rajdhani for Cyberpunk chrome, with IBM Plex Sans for the readable body and
IBM Plex Mono for data. `src/styles/fonts.css` and `src/fonts/` carry the faces and licences;
`check:dist` asks for every font URL and its MIME type.

**One command model** (`src/commands/model.ts`) feeds the toolbar, the menu bar and the palette.
The toolbar buttons stay visible and keep their accessible names; menu items are `menuitem`s and
palette rows are `option`s, so a second surface never duplicates a `getByRole("button")`.

`src/styles/tokens.css` owns every value: the nine tone cells, the three layout shape blocks and
the derived paints. `src/styles/layouts/*.css` carry the boxes, captions and windows but no
colour or duration literal; `test/tokens.test.ts` enforces the contract, the nine cells and the
shell invariants. `test/contrast.test.ts` measures the text pairs; `test/shell.test.tsx` moves,
resizes and resets a Macintosh window. `src/state/appearance.ts` reads and writes
`azoth.appearance`; `index.html` carries a pre-paint copy of the same key and defaults.

To add a layout: add its id to the union and the chooser, add the full shape set in a
`[data-layout]` block plus one tone cell per colour in `tokens.css`, and add a
`layouts/<id>.css` imported before `report.css`. The token and contrast tests name what is
missing.

**Results are two layers, from one envelope.** The model spec's `[outputs]` block (unit, kind,
description and `optional`) is generated into `model_outputs_gen.rs` and carried on each
unit-op form, so `state/results.ts`'s `outputRowsOf` can label a value, convert a bare float from
its declared unit and show the library's own sentence — often the equation behind the number —
without a hand-written table. Only fields the run publishes are rows: a `null`, an `undefined` or
an empty array is absent rather than a dash, while `false` and `0` stay. `property/ResultsSheet`
and the report's cover are the human reading; the sheet's *Machine-readable fields* block, the
Provenance tab and the report's appendix keep the model's field names, library units and hashes.
`report/runRecord.ts` builds the downloadable JSON record from the envelope — document, raw
session, provenance, diagnostics — and `state/verification.ts` is the one place the four model
verification statuses are explained; they are provenance, not a badge on an answer.

**A layout hides decoration, never a control.** `test/layoutVisibility.test.ts` scans the layout
stylesheets and refuses a `display: none` that names a switch, splitter, tab, menu trigger or any
other control; the Cyberpunk workbook's splitter and the Macintosh workspace switch are the two
regressions the rule exists for. The Macintosh title bar draws no close or zoom box, because a
control that does nothing is worse than none.

**Connectors are routed, not guessed.** `src/pfd/route.ts` treats each node plate as an obstacle
and finds an orthogonal path on a visibility grid; a recycle gets a return lane below the
drawing, one lane per recycle. `src/components/PfdEdge.tsx` draws that path and keeps the old
temperature/pressure tag. `test/pfdRoute.test.ts` checks that no segment crosses a plate, the
return uses its lane, and no line is retraced in the opposite direction.

## Publish it

`.github/workflows/pages.yml` builds the wasm and `dist/` on every push to `main` that reaches
the editor or the engine, runs the tests and `check:dist`, and deploys to GitHub Pages. The build
uses a relative `base`, so `dist/` is servable under `/azoth/` or at a site root. Pages must be
enabled for the repository with **Source: GitHub Actions**.
