// The wasm module, driven the way a browser drives it.
//
//     wasm-pack build crates/azoth-wasm --target nodejs --out-dir pkg --release \
//         --no-pack --no-typescript
//     node crates/azoth-wasm/test/smoke.mjs crates/azoth-wasm/pkg/azoth_wasm.js
//
// **What this covers that the Rust tests cannot.** Every function in the crate is a call into the
// wire layer and a string out, so what is untested there is the *boundary*: that a wasm-bindgen
// constructor works, that a `Result` becomes a JS exception rather than a trap, and that the
// module carries its own palette with no filesystem and no fetch. This is that, and nothing else.
//
// It reads the shipped document, so it is run from the repository root.

import { readFileSync } from "node:fs";
import { pathToFileURL } from "node:url";

const modulePath = process.argv[2];
if (!modulePath) {
  console.error("usage: node smoke.mjs <path to azoth_wasm.js>");
  process.exit(2);
}
const azoth = await import(pathToFileURL(modulePath).href);

const failures = [];
const check = (what, condition, detail = "") => {
  if (!condition) failures.push(`${what}${detail ? `: ${detail}` : ""}`);
};

// 1. The palette the module carries, with no filesystem anywhere in sight.
const catalogue = JSON.parse(azoth.palette_json(false));
check("29 palette entries", catalogue.unit_ops.length === 29, `${catalogue.unit_ops.length}`);
check(
  "27 with a model",
  catalogue.unit_ops.filter((entry) => entry.model !== null).length === 27,
);
check(
  "26 runnable",
  catalogue.unit_ops.filter((entry) => entry.runnable).length === 26,
);
const withTools = JSON.parse(azoth.palette_json(true));
check("15 tools", withTools.tools.length === 15, `${withTools.tools.length}`);

// 2. Opening the shipped document.
const document = readFileSync("specs/flowsheets/demo.toml", "utf8");
const editor = new azoth.Editor(document);
check("a fresh document is dirty", editor.dirty === true);
check("and sound", editor.ok === true);

// 3. A command, and the envelope it answers with.
const after = JSON.parse(
  editor.apply('{"command":"set_position","node":"instance:sep1","x":1234.0,"y":56.0}'),
);
check("the envelope says OK", after.ok === true, JSON.stringify(after.diagnostics));
check("6 nodes", after.flowsheet.graph.nodes.length === 6);
check("6 edges", after.flowsheet.graph.edges.length === 6);
const sep1 = after.flowsheet.graph.nodes.find((node) => node.id === "instance:sep1");
check("the position landed", sep1.position.x === 1234.0 && sep1.position.y === 56.0);
check("the edit is in the document", after.flowsheet.document.includes("[layout.instances]"));
check("and the values are stale", after.dirty === true);
check("with no session", after.session === null);

// 4. A run, and a value by path.
const ran = JSON.parse(editor.run());
check("it converged", ran.session.converged === true);
check("in two passes", ran.session.iterations === 2, `${ran.session.iterations}`);
check("with the tear", ran.session.tears[0].stream === "recycle_1");
check("and 43 values", ran.paths.length === 43, `${ran.paths.length}`);
check("a value by path", editor.value("p1.outlet.P") === 2000000.0);
check("the values are current now", editor.dirty === false);

// 5. A command that is not one, and a document that is not one: exceptions, not traps.
let refused = false;
try {
  editor.apply('{"command":"delete_everything"}');
} catch (error) {
  refused = String(error).includes("delete_everything");
}
check("an unknown command raises with the name in it", refused);

refused = false;
try {
  new azoth.Editor("this is not TOML");
} catch (error) {
  refused = String(error).length > 0;
}
check("a document that does not read raises", refused);

// 6. The execution order, which is a property of the session and not of a widget.
const reordered = JSON.parse(editor.set_order("topological"));
check("the order is in the envelope", reordered.execution_order === "topological");
check("and setting it is not a document change", reordered.dirty === false);
check("the document is untouched", reordered.flowsheet.document.includes("flowsheets.demo"));

let refusedOrder = false;
try {
  editor.set_order("kahn");
} catch (error) {
  refusedOrder = String(error).includes("kahn") && String(error).includes("insertion");
}
check("a name that is neither order raises with both in it", refusedOrder);

// 7. And the module is still alive after every refusal - a trap would have taken it with them.
check("the module survives its refusals", editor.value("p1.outlet.T") > 300.0);

if (failures.length > 0) {
  console.error(`wasm smoke: ${failures.length} failure(s)`);
  for (const failure of failures) console.error(`  ${failure}`);
  process.exit(1);
}
console.log("wasm smoke: OK (29 palette entries, 15 tools, 6 nodes, 43 values)");
