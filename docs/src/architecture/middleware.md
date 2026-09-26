# The middleware

This page states what sits between the [process schema](./process-schema.md) and its
front-ends — a flowsheet editor (web front-end, [the stack](#the-stack)), a notebook
(Jupyter/conda), and an agent
— and it states it as a use of [reflection](../calculus/rho.md), not a parallel layer. A
flowsheet is a process; quoting it and dropping it back is the identity (`*@P ≅ P`, proved
as structural congruence), and that round trip is what an editor, a file and a notebook
all do to the same data. The middleware is that operator made concrete.

**It is built, and every front-end binds one surface rather than its own.** The layers are
`crates/azoth-process`'s `middleware` module; the bindings are `azoth-wasm` for a browser,
`azoth.process` for a notebook, `azoth edit`/`azoth forms` for a shell and `azoth mcp` for an agent
— each a thin wrapper over the same calls, each answering with the same document. What is not built
is named at the foot of this page.

## The gap

| a front-end needs | exists |
|---|---|
| run a flowsheet | **yes** — `azoth_process::executor`, with the dispatch table, the tear's fixed point and a NeqSim `ProcessSystem` capture behind it; a flowsheet declares its own inputs, so it runs from the library, `azoth run --flowsheet` and `azoth.process.run_flowsheet` with no argument but the file |
| a session holding named results | **yes** — `middleware::session::Workspace`, where every value is `<endpoint>.<field>`, `paths()` enumerates them, and `dirty` says whether they are older than the document |
| read *and* write a flowsheet | **yes** — `Flowsheet::to_toml`, with a test over `specs/flowsheets/` that the value survives and that writing is a fixed point of itself |
| a form per unit op (ports + parameters) | **yes** — `middleware::form`, one object per palette entry: the ports, the parameters with their units and dimensions, the *kind* that chooses a widget (from the model's input declaration, via `model_inputs_gen`), and the model's own bounds with the sentence that explains each. The bounds the model states on inputs no parameter carries are listed rather than dropped |
| structured, locatable diagnostics | **yes** — every `Diagnostic` carries a `Severity`, a `Location` **and a `Target`**: a stable code and the node, handle, parameter or edge a front-end puts a mark on, both derived from the variant's own fields |
| dispatch an instance to its kernel | **yes** — the `DISPATCH` table, 27 entries and 2 refusals by name |
| an edit as one typed command | **yes** — `middleware::command`, fifteen variants, each a structure edit that leaves the checker to rule on the result |
| a graph a canvas draws | **yes** — `middleware::graph`, in xyflow's own node/edge shape, with the layout in one `[layout]` table the document carries |
| the whole of it as one document | **yes** — `middleware::envelope`, which is what every binding answers with |
| a result as JSON | **yes** — a *session's* result writes JSON from Rust (`executor::json`, which the envelope embeds) and a *calculation's* writes it from Python (`azoth.core.serialise`). That second one is one walk over `dataclasses.fields` which every result inherits through `_HasWarnings`, so **there is no codec per model** — a result is serialisable by being a frozen dataclass. The two writers name the *unit* in the same key and the magnitude in different ones, deliberately: a stream record's five fields are SI by construction and a spec's unit need not be. **Every registered calculation and every registered model is run and its result serialised** by `python/tests/test_result_json.py` — the twenty-eight `process.*` models among them, so a `PumpResult` reads field by field with no codec of its own |
| a tool schema for an agent | **yes** — `middleware::tools`, one tool per command, projected from the command model rather than written beside it, and **served over MCP** by `azoth mcp`, which is that schema's projection and not a second one |

## The layers

Each is a consequence of the calculus, named here because a front-end has to address it.

- **session/document** — `Workspace`: the live flowsheet, the palette it is checked against, the
  last check, the last run, and a dirty flag. Every value has a stable path
  (`p1.outlet.P`, `recycle_1`), the path is `<endpoint>.<field>`, and the field is one the port
  declaration names — `n`, `z`, `P`, `T`, `h` — so a front-end derives the path from the
  declaration rather than inventing a spelling for it. **A check is cheap and happens on every
  edit; a run is the physics and happens when asked**, which is why the diagnostics are always
  current and the values are not.
- **command model** — every edit is a typed command (add/remove an instance, connect two ports,
  set a parameter, place a node) that re-runs the checker. This is
  [linearity](../calculus/process.md) made incremental: a red arrow is a structured `OverfedPort`,
  not a string. **The checker is the only rule set** — a command that names something absent is a
  no-op and the checker's own diagnostic is the answer, and a command is refused only where the
  *command* is malformed.
- **quote/drop** — the round trip in both directions: `toml::to_string` for the flowsheet, and the
  envelope's JSON for a session. The round trip is executed rather than claimed: every edit passes
  the document through `from_toml`/`to_toml`, and a test holds `from_toml(to_toml(f)) == f` and
  that writing is a fixed point of itself, per command.
- **form schema** — one uniform object per unit op: its declared free variables with units,
  dimensions, kinds and valid ranges. This is the [adequacy claim](../calculus/process.md) made
  inspectable — a unit-op window *is* its port declaration.
- **structured diagnostics** — a severity, a location and a target, not a `Debug` line.
- **dispatch** — the `unit_ops.*` id → kernel mapping, with typed parameter coercion.
- **the graph** — a flowsheet's connection graph, in xyflow's node/edge shape. A node's id is
  `{role}:{name}` and **a handle's id is the stream's own path**, so an edge's handle, a session
  key and an edge's `data.path` are one string rather than three that have to be reconciled.

## The bindings

One surface, five doors. None of them holds a second implementation, and a test on the host
covers each: the wire layer's tests are the behaviour, and a binding is exercised for the boundary
only.

| door | what it is | what covers it |
|---|---|---|
| a browser | `crates/azoth-wasm`, compiled with `wasm-bindgen` and carrying its own palette and databank — no fetch, no filesystem | a Node driver over the built module, in the crate |
| a notebook | `azoth.process.Session`, which is the same `Workspace` as a Python object, with `forms()` and `tools()` beside it | `python/tests/test_process.py` |
| a shell | `azoth forms` and `azoth edit --flowsheet F --command JSON [--run] [--json]` | `crates/azoth-cli/tests/wire.rs` |
| an agent | `azoth mcp --flowsheet F [--no-run]`: the same tools over stdio, one session held across calls, every answer the envelope — and the same messages again at `POST /mcp` on `azoth serve`, for a client that cannot start a process | `crates/azoth-cli/tests/mcp.rs`, which asserts the served tools equal `middleware::tools` field by field, and `crates/azoth-cli/tests/mcp_http.rs` for the HTTP lane |
| a client that cannot run the kernels | `azoth serve --flowsheet F [--port N] [--allow-origin ORIGIN]`: one document over HTTP, the same calls, the same envelope — it hosts the editor's `POST /rpc` and the agent's `POST /mcp` over **one** `Session` | `crates/azoth-cli/tests/serve.rs`, and `crates/azoth-cli/tests/mcp_http.rs` for the MCP endpoint |

`--json` prints the envelope, which is byte for byte the object a browser is handed — so the CLI
is a way to look at the wire without a front-end, and a way to capture a fixture for one.

**The network transports call the same code.** `azoth mcp`, `azoth serve`'s `/rpc` and its `/mcp`
all hold a `crate::session::Session` and hand it a tool name and its arguments; none knows what an
edit is, and what differs between them is how a refusal is *expressed* — a JSON-RPC error code, an
HTTP status, or a status the MCP transport's own specification fixes. **`/rpc` and `/mcp` are one
`Session` in one process**, so an edit through either is visible through the other and an agent and
a person are looking at the same document. A second dispatch table, a second place to read a
command, is what that shared session exists to make impossible.

## The GUI

`ui/` is the editor: vite, React and xyflow, a PFD between three docks and a status bar. The widgets
of a HYSYS/UniSim-style editor are these layers viewed one way:

| widget | calculus | fed by |
|---|---|---|
| flowsheet view | the quoted connection graph | session + command model |
| the symbols | the machine each palette entry declares, drawn from its own leaf | the palette |
| unit-op window | the adequacy declaration (ports + parameters) | form schema |
| its result sheets | the operation's registered result, beside its outlet streams | `session.results` |
| input/output fields | the field record, one field name to a dimension | form schema + units |
| workbook | every stream, by the endpoint that produced it | session + graph |
| messages | the checker's diagnostics | the checker |
| palette dropdown | the `UnitOpSpec` registry, grouped by the family its file sits in | the palette |
| top-bar menu | new / open / save / solve | command model + quote/drop |

**It holds no copy of the flowsheet.** Every gesture goes out as a command and comes back as the
whole envelope, so the canvas, the form and the diagnostics are three readings of one object. The
drag is the one piece of local state: xyflow applies a position change every frame and the library
is told once, on release — a `set_position` per frame would be a command and a re-projection for a
figure that has not moved yet. `ui/test/app.test.tsx` drags a node through the real module and
asserts the position the canvas draws is the position the document carries, which is the difference
between a gesture that sent a command and one that only moved a figure.

**A node is drawn as the machine it is.** `ui/src/state/glyphs.ts` maps each palette entry's leaf to
a symbol — a pump is a circle with a discharge triangle, a compressor and an expander are opposed
trapezoids, a column carries tray lines and a packed column a hatch — and to the side each declared
port leaves from, because a column's distillate leaves the top and its bottoms the bottom whatever
its direction says. Both tables are held to the catalogue in both directions by
`ui/test/glyphs.test.ts`, which is the test that caught four entries whose port names in the side
table did not exist: a handle on the wrong side of a vessel is a drawing of a different machine.

**The workbook's rows come from the document and its values from the run.** They are the *producing*
endpoints the graph declares — a unit op's `outlet` handle and the next one's `inlet` are one stream,
so only the producer gets a row — which is what makes the grid the same shape before and after Solve
and a stream the run did not reach a row of dashes rather than a missing row.

**What a run computed is published beside its streams.** A kernel that solved a column or integrated
a reactor reached numbers that are on no outlet stream, and `session.results` carries them per
instance under the field names that operation's *model* declares — the same names a case uses and
the same names the Python half publishes, since a result that differed between two doors would be
two answers to one question. The window's result sheets read those names and enumerate nothing, so a
field added to a model appears without a line of front-end code; a unit op whose whole answer is its
streams publishes no entry at all, which is a statement about the arithmetic rather than a gap.

**The editor is one front-end over two of those doors.** `ui/src/wire/session.ts` is the seam — a
palette and three calls, every one of them a promise, because a `fetch` cannot answer a synchronous
call and one interface both doors implement is the only version where a panel does not have to know
which it is behind. `?serve=http://127.0.0.1:4000` on the page's URL is the served door and no
parameter is the wasm module; the server has to be told the editor's origin
(`--allow-origin`), because a `Content-Type: application/json` POST is not a simple request and the
answer is unreadable without the CORS header it earns.

**What the served door cannot do is the server's shape, not an omission.** `azoth serve` holds one
document for the life of its process and the route has no call that replaces it, so New, Open and
Demo are refused there with the reason on the control, and a save still works because every envelope
carries the document. The palette is the exception in the other direction: a form per unit op is on
no envelope, so `/rpc` answers `{"catalogue": …}` as a *read* of the palette the process loaded —
enough for a palette panel and a unit-op window, and not an opening of anything.

### The stack

- **graph** — xyflow (formerly React Flow). A flowsheet's connection graph is already a
  node/edge set, so the editor's model is the schema's rather than a translation of it.
- **language** — JS/TS, the working choice: the front-end ecosystem is far more developed
  than Rust's. Rust/WASM takes the parts where performance or code safety warrants it, and
  is the longer-term direction for more of the front-end.
- **the kernels** — remain Rust. WASM is how the front-end reaches them, which is what keeps
  the front-end a consumer of the middleware rather than a second implementation of it.

## The agent

An agent's tool is `@calc` and its call is `*tool` — the same reflection surface as the
GUI, not a third one ([orchestration](../agentic/orchestration.md) states a team as a
process the same way). The agentic middleware is a **tool schema and a session**: JSON Schema
descriptions for the tools, derived from the command model (`list_unit_ops`, `run_calc`,
`read_stream`, `set_parameter`, `add_instance`, `validate`, `load`, `save` in the page's own
words — one tool per command, in the code), and an in-process runner that executes them against
the same session a human edits. It is hosted twice — a Jupyter/Colab cell and a GUI side panel —
and an MCP server is a projection of the same schema — `azoth mcp`, over stdio, for one document
at a time. **The same projection is served over HTTP** at `POST /mcp` on `azoth serve`, which
speaks the current protocol revision only: the handshake revisions are the stdio door's, one
process away, and an HTTP client that asks for one is told so rather than left guessing.

**Only edits are tools.** Every call returns the whole envelope, so reading a stream is reading the
answer to the call that changed it, and loading, saving and validating are the session's own
openings rather than edits to it. A read tool would return the same document as every other tool.

## What is not built

Named with the class that would close each, because a page that lists only what works is a page
that reads as finished.

- **MCP resources.** The tools are served; a *resource* is a second surface for a document every
  call already returns, with a URI grammar to invent and a cache a client may serve stale — which is
  what `dirty` exists to prevent. The palette is already `add_instance`'s `unit` enum. **The
  messages are checked against the specification now**: its JSON Schema for each revision is
  vendored under `python/tests/validation/vendor/mcp/`, and `test_mcp_conformance.py` holds every
  message both transports write to the definition for the revision that message belongs to — two
  schemas because the protocol has two lanes, since `2026-07-28` removed the handshake and
  `2025-11-25` has no discovery, and two transports because `POST /mcp` writes the same messages as
  stdio. What that cannot see is a *new* revision: the vendored file is the pin, so upstream moving
  is invisible until somebody re-fetches it.
- **Streaming, and the things that only exist to carry a stream.** `POST /mcp` answers with one
  JSON object, which the specification allows the server to choose; the other option is a
  request-scoped SSE stream, and its reason to exist is a `notifications/progress` while a long call
  runs. Every tool here is synchronous and emits none, so a stream would carry exactly one event, and
  cancellation-by-disconnect — which closing that stream *is* — has nothing to cancel.
  `subscriptions/listen`, whose whole response is an open stream, needs a notification source and
  there is none: no resources, no list changes. `x-mcp-header` annotations are a server's option and
  this one publishes none, so there is nothing to mirror into `Mcp-Param-*` and nothing to validate.
  The legacy lane over HTTP — a session id, a standalone GET stream, `Last-Event-ID` resumption — is
  deliberately not built either: this endpoint is modern-only and says so, and the handshake is at
  the stdio door.
- **Authentication, and who edited what.** `azoth serve` hosts one document over HTTP, and what it
  does not have is named in its own module: no identity, no conflict resolution beyond the order
  requests arrive in, no TLS. It binds `127.0.0.1` and refuses every origin it was not told to
  allow, which is the posture a single-user local tool can defend — anything past that is a
  deployment, and a deployment is not built. **The seam a tokenized deployment needs exists and is
  dormant**: `mcp_http::payment_refusal` is called on every request and refuses nothing, because
  this server has no token, no issuer and no quota — and `402 Payment Required` is in the status
  map beside it, because the transport must be able to answer *payment required* rather than *bad
  request*, which are different instructions to a client. A test drives both halves.
- **A published distribution, and the two things only a person can do.** The wheel and the sdist
  are built, signed, verified and — by the `pypi` job, on a `v*` tag — uploaded to the index under
  the name **`azoth-engine`**, which Trusted Publishing exchanges this workflow's OIDC token for.
  `import azoth` does not change: the distribution's name and the module's are allowed to differ,
  which is what makes the taken name `azoth` a non-problem rather than a rebranding. What is not
  done is the two things no CI job can do for you: a PyPI account with 2FA, and a *pending
  publisher* for `azoth-engine` naming this repository, `release.yml` and no environment. Until
  that exists, the job fails at the exchange — with PyPI's own message about the claim, not a
  silent no-op. The `sdist` and the wheel are also walked by `check_wheel_data.py` for what they
  must not carry (`ui/`, `specs/`, the built browser module), so "a distribution is not a copy of
  the repository" is a measurement rather than an assumption about maturin's defaults.

## What it gates on

The middleware waited on the backend, in the order
[the specification's tranche table](./specification.md) fixes — the one data path first,
then transport properties, the physics, a kernel for every unit op, and the executor.
That order is stated there and is not repeated here.
