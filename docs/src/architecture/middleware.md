# The middleware

This page states what sits between the [process schema](./process-schema.md) and its
front-ends — a flowsheet editor (web front-end, [the stack](#the-stack)), a notebook
(Jupyter/conda), and an agent
— and it states it as a use of [reflection](../calculus/rho.md), not a parallel layer. A
flowsheet is a process; quoting it and dropping it back is the identity (`*@P ≅ P`, proved
as structural congruence), and that round trip is what an editor, a file and a notebook
all do to the same data. The middleware is that operator made concrete.

It is **not built**. It is a requirement, deferred until the backend below it closes
through P12, and written here so it is gated on the roadmap rather than lost.

## The gap

| a front-end needs | exists |
|---|---|
| run a flowsheet | no executor — schema, checker and 26 kernels (26 of the 29 palette entries) |
| a session holding named results | no solve state; the runtime `Stream` has no name and is not serialisable |
| read *and* write a flowsheet | TOML read only; no `toml::to_string` |
| a form per unit op (ports + parameters) | the spec exists; nothing turns it into a form |
| structured, locatable diagnostics | `Debug` strings |
| dispatch an instance to its kernel | no `unit_ops.*` → kernel table |
| a result as JSON | result dataclasses are not JSON-serialisable |
| a tool schema for an agent | skills and one runtime stub; no tool schema, no MCP |

## The layers

Each is a consequence of the calculus, named here because a front-end has to address it.

- **session/document** — the live flowsheet plus its named stream results and a dirty/clean
  flag. Every value has a stable path (`p1.outlet_pressure`, `recycle_1`), which is what a
  widget and an agent both point at.
- **command model** — every edit is a typed command (add/remove an instance, connect two
  ports, set a parameter) that re-runs the checker. This is
  [linearity](../calculus/process.md) made incremental: a red arrow is a structured
  `OverfedPort`, not a string.
- **quote/drop** — the round trip in both directions: `toml::to_string` for the flowsheet,
  a JSON codec for a stream result.
- **form schema** — one uniform object per unit op and calculation: its declared free
  variables with units, dimensions and valid ranges. This is the
  [adequacy claim](../calculus/process.md) made inspectable — a unit-op window *is* its port
  declaration.
- **structured diagnostics** — a severity and a location, not a `Debug` line.
- **dispatch** — the `unit_ops.*` id → kernel mapping, with typed parameter coercion (a
  palette parameter is a `toml::Value` today).

## The GUI

The widgets of a HYSYS/UniSim-style editor are these layers viewed one way:

| widget | calculus | fed by |
|---|---|---|
| flowsheet view | the quoted connection graph | session + command model |
| unit-op window | the adequacy declaration (ports + parameters) | form schema |
| input/output fields | the field record, one field name to a dimension | form schema + units |
| palette dropdown | the `UnitOpSpec` registry | the palette |
| top-bar menu | new / open / save / solve | command model + quote/drop |

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
process the same way). The agentic middleware is a **tool schema and a session**: JSON
Schema descriptions for the tools, derived from the form schema and the command model
(`list_unit_ops`, `run_calc`, `read_stream`, `set_parameter`, `add_instance`, `validate`,
`load`, `save`), and an in-process runner that executes them against the same session a
human edits. It is hosted twice — a Jupyter/Colab cell and a GUI side panel — and an MCP
server is a projection of the same schema, already [beyond P12](./specification.md).

## What it gates on

The middleware waits on the backend, in the order
[the specification's tranche table](./specification.md) fixes — the one data path first,
then transport properties, the physics, a kernel for every unit op, and the executor.
That order is stated there and is not repeated here.
