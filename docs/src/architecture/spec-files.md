# Spec files

<!-- Hand-written, unlike the pages under the namespaces, which are generated from
     the specs. It is listed in SUMMARY.md by tools/gen_docs.py, and
     tools/check_links.py fails the build if a page under docs/src is missing
     from the summary. -->

**A file under `specs/` is a TOML data sheet.** It declares what a calculation or a model
*is*: what it is called, what goes in, what comes out, the parameters the solver runs
with, the bounds a caller is warned about, and the cases that pin it. Everything on
this page is about that format.

TOML is the repository's one data format, read by `tomllib` in the standard library on
the Python side and by `toml` + `serde` on the Rust side. A number is a number in every
spelling of its exponent, so a spec cannot state a float that arrives as a string.

A spec file is **code**, not documentation. `tools/gen_registry.py` and
`tools/gen_models.py` compile it into the Rust tables and the Python modules both
languages read at run time; `tools/gen_docs.py` and `tools/gen_stub.py` render the
reference pages and the type stubs from the same source. Editing one without
regenerating fails the build.

The *specification* of what a calculation means — its contract, its domain, why a
bound takes the value it does, what has not been confirmed — is a document, and it
lives in this book. [The keycard is where responsibility
sits](./specification.md#the-keycard-is-where-responsibility-sits) is where that boundary is
stated: **the library implements; the engineer decides.**

## The rule

> A spec file declares. It does not explain, justify, narrate, or record status.

Three consequences, and they are the reason this page exists:

- **No prose fields.** A spec carries values and short labels. `tools/spec_lint.py`
  enforces the limit: 300 characters for any string in a spec, and 2,000 for a worked
  substitution, which is arithmetic a reader retraces rather than prose. A sentence
  that needs a paragraph is a page in this book, not a spec field.
- **No status.** There is no `verified`, `unverified` or `source_needed` field on a
  spec, and there must not be. A spec declaring whether it is verified would be the
  library grading an engineer's judgement. What a calculation has not had confirmed is
  stated in its own `assumptions` and `source`, next to the thing it is about.
- **No registration.** Nothing lists what exists. A calculation's id *is* its address —
  the module path and function name follow from it — so adding one adds files and
  registers nothing.

`additionalProperties: false` is set at every level of both schemas. That is the
mechanism: a status, a citation-status or a provenance field cannot appear in a spec
without a deliberate edit to `specs/schema/`, which is a change a reviewer sees.

## The two kinds

| | `specs/calcs/**` | `specs/models/**` |
|---|---|---|
| Schema | `calc.schema.json` | `model.schema.json` |
| Is | one equation, scalars in, a result out | a procedure, or a computation over vectors |
| Cases are | `tests` | `cases` |
| Solver | optional `solver` block | `algorithm`, required unless `kind: direct` |
| Generated into | `spec_gen.rs`, `_registry_gen.py` | `model_gen.rs`, `_models_gen.py` |

## Fields

Types are as declared in the schemas, which are the authority if this table and they
disagree. "Reads" names what consumes the field.

### Identity

| Field | Type | Req | Read by |
|---|---|---|---|
| `id` | string, `^[a-z][a-z0-9_]*(\.[a-z][a-z0-9_]*)+$` | yes | everything; it is the address |
| `name` | string | yes | the page title, the registry |
| `kind` | enum `procedure` \| `direct` | models | `gen_models`; selects whether `algorithm` is required |
| `equation` | string | calcs | `gen_docs` (as code), the contract tests |
| `latex` | string | calcs | `gen_docs`, via KaTeX |

### Interface

| Field | Type | Req | Read by |
|---|---|---|---|
| `inputs` | map of name → quantity | yes | signatures, the stub, `gen_docs` tables, range-check `on_input` |
| `outputs` | map of name → quantity or vector or matrix | yes | the result types, `gen_docs` tables |

A quantity carries `type`, `unit`, `description`, and optionally `optional`,
`default`, `interval` and, for vectors and matrices, `length` or `shape`. `interval:
true` marks a temperature *difference* rather than an absolute temperature; a flag the
spec holds and the call site ignores is a check that exists on paper only, so
`input_to_si` reads it.

### Behaviour

| Field | Type | Req | Read by |
|---|---|---|---|
| `algorithm` | scheme, convergence, tolerance, max_iterations, initialisation, `bracket`, `inner` | models, unless `direct` | `gen_models`, both implementations |
| `solver` | kind, tolerance, max_iterations, convergence, initial_guess | calcs, optional | `gen_registry`, `gen_docs` |
| `valid_range` | list of range checks | no | the runtime, via `azoth_core::spec` |
| `assumptions` | list of short strings | no | `gen_docs`; each is one sentence |
| `data` | map of name → path | no | `gen_docs` |

A range check carries `quantity`, `severity` (`warning` or `error`), at least one of
`min`/`max`/`equals`/`enum`, the `min_inclusive`/`max_inclusive` flags, and optionally
`when`, `code` and `computed_from`. Its `rationale` is **rendered into the runtime
warning message**, so it is part of the format rather than commentary: one sentence,
and it says why the bound exists, not how it was decided.

**Errors are typed; range violations are warnings.** A value outside its validated range
is still a value — the library reports it and does not refuse it, unless the check
declares `severity: error`.

### Verification

| Field | Type | Req | Read by |
|---|---|---|---|
| `tests` | list of tests | calcs | `gen_registry`, the test files |
| `cases` | list of cases | models | `gen_models`, the test files |
| `worked_example` | source, inputs, expected, tolerance, `derivation` | calcs | `gen_registry`, `gen_docs` |

A test carries `id`, `type` (`worked_example`, `reference`, `property`), `status`
(`active` or `skipped`), and — when skipped — `skip_reason`. A skipped test is allowed;
skipping one silently is not.

A **worked example must be retraceable by hand**, and `derivation` is where the
substitution is written out. It is the one long string the format keeps, because a
number nobody can follow is a number somebody typed.

### Attribution

| Field | Type | Req | Read by |
|---|---|---|---|
| `source` | `standard` (required if present), `edition`, `equation`, `doi`, `url` | no | `gen_docs`, the `## Source` section |
| `references` | list of strings | no | `gen_docs` |

**A citation, and nothing beyond it.** A bibliographic reference, or — for a port — the
upstream class and the line ranges taken and not taken. That second form is a licence
obligation and a citation rather than an argument, which is why the line numbers stay
here. Attribution for anything ported lives in
[`NOTICE`](https://github.com/PatrickMockridge/Azoth/blob/main/NOTICE), once, rather
than in a per-calc block — a block repeated in twenty specs is a block nobody reads.

**What a port *changed*, and why, belongs in the spec's `assumptions`**, not in a
citation field. The line between the two is whether the sentence identifies a source or
argues about it.

## Fields that are not part of the format

The schema rejects each of these, because `additionalProperties: false` is set at every
level — a spec carrying one fails to validate rather than being quietly ignored.

| Field | Where its content belongs instead |
|---|---|
| `notes` | Required improvements, for the parts that report a defect or a limitation. The rest is the code's job, or nothing |
| `description` | the calculation's page, or this book |
| `cases[].note`, `cases[].source`, `tests[].note` | Required improvements, or deleted |
| any `status`, citation-status or provenance field | nowhere — see [the keycard's
rules](./specification.md#the-rules-a-keycard-follows) |

**One field is only partly clean.** `source` is a citation, and for a port it keeps the
upstream class and the line ranges taken and not taken, because that is what a licence asks
for. But several `source.edition` values *argue* about the source rather than naming it, and
an argument is a limitation report. That is the last of it, and the sentence to keep is the
one that identifies the source.

The reason the removed fields could not stay is not taste. A spec is data, and prose that
convinces is prose that goes on convincing after the value beside it has changed. What a
spec declares is a value and a short label; an argument belongs in a page of this book,
where it is written as markdown in a markdown file.

## What a valid sheet looks like

```toml
id = "hydraulics.reynolds_number"
name = "Reynolds number"
equation = "rho * v * d / mu"
latex = "Re = \\frac{\\rho v d}{\\mu}"

[inputs.rho]
unit = "kg/m**3"
description = "density"

[inputs.v]
unit = "m/s"
description = "bulk velocity"

[inputs.d]
unit = "m"
description = "internal diameter"

[inputs.mu]
unit = "Pa*s"
description = "dynamic viscosity"

[outputs.Re]
unit = "dimensionless"
description = "the Reynolds number"

[[valid_range]]
quantity = "d"
min = 0.0
min_inclusive = false
severity = "error"
rationale = "a diameter of zero has no flow in it"

[implementations]
python = "azoth.hydraulics.reference.reynolds_number"
rust = "azoth_hydraulics::reynolds_number"

[worked_example]
source = "retraced by hand from the definition"
inputs = {rho = 998.0, v = 2.0, d = 0.05, mu = 1.002e-3}
expected = {Re = 99600.79840319361}
tolerance = 1.0e-9

[[tests]]
id = "the_worked_example"
type = "worked_example"
status = "active"
inputs = {rho = 998.0, v = 2.0, d = 0.05, mu = 1.002e-3}
expected = {Re = 99600.79840319361}
tolerance = 1.0e-9
```

## How it is checked

| Gate | What it catches |
|---|---|
| `tools/spec_lint.py` | schema conformance, the semantic rules a schema cannot express, and every string's length and shape |
| `tools/check_json_keys.py` | a duplicated key, which `json.loads` silently resolves to the last one |
| `tools/check_numerics.py` | an integral exponent written as a float, and a decimal exponent that approximates a rational |
| `test_registry_contract.py` | spec ↔ registry ↔ code agreement, in both languages |
| `test_model_contract.py` | the same for models, plus the algorithm-scheme vocabulary |
| the `docs-drift` job | a spec edited without regenerating its output |
