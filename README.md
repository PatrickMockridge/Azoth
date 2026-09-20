# azoth

Chemical-engineering calculations — equations of state, flashes, property models, and the
hydraulics and heat transfer around them — in **Python for the ecosystem and Rust for the
engine**, with an agentic layer on top. The algorithms are ported from
[NeqSim](https://github.com/equinor/neqsim) under Apache-2.0 and credited in
[`NOTICE`](NOTICE).

## Why azoth

- **Two languages, two jobs.** Python is the surface: basic calculations, notebooks,
  data science, AI SDKs. Rust is the engine: compile-time safe, it expresses the process
  calculus natively and composes into large or many concurrent simulations that run in
  minutes where Python takes days.
- **A spec file is the source.** One TOML per calculation declares its inputs, outputs,
  ranges, assumptions and tests; generators compile it into Rust, Python, the type stubs,
  the docs and the tests, and CI fails on any drift.
- **Provenance is split three ways.** A validation case is `verified`, `unverified` or
  `source_needed`; a shipped value carries a `verify_status`; a spec refuses to carry a
  status field at all.
- **NeqSim is what azoth is a port of, so a disagreement with it is the interesting
  failure.** Its numbers are committed beside the Java driver that prints them rather than
  copied on trust, and a divergence is a finding to *explain* - starting from the
  assumption that azoth is the one that is wrong. "A divergence is a finding, never a
  failure" reads the other way round, and that reading cost a session: a probe that
  subtracted the wrong quantity reported a NeqSim defect that was azoth's own.
- **An agentic layer.** Ninety skills under `skills/` teach an agent to call the library,
  and a four-role HAZOP team chains them.
- **A formal layer in Lean.** The dimension group is proved (`lean/Azoth/Dim.lean`); the
  33-unit vocabulary compiles into theorems (`lean/Azoth/Vocabulary.lean`).

## Install

```bash
pip install azoth          # the Python package, with the compiled Rust core
```

## Quick start

```python
import azoth

q = azoth.ureg.Quantity

r = azoth.hydraulics.darcy_weisbach(
    0.02, q(100.0, "m"), q(0.1, "m"), q(998.0, "kg/m**3"), q(1.5, "m/s")
)
r.dp  # 22455.0 pascal
```

Units cross the API as `pint` quantities; a bare number where a length is expected raises
`UnitMismatchError`.

## The book

The full catalog — 61 calculations and 65 models, each with its equation, source, valid
range and a worked example — is in [the book](docs/src/index.md), generated from the same
spec files as the code.

## Entry points

- **Use the library** — [the book front door](docs/src/index.md).
- **Build on the agentic layer** — [`CLAUDE.md`](CLAUDE.md) and [`AGENTS.md`](AGENTS.md),
  and [`agents/README.md`](agents/README.md) for the HAZOP team and its runtime.

## Licence

Code is **AGPL-3.0-or-later**; documentation and data are **CC-BY-4.0**. See
[`LICENSE`](LICENSE) and [`LICENSE-CC-BY-4.0`](LICENSE-CC-BY-4.0).
