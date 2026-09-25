# azoth

Chemical-engineering calculations — equations of state, flashes, property models, and the
hydraulics and heat transfer around them — in **Python for the ecosystem and Rust for the
engine**, with an agentic layer on top. The algorithms are ported from
[NeqSim](https://github.com/equinor/neqsim) under Apache-2.0 and credited in
[`NOTICE`](NOTICE).

## Why azoth

- **Two implementations, mirrored.** Every calculation is written twice by hand — one Rust
  kernel, one Python kernel — from the same declaration, and neither is a wrapper over the
  other. What differs between them is what each is good for, not which one computes: Rust is
  compile-time safe, expresses the process calculus natively and composes into large or many
  concurrent simulations that run in minutes where Python takes days, while Python is where
  the ecosystem is — notebooks, data science, AI SDKs. Both are complete, so the two are
  compared case by case, and a disagreement is a finding.
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

## Inherent safety

The worst thing this library can do is not crash. It is **return a wrong number that looks
reasonable** — a value plausible enough to go into a design, with nothing about it saying
it was never computed from anything. That is the hazard, and the design removes it rather
than adding a layer to contain it, which is the preference process safety has for a
chemistry that cannot run away over a relief valve that catches it.

- **Failures are loud, and they name the input.** Every error is typed and carries the
  offending field, because "calculation failed" is not actionable. A two-phase stream has
  no one density, so asking for one is refused rather than answered with phase 0's. A
  flash that converged onto the feed reports `beta` absent rather than zero.
- **A partial function is total, or refused.** `sqrt`, `ln`, `powf` with a non-integral
  exponent and `/` are each undefined somewhere, and `f64` returns a number there anyway.
  A kernel that evaluates one without establishing its domain fails *silently*, so it is
  not allowed to. `NaN`, zero and `null` are values with meanings here, not gaps.
- **Nothing is completed from a similar substance.** A name the data cannot answer is an
  error rather than a default, and a coefficient nobody chose is an error rather than a
  plausible number. Inventing data is the one thing a calculation may not do.
- **A keycard is refused at the door.** It declares which data its holder is entitled to
  use; a section, parameter or unit the reader does not know is refused when the card is
  *loaded* rather than when it is finally used, because a value nothing reads is data that
  looks in use and is not. The type is sealed — its parse state is private, with no default
  constructor — so a card cannot exist without having passed every refusal.
- **Two implementations, mirrored, and compared.** Because both are complete, neither can
  cover for the other, so the cross-implementation tests run with `AZOTH_REQUIRE_RUST=1`:
  a missing Rust core is a hard failure rather than a quiet fallback, because a suite that
  stays green on the Python path while the Rust one is absent is a guarantee that has
  evaporated.
- **A formal layer underneath, and a gate on it.** A quantity's dimension and the
  arithmetic over it are meant to be a proved construction, so a dimension that does not
  compose fails rather than converting. What reaches Lean today is the dimension group, the
  37-unit vocabulary, the fractional-power conventions, the raw-to-normalised derivative,
  the first-order sensitivity of a solution and the keycard's non-amplification — each
  behind an axiom gate that fails the build on a `sorry`, including one inside a vendored
  dependency. **The physics kernels are measured, not proved**, and
  [the chapter](docs/src/safety.md) states the boundary in full.

Porting against an oracle found this class of defect in *the oracle*, which is the argument
for the checks rather than for the port:

| Found in NeqSim | |
|---|---|
| `getHID` multiplies the ideal-gas enthalpy of formation by zero, so the term is absent from both sides of every comparison | [#3991](https://github.com/equinor/neqsim/issues/3991) |
| a zero-filled vapour-pressure row is given a synthetic curve from its normal boiling point, in the wrong unit | [#3994](https://github.com/equinor/neqsim/issues/3994) |
| `performGibbsMinimization` performs no minimization, and logs that it completed | [#3993](https://github.com/equinor/neqsim/issues/3993) |
| the parachor dispatch tests exact class names, so `ComponentSrkCPAMM` reads the plain `PARACHOR` column | [#3990](https://github.com/equinor/neqsim/issues/3990) |

Each is filed upstream and found by azoth's own tests rather than by trusting provenance.
The rule that keeps that honest is that a divergence is a finding to *explain*, starting
from the assumption that azoth is the one that is wrong — and one of these turned out not
to reproduce in situ, with the correction on the issue.

## Standalone

The Rust implementation has no Python dependency. The keycard is a value a caller passes
rather than a global in force, so a Rust-native caller supplies their own authority instead
of inheriting a Python session's, and the unit-operation and executor layers compose into
large or many concurrent simulations. There is a CLI. Being a standalone simulator is the aim
rather than the status: port coverage is incomplete, everything ported ships `unverified`,
and the book's front door says [not for design work
yet](docs/src/index.md#not-for-design-work-yet) about the placeholder fitting coefficients.

## Install

```bash
pip install azoth-engine   # the distribution; `import azoth` is unchanged
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

The full catalog — 67 calculations and 117 models, each with its equation, source, valid
range and a worked example — is in [the book](docs/src/index.md), generated from the same
spec files as the code.

## Entry points

- **Use the library** — [the book front door](docs/src/index.md).
- **Build on the agentic layer** — [`CLAUDE.md`](CLAUDE.md) and [`AGENTS.md`](AGENTS.md),
  and [`agents/README.md`](agents/README.md) for the HAZOP team and its runtime.

## Licence

Code is **AGPL-3.0-or-later**; documentation and data are **CC-BY-4.0**. See
[`LICENSE`](LICENSE) and [`LICENSE-CC-BY-4.0`](LICENSE-CC-BY-4.0).
