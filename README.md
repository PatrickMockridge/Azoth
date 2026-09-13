# azoth

Open, validated, citable chemical engineering calculations.

> azoth is an open library of standard chemical engineering calculations where
> every one carries its equation, the standard it came from, the range over which
> it is validated, its assumptions, a worked example and automated tests. Python
> reference implementation, Rust performance core, PyO3 binding them. Docs are
> generated from the same machine-readable specs, so they cannot drift from the
> code.

**Status: early.** One vertical slice is implemented - the hydraulics kernel
through Darcy-Weisbach pressure drop. The fitting coefficients it uses are
**placeholders, not engineering data**. See
[Not for design work yet](#not-for-design-work-yet).

## Install

```bash
pip install azoth          # the Python package, with the compiled Rust core
```

From source:

```bash
git clone https://github.com/PatrickMockridge/Azoth
cd Azoth

uv venv && uv pip install maturin pytest ruff mypy
maturin develop            # builds the Rust extension into the venv

cargo test                 # Rust core
pytest                     # Python package
mdbook build docs          # documentation
```

## Quick start

```python
import azoth

q = azoth.ureg.Quantity

r = azoth.hydraulics.darcy_weisbach(
    0.02, q(100.0, "m"), q(0.1, "m"), q(998.0, "kg/m**3"), q(1.5, "m/s")
)

r.dp  # 22455.0 pascal
r.warnings  # (Warning(RANGE_CHECK_SKIPPED, ...),)
r.is_clean  # False
```

Units cross the API as `pint` quantities, so a bare number where a length is
expected is a `UnitMismatchError` rather than a silent thousand-fold error:

```python
>>> azoth.hydraulics.darcy_weisbach(0.02, 100.0, q(0.1, "m"), ...)
UnitMismatchError: `L` must be a quantity convertible to 'meter', got 'float 100.0'
```

### From the command line

```bash
azoth pipe --fluid water --flow 10 --diameter 0.1 --length 100 \
    --fittings "90_elbow,gate_valve_open"
```

```
  reynolds number      35233.6  (turbulent)
  friction factor      0.02391984  (colebrook, 13 iterations)
    straight pipe           1493.3479 Pa
    fittings                  56.7472 Pa   (3.7% of total)
    total                   1550.0951 Pa   (1.550095 kPa, 0.015501 bar)

  1 warning(s)
    [ESTIMATED_DATA] 2 of 2 fitting(s) use ESTIMATED DUMMY coefficients...
```

## The two ideas this library is built around

**Warnings are not errors.** A value outside the range in which a correlation was
validated is still a value. Refusing to return it would be less useful than
returning it with a warning - but the library must never return it *silently*.
Check `result.warnings`, or `result.is_clean` if you are willing to see every
caveat at once.

**A check that could not run is not a check that passed.** When an optional input
is missing, the range check that depends on it reports `RANGE_CHECK_SKIPPED`
rather than quietly succeeding. In `darcy_weisbach`, omitting the viscosity
leaves the flow regime unchecked, and the result says so.

## What is implemented

| Calc | What it does |
|---|---|
| `hydraulics.reynolds_number` | Reynolds number and flow regime |
| `hydraulics.friction_factor_colebrook` | Implicit Colebrook-White friction factor |
| `hydraulics.friction_factor_swamee_jain` | Explicit approximation to it |
| `hydraulics.friction_factor_haaland` | A second explicit approximation, fitted differently |
| `hydraulics.crane_k_factors` | Fitting losses by the equivalent-length method |
| `hydraulics.darcy_weisbach` | Pressure drop over a straight pipe |
| `hydraulics.pump_power` | Shaft power from flow, head and efficiency |
| `hydraulics.orifice_flow` | Flow through an orifice from its pressure difference |
| `hydraulics.control_valve_cv` | Liquid flow through a control valve |
| `hydraulics.choked_flow_area` | Throat area for a choked gas flow |
| `thermal.conduction_plane_wall` | Steady conduction through a slab |

Pipe *with* fittings is a composition of the last two hydraulics calcs, done by the
`azoth pipe` CLI rather than by a calc of its own, because the two losses use
different methods and adding them is a modelling decision worth seeing explicitly.

Relief valve *sizing* to a standard is not implemented. `hydraulics.choked_flow_area`
is the isentropic basis - the throat area a given choked mass flow needs - and the
de-rating coefficients a standard applies are the caller's to compose.

The calc ids are namespaced by **domain** (`hydraulics.*`, `thermal.*`), not by
project. They appear in provenance records and citations, so renaming the project
does not - and should not - invalidate them. Two namespaces exist deliberately:
the second is what proves the spec pipeline is domain-agnostic rather than shaped
around pipe flow, since it runs through the same specs, generators, tests and
documentation with no special case anywhere.

## Architecture

```
specs/calcs/**/*.yaml     the registry: one file per calculation
        │
        ├─► tools/gen_registry.py ─► crates/azoth-hydraulics/src/spec_gen.rs
        │                          └► python/src/azoth/_registry_gen.py
        ├─► tools/gen_docs.py     ─► docs/src/**
        └─► tools/provenance.py   ─► provenance.json
```

Python is the reference implementation and Rust the core, with PyO3 binding them.
Units cross the public API as `pint` quantities and become plain floats inside;
the Rust side does the same with `uom`. Both implementations therefore run the
same arithmetic on the same numbers rather than each trusting its units library
to arrive there by a different route, and they agree bit-for-bit on the worked
examples. The test suite runs every spec case through both.

### Why Rust, if not for speed

For a single call, the Rust core is **not** faster than Python - PyO3 call
overhead exceeds the cost of the arithmetic. Its value here is being a *second
independent implementation* to cross-check the first against, plus a path to
performance that does not require rewriting the maths.

If performance is your goal, the win is the batch API, which **is** implemented:
`azoth.batch.hydraulics` and `azoth.batch.thermal` evaluate one calculation over
arrays, crossing the language boundary once instead of N times. It is a loop over
the same scalar kernels, not a vectorised second implementation, so it does not
weaken the two-implementation claim - see `docs/src/batch.md`.

The single-call path is still not the fast one, and this section exists so nobody
adopts the current design for a reason it does not support.

## The specs are the source of truth

`specs/calcs/**/*.yaml` defines each calculation: equation, source, valid range,
assumptions, worked example, tests. Each calc reads its range checks from the
generated table, so a bound edited in a YAML file changes behaviour in both
languages with no second edit.

CI regenerates the docs and the registries and fails on any diff, which is what
makes "the docs cannot drift from the code" a property of the build rather than a
claim in this file.

Adding a calculation means writing one spec, one Python function, one Rust
function, and declaring the tests. The docs, the range checks and the test cases
follow.

## Not for design work yet

`data/fittings/crane_k_factors.csv` holds **estimated dummy values** - plausible
magnitudes chosen so the software has something to run against. They are not from
Crane TP-410 or any other standard. A pressure drop computed from them can be
wrong by a factor of two and look entirely reasonable. Supply your own with
`azoth-data.example.yaml`; see
[Copyright and licensed data](docs/src/copyright.md).

Every affected result carries an `ESTIMATED_DATA` warning, `tools/spec_lint.py`
prints the count on every run, and a test fails the day someone populates the
file properly. Values that are cited but not confirmed by a named verifier carry
`UNVERIFIED_SOURCE` instead; only confirmed values are silent.

Water and air properties under `data/fluids/` are a different case: real
published values, marked `unverified` because they have not been checked against
a primary formulation.

## Verifying a result

See [TRUST.md](TRUST.md): how to verify a release tag, check a wheel with cosign,
reproduce a calculation by hand, and check a `provenance.json` record. It also
states plainly what none of that proves - a verified artifact means the code is
what it claims to be, not that the correlation is right for your fluid,
roughness or Reynolds number.

## Documentation

```bash
mdbook build docs && xdg-open docs/book/index.html
```

Every calc page carries the equation in LaTeX and in the form the library
evaluates, its source and verification status, inputs and outputs, the validated
range with the reason for each bound, the assumptions that are *not* checked, a
worked example, and the tests - including which are deliberately skipped and why.

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md) and the
[PR template](.github/PULL_REQUEST_TEMPLATE.md). Commits are signed.

Participation is covered by the [Code of Conduct](CODE_OF_CONDUCT.md), and
anything that could produce a wrong number should go through
[SECURITY.md](SECURITY.md) rather than the public issue tracker.

## Licence

**Code is AGPL-3.0-or-later.** If you run a modified version of this library as
a network service, the AGPL requires you to offer the modified source to that
service's users. That is deliberate: a validated calculation library is only
worth what its validation is worth, and validation that cannot be read cannot be
checked.

**Documentation and data are CC-BY-4.0** - the reference data under `data/` and
the prose under `docs/`. Attribution only, no copyleft. The equations and the
coefficients are the part most people want to reuse or cite, and they should not
require adopting a copyleft obligation to do it.

| What | Licence | File |
|---|---|---|
| Rust crates, Python package, CLI, tools | AGPL-3.0-or-later | [LICENSE](LICENSE) |
| `docs/`, `data/` | CC-BY-4.0 | [LICENSE-CC-BY-4.0](LICENSE-CC-BY-4.0) |

If you need a permissive licence for commercial use, the calculations themselves
are standard published equations and the citations are in each spec - you are
free to reimplement from the sources we cite.
