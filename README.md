# chemeng

Open, validated, citable chemical engineering calculations.

Every calculation ships with its equation, the source it came from, the range in
which it is validated, its assumptions, a worked example, and tests. The
documentation is generated from the same machine-readable specs the code is
generated from, so the two cannot drift apart.

> **Status: early.** This repository currently implements one vertical slice - the
> hydraulics kernel through Darcy-Weisbach pressure drop. The fitting
> coefficients it uses are **placeholders, not engineering data**. See
> [Not for design work yet](#not-for-design-work-yet).

## Install and run

```bash
cargo test                       # Rust core
pytest                           # Python package (from python/)
mdbook build docs                # documentation
cargo run -p chemeng-cli -- pipe --help
```

The end-to-end check:

```bash
chemeng pipe --fluid water --flow 10 --diameter 0.1 --length 100 \
    --fittings "90_elbow,gate_valve_open"
```

## What's implemented

| Calc | What it does |
|---|---|
| `hydraulics.reynolds_number` | Reynolds number and flow regime |
| `hydraulics.friction_factor_colebrook` | Implicit Colebrook-White friction factor |
| `hydraulics.friction_factor_swamee_jain` | Explicit approximation to it |
| `hydraulics.crane_k_factors` | Fitting losses by the equivalent-length method |
| `hydraulics.darcy_weisbach` | Pressure drop over a straight pipe |

Pipe *with* fittings is a composition of the last two, done by the `chemeng pipe`
CLI rather than by a calc of its own, because the two losses use different
methods and adding them is a modelling decision worth seeing explicitly.

## Two ideas the library is built around

**Warnings are not errors.** A value outside the range in which a correlation was
validated is still a value. Refusing to return it would be less useful than
returning it with a warning - but the library must never return it *silently*.
Check `result.warnings`, or `result.is_clean` if you are willing to see every
caveat at once.

**A check that could not run is not a check that passed.** When an optional input
is missing, the range check that depends on it reports `RANGE_CHECK_SKIPPED`
rather than quietly succeeding. In `darcy_weisbach`, omitting the viscosity
leaves the flow regime unchecked, and the result says so.

## Architecture

Python for the reference implementation, Rust for the core, PyO3 binding them.

```python
import chemeng

q = chemeng.ureg.Quantity

r = chemeng.hydraulics.darcy_weisbach(
    0.02, q(100.0, "m"), q(0.1, "m"), q(998.0, "kg/m**3"), q(1.5, "m/s")
)
r.dp  # 22455.0 pascal
r.warnings  # (Warning(RANGE_CHECK_SKIPPED, ...),)
```

Units cross the public API as `pint` quantities and become plain floats inside;
the Rust side does the same with `uom`. Both implementations therefore run the
same arithmetic on the same numbers rather than each trusting its units library
to arrive there by a different route, and they agree bit-for-bit on the worked
examples.

### Why Rust, if not for speed

For a single call, the Rust core is **not** faster than Python - PyO3 call
overhead exceeds the cost of the arithmetic. Its value here is being a *second
independent implementation* to cross-check the first against, plus a path to
performance that does not require rewriting the maths. The cross-language test
suite runs every spec case through both.

If performance is your goal, the win is a batch API operating on arrays, which is
not implemented. This section exists so nobody adopts the current design for a
reason it does not support.

## The specs are the source of truth

`specs/calcs/**/*.yaml` defines each calculation: equation, source, valid range,
assumptions, worked example, tests. From those, `tools/` generates:

- `crates/chemeng-hydraulics/src/spec_gen.rs` - range checks and test cases
- `python/src/chemeng/_registry_gen.py` - the same, as Python data
- `docs/src/**` - the documentation

CI regenerates and fails on any diff, so the docs cannot drift from the code.
Each calc reads its range checks from the generated table, meaning a bound edited
in a YAML file changes behaviour with no second edit.

## Not for design work yet

`data/fittings/crane_k_factors.csv` holds **estimated dummy values** - plausible
magnitudes chosen so the software has something to run against. They are not from
Crane TP-410 or any other standard. A pressure drop computed from them can be
wrong by a factor of two and look entirely reasonable. Every affected result
carries an `ESTIMATED_DATA` warning, `tools/spec_lint.py` prints the count on
every run, and a test fails the day someone populates the file properly.

Water and air properties under `data/fluids/` are a different case: real
published values, marked `unverified` because they have not been checked against
a primary formulation.

## Verifying a result

See [TRUST.md](TRUST.md) for verifying a release tag, checking a wheel with
cosign, reproducing a calculation, and reading `provenance.json`.

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md). Adding a calc means writing one YAML
spec, one Python function, one Rust function, and one test - the docs and the
generated registries follow.

## Licence

Code under dual MIT / Apache-2.0. Documentation and data under CC-BY-4.0. See
`LICENSE-MIT`, `LICENSE-APACHE` and `LICENSE-CC-BY`.
