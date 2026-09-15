# azoth

**An opinionated port of NeqSim to Rust, with every calculation mirrored in Python.**

azoth is a thermodynamic and process library: equations of state, flashes, property
models and unit operations. The algorithms are carried across from
[NeqSim](https://github.com/equinor/neqsim), Equinor's open-source Java process
simulator, under Apache-2.0 and credited in [`NOTICE`](NOTICE) rather than re-derived.

What is azoth's own is the structure around them, and that is where the opinions are:

- **Composable by construction.** azoth's formal layer is a development in the
  calculus of constructions — the dependent type theory Lean implements — in which the
  dimension group is proved (`lean/Azoth/Dim.lean`) and the 24-unit vocabulary is data
  compiled into proved theorems (`lean/Azoth/Vocabulary.lean`). Above them, pi and rho
  specify a unit operation as a pure function on typed channels, so conservation is
  linearity and the balances are lemmas rather than checks. The crates are a small core
  with no cross-dependency, so a new domain is a new crate. What
  [the calculus of thermodynamic dimensionality](docs/src/calculus/index.md) states
  once, every calculation is an instance of.
- **A standard, not an encyclopaedia.** No single project ships every fluid and every
  correlation. azoth ships a *standard* others author against: the keycard (TOML data,
  checked on load), a new equation as one spec plus one Rust file and one Python file,
  and a `PropertyProvider` for fluid data it cannot ship. There is no runtime plugin
  registry - a calculation exists twice, once per language, or not at all.
- **The keycard is where responsibility sits.** The library implements; the engineer
  decides. A name outside the vocabulary is refused when the card is loaded; a
  coefficient must declare a unit, checked against the spec's declared unit; errors are
  typed and range violations are warnings. What azoth ships is NeqSim's, vendored; what
  you add is yours, and the responsibility for it is yours too.
- **Python is a real second implementation**, not a binding to a black box. Every
  calculation exists twice and the two are compared case by case, so azoth works in a
  notebook or a conda environment with no Rust toolchain at all.

[the specification](docs/src/architecture/specification.md) is the normative statement
of all four, and [`SPEC.md`](SPEC.md) points at it - why Rust rather than Java, what is
in scope and what deliberately is not, and what it costs to add a calculation. Where
another document disagrees with the specification, that document is wrong.

**Status: early.** What is implemented is the table below, and it is generated from the
specs rather than maintained by hand: a hydraulics kernel through Darcy-Weisbach
pressure drop, steady conduction, the Peng-Robinson equation of state through a
two-phase flash, stability testing and mixture critical points. The unit-operation tier
and the **flowsheet** above it are designed in
[the specification](docs/src/architecture/specification.md) and not built. The fitting coefficients azoth
ships are **placeholders, not engineering data**;
see [Not for design work yet](#not-for-design-work-yet), and
[the specification](docs/src/architecture/specification.md) has the programme. The
order in which the rest of NeqSim is ported is [ROADMAP.md](ROADMAP.md).

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

  no warnings: every range check passed
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

<!-- BEGIN GENERATED: implemented -->
| Calculation | What it does |
|---|---|
| `eos.antoine_vapor_pressure` | Antoine vapour pressure from NeqSim's correlation |
| `eos.bubble_pressure` | Bubble-point pressure — a *model* |
| `eos.chung_conductivity` | Gas thermal conductivity from the Chung correlation |
| `eos.chung_viscosity` | Gas viscosity from the Chung correlation |
| `eos.costald_molar_volume` | Saturated liquid molar volume from the COSTALD equation |
| `eos.critical_point` | Mixture critical point — a *model* |
| `eos.dew_pressure` | Dew-point pressure — a *model* |
| `eos.hayduk_minhas_diffusivity` | Liquid binary diffusivity from the Hayduk-Minhas correlation |
| `eos.heat_of_vaporization` | Heat of vaporisation from NeqSim's correlation |
| `eos.ideal_gas_cp` | Ideal-gas heat capacity from a polynomial |
| `eos.liquid_heat_capacity` | Liquid heat capacity from NeqSim's polynomial |
| `eos.mason_saxena_conductivity` | Gas mixture conductivity by Mason-Saxena mixing over Chung pure-component conductivities — a *model* |
| `eos.molar_enthalpy_entropy` | Molar enthalpy and entropy of a mixture — a *model* |
| `eos.ph_flash` | Pressure-enthalpy flash — a *model* |
| `eos.pr78_kappa` | Peng-Robinson (1978) attraction-parameter coefficient |
| `eos.pr_alpha_ab` | Peng-Robinson alpha function and reduced attraction parameters |
| `eos.pr_departure` | Peng-Robinson fugacity coefficient and departure functions |
| `eos.pr_kappa` | Peng-Robinson attraction-parameter coefficient |
| `eos.pr_mass_density` | Mass density from a molar volume |
| `eos.pr_molar_volume` | Molar volume from a compressibility factor |
| `eos.pr_peneloux_shift` | Peng-Robinson Peneloux volume-translation parameter |
| `eos.pr_z_factor` | Peng-Robinson compressibility factor |
| `eos.prsv_kappa` | Peng-Robinson-Stryjek-Vera alpha-function coefficient |
| `eos.ps_flash` | Pressure-entropy flash — a *model* |
| `eos.pt_flash` | Pressure-temperature flash — a *model* |
| `eos.pure_saturation` | Pure-component saturation pressure — a *model* |
| `eos.rachford_rice_binary` | Rachford-Rice vapour fraction, for a binary |
| `eos.rackett_molar_volume` | Saturated liquid molar volume from the Rackett equation |
| `eos.rk_alpha_ab` | Redlich-Kwong alpha function and reduced attraction parameters |
| `eos.rk_departure` | Redlich-Kwong fugacity coefficient and departure functions |
| `eos.siddiqi_lucas_diffusivity` | Liquid binary diffusivity from the Siddiqi-Lucas correlation |
| `eos.srk_alpha_ab` | Soave-Redlich-Kwong alpha function and reduced attraction parameters |
| `eos.srk_departure` | Soave-Redlich-Kwong fugacity coefficient and departure functions |
| `eos.srk_kappa` | Soave-Redlich-Kwong attraction-parameter coefficient |
| `eos.srk_peneloux_shift` | Soave-Redlich-Kwong Peneloux volume-translation parameter |
| `eos.srk_z_factor` | Soave-Redlich-Kwong compressibility factor |
| `eos.stability_test` | Tangent-plane stability test — a *model* |
| `eos.twu_kappa` | Twu attraction-parameter coefficient |
| `eos.tyn_calus_diffusivity` | Liquid binary diffusivity from the Tyn-Calus correlation |
| `eos.vdw1f_mix_binary` | van der Waals one-fluid mixing, for a binary |
| `eos.wilke_chang_diffusivity` | Liquid binary diffusivity from the Wilke-Chang correlation |
| `eos.wilke_viscosity` | Gas mixture viscosity by Wilke's rule over Chung pure-component viscosities — a *model* |
| `hydraulics.choked_flow_area` | Choked-flow throat area for an ideal gas |
| `hydraulics.control_valve_cv` | Liquid flow through a control valve from its flow coefficient |
| `hydraulics.crane_k_factors` | Fitting resistance coefficients by the equivalent-length method |
| `hydraulics.darcy_weisbach` | Darcy-Weisbach pressure drop |
| `hydraulics.friction_factor_colebrook` | Colebrook-White friction factor (implicit) |
| `hydraulics.friction_factor_haaland` | Haaland friction factor (explicit) |
| `hydraulics.friction_factor_swamee_jain` | Swamee-Jain friction factor (explicit) |
| `hydraulics.orifice_flow` | Flow through an orifice from the pressure difference across it |
| `hydraulics.pump_power` | Pump shaft power from flow, head and efficiency |
| `hydraulics.reynolds_number` | Reynolds number for pipe flow |
| `thermal.conduction_plane_wall` | Steady conduction through a plane wall |
<!-- END GENERATED: implemented -->

Pipe *with* fittings is a composition of the last two hydraulics calcs, done by the
`azoth pipe` CLI rather than by a calc of its own, because the two losses use
different methods and adding them is a modelling decision worth seeing explicitly.

Relief valve *sizing* to a standard is not implemented. `hydraulics.choked_flow_area`
is the isentropic basis - the throat area a given choked mass flow needs - and the
de-rating coefficients a standard applies are the caller's to compose.

A mixture **critical point** (`eos.critical_point`) is implemented, and the way it is
worth describing is by the route it does *not* take. The obvious one - solving
`dP/dV = d2P/dV2 = 0` at fixed composition - is exact for a pure component and predicts
the *same* `Z_c` for every mixture, because in reduced variables those two conditions
have a single universal root. It is a plausible-looking wrong number, and a
pure-component test cannot tell it apart from a right one. This uses Heidemann & Khalil
(1980), *AIChE Journal* 26(5), 769-779, whose two conditions do generalise: a mixture's
`Z_c` moves with composition, and a test asserts that it moves.

The method is ported from NeqSim's `CriticalPointFlash` (Apache-2.0) - see
[`NOTICE`](NOTICE) - which supplies the Q matrix as
the scaled Helmholtz Hessian at constant temperature and volume, and the nested Newton
that drives its smallest eigenvalue to zero. **That implementation validates its result
nowhere**, so it is a source for the method and not for the answer; everything this
model is checked against is this project's, including the closed-form `Z_c` a pure
Peng-Robinson fluid has. The spec's notes say plainly what has and has not been
confirmed, which is that nobody has read the paper.
[the specification](docs/src/architecture/specification.md) compares the two libraries in full.

The ids are namespaced by **domain** - `hydraulics.*`, `thermal.*`, `eos.*` - not by
project. They appear in provenance records and citations, so renaming the project does
not, and should not, invalidate them. A unit-operation or flowsheet namespace above
them will not be a fourth domain sitting beside the others but a composition tier,
and [the specification](docs/src/architecture/specification.md) sets out the programme for it.

The namespaces are also what proved the pipeline is domain-agnostic rather than shaped
around pipe flow. `eos` is where the shapes stop matching: an equation of state is
written in reduced variables, so it is dimensionless end to end and carries no unit at
all - and the same specs, generators, tests and documentation absorbed that with no
special case anywhere.

## Architecture

```
specs/**/*.toml        one file per calculation, model, case, and the unit vocabulary
        │
        ├─► tools/gen_registry.py   ─► crates/azoth-*/src/spec_gen.rs
        │                              python/src/azoth/_registry_gen.py
        ├─► tools/gen_models.py     ─► crates/azoth-eos/src/model_gen.rs
        │                              python/src/azoth/_models_gen.py
        ├─► tools/gen_vocabulary.py ─► crates/azoth-core/src/unit_vocab_gen.rs
        │                              python/src/azoth/core/_units_gen.py
        │                              specs/schema/unit.schema.json
        │                              lean/Azoth/Vocabulary.lean, Gate.lean
        ├─► tools/gen_docs.py       ─► docs/src/**, including SUMMARY.md
        ├─► tools/gen_stub.py       ─► python/src/azoth/_core.pyi
        └─► tools/provenance.py     ─► provenance.json
```

Every generated file comes from a file in this repository, and CI regenerates all of
them and fails on a diff - so a spec is not a document that is supposed to match the
code, it is the thing the code was made from. Python is the reference implementation
and Rust the second one, with PyO3 binding them; both run the same arithmetic on the
same numbers, and the test suite compares them by tolerance.

The pipeline in full, what is true and where it is written, and the three rules a port
follows are on **[Architecture](docs/src/architecture/index.md)**.

## Extending it

There are three ways in. There is no runtime plugin registry - no `register()` call and
no loading at run time - which means a calculation always exists twice, once in each
language, as the two-implementation rule requires.
[the specification](docs/src/architecture/specification.md) has the three.

Most of what you would want to change is **data, not arithmetic**. A **keycard** is one
TOML file that overrides or extends what the library ships — a component's critical
constants, a binary interaction parameter, a fluid's property table, a fitting's
equivalent length, a coefficient a calculation takes, or a named model variant:

```python
import azoth

card = azoth.keycard.load("keycard.toml")
azoth.eos.component("methane", card=card)  # your values, not the databank's
```

**The card is a value you pass, not a setting you make.** `load` reads a file and
returns what it says; it stores nothing, and a call handed no card reads what the
library ships. A process-wide card would give two answers for one calculation as soon
as somebody loaded a file between the two calls, and which is which would depend on
when they did it.

Nothing needs registering to make it apply, in either language. `keycard.example.toml`
is the template, `python tools/check_user_data.py` checks yours, and
[the specification](docs/src/architecture/specification.md) documents every section - including the two sections
that are compiled into the shipped data files by `tools/gen_user_data.py` rather than
read at run time, and so need a rebuild. What the library ships and where it came from is
in [`NOTICE`](NOTICE) and `databank/manifest.toml`.

**A new equation or procedure is code**, and it is a spec plus one Rust file and one
Python file — the section below is the whole contract. **A new source of fluid
properties** is a `PropertyProvider` on the Python side — a name, a density and a
viscosity — for the one case where the right answer depends on data this library cannot
ship.

## Not for design work yet

`data/fittings/crane_k_factors.csv` holds **estimated dummy values** - plausible
magnitudes chosen so the software has something to run against. They are not from
Crane TP-410 or any other standard. A pressure drop computed from them can be
wrong by a factor of two and look entirely reasonable. Supply your own with
`keycard.example.toml`; the attribution obligations are in [`NOTICE`](NOTICE).

The `verify_status` column records that, and a test fails the day someone populates
the file properly. Water and air under `data/fluids/` are a different case: real
published values, marked `unverified` because nobody has checked them against a
primary formulation.

That column exists on the data *this repository ships*, and not on the rows of a
keycard, which is a deliberate asymmetry rather than a leftover —
[the specification](docs/src/architecture/specification.md) is where the reasoning lives.

## Verifying a result

A release tag is signed and its wheels are signed with cosign, and every build emits a
`provenance.json` recording the commit and a SHA-256 of every spec, implementation, test
and data file. None of that proves the correlation is right for your fluid, roughness or
Reynolds number: a verified artifact means the code is what it claims to be.

## Documentation

```bash
mdbook build docs && xdg-open docs/book/index.html
```

Every calc page carries the equation in LaTeX and in the form the library
evaluates, its source, its notes, inputs and outputs, the validated range with the
reason for each bound, the assumptions that are *not* checked, a worked example, and
the tests - including which are deliberately skipped and why. A model page carries
its algorithm where a calc page carries its solver, because a model's spec fixes a
*loop* rather than an equation.

## Contributing

See the [PR template](.github/PULL_REQUEST_TEMPLATE.md). Commits are signed.

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
