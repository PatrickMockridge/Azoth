# Python scripting and Jupyter

The library is a Python package with a Rust core. A script and a notebook use the same
surface: `import azoth`, pass `pint` quantities in, and read a frozen result out. This page is
the route for both, and the five worked notebooks under
[`notebooks/`](https://github.com/PatrickMockridge/azoth/tree/main/notebooks) are the long form.

## Install and import

Python 3.12 or later:

```bash
pip install azoth-engine
```

The distribution is `azoth-engine` because `azoth` is taken on PyPI; the import name is
`azoth`.

```python
import azoth

q = azoth.ureg.Quantity
print(azoth.__version__)
```

For the notebooks and the plotting they do, install the tools beside it:

```bash
pip install jupyterlab matplotlib pandas
```

The repository also declares a `notebook` extra naming those tools; it ships with the next
release. Until then, install them explicitly.

## The shape of a call

Every public function lives in a namespace — `hydraulics`, `thermal`, `eos`, `process`,
`batch`, `reactions`, `standards`, `properties`, `keycard` — and returns a frozen dataclass.

- **Dimensioned inputs are `pint` quantities.** A bare number where a length is expected raises
  `UnitMismatchError`.
- **Dimensionless inputs are plain numbers** — Reynolds number, friction factor, mole fractions.
- **The result carries its warnings**: `result.warnings`, `result.is_clean`, and
  `result.provenance`. `to_dict()` and `to_json()` include the provenance block.
- **A warning is not an error.** A value outside its validated range is returned with
  `OUT_OF_VALID_RANGE`; a value that could not be checked is returned with
  `RANGE_CHECK_SKIPPED`. The two states stay different.

[Reading a calculation page](./reading-a-calculation.md) explains the sections behind each
generated reference page; [what a claim is worth](../architecture/claims.md) explains the
verification status a result reports.

## Units, and the one rule that surprises people

An absolute temperature and a temperature difference are different quantities. A `dT` argument
refuses `degC` and accepts `delta_degC`, `delta_degF` or `K`:

```python
azoth.thermal.conduction_plane_wall(
    q(0.6, "W/(m*K)"), q(12.0, "m**2"), q(30.0, "degC"), q(0.2, "m")
)
# UnitMismatchError: `dT` ... is a temperature difference, so pass `delta_degC`/`delta_degF`

azoth.thermal.conduction_plane_wall(
    q(0.6, "W/(m*K)"), q(12.0, "m**2"), q(30.0, "delta_degC"), q(0.2, "m")
).q  # 1080 W
```

[The calculus of dimensions](../calculus/dimensions.md) is why the refusal is at the boundary
rather than in the answer.

## A script in ten lines

A pipe, its Reynolds number, a friction factor and the pump that moves it:

```python
import azoth

q = azoth.ureg.Quantity

rho, mu = q(998.0, "kg/m**3"), q(1.002e-3, "Pa*s")
D, v = q(0.1, "m"), q(1.5, "m/s")

re = azoth.hydraulics.reynolds_number(rho, v, D, mu)
f = azoth.hydraulics.friction_factor_colebrook(re.re, 4.5e-4).f
dp = azoth.hydraulics.darcy_weisbach(f, q(100.0, "m"), D, rho, v, mu)
power = azoth.hydraulics.pump_power(rho, q(0.02, "m**3/s"), q(30.0, "m"), 0.72)

print(re.re, re.regime, dp.dp.to("kPa"), power.power.to("kW"))
```

`darcy_weisbach` also takes the velocity and the viscosity directly; passing `mu=None` is what
makes the Reynolds bound unchecked, and the result then reports `RANGE_CHECK_SKIPPED` with
`skipped_checks == ('re',)`.

The generated pages for
[`reynolds_number`](../hydraulics/reynolds_number.md),
[`darcy_weisbach`](../hydraulics/darcy_weisbach.md) and
[`pump_power`](../hydraulics/pump_power.md) carry the equations and the checked ranges.

## Thermodynamics

A mixture is built once, then flashed at a state:

```python
from azoth.eos.components import mixture_of

mix, ideal = mixture_of(["methane", "n-butane"])
z = [0.9, 0.1]

flash = azoth.eos.pt_flash(mix, q(250.0, "K"), q(20.0, "bar"), z)
flash.phase, flash.vapour_fraction, flash.x, flash.y
```

`mixture_of` lives in `azoth.eos.components` and returns the `Mixture` and the ideal-gas model
the enthalpy calculation needs. The mixture has no single molar-mass property; each component
carries its own `molar_mass`.

- `pt_flash` returns `phase` and `vapour_fraction`. Outside a two-phase split the vapour
  fraction is the negative-flash value and can lie outside `[0, 1]`, so `phase` is the field to
  read first.
- `pt_phase_envelope` returns the dew and bubble branches as parallel temperature and pressure
  tuples.
- `molar_enthalpy_entropy` takes the composition and the compressibility root the flash chose.
- A flash of a mixture whose components are not in the databank is an error, not a default.

See [`pt_flash`](../eos/pt_flash.md), [`pt_phase_envelope`](../eos/pt_phase_envelope.md) and
[`molar_enthalpy_entropy`](../eos/molar_enthalpy_entropy.md).

## Process simulation

A flowsheet is a TOML document. `Session` holds it, checks each edit, and returns one envelope:

```python
import azoth.process

session = azoth.process.Session(document)  # document is TOML text
env = session.run()

env["ok"], env["dirty"], env["run_error"]
env["session"]["converged"], env["session"]["iterations"]
env["session"]["streams"]["p1.outlet"]["T"]
```

The envelope carries the document, the graph, the diagnostics, the paths, the run's report and
the provenance of the units it ran. A scalar inside a stream is
`{"magnitude_si": ..., "unit": ...}`.

An edit returns a new envelope and marks the values stale; it does not pay the physics:

```python
env = session.apply({
    "command": "set_parameter",
    "instance": "p1",
    "name": "outlet_pressure",
    "value": 4.0e6,
})
env["dirty"]  # True
env = session.run()
session.value("p1.outlet.P")  # 4000000.0, SI
```

The unit operations are also functions, so one machine can be checked without a document:

```python
azoth.process.pump(
    ["methane", "n-butane"], q(1.0, "mol/s"), [0.9, 0.1],
    q(5.0e5, "Pa"), q(300.0, "K"), 0.75, q(4.0e6, "Pa"),
).power
```

**From an installed wheel, use `Session`.** `run_flowsheet` and `validate` take a real
`palette_dir` string that defaults to `specs/unit_ops` in a checkout; the wheel does not ship
`specs/`. `Session` (and `forms`, `tools`, `catalogue`) default to the palette compiled into
the extension.

[The middleware](../architecture/middleware.md) describes the envelope and the command model;
[the process schema](../architecture/process-schema.md) is the document format;
[session](../calculus/session.md) is the argument for holding the document rather than the
values.

## Notebooks

The worked examples are committed with their outputs, so GitHub renders their tables and
figures without running anything:

- [`01-simple-calculations.ipynb`](https://github.com/PatrickMockridge/azoth/blob/main/notebooks/01-simple-calculations.ipynb)
  — units, Reynolds, friction factors, Darcy-Weisbach, a pump and a wall, and the provenance
  block.
- [`02-intermediate-thermodynamics.ipynb`](https://github.com/PatrickMockridge/azoth/blob/main/notebooks/02-intermediate-thermodynamics.ipynb)
  — a mixture, a flash, the P-T envelope, and properties at the two roots.
- [`03-advanced-calculations.ipynb`](https://github.com/PatrickMockridge/azoth/blob/main/notebooks/03-advanced-calculations.ipynb)
  — an LMTD correction factor, a packed column's operating window, and a check that did not run.
- [`04-process-simulation.ipynb`](https://github.com/PatrickMockridge/azoth/blob/main/notebooks/04-process-simulation.ipynb)
  — the palette, a flowsheet with an active recycle, an edit-and-run, and a unit operation by
  hand.
- [`05-data-science.ipynb`](https://github.com/PatrickMockridge/azoth/blob/main/notebooks/05-data-science.ipynb)
  — the batch API, pandas tables, compressibility and conduction sweeps, and ISO 6976.

They read no repository file; the flowsheet in notebook 04 is embedded as text. CI executes
them, so an example that stops running fails the push rather than the reader.

## Data science

`azoth.batch` evaluates the same scalar kernels over arrays with one crossing into the core.
It is a loop over the same kernels, not a second implementation:

```python
import azoth.batch.hydraulics

result = azoth.batch.hydraulics.darcy_weisbach(
    f=[0.02, 0.02], L=[100.0, 100.0], D=[0.1, 0.1],
    rho=[998.0, 998.0], v=[1.0, 2.0], mu=[1.002e-3, 1.002e-3],
)
list(result.dp), result.units, result.quantities("dp")
```

Inputs are plain numbers in the spec's canonical unit, so batch gives up unit checking; outputs
are SI base magnitudes with a unit map. There is no batch flash. A DataFrame is one dict away,
and [the batch API](../index.md#the-batch-api) states the rest.

## Provenance

Every result answers what it rests on:

```python
result.provenance.calc_id          # the spec id
result.provenance.spec_sha256      # the spec the bounds came from
result.provenance.python_sha256    # this implementation's bytes
result.provenance.rust_sha256      # the other implementation's bytes
result.provenance.verification     # verified / partially_verified / source_needed / unverified
result.provenance.skipped_checks   # declared checks this call could not evaluate
```

`to_json()` carries the same block, so a result written to a file, served over MCP, or handed
to an agent carries its provenance with it. `checked and fine` and `never checked` are
different states, and this is where they stay different.

## Troubleshooting

- **`AttributeError: module 'azoth' has no attribute 'standards'`** — `standards` is not bound
  by `import azoth`; use `import azoth.standards` or `from azoth import standards`.
- **`UnitMismatchError` on a temperature difference** — pass `delta_degC`, `delta_degF` or `K`.
- **A palette error from an installed package** — `run_flowsheet`/`validate` need a real
  `palette_dir`; use `Session`, whose palette is embedded.
- **`ImportError` on `azoth.process` from a source tree** — the compiled extension is not
  built; install the wheel or run `maturin develop`.
- **A warning instead of a refusal** — read it. It names the field and the bound; the value is
  still returned because a warning is not an error.

## Before you rely on a number

The port is incomplete and every shipped model is `unverified` or `partially_verified`.
**These are screening values, not design deliverables** — [the front page](../index.md#not-for-design-work-yet)
says what is missing.
