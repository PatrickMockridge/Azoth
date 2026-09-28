# Probe captures

A **capture** is a probe's standard output, committed verbatim. The probes live one
directory up, are compiled and run by hand against the pinned `neqsim-f0c7436.jar`, and what
they print is what lands here.

```bash
cd validation/neqsim
javac -proc:none -cp neqsim-f0c7436.jar CpaSweep.java
java -cp .:neqsim-f0c7436.jar CpaSweep > captures/cpa_sweep.tsv
```

**A probe that declares NeqSim's own package compiles with `-d .` and runs under that
package's name**, because what it needs is package-private; `WaxReferenceProbe` and
`PcsaftCompositionProbe` are the two, and each says in its header why.

```bash
javac -proc:none -cp neqsim-f0c7436.jar -d . PcsaftCompositionProbe.java
java -cp .:neqsim-f0c7436.jar neqsim.thermo.component.PcsaftCompositionProbe \
  > captures/pcsaft_composition_probe.tsv
```

**A state is passed as `T_K P_bara name z name z ...`** - the name and its mole fraction are
two arguments, not one `name:z` token. Seven probes take that form, and each refuses an odd
argument count rather than silently dropping a name whose fraction is missing. The other
three that read arguments take `T_K P_bara n_water` triplets.

**A state that a probe only prints when it is invoked for one needs its own recipe here**,
because this file is the registry `tools/oracle_sweep.py` reads: a run with no recipe has
no committed capture to compare against, so the next pin move cannot say whether it moved.
`PcsaftProbe`'s three aside from its two defaults are the ones a test asserts.

```bash
java -cp .:neqsim-f0c7436.jar PcsaftProbe 150 50 methane 1 > captures/pcsaft_probe_methane_150_50.tsv
java -cp .:neqsim-f0c7436.jar PcsaftProbe 300 100 propane 1 > captures/pcsaft_probe_propane_300_100.tsv
java -cp .:neqsim-f0c7436.jar PcsaftProbe 400 50 methane 1 > captures/pcsaft_probe_methane_400_50.tsv
```

**`ProcessProbe` takes a unit operation's name and prints that unit operation's inlet and
outlet on the palette's record**, so one recipe per unit operation and one capture each. It
is the only probe here that drives `neqsim.process` rather than `neqsim.thermo`.

```bash
java -cp .:neqsim-f0c7436.jar ProcessProbe pump > captures/process_pump.tsv
java -cp .:neqsim-f0c7436.jar ProcessProbe splitter > captures/process_splitter.tsv
java -cp .:neqsim-f0c7436.jar ProcessProbe mixer > captures/process_mixer.tsv
java -cp .:neqsim-f0c7436.jar ProcessProbe separator > captures/process_separator.tsv
java -cp .:neqsim-f0c7436.jar ProcessProbe throttling_valve > captures/process_throttling_valve.tsv
java -cp .:neqsim-f0c7436.jar ProcessProbe heat_exchanger > captures/process_heat_exchanger.tsv
java -cp .:neqsim-f0c7436.jar ProcessProbe stream > captures/process_stream_properties.tsv
java -cp .:neqsim-f0c7436.jar ProcessProbe column > captures/process_column.tsv
java -cp .:neqsim-f0c7436.jar ProcessProbe column_solvers > captures/process_column_solvers.tsv
java -cp .:neqsim-f0c7436.jar ProcessProbe column_efficiency \
  > captures/process_column_efficiency.tsv
java -cp .:neqsim-f0c7436.jar ProcessProbe column_divergence \
  > captures/process_column_divergence.tsv
```

**`ProcessProbe absorber_efficiency` is the absorber's Murphree override, and NeqSim's own solve
does not converge with it.** It runs the lean-oil absorber at `0.6` column-wide and at `1.0` with
methane overridden to `0.6`: **both take the 80-iteration cap**, reporting `FAILED` and
`FALLBACK_PRODUCTS` with temperature residuals of `3.25` and `22.66` K - where the same state with
no correction converges in 17 passes at `9.4e-5` and `RIGOROUS_CONVERGED`. So the correction is
what breaks it, and this capture is **evidence of a divergence rather than an oracle for a state**.
What it does pin to its own numbers is the resolution the override reads, printed per stage as
`tray_murphree_efficiency=` and `component_murphree_efficiency=`.

```bash
java -cp .:neqsim-f0c7436.jar ProcessProbe absorber_efficiency \
  > captures/process_absorber_efficiency.tsv
```

**`ProcessProbe column_tear` is the coordinated tear, and NeqSim does not converge it.** The
column carries a side-draw flow specification *and* a pumparound, which puts it on
`solveWithColumnTearVariables`' outer loop over every active tear variable rather than either of
its two fast paths. **Neither row converges, and the pair is a thousand times worse**: the side
draw alone ends at a residual of `4.4e-4` after 30 iterations with **18 of 30 candidates
rejected**, and with the pumparound beside it at `0.447` with an **empty candidate history** - so
the specification's search never runs and the loop is not the fast path with one more variable.
This capture is the *measurement* behind this port's refusal to carry the pair, not an oracle for
a state.

```bash
java -cp .:neqsim-f0c7436.jar ProcessProbe column_tear \
  > captures/process_column_tear.tsv
```

**`ProcessProbe column_divergence` is the Murphree divergence instrument.** It runs the same
binary column at `0.6` and `0.85` and prints two things `process_column.tsv` does not carry: the
class's per-pass `getConvergenceHistory()`, and per tray the *flashed* system's two phase
compositions beside the vapour the tray hands up with its phase count and phase-0 type. Without
those, the `7` K between the two implementations on tray 3 is a number with no reading attached.
`0.6` is the control: this port reproduces that row to `1e-4` K, so installing its endpoint and
running one pass of the port's map must move no tray - and it does not.

**`ProcessProbe column_efficiency` measures the per-stage efficiency rule.** It runs the
Murphree binary column three times: at `0.6` column-wide, at `0.6` with stage 3 overridden to
`0.85` through `setMurphreeEfficiency(3, 0.85)`, and at `0.6` with the same override spelled as
the array `setMurphreeEfficiencies` takes with `NaN` everywhere else. **The second and third
blocks are identical apart from their label**, which is what makes `NaN` the fall-through
`getEffectiveMurphreeEfficiency` documents rather than an absence. Each override block carries a
`tray_murphree_efficiency=` line, which is the class's *resolved* vector over all six stages -
and the first block has none, because it has no override to resolve.

**`ProcessProbe column_solvers` takes a solver's name after the subcommand** to run one
strategy alone, which is how a slow rung is measured without paying for the others. The
capture is the whole ladder: the binary column under all ten of `ColumnSolverFactory`'s
strategies and the deethanizer under `NAPHTALI_SANDHOLM`, which is the only one that converges
it. Two of the ten are left out of the deethanizer rows on purpose, and the probe says which:
`MESH_RESIDUAL` takes 594 s for the pair and `AUTO` had not returned after 900 s.

**`ProcessProbe packed_column` is the equilibrium packed column, and it prints the report's
inputs as well as its outputs.** Sixteen rows: eleven stage-count rows, where the height's
`ceil(h / 0.5)` rule and its two-tray floor are visible, and five solved ones. **Every solved row
prints, per tray, the four properties `ColumnInternalsDesigner.getTrayProperties` reads** - the
vapour's and the liquid's density, the liquid's viscosity, and the interphase surface tension -
plus the calculator's own verdicts beside the column's.

```bash
java -cp .:neqsim-f0c7436.jar ProcessProbe packed_column \
  > captures/process_packed_column.tsv
```

**Four things those keys pin, and each was unread before.** `isDesignOk()`, `isWettingOk()` and
the column's `isHydraulicsOk()` are three predicates that all read as "the hydraulics are fine",
and on these rows **all three are false**, so the row alone cannot tell them apart - which is
what makes publishing one under another's name a way to reproduce every row while answering a
different question. `getPressureDropPerMeter()` beside the column's `getPackingPressureDrop()`
settles a name that does not say which it is: on the 2.3 m row they are `0.27862153390634686` and
`0.6408295279845977`, **a ratio of exactly the bed height**, so the column's is the total.
`resolved_packing` and `resolved_packing_factor` (`Pall-Ring-50`, `180.0`) and
`calculated_diameter_m` (`0.3`) are the geometry and the sized diameter the report used.

**And `interphase_surface_tension_N_per_m` is `0.0` on every tray**, which is the measurement
behind the surface tension the report takes: `getTrayProperties` asks
`fluid.getInterphaseProperties().getSurfaceTension(0, 1)` and answers its own `0.02` fallback when
that is not finite and positive. So the route is the fallback, shown at the call rather than
inferred from the wetted area - and it is a **hydrocarbon** pair here, not the gas-and-aqueous one
the rate-based column's `estimateSurfaceTension` meets, so the `0.0` is not specific to water.

**`PackingProbe` drives `PackingHydraulicsCalculator` directly, on states of its own.** The
equilibrium packed column's capture prints the calculator's *outputs* and not the inputs it
reads - the two mass flows, the four transport properties, the two diffusivities and the
packing geometry - so this probe states them and prints both sides.

```bash
java -cp .:neqsim-f0c7436.jar PackingProbe > captures/packing_probe.tsv
```

**`PackingProbe size` is the second capture it writes, and the one the sizing path needs.** The
four rows above *set* a diameter, and `PackedColumn`'s own capture reaches sizing only through a
solved middle tray, so `sizeColumnDiameter` - the trial flood at **1.0 m**, the design velocity
from `designFloodFraction`, the area from the vapour's volumetric flow, and the round up to the
class's standard-diameter table - was measured by nothing. Five states: the absorber, structured
packing, a high liquid load, a light load that reaches the table's **0.3 m floor**, and the
absorber again at a `0.5` flood fraction. **That last pair is the measurement of the one claim** -
the flood fraction moves the sized diameter from `0.4` to `0.5` m and leaves the calculation
untouched - and **the light and high-liquid rows share a flooding velocity to the last digit**,
which is the Eckert fit's `FLV` clamped at its ceiling of `5` on both.

```bash
java -cp .:neqsim-f0c7436.jar PackingProbe size \
  > captures/packing_sizing_probe.tsv
```

**`WaterCpSentinel` asks a question about the *data* rather than about a model.** `COMP.csv`
gives 131 of its 389 rows the whole of water's ideal-gas Cp polynomial - the same five numbers
`devtools/generate_water_caloric_alpha_reference.py` fits for water - and 130 of those rows are
substances other than water. This drives one single-component `SystemSrkEos` per name in the
table and compares `getCp0` against water's **exactly**, because two different molecules cannot
share an ideal-gas heat capacity: an exact match is not a tolerance question, and a near one
would be evidence of nothing. It parses the `NAME` column quote-aware, since names carry commas
inside them, and its capture is the sorted list of the names that match.

```bash
java -cp .:neqsim-f0c7436.jar WaterCpSentinel > captures/water_cp_sentinel.tsv
```

**`ProcessProbe rate_based` runs every state on two cubics.** `RateBasedPackedColumnTest` uses
SRK and this library's process layer resolves PR, so each state is run on both and the pair
measures what the cubic moved; the PR rows are the ones the port is held to. `teg_dehydration`
and `teg_circulation` are `SystemSrkCPAstatoil` with mixing rule 10, which PR cannot
represent, so they are evidence for a refusal rather than an oracle for a state.

```bash
java -cp .:neqsim-f0c7436.jar ProcessProbe rate_based \
  > captures/process_rate_based_packed_column.tsv
```

**`ProcessProbe rate_based_billet` is the same absorber state run twice**, once at each
`MassTransferCorrelation`, because that parameter is not a correlation: it scales the two film
coefficients by constants the packing row carries, so one state at each value is the whole
measurement. `Pall-Ring-50` resolves to the file's plastic row, `cp = 0.698` and `ch = 2.725`,
so the pair is `1.1354166...` on `kGa` and `0.698` on `kLa`; the capture's own two ratios are
what the port's test reads, and the liquid one is `0.6978956` rather than `0.698` because the
base `kLa` moved `1.495e-4` between the two converged states.

```bash
java -cp .:neqsim-f0c7436.jar ProcessProbe rate_based_billet \
  > captures/process_rate_based_billet.tsv
```

**`ProcessProbe rate_based_solvers` measures the two solvers this port refuses.** It runs the
class's own absorber state - the PR row `process_rate_based_packed_column.tsv` holds and the port
reproduces to `1e-4` - five ways: the pair this port carries as the control, then
`SegmentSolver.SIMULTANEOUS_RESIDUAL` and `ColumnSolver.EQUATION_ORIENTED` alone, then both, then
`EQUATION_ORIENTED` at the settings of the class's own test for it. **No refused row converges**,
and the control reproduces the baseline exactly (fifteen passes at `6.009670053264138e-10` mol/s),
which is what makes the rest a reading rather than a mis-set-up state. The segment solve takes the
full 20-pass cap at a `1.76` mol/s outlet residual against a `1e-9` gate and publishes a vapour at
`2.2e-29` K; the column-wide Newton stalls after two passes at a residual of `65.9` against its own
`1e-6` gate and publishes a `153` K vapour beside a `454` K liquid, and at its own test's settings
the norm is `30.5` against `1e-5`. **That test asserts finiteness and never convergence** - it caps
the Newton at two iterations and removes the homotopy - and its companion caps them at one.

```bash
java -cp .:neqsim-f0c7436.jar ProcessProbe rate_based_solvers \
  > captures/process_rate_based_solvers.tsv
```

They are committed rather than regenerated in CI because the jar is gitignored, so a gate
that ran the JVM could not run on a runner that has no NeqSim checkout. Committing the
output is what lets `tools/gen_neqsim_cases.py --check` be a build gate: it reads a capture
and emits `validation/eos/*.json`, so a case and the probe output it came from move
together or the gate fails.

**`tools/oracle_sweep.py` runs them all at once and is the step a pin move starts with.**
It re-runs every probe and every recipe below, compares each against the capture committed
here, and reports the keys that moved - and then any number in a Rust test that still holds
the *previous* value, which is what a commit regenerating these cannot see on its own. It
also names the probes with no capture and the captures no run reproduces; both are gaps of
the same kind, one on each side.

**A capture is a pinned artifact, not a live one.** Regenerating one means running the probe
again, and the cases derived from it move with it — which is the point, and also the reason
they are here rather than fetched.

**One capture cannot be regenerated at all.** `multi_salt_probe.tsv` is the multi-mineral
complementarity loop's only oracle — a brine of `Ca++`, `SO4--`, `CO3--` and `HCO3-` whose
three minerals take three rounds — and it was written before `PitzerParameterCoverage`
existed. Both pinned revisions refuse that brine, because the legacy dataset carries no theta
for `Ca++|Na+` or `Cl-|SO4--`, so the capture is a *record* of a state no current NeqSim will
run rather than something to re-measure. `tools/oracle_sweep.py` reports it as a capture with
no probe, which is the honest form of that.

Each capture's own header says what it is and, where the probe emits a key table, what each
key came from.

## Which captures a layer diff can read

`tools/neqsim_layer_diff.py` compares a model's intermediates against a capture key by key,
and it can only read `key = value` rows. **Two of the electrolyte captures are not in that
shape at all**: `pitzer_arithmetic.tsv` and `soreide_alpha_derivatives.tsv` print aligned
columns of numbers with no keys - `T`, `alpha`, `dalpha/dT`, and the Debye-Hückel parameter
against temperature - so they are read by hand, and the models built against them are pinned
by their cases rather than by a layer diff.

`pitzer_probe.tsv`, `soreide_whitson_probe.tsv` and `ge_electrolyte_probe.tsv` *are* keyable,
and still have none, because what they print is **which parameters a model selected** - the
dataset id, the active ions, the missing binary terms, the `kij` table - rather than a
model's own arithmetic. A layer diff exists to say which layer of a calculation moved, and
that is not a question those captures answer.

`furst_probe.tsv` was written in the readable shape and in the order the phase builds its
layers, for exactly this reason, and `python/src/azoth/eos/layers.py` carries the dumper
that meets it.
