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

**`TrayHydraulicsProbe` drives `TrayHydraulicsCalculator` directly, on states of its own.** This is
the tray hydraulics of `distillation_column`, `absorption_column` and `stripping_column`, and **none
of their captures holds its output**: all three declare `max_allowable_gas_load_factor` and refuse it
rather than reading it, and `getFsFactor()` is the quantity that would answer it. So this probe states
the geometry and the load and prints both sides - sixteen inputs, twenty-four answers and one
**sizing** - which is what makes the port's case set an oracle rather than a self-consistency
check.

**The twenty-fifth line is `sized_column_diameter_m`, and it is printed last because it is the
one call that changes the object it is read from**: `sizeColumnDiameter` writes the trial `1.0` m
into `columnDiameter`, re-derives the areas and the flooding velocity there, and leaves the sized
value behind - so every line above a row is read at the diameter that row stated. The eight
answers are `0.8`, `0.8`, `0.9`, `1.1`, `0.5`, `0.5`, `0.8` and `1.4` m, **every one of them an
entry of the class's own thirty-one-size table**: the two `0.5`s are the two low-vapour rows and
the two `0.8`s include the stated-weir row, because a stated weir length moves nothing. The spec
declares no output for it - it is a companion of `hydraulics.tray_hydraulics` rather than one of
its results - so the two languages hold it in their own tests instead.

Eight states: the class's own defaults, the other two tray types it names (`valve` and `bubble-cap`)
at that same state, a load above flood and a load where weeping trips, the weeping state again on a
`valve` tray, a stated weir length where the class's own is `-1.0`, and a wider column so the areas
are not all one diameter's.

```bash
java -cp .:neqsim-f0c7436.jar TrayHydraulicsProbe > captures/tray_hydraulics_probe.tsv
```

**Two behaviours the rows carry that a single state would not show.** A `bubble-cap` tray's minimum
vapour velocity is `0.0`, so `isWeepingOk` can never go false for it and its turndown ratio is `0.0`
at every load - the verdict is the *type's*. And a stated weir length is read back
(`calculated_weir_length_m` answers the stated `0.7` rather than the derived `0.73`) while moving
none of the flood, the pressure drop or the turndown, so the class reads the field and the
arithmetic does not.

**`ProcessProbe capacity` is the two capacity-limit families, on solved columns.**
`DistillationColumn` computes `Fs = u·sqrt(rho_g)` over the **total** cross-section against
`maxAllowableFsFactor` (`2.5`), with a utilization, a within-limit verdict and a minimum-diameter
inverse; `AbsorptionColumn` adds the Souders-Brown `Ks = u·sqrt(rho_g/(rho_l - rho_g))` against
`maxAllowableGasLoadFactor` (`0.15`) and overrides the Fs limit to `3.0`. **No capture in this
directory measured either family before**: `distillation_column`, `absorption_column` and
`stripping_column` all declare `max_allowable_gas_load_factor` and refuse it rather than reading
it, and the Fs family is declared nowhere.

Seventeen rows: the port's own `binary_methane_butane_4_stages` at three diameters and two stated Fs
limits, the lean-oil absorber and the hydrocarbon stripper, three stated limits on the absorber,
three rows whose *outlet* is reflashed - because the outlet is what both families read - and three
on a solved `PackedColumn`, whose two diameters are not one.

```bash
java -cp .:neqsim-f0c7436.jar ProcessProbe capacity > captures/process_column_capacity.tsv
```

**Five things the rows carry that one row would not.** `getPhase(0)` of an outlet is not the phase
a reader expects: on the solved absorber the **liquid** outlet is two-phase and its phase 0 is the
*gas* it carries (`11.576297622897625` against the system's `645.1739799825267`), so `rho_l - rho_g`
is `0.30`, the class's `10.0` floor is met and its `DEFAULT_LIQUID_DENSITY = 1000.0` is substituted
on the ordinary row. The reflashed-to-one-liquid row is the control: `K` is
`0.00826846105398939` against the fallback's `0.006694885232955643`, a factor of `1.235`, so the
substitution and not the field is what the answer rests on - and reflashing the same outlet to a
*gas* reproduces the fallback row to the last digit. Second, the two families read different
densities, which the two-phase gas-out row separates: the Fs family takes the system's
`15.34351809893286` and the K family phase 0's `14.4273646532561`. Third, `getFlowRate("m3/sec")`
is `n·M/rho` on every row - `_mass_over_density_m3_per_s` equals `_volumetric_flow_m3_per_s` to the
last digit - which is the identity the port spells a volumetric flow with. And fourth, both
verdicts are in the capture: `fs_factor_within_design_limit=false` at a stated `0.01` and
`gas_load_factor_within_design_limit=false` at a stated `0.005`.

**And the stripper row is what settles which vapour that is.** `getLiquidOutStream` answers a
stream carrying the *tray's* two-phase system, so its phase 0 is the vapour that tray leaves with -
`9.221459471348782` against the gas outlet's `11.267630874486859` on `hydrocarbon_stripper`, a
*negative* difference the class's `10.0` floor turns into the substitution. The liquid product's
own density would answer `0.0011870273688000244` against the class's `0.0009402581113639056`: a
port that read the product is exact on the absorber row and 26 per cent out on the stripper's,
which is why one row cannot settle it and two can.

The three minimum diameters are invariant under the diameter they are compared against
(`0.09032514302269502` on all three binary rows), since they are functions of the flow, the
density and the limit alone.

**A `PackedColumn` has two diameters and the Fs family reads only one of them.**
`setColumnDiameter` writes the class's own `columnDiameter`; `calcPackingHydraulics` *ends* by
writing the inherited `internalDiameter` from it, or from the designer's `getRequiredDiameter` when
it is not stated - so an unstated packed column carries `0.3` m, the diameter the sizing resolved,
and `fs_factor=0.22662865172417535` is measured against that and not against a caller's number. The
stated `0.5` row reads back `0.5` and moves the factor to `0.08158631462070312`, which is the
*distillation* column's own `diameter_0_5` row exactly: the packing's height is what makes that
column's stage count, and the capacity family is indifferent to the packing. Both packed rows also
read `max_allowable_fs_factor=2.5`, so `PackedColumn` overrides neither default.

**`ProcessProbe designer` is `ColumnInternalsDesigner`'s trayed half - the class the tray
calculator is driven *by*.** It walks the column's trays, picks the **controlling** one by the
largest vapour mass flow, sizes a diameter off that one and runs a calculator per tray at the
diameter it resolved, summing them. **No capture in this directory held it before**: it is the
class `DistillationColumn.updatePressureProfileFromHydraulics` builds inside its tear loop, and
`hydraulics.tray_hydraulics` is its arithmetic rather than its report.

```bash
java -cp .:neqsim-f0c7436.jar ProcessProbe designer > captures/process_internals_designer.tsv
```

Seven rows: the binary column at the designer's own defaults, at the other two tray types, at a
stated `columnDiameterOverride`, at a looser flood fraction, and at a different geometry - then
the lean-oil absorber, which is the no-ends case.

**Four things the rows carry.** `getTrays()` **includes the ends**: the binary column's four
stages walk as six trays and the absorber's five as five, so a port that read the middle trays
would miss the reboiler and the condenser. **The controlling tray's calculator is discarded** -
it exists only to size the diameter - and it is built **without a relative volatility** while
every tray in the summation loop carries one, so the sizing branch runs at the class's own
default `2.0` on a state whose trays read `9.25` and `20.0`. **The relative volatility is a
spread and not a K-value**: `getTrayProperties` takes `getPhase(0)`'s mole fraction over the
liquid phase's for every component, divides the largest ratio by the smallest and **caps the
result at `20.0`**, which two of the binary column's six trays reach. And `design_ok` is false on
all seven rows while the sized diameter is `0.5` m on every binary row - the standard-diameter
table's granularity, which a stated `0.5` override reproduces to the last digit.

**`ProcessProbe coupling` is the pressure-drop coupling, and it refutes the premise it was built
to test.** `setHydraulicPressureDropCouplingEnabled(true)` makes `updatePressureProfileFromHydraulics`
run inside `updateColumnTearVariables`: it builds a designer at the stated internals type, takes
its `getTotalPressureDrop`, and rewrites one end so the *difference* between the two equals that
drop - `bottomTrayPressure = topTrayPressure + drop/1e5` where the top is positive,
`topTrayPressure = max(1e-6, bottomTrayPressure - drop/1e5)` otherwise - then rebuilds the
profile. Six rows: the captured binary column with the coupling off and on, on at a looser tear
tolerance and under another internals type, and the lean-oil absorber with it off and on.

```bash
java -cp .:neqsim-f0c7436.jar ProcessProbe coupling > captures/process_hydraulic_coupling.tsv
```

**The flag does not gate the coupling alone.** `hasActiveColumnTearVariables` is
`!sideDrawSpecifications.isEmpty() || !pumparounds.isEmpty() || hydraulicPressureDropCouplingEnabled`,
so **turning the coupling on is what puts the column on `solveWithColumnTearVariables`** - the
coordinated loop `process_column_tear.tsv` measured as not converging on a side-draw flow
specification and a pumparound, which is why this capture was written before anything was ported.

**It converges on every row, and that is the finding.** `binary_coupling_on` reaches
`RIGOROUS_CONVERGED` in **2 tear iterations with 0 rejected candidates and 0 rollbacks**, a tear
residual of `7.564364502033816E-7`, and the same for the valve type, a looser tolerance and the
absorber (`2`, `51` inner iterations, residual `3.285096581696264E-7`). So this tear is **not**
the coordinated one of the side-draw row: that one searches a *candidate list* over several
coupled variables and rejects 18 of 30, while this one is a single scalar - the end pressure -
updated by a contraction, so two passes suffice.

**And it reproduces `applyHydraulicPressureDrop`'s own arithmetic to the last digit**: the
binary column's stated `20.0` bara becomes `19.017815747452257`, which is `19.0` plus the
`1781.574745225591` Pa the designer summed over `1e5`. The inner-iteration count is the sum over
the loop's three solves - `15 + 22 + 22 = 59` - which is what makes the tear's two iterations
visible in the sweep count rather than in the answer.

**`ProcessProbe mechanical` is `DistillationColumnMechanicalDesign`, and it found one defect
upstream.** `calcDesign` reads **tray 0** and nothing else - the vapour outlet's molar flow and
its fluid's density and molar mass, then the liquid outlet's - and from those sizes the diameter,
the height and the wall thickness. Two rows: the binary column this port already reproduces, and
the lean-oil absorber. **No capture in this directory held the class**, which is 1,348 lines with
a public surface of forty-odd getters.

```bash
java -cp .:neqsim-f0c7436.jar ProcessProbe mechanical > captures/process_column_mechanical_design.tsv
```

**The two density reads never ask the fluid to initialise, and on the absorber the liquid one
answers `0.0`.** `getWeirLoading()` is then **`Infinity`** - the volume flow is a mass flow over
that zero - and `getTrayPressureDrop()` collapses to its own `5.0` mbar constant, dropping the
liquid head without a word. **The library's own `ColumnInternalsDesigner.getTrayProperties` reads
the same tray and does call `fluid.initProperties()`**, which is what makes it an omission rather
than a convention. The capture carries both sides: the binary row's `435.74537901070215` against
the absorber's `0.0`, and the same read after one `initProperties()` - `645.1756172849634`, with
`2.664675513235188` and `8.164586402782746` after a second `calcDesign()`. **Filed as
[equinor/neqsim#4140](https://github.com/equinor/neqsim/issues/4140)**, with the minimal
reproduction and the class's own test as its starting point.

**Three quantities of the same class are *not* defects, and an earlier version of this entry said
they were.** `getTotalPressureDrop()` is in **bar**, `getReboilerDuty()` in **kW** and
`getColumnWallThickness()` in **mm**, each as its javadoc says; the `_Pa`, `_W` and metre labels
beside them in this capture are **the probe's**, and reading them as the class's units is what
produced the false claim. **And the `216.19496855345912` mm wall is the probe's own artefact**:
`DistillationColumn` sets the design pressure from `getEconomicDesignPressure()` before it calls
`calcDesign`, and this row constructs the design directly and so sizes against the base class's
`100.0` bara initialiser - a property of the measurement, not of the class. Both corrections are
in the issue's own "what I am not claiming" section.

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
