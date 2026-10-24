# Probe captures

A **capture** is a probe's standard output, committed verbatim. The probes live one
directory up, are compiled and run by hand against the pinned `neqsim-f0c7436.jar`, and what
they print is what lands here.

```bash
cd validation/neqsim
../../tools/gated.sh javac -proc:none -cp neqsim-f0c7436.jar CpaSweep.java
../../tools/gated.sh java -cp .:neqsim-f0c7436.jar CpaSweep > captures/cpa_sweep.tsv
```

**Eight captures are run against the pin plus the fixes that moved them.**
`packing_probe`, `packing_sizing_probe`, `process_packed_column`,
`process_rate_based_packed_column`, `process_rate_based_billet` and
`process_rate_based_solvers` were taken from a jar built at `f0c7436c` with only the packing fix
`834b5c4f` applied: its `Packing.csv`, `PackingSpecificationLibrary.java`,
`PackingHydraulicsCalculator.java`, and the single `getOrDefault` -> `require` line of
`RateBasedPackedColumn.java`. `process_shortcut_distillation_column` adds the Kirkbride fix
`1f7581f9`'s `ShortcutDistillationColumn.java`.

**The same jar carries the re-vendored `COMP.csv`** (upstream `8bb7bace`), because the port's
component data is that file: `ge_henry_probe`, `ge_electrolyte_probe`, `pitzer_fugacity_probe`,
`pitzer_hydrate_probe` and `aqueous_viscosity_probe` are taken with it.
`pseudo_component_combiner_probe` adds the characterization fix `c93c8153`'s one hunk: the
mass-quantile candidate is projected into `nextUp`/`nextDown` of the reference gap rather than
discarded for the midpoint (upstream PR #4330).

**The jar also carries the re-vendored `REACTIONDATA.csv` and `REACTIONDATAPITZER.csv`**
(upstream `1acf0652`), so `reactive_flash_probe`, `reactive_ph_flash_probe`,
`hybrid_eos_ge_reactive_probe` and `hybrid_eos_ge_reactive_loop_probe` are taken with them. Everything else is the pin, so each capture's delta is those changes rather than the
five hundred commits between the revisions. The built-in packing rows are
the reference after the fix: `Pall-Ring-50` is `120.0`, `0.96`, factor `66`, and the file's
plastic ring is reachable as `Pall-Ring-50-plastic` (`111.1`, `0.919`, `180`).

**A probe that declares NeqSim's own package compiles with `-d .` and runs under that
package's name**, because what it needs is package-private; `WaxReferenceProbe`,
`PcsaftCompositionProbe` and `WaterDewPointSeedProbe` are the three, and each says in its
header why - the last one for `TPmultiflash.multiPhaseTest`, which is the flag its whole
chain of seeding hinges on and which no probe outside the package can read.

```bash
../../tools/gated.sh javac -proc:none -cp neqsim-f0c7436.jar -d . PcsaftCompositionProbe.java
../../tools/gated.sh java -cp .:neqsim-f0c7436.jar neqsim.thermo.component.PcsaftCompositionProbe \
  > captures/pcsaft_composition_probe.tsv
```

**`TpMultiflashSeedProbe` attributes every phase addition in the multiphase flash.** It
declares NeqSim's package because `multiPhaseTest` and the two ablation flags are
package-private and the gas/oil seeds are private methods, so a probe outside would have to
read the chain rather than measure it. Each construct runs on its own fresh fluid, and the
full chain is run twice more with the post-flash enhanced check and the inline aqueous seed
ablated: a construct whose ablation changes the phase inventory fired. **Every measurement
builds a new fluid**, because `setMultiPhaseCheck(true)` is a system flag and a reused fluid
makes the next `TPflash` a multiflash rather than the two-phase seed the probe starts from.

```bash
../../tools/gated.sh javac -proc:none -cp neqsim-f0c7436.jar -d . TpMultiflashSeedProbe.java
../../tools/gated.sh java -cp .:neqsim-f0c7436.jar neqsim.thermodynamicoperations.flashops.TpMultiflashSeedProbe \
  > captures/tp_multiflash_seed_probe.tsv
```

**A state is passed as `T_K P_bara name z name z ...`** - the name and its mole fraction are
two arguments, not one `name:z` token. Seven probes take that form, and each refuses an odd
argument count rather than silently dropping a name whose fraction is missing. The other
three that read arguments take `T_K P_bara n_water` triplets.

**A state that a probe only prints when it is invoked for one needs its own recipe here**,
because this file is the registry `tools/oracle_sweep.py` reads: a run with no recipe has
no committed capture to compare against, so the next pin move cannot say whether it moved.
`PcsaftProbe`'s three aside from its two defaults are the ones a test asserts. **A probe that
takes no argument needs one too**, because the recipe is also what names the capture a bare run
writes - `IapwsIf97Probe` prints its whole state grid and reads nothing.

```bash
../../tools/gated.sh java -cp .:neqsim-f0c7436.jar IapwsIf97Probe > captures/iapws_if97_probe.tsv
```

**`SurfaceTensionProbe` records model 0 and the two calculators it installs.** With no
argument it prints the default pair - Parachor for gas/oil and gas/aqueous, Firozabadi-Ramley
for liquid/liquid - and calls both calculators directly on the same states. The `solvers`
invocation adds models 1, 2 and 4, which build the GT and LGT density-profile solvers; those
rows are the next tranche's expected values and are kept out of the bare run so the weekly
sweep does not pay for a solver no test asserts yet. `TransportMethodsProbe` selects
`Filippov` conductivity and `friction theory` viscosity on the densest phase and prints every
intermediate, because both correlations carry bar/Pa/micropoise conversions that a final
number alone would hide.

```bash
cd validation/neqsim
../../tools/gated.sh javac -proc:none -cp neqsim-f0c7436.jar SurfaceTensionProbe.java TransportMethodsProbe.java
../../tools/gated.sh java -cp .:neqsim-f0c7436.jar SurfaceTensionProbe > captures/surface_tension_probe.tsv
../../tools/gated.sh java -cp .:neqsim-f0c7436.jar SurfaceTensionProbe solvers > captures/surface_tension_solvers_probe.tsv
../../tools/gated.sh java -cp .:neqsim-f0c7436.jar TransportMethodsProbe > captures/transport_methods_probe.tsv
```

**`LgtProfileProbe` records the gradient-theory profile behind model 2, not only its final
number.** It prints the influence parameter, aT and b of every component, then the 200-step
molar-density path, the pressure at each step (negative from step 53 to 198, where the cubic's
homogeneous path crosses its loop) and the accumulated `z`. `LgtLocalProbe` rebuilds the
class's one-phase local system outside the class: at the negative-pressure steps
`getGibbsEnergy(T, P)` is NaN while `getChemicalPotential(phase)` is finite, because
`ComponentEos` overrides the potential with the T,V form
`H_id - T S_id + R T ln(n R T/V / P_ref) + R T dF/dN`. That is the primitive a port needs, and
it is why the profile does not stop at the loop.

```bash
cd validation/neqsim
../../tools/gated.sh javac -proc:none -cp neqsim-f0c7436.jar LgtProfileProbe.java LgtLocalProbe.java
../../tools/gated.sh java -cp .:neqsim-f0c7436.jar LgtProfileProbe > captures/lgt_profile_probe.tsv
../../tools/gated.sh java -cp .:neqsim-f0c7436.jar LgtLocalProbe > captures/lgt_local_probe.tsv
```

**`GtSimpleProfileProbe` records the profile `GTSurfaceTensionSimple`'s inner Newton
produces**, and `GtSimpleInitialProbe` the step-0 chemical potentials and `dmu/dn` it starts
from. The CO2/water block is the one to read twice: the class's densities are **zero from
step 3** and it returns the two-step partial integral `0.0145737553`, which is why the port
refuses that state.

```bash
cd validation/neqsim
../../tools/gated.sh javac -proc:none -cp neqsim-f0c7436.jar GtSimpleProfileProbe.java GtSimpleInitialProbe.java
../../tools/gated.sh java -cp .:neqsim-f0c7436.jar GtSimpleProfileProbe > captures/gt_simple_profile_probe.tsv
../../tools/gated.sh java -cp .:neqsim-f0c7436.jar GtSimpleInitialProbe > captures/gt_simple_initial_probe.tsv
```

**`GeHenryProbe` measures `ComponentGE.fugcoef`'s two reference states.** It reads nothing and
prints its whole grid, so it needs the recipe for the same reason. The states cover the IAPWS
arm in water, the databank arm with and without water, the `900` sentinel, the four-zero rows
(which NeqSim evaluates to 1.802 bar), the three solvent-tagged IAPWS species whose
one-component reference phase throws, the ions' 1e12 cap, and the UNIFAC group refusals; the
last block is a solute-bearing gamma-phi flash transcribed from `NrtlGammaPhiFlash.java`.
`gamma_inf_water_synced` is the UNIFAC reference phase rebuilt with its groups and `calcaij`,
so the raw path can be checked against the intended one.

```bash
cd validation/neqsim
../../tools/gated.sh javac -proc:none -cp neqsim-f0c7436.jar GeHenryProbe.java
../../tools/gated.sh java -cp .:neqsim-f0c7436.jar GeHenryProbe > captures/ge_henry_probe.tsv
```

**`BrineProbe` prints the electrolyte-CPA phase's layers and then runs the four brine
operations.** It reads nothing and covers its whole grid, so it needs the recipe for the same
reason. Two things in it are measured rather than assumed: the two halves of
`PhaseElectrolyteCPA`'s sum are printed separately (`FCPA` and `F_super`), and the statoil RDF
is printed beside the interface default's value. **Its reactive fluids call
`chemicalReactionInit()` before `setMixingRule`** - the reaction tables grow the system from
five components to nine, and the electrolyte mixing rule's `wijCalcOrFitted` is built at the
count of the moment, so the other order throws `Index 4 out of bounds for length 4` from
`calcWij`.

```bash
../../tools/gated.sh java -cp .:neqsim-f0c7436.jar BrineProbe > captures/brine_probe.tsv
```

**`ReactiveBrineProbe` carries what that capture does not: the nine species and the residuals.**
`BrineProbe`'s reactive block prints the component *count* and never the names - and a
`removeJunkReactions` reconstruction from the reaction tables gives eight, where
`chemicalReactionInit()` grows the brine to nine. It prints the residuals the coupling is
accepted on as well, because `ReactiveCO2BrinePhaseEquilibrium.run` reads its reaction residual
from the pass *before*, and on the first pass that accessor answers `NaN` - which is why a pass
that has solved nothing cannot be accepted.

```bash
../../tools/gated.sh java -cp .:neqsim-f0c7436.jar ReactiveBrineProbe > captures/reactive_brine_probe.tsv
```

```bash
../../tools/gated.sh java -cp .:neqsim-f0c7436.jar PcsaftProbe 150 50 methane 1 > captures/pcsaft_probe_methane_150_50.tsv
../../tools/gated.sh java -cp .:neqsim-f0c7436.jar PcsaftProbe 300 100 propane 1 > captures/pcsaft_probe_propane_300_100.tsv
../../tools/gated.sh java -cp .:neqsim-f0c7436.jar PcsaftProbe 400 50 methane 1 > captures/pcsaft_probe_methane_400_50.tsv
```

**`ProcessProbe` takes a unit operation's name and prints that unit operation's inlet and
outlet on the palette's record**, so one recipe per unit operation and one capture each. It
is the only probe here that drives `neqsim.process` rather than `neqsim.thermo`.

```bash
../../tools/gated.sh java -cp .:neqsim-f0c7436.jar ProcessProbe pump > captures/process_pump.tsv
../../tools/gated.sh java -cp .:neqsim-f0c7436.jar ProcessProbe splitter > captures/process_splitter.tsv
../../tools/gated.sh java -cp .:neqsim-f0c7436.jar ProcessProbe mixer > captures/process_mixer.tsv
../../tools/gated.sh java -cp .:neqsim-f0c7436.jar ProcessProbe separator > captures/process_separator.tsv
../../tools/gated.sh java -cp .:neqsim-f0c7436.jar ProcessProbe throttling_valve > captures/process_throttling_valve.tsv
../../tools/gated.sh java -cp .:neqsim-f0c7436.jar ProcessProbe heat_exchanger > captures/process_heat_exchanger.tsv
../../tools/gated.sh java -cp .:neqsim-f0c7436.jar ProcessProbe stream > captures/process_stream_properties.tsv
../../tools/gated.sh java -cp .:neqsim-f0c7436.jar ProcessProbe column > captures/process_column.tsv
../../tools/gated.sh java -cp .:neqsim-f0c7436.jar ProcessProbe column_solvers > captures/process_column_solvers.tsv
../../tools/gated.sh java -cp .:neqsim-f0c7436.jar ProcessProbe column_efficiency \
  > captures/process_column_efficiency.tsv
../../tools/gated.sh java -cp .:neqsim-f0c7436.jar ProcessProbe column_divergence \
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
../../tools/gated.sh java -cp .:neqsim-f0c7436.jar ProcessProbe absorber_efficiency \
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
../../tools/gated.sh java -cp .:neqsim-f0c7436.jar ProcessProbe column_tear \
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
../../tools/gated.sh java -cp .:neqsim-f0c7436.jar ProcessProbe packed_column \
  > captures/process_packed_column.tsv
```

**Four things those keys pin, and each was unread before.** `isDesignOk()`, `isWettingOk()` and
the column's `isHydraulicsOk()` are three predicates that all read as "the hydraulics are fine",
and on these rows **all three are false**, so the row alone cannot tell them apart - which is
what makes publishing one under another's name a way to reproduce every row while answering a
different question. `getPressureDropPerMeter()` beside the column's `getPackingPressureDrop()`
settles a name that does not say which it is: on the 2.3 m row they are `0.08962289003108304` and
`0.20613264707149095`, **a ratio of exactly the bed height**, so the column's is the total.
`resolved_packing` and `resolved_packing_factor` (`Pall-Ring-50`, `66.0`) and
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
../../tools/gated.sh java -cp .:neqsim-f0c7436.jar PackingProbe > captures/packing_probe.tsv
```

**`PackingProbe size` is the second capture it writes, and the one the sizing path needs.** The
four rows above *set* a diameter, and `PackedColumn`'s own capture reaches sizing only through a
solved middle tray, so `sizeColumnDiameter` - the trial flood at **1.0 m**, the design velocity
from `designFloodFraction`, the area from the vapour's volumetric flow, and the round up to the
class's standard-diameter table - was measured by nothing. Five states: the absorber, structured
packing, a high liquid load, a light load that reaches the table's **0.3 m floor**, and the
absorber again at a `0.5` flood fraction. **That last pair is the measurement of the one claim** -
the flood fraction moves the sized diameter from `0.3` to `0.4` m and leaves the calculation
untouched - and **the light and high-liquid rows share a flooding velocity to the last digit**,
which is the Eckert fit's `FLV` clamped at its ceiling of `5` on both.

```bash
../../tools/gated.sh java -cp .:neqsim-f0c7436.jar PackingProbe size \
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
../../tools/gated.sh java -cp .:neqsim-f0c7436.jar TrayHydraulicsProbe > captures/tray_hydraulics_probe.tsv
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
../../tools/gated.sh java -cp .:neqsim-f0c7436.jar ProcessProbe capacity > captures/process_column_capacity.tsv
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
../../tools/gated.sh java -cp .:neqsim-f0c7436.jar ProcessProbe designer > captures/process_internals_designer.tsv
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
../../tools/gated.sh java -cp .:neqsim-f0c7436.jar ProcessProbe coupling > captures/process_hydraulic_coupling.tsv
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
the height and the wall thickness. Four rows: the binary column this port already reproduces, the
same column at a stated `2.0` m diameter override, the same column with a `valve` tray, and the
lean-oil absorber. **No capture in this directory held the class**, which is 1,348 lines with a
public surface of forty-odd getters.

```bash
../../tools/gated.sh java -cp .:neqsim-f0c7436.jar ProcessProbe mechanical > captures/process_column_mechanical_design.tsv
../../tools/gated.sh java -cp .:neqsim-f0c7436.jar ProcessProbe reactor_mechanical \
  > captures/process_reactor_mechanical_design.tsv
../../tools/gated.sh java -cp .:neqsim-f0c7436.jar ProcessProbe gibbs_algorithms \
  > captures/process_gibbs_algorithms.tsv
```

**`calcDesign` runs two sizings and publishes the second one's diameter beside the first one's
numbers, and the override row is what makes that a measurement rather than a reading of the
source.** The Souders-Brown pass fixes the diameter, recomputes a flooding factor at it and takes
a weir loading from `0.7 D`; then `calculateContactorCapacity` builds the internals designer,
hands it `resolveRatingDiameter` as an override, and **replaces the diameter with the designer's**
- while the flooding factor and the weir loading keep the value they were computed at. So the
`2.0` m row reports `column_diameter_m=2.0` beside a `flooding_factor` of
`0.029981966003039984` and a `weir_loading` of `4.992709047038033`, **which are the defaults
row's, to the last digit**, while `column_wall_thickness_m` moves to `864.7798742138365`. The
`valve` row moves the Souders-Brown `kFactor` from `0.1` to `0.12` and the designer's internals
type with it, which is one tray type answering in two places.

**Four keys were added to every row so the class's own units are in the capture rather than in
this entry.** The `_m`, `_Pa` and `_W` labels on the block *below* them are the probe's and are
wrong for three of the quantities; `column_wall_thickness_mm`, `total_pressure_drop_bar`,
`tray_pressure_drop_mbar_per_tray` and the two `_kw` keys state the same numbers in the unit each
getter's javadoc gives. **And three keys name what the designer inside `calcDesign` answered**
(`internals_required_diameter_m`, `internals_max_percent_flood`,
`internals_total_pressure_drop_pa`), so the two halves of the call are separable in the record
rather than only in the source.

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

**The plan's second item, `AbsorberMechanicalDesign`, is not a port target, and that is measured
rather than assumed.** Its only constructor call site in `src/main` is
`SimpleAbsorber.getMechanicalDesign()` (`SimpleAbsorber.java:383`) - and `SimpleAbsorber extends
Separator`, so driving the class on an `AbsorptionColumn` (a `DistillationColumn`) throws
`ClassCastException` from `PressureVesselDesignStandard.java:36` through
`SeparatorMechanicalDesign.readDesignSpecifications:271`. It is the separator family's sizing for
a machine this port already refuses, which is why `process.absorption_column`'s mechanical design
is the *column's* own design and nothing here has a row for it.

**Three quantities of the same class are *not* defects, and an earlier version of this entry said
they were.** `getTotalPressureDrop()` is in **bar**, `getReboilerDuty()` in **kW** and
`getColumnWallThickness()` in **mm**, each as its javadoc says; the `_Pa`, `_W` and metre labels
beside them in this capture are **the probe's**, and reading them as the class's units is what
produced the false claim. **And the `216.19496855345912` mm wall is the probe's own artefact**:
`DistillationColumn` sets the design pressure from `getEconomicDesignPressure()` before it calls
`calcDesign`, and this row constructs the design directly and so sizes against the base class's
`100.0` bara initialiser - a property of the measurement, not of the class. Both corrections are
in the issue's own "what I am not claiming" section.

**A fifth row is the stripping column's own**, and it is the same class a second time.
`StrippingColumn extends AbsorptionColumn extends DistillationColumn`, so it inherits
`DistillationColumnMechanicalDesign` and overrides no getter of it. The row is built from the same
recipe as `process_absorber.tsv`'s `hydrocarbon_stripper` row - five trays at 12 bara, 150 kg/hr of
stripping gas against 900 kg/hr of rich liquid, both at 70 C, the `1e-4` K gate and direct
substitution - so the two captures describe one state rather than two copies of one. It carries
the same `#4140` read: `tray0_liquid_density_kg_per_m3=0.0`, `weir_loading=Infinity` and
`tray_pressure_drop_Pa=5.0`, which one `initProperties()` moves to `631.6396155323132`,
`3.480252186825798` and `8.098192314185997`. The vessel is `0.5` m by `8.8` m at
`216.19496855345912` mm - the base class's `100.0` bara initialiser, the same artefact the binary
rows carry - with `8` actual trays at the `0.65` efficiency and both duties zero, because a
stripper has no ends.

**`ProcessProbe reactor_mechanical` is `ReactorMechanicalDesign`, which is `GibbsReactor`'s own
vessel.** `GibbsReactor` overrides `getMechanicalDesign()` (`GibbsReactor.java:700`) and lazily
builds it; `PlugFlowReactor` overrides nothing, so its `getMechanicalDesign()` is
`ProcessEquipmentBaseClass`'s, whose `calcDesign()` only reads design specifications and computes
no vessel - **there is nothing there to capture.** The reactor class reads the **inlet** stream of
a `TwoPortInterface` and nothing else: its pressure, temperature, volumetric flow
(`getFlowRate("m3/hr")`) and the fluid's system density (`getDensity("kg/m3")`), then a catalyst
volume from a space velocity, a vessel from a bed fill fraction and a target L/D, an ASME VIII
Div 1 shell and head, an Ergun bed drop and a weight. **No capture in this directory held the
class.**

**The three rows separate every branch that changes an answer.** The defaults row is the ammonia
feed at 300 bara, where `flow/GHSV` floors the catalyst volume at `0.01` m3 and the shell is
`69.85605419136324` mm at a `330.0` bara design pressure. The CSTR row takes the LHSV branch,
which is the only place `reactorType` is read, and its `1.5` 1/hr keeps the vessel **above** the
`0.5` m diameter floor - at every default-GHSV flow the floor answers, and a capture that only
ever reads `0.5` does not measure the cube root: this row reports
`vessel_diameter_m=0.5385511691956425` and `vessel_length_m=2.15420467678257`. The high-pressure
row feeds 1 bara with `maxOperationPressure` stated as `150.0`, so `max(inlet, maxOperation)`
takes the max branch (`165.0` bara design) and its 299 m3/hr over the default `5000` keeps the
GHSV division above the floor on a gas whose density is `0.2767` kg/m3, where the Ergun inertial
term is a third of the answer.

**Two branches cannot answer, and the capture is what makes that a measurement.**

- **The head thickness is always the shell's.** The 2:1 ellipsoidal form divides by
  `2SE - 0.2p` and the shell by `2SE - 1.2p`, so at every positive design pressure the head is
  thinner and `Math.max(head, shell)` reports the shell. `head_thickness_mm` equals
  `shell_thickness_mm` on all three rows, to the last digit.
- **The L/D clamps are unreachable.** The diameter is built at the target `4.0` and the length is
  `4.0 * D`, so the ratio is exactly `4.0` and neither the `2.0` floor nor the `8.0` ceiling can
  fire on any state the class can reach.

**Five configured parameters are setters with no getters.** `ghsv`, `lhsv`,
`catalystBulkDensity`, `catalystParticleDiameter` and `bedVoidFraction` therefore appear as
`requested_` keys - the probe's statement of what it configured, with the class's own
initialisers on the defaults row - and are not a reading of the class. The six that do have
getters are printed as `class_` keys, and `calcDesign`'s own `super.calcDesign()` leaves
`joint_efficiency=1.0` and `corrosion_allowance=0.0` at the pin, which is why those two are class
defaults and not design-standard output.

**`ProcessProbe gibbs_algorithms` is the Gibbs solver's two step controls, on the states
`GibbsReactorAlgorithmTest` solves.** `useArmijoLineSearch` and `useAdaptiveStepSize` are both
`false` by default, so no row of `process_gibbs_reactor.tsv` reaches either; and
**`useAdaptiveStepSize` is named nowhere in this tree's record** - the rules below are its first
statement. Five rows, each printing the whole path rather than the class's own "it converged":
the iteration count, the final error, the outlet, and both diagnostic histories.

**The Armijo row does not backtrack, and printing the step history is how that is known.** The
class's own Armijo test caps the step at its `damping` of `0.01`, and every one of the row's 100
steps is exactly `0.01` - the line search accepts `alpha_max` on every pass. Its baseline
(`fixed_methane_adiabatic`, the same feed, damping and gate with the search off) agrees to `1e-15`
on the whole Gibbs history, which is the class's own assertion made measurable rather than assumed.

**The adaptive pair separates the two branches, and on this state Armijo adds nothing.**
`adaptive_hydrogen` and `armijo_adaptive_hydrogen` carry step histories that differ in the
fifteenth digit: the adaptive limit answers `1.0, 1.0, 0.0193, 0.0247, 0.0127, ...` and the line
search accepts each of them.

**The row that does backtrack is the undamped one.** `armijo_ammonia_undamped` starts from a full
Newton step (`damping = 1.0`) and its history reads `1.0, 1.0, 1.0, 1.0, 1.0, 0.5, ... 0.001953125
... 0.03125, ...` - contractions by `armijoRho = 0.5`, nine of them in one pass. It converges at
the `min_iterations` floor with a final error of `7.0e-16`, **and the class's own
`getMassBalanceConverged()` reports false on the row** - a second stopping flag this port does not
carry, printed so the capture states it.

**`WaterCpSentinel` asks a question about the *data* rather than about a model.** `COMP.csv`
gives 131 of its 389 rows the whole of water's ideal-gas Cp polynomial - the same five numbers
`devtools/generate_water_caloric_alpha_reference.py` fits for water - and 130 of those rows are
substances other than water. This drives one single-component `SystemSrkEos` per name in the
table and compares `getCp0` against water's **exactly**, because two different molecules cannot
share an ideal-gas heat capacity: an exact match is not a tolerance question, and a near one
would be evidence of nothing. It parses the `NAME` column quote-aware, since names carry commas
inside them, and its capture is the sorted list of the names that match.

```bash
../../tools/gated.sh java -cp .:neqsim-f0c7436.jar WaterCpSentinel > captures/water_cp_sentinel.tsv
```

**`ProcessProbe rate_based` runs every state on two cubics.** The probe states an iteration cap
of 60 rather than the class's 20 because the packing fix moves the stripper to 21 passes and the
TEG profile to 33; every row's converged state is independent of the cap. `RateBasedPackedColumnTest` uses
SRK and this library's process layer resolves PR, so each state is run on both and the pair
measures what the cubic moved; the PR rows are the ones the port is held to. `teg_dehydration`
and `teg_circulation` are `SystemSrkCPAstatoil` with mixing rule 10, which PR cannot
represent, so they are evidence for a refusal rather than an oracle for a state.

```bash
../../tools/gated.sh java -cp .:neqsim-f0c7436.jar ProcessProbe rate_based \
  > captures/process_rate_based_packed_column.tsv
```

**`ProcessProbe rate_based_billet` runs the same absorber state twice at each
`MassTransferCorrelation`, on two packings.** That parameter is not a correlation: it scales the
two film coefficients by constants the packing row carries. On the built-in `Pall-Ring-50` those
constants are `1.0` and `0.40`, so the pair is the identity; on `Pall-Ring-50-plastic` (`cp =
0.698`, `ch = 2.725`) it is `1.1354166...` on `kGa` and `0.698` on `kLa` - the state the pin's
`Pall-Ring-50` resolved to. The capture's own two ratios are what the port's test reads, and the
liquid one is `0.6978956` rather than `0.698` because the base `kLa` moved `1.495e-4` between the
two converged states.

```bash
../../tools/gated.sh java -cp .:neqsim-f0c7436.jar ProcessProbe rate_based_billet \
  > captures/process_rate_based_billet.tsv
```

**`ProcessProbe rate_based_solvers` measures the two solvers this port refuses.** It runs the
class's own absorber state - the PR row `process_rate_based_packed_column.tsv` holds and the port
reproduces to `1e-4` - five ways: the pair this port carries as the control, then
`SegmentSolver.SIMULTANEOUS_RESIDUAL` and `ColumnSolver.EQUATION_ORIENTED` alone, then both, then
`EQUATION_ORIENTED` at the settings of the class's own test for it. **No refused row converges**,
and the control reproduces the baseline exactly (fifteen passes at `6.770219315721526e-10` mol/s),
which is what makes the rest a reading rather than a mis-set-up state. The segment solve takes the
full 20-pass cap at a `0.25` mol/s outlet residual against a `1e-9` gate and publishes a vapour at
`4.3e-31` K; the column-wide Newton stalls after two passes at a residual of `118.6` against its own
`1e-6` gate and publishes a `300` K vapour beside a `463` K liquid, and at its own test's settings
the norm is `30.5` against `1e-5`. **That test asserts finiteness and never convergence** - it caps
the Newton at two iterations and removes the homotopy - and its companion caps them at one.

```bash
../../tools/gated.sh java -cp .:neqsim-f0c7436.jar ProcessProbe rate_based_solvers \
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

**The scale family's three probes** record the water-analysis scale calculations: the
saturation-index calculator with every intermediate it exposes (Ksp, ionic strength,
Debye-Huckel A, the free ions, the five SIs) across the upstream tests' brines and both pH
paths; the binary Pitzer model with its raw eight-salt parameter table read by reflection and
the 29-point Hamer-Wu NaCl curve; and the solid solution, mass calculator, multi-mineral
equilibrium in three activity modes, flowline profile and water-compatibility screener. The
last two compose the first, and both inert pH arguments are proven inert by running the same
state twice.

```bash
../../tools/gated.sh java -cp .:neqsim-f0c7436.jar ScalePredictionProbe > captures/scale_prediction_probe.tsv
../../tools/gated.sh java -cp .:neqsim-f0c7436.jar PitzerScaleActivityProbe > captures/pitzer_scale_activity_probe.tsv
../../tools/gated.sh java -cp .:neqsim-f0c7436.jar ScaleFamilyProbe > captures/scale_family_probe.tsv
```

**The hydrate group's two probes** record the Pitzer-brine saturation the hybrid EoS/GE
flash exists for: the temperature solve and the pressure solve on CO2/water, CO2 with NaCl,
MgCl2 and CaCl2 brines, methane and a CO2/methane mixture; the structure pinned to auto, sI
and sII; the public `ThermodynamicOperations` round trip; and the two refusals a feed without
a guest and one without water produce. The second probe walks the salt suppression curve over
the same three salts. The fluid is the upstream test's one-kilogram water basis, with the
PHREEQC CO2-chloride parameters at its zero zeta map.

```bash
../../tools/gated.sh java -cp .:neqsim-f0c7436.jar PitzerHydrateProbe > captures/pitzer_hydrate_probe.tsv
../../tools/gated.sh java -cp .:neqsim-f0c7436.jar PitzerHydrateSalinityProbe > captures/pitzer_hydrate_salinity_probe.tsv
```

**The process-scale group's probe** records the five process-chemistry scale calculators:
the activity-corrected `ElectrolyteScaleCalculator` over the upstream brucite seawater and
advanced-models brines, the pH and temperature grid, the Debye-Hueckel floor, the absent-ion
state and the 50-bar pressure correction; the `BrineMixingScaleEvaluator` over the two
upstream formation/seawater pairs and its steps floor; the `ScaleInhibitorPerformance` MIC
grid with every warning and the efficiency on both sides of the MIC; the `ScaleKinetics`
supersaturation grid with the transport and reaction branches and the clamped setters; and
the `ScaleControlAssessor` over the RCA and skill states, the absent-scale NaN, the zero-dose
identity and the 0.999 efficiency cap.

```bash
../../tools/gated.sh java -cp .:neqsim-f0c7436.jar ProcessScaleProbe > captures/process_scale_probe.tsv
```

**The stream analyser probe** records `StreamScaleAnalyzer` on the equilibrium capture's own
brines, so the two layers can be compared directly: the sulphate-limited water in Davies and
B-dot, the co-precipitating water, and the corners - zero flow, no minerals, an absent
barium, and the auto-pH path, where the class's broken Henry constant removes the carbonate
row rather than moving it. Every row prints the per-mineral initial and final indices, the
precipitated mg/L and the kg/day rate, with the total and the dominant mineral.

```bash
../../tools/gated.sh java -cp .:neqsim-f0c7436.jar StreamScaleAnalyzerProbe > captures/stream_scale_analyzer_probe.tsv
```

**The valve drift probe** records `ValveScaleDrift`: a uniform radial deposit on a trim
port, its open-area fraction, fouling fraction and effective Kv, and the time to plug and
to pin the opening for a stated clean opening. The valve is constructed without a stream
and given a Kv directly - nothing in the drift arithmetic reads a fluid. The states are
the not-advanced, thirty-day, one-year and pinned deposits, zero growth, an unset port, a
clean opening already at 100 percent, a negative advance the class ignores, a reset, and a
growth rate taken from the evaluated `ScaleKinetics`.

```bash
../../tools/gated.sh java -cp .:neqsim-f0c7436.jar ValveScaleDriftProbe > captures/valve_scale_drift_probe.tsv
```

**The deposition probe** records `ScaleDepositionAccumulator` over a real Beggs-and-Brills
pipe: the pipe's own segment temperatures and pressures (which are the state the port's model
takes directly), then the accumulator's per-segment CaCO3 index, deposit mass and wall
thickness, its total, its worst thickness and its time to half-block. Two profiles are run -
an isothermal one whose temperature list is a single value, so the class evaluates one
segment over the whole length, and a cooled one with eleven segments - across the inhibitor,
five-year service, doubled base rate and no-calcite corners.

```bash
../../tools/gated.sh java -cp .:neqsim-f0c7436.jar ScaleDepositionProbe > captures/scale_deposition_probe.tsv
```

**The treatment-scenario probe** records `ProductionChemicalScaleScenario`: caustic, MDEA and
HCl pH adjusters, two triazine scavenger doses, a scale and a corrosion inhibitor, no
bicarbonate, a large caustic dose that trips the chemical-induced calcite warning, the B-dot
and Pitzer activity models, and an unsupported active. Each row prints the treated pH and
bicarbonate, the dose equivalents, the scavenger capacity and residual H2S, the baseline and
treated indices and their change, and the warning map with its messages.

```bash
../../tools/gated.sh java -cp .:neqsim-f0c7436.jar ProductionChemicalScenarioProbe > captures/production_chemical_scenario_probe.tsv
```

**The Beggs-and-Brills probe** records the pipe's per-segment state: the pressure, the flow
regime, the liquid holdup, the three superficial velocities, the no-slip density and
viscosity, the Reynolds number and the pressure drop, plus the outlet and the total drop.
The states are a horizontal wet-gas line at two diameters, an uphill and a downhill leg, a
liquid-rich line, a three-phase oil/water/gas line, single-phase gas and liquid, a rough
wall, and the three liquid-rich or high-flow states that reach the INTERMITTENT and
DISTRIBUTED branches the others do not; one high-flow state is refused upstream for a
negative outlet pressure and is recorded as such.

```bash
../../tools/gated.sh java -cp .:neqsim-f0c7436.jar PipeBeggsAndBrillsProbe > captures/pipe_beggs_and_brills_probe.tsv
```

**The closed-loop deposition probe** records `ClosedLoopDepositionSolver` over the same pipe:
the iteration count, convergence, the final effective diameter, the last accumulator pass and
the diameter, thickness and mean-velocity histories. The states are a one-year service that
meets the tolerance on its first step, a thousand-year service that iterates three times, the
iteration cap reached before convergence, a blocked line, a brine that deposits nothing, an
inhibited line, a tolerance wider than the step, and a Peng-Robinson fluid.

```bash
../../tools/gated.sh java -cp .:neqsim-f0c7436.jar ClosedLoopDepositionProbe > captures/closed_loop_deposition_probe.tsv
```

**The flow-assurance remainder probe** records the six reachable calculators of
`pvtsimulation/flowassurance` outside the scale set: de Waard–Milliams CO2/H2S/S8 corrosion,
the API RP 14E / DNV RP O501 erosion screen, the six emulsion-viscosity models and the inversion
point, the lumped pipeline cooldown march (supplied U and layer resistances, a step that does not
divide the horizon, target crossings), the hydrate risk map over four profile points (including a
dry fluid and a zero-water fluid), and the surface cooldown no-touch analyser (required-time
bands, an operating-temperature override, a seabed above the hydrate temperature, and the
no-water refusal).

```bash
../../tools/gated.sh java -cp .:neqsim-f0c7436.jar FlowAssuranceRemainderProbe > captures/flow_assurance_remainder_probe.tsv
```

**The characterization fluid-rebuild probe** records `PseudoComponentCombiner`'s two fluid
builders and its binary-interaction transfer: two and three unequal reservoir fluids combined at
a target pseudo-component count, a target beyond the pool, base-only fluids and a base-only fluid
beside a pseudo one; the common slate at equal and weighted settings, an inferred and an explicit
count, a zero-weight fluid and the no-pseudo clone path; and the BIP transfer from a reference to
a target with different pseudo names, fewer pseudo components, and a missing base name. Each
component is one wide row (`fluid|name|pseudo|` then the nineteen pseudo fields), and the BIP
rows print the reference's and the target's initial matrices beside the transferred one. Every
source is a real TBP fraction, so the properties are finite and the TBP model's own defaults
cannot leak into a row.

```bash
../../tools/gated.sh java -cp .:neqsim-f0c7436.jar CharacterizationFluidRebuildProbe > captures/characterization_fluid_rebuild_probe.tsv
```

**The TBP fraction component probe** records what `SystemThermo.addTBPfraction` makes of one
cut: the component's own getters after the call, and the intermediates that are not on it - the
model's Edmister acentric factor and the fitted `m`, and the molar volume of the one-component
SRK reference phase `calcRacketZ` flashes at 15 C and one atmosphere. Six cuts span the
boiling-point switch at 540 g/mol and the heavy coefficient switch at 1120, and the last one
pins the class's clamp on a fitted `m` whose root runs past -0.99.

```bash
../../tools/gated.sh java -cp .:neqsim-f0c7436.jar TbpFractionComponentProbe > captures/tbp_fraction_component_probe.tsv
```

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
