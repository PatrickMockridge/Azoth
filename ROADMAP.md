# NeqSim port roadmap

What azoth has ported of [NeqSim](https://github.com/equinor/neqsim), what it has not,
and the order the rest is taken in — the physics the most callers need before the
physics the fewest do, so a property ordinary oil, gas and water work uses comes
before one only a specialist reaches for.

The mapping is one to one. Every class under NeqSim's `thermo/`, `thermodynamicoperations/`
and `physicalproperties/` trees is assigned to a tier below, each bullet naming the azoth
id that carries it where it is ported, the NeqSim class that would close it where it is
not, and the measurement where neither applies. The trees that are not physics are named
at the end rather than silently dropped, and they are of two kinds: **`process/`** — unit
operations and the flowsheet — is not physics but is a port target in its own right, the
specification putting it at tranche P11 and P12; while **`fluidmechanics/`** (azoth has
its own hydraulics) and the support packages are not ported at all. The column-by-column
record is `databank/manifest.toml`, which this file orders; it does not re-inventory it.

## The rules every port follows

[The specification](docs/src/architecture/specification.md) is normative. The three
rules a port follows, and the tie-breaker, are there; in short:

- a port names the paper for the method, the implementation for the port, and its own
  notes for the changes;
- a port never upgrades `verification.status` — reading NeqSim's Java is not reading
  the paper, so what is ported ships `unverified`;
- a port is accepted on azoth's tests, never on its provenance.

Where azoth and NeqSim would otherwise disagree, NeqSim wins. And nothing here is
"out of scope": a physics family NeqSim implements is *not ported*, with the class that
would close it named — never "out of scope".

## The tiers

The tiers are the order the port is taken in, and each carries what is ported and what
is not, so that "azoth does not do this" and "azoth has not reached this" stay different
statements.

### Tier 0 — the one data path

The cubic Peng-Robinson core is ported: `ComponentPR`, `AttractiveTermPr`,
`PhasePrEos`, `SystemPrEos`, `TPflash`, `PHflash`, `PSFlash`, the tangent-plane
stability test, `CriticalPointFlash`, the phase envelope, `RachfordRice` and the
classical vdW1f mixing rule, reading 14 of the databank's columns. What is still open on
it is what the rest of the port stands on:

- the single `databank → keycard → every calculation` path, so a calculation names its
  components and every constant it reads comes through one route;
- the `not-yet` columns the manifest names — the collision and liquid-viscosity
  constants `eos.viscosity` and `eos.thermal_conductivity` read. `eos.molar_enthalpy_entropy`
  already reads the databank, and `eos.ideal_gas_cp` stays a scalar calc by design.

### Tier 1 — oil and gas

The cubic family in breadth, the activity-coefficient models, the standard component
correlations, petroleum-fraction characterisation, and the transport properties every
process model needs.

- **Cubic EOS breadth.** `ComponentSrk`, `ComponentRK`, `ComponentPRvolcor`,
  `ComponentSrkvolcor`, `ComponentSrkPeneloux`, `ComponentCSPsrk`, `ComponentTST`,
  `ComponentBWRS`, `ComponentBNS`, `ComponentAmmoniaEos`, and the alpha terms
  `AttractiveTermSrk`, `AttractiveTermRk`, `AttractiveTermPr1978`, `AttractiveTermTwu`,
  `AttractiveTermTwuCoon`, `AttractiveTermTwuCoonParam`, `AttractiveTermTwuCoonStatoil`,
  `AttractiveTermMatCop`, `AttractiveTermMatCopPR`, `AttractiveTermMatCopPRUMR`,
  `AttractiveTermMatCop5PRUMR`, `AttractiveTermMollerup`, `AttractiveTermPrLeeKesler`,
  `AttractiveTermPrDanesh`, `AttractiveTermPrDelft1998`, `AttractiveTermPrGassem2001`,
  `AttractiveTermSchwartzentruber`, with their `Phase*` and `System*` counterparts
  (`SystemSrkEos`, `SystemSrkPenelouxEos`, `SystemSrkTwuCoonEos`,
  `SystemSrkMathiasCopeman`, `SystemPrMathiasCopeman`, `SystemPrEos1978`,
  `SystemPrLeeKeslerEos`, `SystemPrDanesh`, `SystemPrGassemEos`, `SystemPrEosDelft1998`,
  `SystemSrkSchwartzentruberEos`, `SystemCSPsrkEos`, `SystemRKEos`, `SystemBWRSEos`,
  `SystemBnsEos`, `SystemAmmoniaEos`, `SystemTSTEos`).
- **Mixing rules.** `EosMixingRuleHandler` (classic, Huron-Vidal, Wong-Sandler),
  `HVMixingRulesInterface`, `MixingRuleHandler`.
- **Activity-coefficient / GE models.** `ComponentGeNRTL`, `ComponentGENRTLmodifiedHV`,
  `ComponentGENRTLmodifiedWS`, `ComponentGEUniquac`, `ComponentGEUniquacmodifiedHV`,
  `ComponentGEWilson`, `ComponentGEUnifac`, `ComponentGEUnifacPSRK`,
  `ComponentGEUnifacUMRPRU`, `ComponentGEVanLaarAcid`, with `PhaseGENRTL`,
  `PhaseGEUniquac`, `PhaseGEWilson`, `PhaseGEUnifac`, `SystemNRTL`, `SystemUNIFAC`,
  `SystemUNIFACpsrk`, `SystemGEWilson`.
- **Component correlations** read from `thermo/component/` (normal boiling point,
  Antoine vapour pressure) and `physicalproperties/` (standard liquid density, heat of
  vaporisation, liquid heat capacity) — the columns the manifest lists under those
  packages.
- **Petroleum-fraction characterisation.** `thermo/characterization/`: `Characterise`,
  `PlusCharacterize`, `TBPCharacterize`, `TBPfractionModel`, `PlusFractionModel`,
  `PedersenPlusModelSolver`, `LumpingModel`, `LumpingConfigBuilder`,
  `PseudoComponentCombiner`, `Recombine`, `OilAssayCharacterisation`, `TbpClosure`,
  `RefineryAssayBlend`, `BioFeedstock`, `BiomassCharacterization`.
- **Transport properties.** `physicalproperties/methods/`: gas and liquid viscosity
  (`ChungViscosityMethod`, `Viscosity`), conductivity (`ChungConductivityMethod`,
  `FilippovConductivityMethod`, `Conductivity`), density (`Rackett`, `Costald`,
  `Density`), diffusivity (`WilkeChangDiffusivity`, `HaydukMinhasDiffusivity`,
  `SiddiqiLucasMethod`, `TynCalusDiffusivity`, `WilkeLeeDiffusivity`,
  `FullerSchettlerGiddingsDiffusivity`), and surface tension (`ParachorSurfaceTension`,
  `GTSurfaceTension`, `FirozabadiRamleyInterfaceTension`, `LGTSurfaceTension`,
  `CDFTSurfaceTension`).
- **Flash breadth.** `TVflash`, `PVflash`, `VUflash`, `THflash`, `TSflash`, `TUflash`,
  `PUflash`, `PVFflash`, `TVfractionFlash`, `QfuncFlash`, `ImprovedVUflashQfunc`.
- **Phase envelopes.** `PTphaseEnvelope`, `HPTphaseEnvelope`, `PTPhaseEnvelopeMichelsen`,
  `CricondenBarFlash`, `CricondenThermFlash`, `SysNewtonRhapsonPhaseEnvelope`.
- **Flow assurance.** `pvtsimulation/flowassurance/` is carried as **one named group — 21
  classes and 12,486 lines** — because it is not one family: seven scale calculations
  (`MultiMineralScaleEquilibrium`, `ScalePredictionCalculator`, `ScaleMassCalculator`,
  `PitzerScaleActivityModel`, `BariteCelestiteSolidSolution`, `FlowlineScaleProfile`,
  `WaterCompatibilityScreener`), six asphaltene screens, one wax (`WaxCurveCalculator`), two
  corrosion (`CO2CorrosionAnalyzer`, `DeWaardMilliamsCorrosion`) and five generic
  flow-assurance calculators (`HydrateRiskMapper`, `PipelineCooldownCalculator`,
  `SurfCooldownAnalyzer`, `ErosionPredictionCalculator`, `EmulsionViscosityCalculator`).
  None of it is ported. **The rest of `pvtsimulation/` is not here**: the simulation suite,
  the model tuning and the reservoir properties are engineering deliverables, and
  [the specification](docs/src/architecture/specification.md) puts them beyond the port.

### Tier 2 — water and gas-water

Water and the aqueous models: the water EOS, water content and dehydration, acid-gas
solubility, electrolytes, salts and scale, freezing, and the reactions an aqueous phase
reaches equilibrium through.

- **The water EOS.** `ComponentWater`, `PhaseWaterIAPWS`, `SystemWaterIF97`,
  `thermo/util/steam/Iapws_if97`.
- **Water content and dehydration.** `WATcalc`, `WaterDewPointTemperatureFlash`,
  `WaterDewPointTemperatureMultiphaseFlash`, `WaterDewPointEquilibriumLine`.
- **Acid gas.** `ComponentSoreideWhitson` and `AttractiveTermSoreideWhitson` with
  `PhaseSoreideWhitson` (`eos.soreide_whitson_phase`) and `SystemSoreideWhitson`;
  `ComponentGEVanLaarAcid` with `PhaseGEVanLaarAcid` (`eos.ge_van_laar_acid_phase`).
  **The brine operations are not ported, and the reason this entry gave was wrong.**
  `CO2BrinePhaseEquilibrium`, `ReactiveCO2BrinePhaseEquilibrium`, `SaturateWithWater` and
  `CalcIonicComposition` were said to need the EoS/GE hybrid seam. They do not: their gate is
  `fluid instanceof SystemElectrolyteCPAstatoil`, which extends `SystemFurstElectrolyteEos`
  extends `SystemSrkEos` — so the hybrid route, whose dispatch needs a `SystemEosGE`, is not in
  the chain at all. What they need is the **electrolyte CPA** route and their own 576 lines:
  `CO2BrinePhaseEquilibrium`'s supported composition is CO2, water and ions alone with
  `water > CO2`, and `SaturateWithWater` and `CalcIonicComposition` are operations of their own.
  **The seam those four were recorded as waiting on is ported** — `eos.hybrid_eos_ge_flash` and
  `reactions.reactive_hybrid_eos_ge_flash` — and it was never the blocker here.
- **Electrolytes.** Ported, each with the id that carries it: `ComponentGePitzer`,
  `PhasePitzer` (`eos.pitzer_phase`) and `SystemPitzer`, with the `Pitzer*` machinery
  (`PitzerNeutralInteraction`, `PitzerElectrostaticMixing`, `PitzerTemperatureFunction`,
  and the parameter catalogs);
  `ComponentKentEisenberg`, `PhaseKentEisenberg` (`eos.kent_eisenberg_phase`) and
  `SystemKentEisenberg`; `ComponentDesmukhMather`, `PhaseDesmukhMather`
  (`eos.desmukh_mather_phase`) and `SystemDesmukhMather`;
  `ComponentModifiedFurstElectrolyteEos`, `PhaseModifiedFurstElectrolyteEos`
  (`eos.furst_electrolyte_phase`) and `SystemFurstElectrolyteEos`, with the same three
  again for Mod2004 (`eos.furst_electrolyte_mod2004_phase`) over the parameters
  `FurstElectrolyteConstants` hardcodes.
- **Reactions.** Ported, each with the id that carries it: `ChemicalReaction`
  (`reactions.equilibrium_constant` - `ln K = K1 + K2/T + K3 ln T + K4 T`, its derivative and
  the heat of reaction, over all three data sources); `ChemicalReactionList`
  (`reactions.reference_potentials` - the greedy independent basis and the potentials
  `sum(nu_i mu_i) = -RT ln K` gives); `ChemicalEquilibrium` (`reactions.chemical_equilibrium`
  - the Smith-Missen Newton solve, with the electroneutrality row among the constraints); and
  `ChemicalReactionOperations` (`reactions.reactive_phase_equilibrium` - that solve as an
  operation on one phase, with the phase search's `-1` skip carried as a result rather than an
  error). The data is NeqSim's own: `element`, `STOCCOEFDATA`, `REACTIONDATA`,
  `REACTIONDATAPITZER` and `REACTIONDATAKENTEISENBERG`, compiled to `data/reactions/`.
  **The element table is what bounds what a fluid can react** - a substance with no row cannot
  enter the balance, and **the glycols have no row at all**, so a P9 inhibitor fluid is not
  reactive. **The three sources are different standard states and are not interchangeable**:
  `CO2water`'s `K1` is `253.235548` in `REACTIONDATA.csv` against `653.705141388` in
  `REACTIONDATAPITZER.csv`, which is why the source is an input rather than a default.
  **Not ported, each with a measured blocker**: the second refinement, which switches
  `useAdaptiveDerivatives` on and **converges on none of the three portable fluids** its oracle
  walks - two iterations, its loop's minimum, and the first refinement's answer left where it
  was; and the **sequential chemical dispatch** inside `TPflash` with the 100-line
  `ChemicalEquilibrium` that loops a system's phases, which are reachable only where
  `isChemicalSystem()` is true - and **every fluid that predicate accepts carries ions**, because
  NeqSim's reaction tables are water chemistry. **That second blocker was measured and it is
  sharper than "it needs ions"**: on a `SystemSrkEos` chemical system the ions are cubic
  components, so those two dispatches take their activity coefficients from a cubic built over
  the ion rows' *filler* critical constants - `IonFillerSensitivityProbe` perturbs them and the
  CO2-water brine's bicarbonate moves by a factor of thirteen at `1.1x`, and to `0.687` mol at
  `2x`. `mixture_of` refuses a cubic over an `ION` by decision, so porting them would rest the
  answer on invented data. **The trace-ion short circuit is ported** (it is a predicate over the
  composition the driver already has), and so is **the RAND solver's ionic branch**, which
  crosses as the three facts only an ionic fluid's caller knows.
  **`Kinetics` is ported as far as it can be pinned**: its rate law (two laws behind a selector)
  and its Krishna-Standart mass-transfer matrix are ids, both oracled on a real fluid, and the
  effective-diffusion assembly the matrix would otherwise need is eight lines whose input - the
  diffusivity model's binary matrix - is not observable from outside its class. `Kinetics` has no
  consumer outside `fluidmechanics/`'s reactive film model, which is not a port target, so its
  two ids have no in-tree caller and that is stated rather than hidden.
- **Duan-Sun, unreachable.** `ComponentGeDuanSun`, `PhaseDuanSun`, `SystemDuanSun` and
  `thermo/util/empiric/DuanSun.java` are **not ported, and no state they accept exists**:
  `SystemDuanSun.addComponent` throws for every name but `CO2`, and the phase it would build
  divides by the moles and molar mass of a component named `water`, which a system that admits
  only `CO2` does not have. **`#3841` closed the last way in**: `SystemThermo`'s model-name
  factory now throws `UnsupportedOperationException` for the conversion rather than building
  the system, with the reason this entry gives, so `new SystemDuanSun(` appears nowhere in
  `src/main` - the carried-as-unreachable case and no longer a reachable one with nothing to
  reproduce (NeqSim issues 3837 and 3839). `validation/neqsim/GeElectrolyteProbe.java` records
  the refusal.
- **Electrolyte-CPA.** `PhaseElectrolyteCPA`, `PhaseElectrolyteCPAAdvanced`,
  `PhaseElectrolyteCPAMM`, `PhaseElectrolyteCPAOld` and their `SystemElectrolyte*` are
  **carried as unreachable upstream**: 5,727 lines with no non-test construction site
  anywhere — each phase is built only by its own `SystemElectrolyte*`, and
  `SystemThermo`'s factory maps `Electrolyte-CPA-EOS` to `CPAstatoil` and never to
  `SystemElectrolyteCPA`. `SystemElectrolyteCPAstatoil` and
  `PhaseElectrolyteCPAstatoil` are **not ported** — 145 lines, eleven `src/main` sites
  that build them, and Fürst plus the Wertheim association this library already has.
- **Salts and scale.** Ported, each with the id that carries it: `CheckScalePotential`
  (`eos.scale_saturation_ratio`, one salt's `IAP/Ksp` against its solubility product) and
  `MultiSaltPrecipitation` over `CalcSaltSatauration` (`eos.salt_precipitation`, one mineral's
  extent, which is the per-mineral half of the complementarity loop). **Not ported**:
  `AddIonToScaleSaturation` and `CalciumSulfatePhaseBoundaryQualification`; the
  `pvtsimulation/flowassurance/` scale set — `MultiMineralScaleEquilibrium`,
  `ScalePredictionCalculator`, `ScaleMassCalculator`, `PitzerScaleActivityModel`,
  `BariteCelestiteSolidSolution`, `FlowlineScaleProfile`, `WaterCompatibilityScreener`; and
  the eleven `process/chemistry/scale/` classes, 3,935 lines, `BrineMixingScaleEvaluator`
  among them.
- **Freezing.** `FreezingPointTemperatureFlash` is ported as `eos.freezing_point`, on **both**
  of its routes. The Helmholtz one needs a `PhaseSolidHelmholtzEos`: `eos.hydrogen_phase` for
  the fluid and `eos.parahydrogen_solid_phase` for the solid, with the solid calibrated to the
  liquid it meets at the triple point. The tabulated one is `ComponentSolid.fugcoef2` over the
  candidate's own melt data, and its residual is the multiphase appearance condition over the
  fluid's phases - a log-sum-exp, where `eos.tp_solid_flash`'s equation is the solid's amount.
  Measured on NeqSim's own LNG test fluid, this reproduces its answer to `3.4e-9` relative at
  5 bara and `8.2e-8` at 50; at 20 bara NeqSim's bracket collapses and this solves anyway, to a
  second root the case records rather than pins. **`FreezingPointTemperatureFlashTR`,
  `FreezingPointTemperatureFlashOld` and `FreezeOut` are unreachable upstream**: none has a
  `new` anywhere in `src/main`, a reference of any kind outside its own file, or a site in
  `src/test`, and the reflection NeqSim does (`this.getClass()...newInstance()`) can only
  reach a class something else constructs. So the legacy pair is not work azoth has skipped,
  and `FreezeOut`'s amount-solve is not a solve that exists.

### Tier 3 — wax, hydrate, hydrogen, asphaltene

The specialist physics.

- **Plus-fraction characterisation.** The wax and asphaltene families' foundation.
  `thermo/characterization/` — 16,204 lines across 39 classes: `PlusFractionModel`
  and `TBPfractionModel`'s eight cut models, `PlusCharacterize`, `PedersenPlusModelSolver`,
  `LumpingModel`, `Recombine`, and the assay and refinery-blend machinery around them. A wax
  or asphaltene fluid is built from TBP and plus fractions, so this is what those models
  stand on. **azoth has the seam and not the subsystem**: `eos.tbp_fraction_properties` is
  Pedersen's `PedersenSRK` cut correlations, so a pseudo-component can be *described* by a
  molar mass and a normal liquid density — but a plus fraction cannot be split into cuts,
  lumped, or recombined into one, and `racketZ` (the Peneloux shift `addTBPfraction` sets from
  a flashed reference system) is carried with it.
- **Hydrate.** Ported, each with the id that carries it: `PhaseHydrate` over
  `ComponentHydratePVTsim`, which is the component model `PhaseHydrate` selects by default,
  through `HydrateFormationTemperatureFlash` (`eos.hydrate_formation_temperature`) and
  `HydrateFormationPressureFlash` (`eos.hydrate_formation_pressure`); `TPHydrateFlash`
  (`eos.hydrate_fraction`), **with the composition taken from both cavities and the material
  balance asserted at the answer**; `HydrateEquilibriumLine` (`eos.hydrate_equilibrium_line`);
  `HydrateInhibitorConcentrationFlash` (`eos.hydrate_inhibitor_concentration`) and
  `HydrateInhibitorwtFlash` (`eos.hydrate_inhibitor_wt`), the dosing pair; and
  `ComponentHydrateGF`, the second fitted component model, behind `hydrate_model = "guo_finch"`
  on all three hydrate ids. The dose pair needed one piece of machinery this library did not
  have: the **phase label** `PhaseEos.init` assigns - a volume ratio against `1.75`, then
  whether a phase's hydrocarbons outweigh its aqueous components - which `eos.hydrate_inhibitor_wt`
  carries because it is the only model here that asks which phase is the aqueous one.
  **Not ported**: `PitzerHydrateFlash` with `ComponentHydratePitzer`, which is reachable through
  the model name. It was recorded as gated on the EoS/GE hybrid seam P8 did not close - **and
  that seam is now closed** (`eos.hybrid_eos_ge_flash`, with the reactive coupling over it in
  `reactions.reactive_hybrid_eos_ge_flash`), so `SystemPitzer`'s two-models-in-one-flash route
  exists here. What is left is the flash itself and the Pitzer hydrate component, which is this
  tier's. `HydrateEquilibriumDiagnostics` is an audit of a
  state rather than a model, and what it asserts belongs in a test. **Carried as unreachable
  upstream**: `ComponentHydrateKluda`, which has no construction site anywhere in `src/main`,
  `ComponentHydrateStatoil` and `ComponentHydrateBallard`, whose only sites are two
  commented-out lines in `PhaseHydrate`, and `OLGAhydrateCurveGenerator`, which only a test
  builds. `process/chemistry/hydrate/` is carried with the process tier.
- **Wax.** Ported, each with the id that carries it: `ComponentWax` (`eos.wax_solid_fugacity`,
  its `fugcoef2`), `TPmultiflashWAX` (`eos.tp_multiflash_wax`) and `PhaseWax`, which is the
  class `SystemThermo` adds when the wax check is on. `eos.tbp_fraction_properties` is
  Pedersen's cut correlations, which `WaxCharacterise` reads. **The three other component
  models `PhaseWax` can be given by name — `ComponentWonWax` (`Won`), `ComponentWaxWilson`
  (`Wilson`) and `ComponentCoutinhoWax` (`Coutinho`) — are reachable and cannot compute**, so
  they are not a port anybody skipped. Measured on the fluid this family's own probe uses
  (`WaxModelProbe`, `captures/wax_model_probe.tsv`): Pedersen's reports a wax fraction of
  `0.146` at 275 K, and all three of the others report **`1.0`** — the entire feed as wax,
  methane included. Two of them do it by returning `NaN`: Won's solubility parameter and
  Wilson's activity coefficient are `NaN` at every temperature a wax exists at, because their
  `sqrt` arguments go negative, and Coutinho's activity coefficient is finite while its
  coefficient comes out at `8.0e+252`. **Not ported** alongside them: `WaxCharacterise`
  itself, `WaxModelInterface`, and `pvtsimulation/flowassurance/WaxCurveCalculator`, which adds
  a sweep and a WAT rather than a model: it clones the fluid down a temperature grid, runs the
  ordinary wax flash at each point and enforces monotonicity on the resulting curve, so the
  flash it drives is `eos.tp_multiflash_wax`'s. **It is constructed, not dead** —
  `neqsim/mcp/runners/FlowAssuranceRunner.java:192` builds it, fully qualified, which is why a
  bare `new WaxCurveCalculator(` grep misses it. `process/chemistry/wax/` is carried with the
  process tier. **Filed upstream as NeqSim `#3914`**, under this finding's own title; the
  solid-vapour-pressure route below is `#3913`.
- **Asphaltene, carried on a defect upstream rather than for want of a port.**
  `AsphalteneCharacterization`, `PedersenAsphalteneCharacterization`,
  `AsphalteneOnsetPressureFlash`, `AsphalteneOnsetTemperatureFlash`, the
  `pvtsimulation/flowassurance/` screens (`DeBoerAsphalteneScreening`,
  `FloryHugginsAsphalteneModel`, `RefractiveIndexAsphalteneScreening`,
  `AsphalteneStabilityAnalyzer`, `AsphalteneMethodComparison`,
  `AsphalteneMultiMethodBenchmark`), the onset fitting pair (`AsphalteneOnsetFitting`,
  `AsphalteneOnsetFunction`) and `process/chemistry/asphaltene/`.
  The two onset flashes are the reachable pair, and what they scan over is a `TPflash` with
  `setSolidPhaseCheck("asphaltene")` — a call that **turns on the multiphase check as a side
  effect** (`SystemThermo.addSolidPhase`), which then skips `Flash`'s own solid hook and hands
  back a two-phase fluid collapsed to the feed. Measured on methane 0.30 / n-heptane 0.70 at
  333.15 K and 60 bara: `GAS 0.092231647152` against `OIL 0.907768352848` without the solid
  check, and one phase at `x = z` with it. So the sweep the onset is found along is not a
  sequence of states, and only its endpoint — reached where a solid forms and `SolidFlash`
  re-derives the phases — is a state at all. The coefficient the family's solid is built on is
  ported (`eos.solid_fugacity`, `eos.tp_solid_flash`); this family is revisited when upstream
  closes the defect, which is **filed upstream as NeqSim `#3913`**; the wax family's three
  models above are `#3914`.
- **Hydrogen and cryogenic.** Ported, each with the id that carries it:
  `thermo/util/leachman/` as `eos.hydrogen_phase`, the Leachman equation of state for the two
  spin isomers with both the dilute and the dense root selectable;
  `thermo/util/solid/ParaHydrogenSolidHelmholtzEquation` as `eos.parahydrogen_solid_phase`;
  and `thermo/util/solid/ArgonSolidHelmholtzEquation` as `eos.argon_solid_phase`.
  **Not ported**: `ComponentGERG2008Eos`, `ComponentGERG2004`, `PhaseGERG2008Eos`,
  `SystemGERG2008Eos`, `PHflashGERG2008`, `PSFlashGERG2008` and `thermo/util/gerg/`;
  `ComponentLeachmanEos`, `PHflashLeachman` and `PSFlashLeachman`, which are the flashes the
  Leachman *component* is used through rather than the equation; and
  `thermo/util/hydrogen/ParaOrthoH2Correction`.

### Tier 4 — the long tail

Port on demand, as a caller needs them:

- **CPA.** `ComponentPrCPA`, `ComponentSrkCPA`, `ComponentSrkCPAs`, `ComponentSrkCPAMM`,
  `ComponentUMRCPA`, `PhasePrCPA`, `PhaseSrkCPA` and its mixing variants,
  `CPAMixingRuleHandler`, `CPAMixingRuleType`, `SystemPrCPA`, `SystemSrkCPA`,
  `SystemUMRCPAEoS`.
- **SAFT.** `ComponentPCSAFT`, `ComponentPCSAFTa`, `ComponentSAFTVRMie`, `PhasePCSAFT`,
  `PhasePCSAFTa`, `PhaseSAFTVRMie`, `SystemPCSAFT`, `SystemPCSAFTa`, `SystemSAFTVRMie`,
  `TPflashSAFT`.
- **Reference equations of state.** `ComponentSpanWagnerEos`, `ComponentVegaEos`,
  `ComponentEOSCGEos`, `PhaseSpanWagnerEos`, `PhaseVegaEos`, `SystemSpanWagnerEos`,
  `SystemVegaEos`, `PHflashVega`, `PSFlashVega`, `EOSCGSaturationVUFlash`, and
  `thermo/util/spanwagner/`, `thermo/util/Vega/`.
- **UMR.** `ComponentUMRCPA`, `ComponentUMRCPAvolcor`, `ComponentGEUnifacUMRPRU`,
  `AttractiveTermUMRPRU`, `SystemUMRPRUEos`, `SystemUMRPRUMCEos`.
- **Flow assurance with no tranche behind it.** The rest of `pvtsimulation/flowassurance/`:
  corrosion (`DeWaardMilliamsCorrosion`, the de Waard-Williams correlation NORSOK M-506
  references, and `CO2CorrosionAnalyzer`, which couples it to an **electrolyte-CPA flash for
  pH** - so the coupled one has no route even upstream, Electrolyte-CPA being carried as
  unreachable), cooldown (`PipelineCooldownCalculator`, a lumped-parameter thermal transient,
  and `SurfCooldownAnalyzer`, which composes it with `eos.hydrate_formation_temperature` to
  give a no-touch time), `ErosionPredictionCalculator` (API RP 14E and DNV RP O501) and
  `EmulsionViscosityCalculator` (oil-water viscosity correlations and the inversion point).
  The specification assigns `flowassurance/` to "the tranche that backs each" and **no tranche
  backs these four** - so they are here, in the long tail, rather than behind a P-number that
  no page defines. The seven scale and six asphaltene screens in the same directory are not
  here: they belong to the salt and asphaltene families above.
- **Solids.** Ported, each with the id that carries it: `ComponentSolid`
  (`eos.solid_fugacity`, the tabulated coefficient) and `SolidFlash`
  (`eos.tp_solid_flash`, the fluid flash carrying one pure solid), over the
  `PhasePureComponentSolid` route. **Not ported**: `ComponentSolidHelmholtzEos`,
  `PhaseSolidComplex`, `PhaseSolidHelmholtzEos`, `SystemSolidHelmholtzEos`, `SolidFlash1`,
  `PHsolidFlash` and `thermo/util/solid/` — `SolidFlash1` because it normalises nothing and
  leaves phases whose compositions do not sum to one (measured on a water/methane feed: `x`
  summing to `0.819`). **Unreachable upstream**: `SolidFlash12`, which has no `new` in
  `src/main` and no reference outside its own file, and `SystemArgonSolidHelmholtzEos`, whose
  only sites are two tests and an inventory file — so the argon solid this library carries as
  `eos.argon_solid_phase` is a phase state with no system that builds it.
- **Sulfur.** `thermo/util/sulfur/SulfurThermodynamics`.
- **Amines.** `thermo/util/amines/` (`AmineSystem`, `AmineKentEisenberg`), and the
  amine viscosity and diffusivity methods.
- **Reactive and equilibrium, the rest.** `thermodynamicoperations/chemicalequilibrium/ChemicalEquilibrium`
  - **a name collision**: two classes are called `ChemicalEquilibrium`, and this one is a 100-line
  operation that loops `solveChemEq` over a system's phases until the composition stops moving,
  while the 1,211-line solver Tier 2 ports is the other. It is **reachable** -
  `ThermodynamicOperations` builds it - and it is **blocked on the same predicate as the flash's
  dispatch**: its whole body sits inside `if (system.isChemicalSystem())`, so on every fluid azoth
  admits it is the empty operation, and on the fluids it acts on the flash has already solved the
  chemistry. Its loop over phases is also a fiction - `solveChemEq` discards the caller's index -
  which its oracle records. `flashops/reactiveflash/` is **no longer here**:
  `ReactiveMultiphaseTPflash` and `ReactiveMultiphasePHflash` are Tier 2's
  `reactions.reactive_tp_flash` and `reactions.reactive_ph_flash`, with the stability analysis,
  the modified-RAND solve and the DIIS accelerator behind them.
- **Black oil.** `blackoil/`.
- **Standards.** `standards/` — standard and regulatory calculations.

## Not a port target

The following NeqSim trees are deliberately not ported, so that the mapping above is
complete rather than silent about them:

**`process/` is not in this list.** It is not physics, but it *is* a port target: the
unit-operation tier, which the specification puts at tranche P11, and the flowsheet executor
at P12. It is founded on the process calculus — `crates/azoth-process` carries the channel
types, the stream record, the palette loader and the checker, and `specs/unit_ops/` declares 29
unit operations of which **fourteen carry kernels**. The rest are owed, with one exception that
is stated rather than outstanding. **`unit_ops.simple_absorber` is refused on measured
evidence**: `SimpleAbsorber` is not the stage-wise absorber its ports describe but a fixed-point
loop over MDEA/CO₂ loading whose `setNumberOfStages` writes a field its `run` never reads, and a
faithful port needs the amine electrolyte chemistry P8 declined — `AmineSystem` and
`AmineKentEisenberg` are the classes that would close it, and P8's tier is where it belongs.
That names what would close it, which is the point: it is not "out of scope".

**The distillation column is no longer parked.** `unit_ops.distillation_column` is the one
entry that is a solver rather than a composition of kernels this library already has — 40,058
lines across 46 files in NeqSim's `process/equipment/distillation/`, `NaphtaliSandholmSolver`
alone 5,153 — and its declaration now carries the feed stage, the two ends and the two product
specifications a rigorous column takes. It runs as a workstream of its own, with five entries
beside it for the classes the rigorous base is the source of, so its numbers are neither owed
nor refused. Both solves are in scope, the sequential-substitution core and
Naphtali-Sandholm, and the oracle is differential because NeqSim carries no absolute reference
numbers for a column.

**P11's palette is closed except for two entries, and both are deferred with a measurement.**
Seven of the nine kernels the column's unparking deferred have landed - `component_splitter`,
`ejector`, `flare`, `gas_scrubber`, `stirred_tank_reactor`, `tank` and
`three_phase_separator` - each a registered `process.*` id with a spec, two implementations, a
NeqSim capture and a case set. What is left is `gibbs_reactor` and `plug_flow_reactor`, and
their palette entries now carry why:

- **`gibbs_reactor` is not the composition the P11 plan expected.** `GibbsReactor.run` does not
  call `ChemicalEquilibrium`, and the class is 3,163 lines carrying its own Lagrange-multiplier
  Newton solve (with an Armijo line search and Tikhonov regularisation) and its own species
  database - `GibbsReactDatabase.csv`, vendored, whose rows carry element vectors *and*
  per-species Gibbs, enthalpy and entropy correlations that are not the databank's formation
  properties. A port composed from P10 would answer with different numbers than the class, so it
  waits for a tier of its own.
- **`plug_flow_reactor` needs an integrator, and the tree has none.** Its `run` is 1,261 lines
  marching the molar flows with a catalyst bed's activity and bulk density, and it *adds* species
  the feed does not carry. The stepper would be new numerical code, and the plan leaves the
  choice open - a spec'd stepper of its own, or the kinetics path declared out and the isothermal
  plug limit ported.

- **`fluidmechanics/`** — azoth has its own hydraulics (`hydraulics.*`); this tree is
  NeqSim's parallel one and is not the port source.
- **`statistics/`, `util/`, `mcp/`, `mathlib/`, `integration/`, `datapresentation/`,
  `api/`** — infrastructure and tooling, not physics.

## Beyond the port

The interoperation surface — the middleware a flowsheet editor, a notebook and an agent
all drive — is built after the port, not beside it. Its shape is
[The middleware](docs/src/architecture/middleware.md), and the order it gates on is
[the specification](docs/src/architecture/specification.md)'s, which states it once.

What the P11/P12 kernels are waiting on is measured, and it is:

- `pipe` waits on density and viscosity assembled from a `Mixture` — the P1 transport
  models that read the collision and liquid-viscosity columns the manifest marks
  `not-yet`. The databank fields it also needs (molar mass, critical volume, dipole)
  are carried.
- `compressor`/`expander` wait on the **molar-entropy field `s`** on the stream record;
  the `entropy_at` helper already exists (`ps_flash`).
