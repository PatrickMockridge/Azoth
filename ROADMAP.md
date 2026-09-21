# NeqSim port roadmap

The order in which azoth ports [NeqSim](https://github.com/equinor/neqsim), from
the physics the most callers need to the physics the fewest do. Impact first: a property
or model that ordinary oil, gas and water work uses is ported before one only a
specialist reaches for.

This maps NeqSim's physics **one to one**. Every class under NeqSim's `thermo/`,
`thermodynamicoperations/` and `physicalproperties/` trees is assigned to a tier below;
the trees that are not physics — `process` (unit operations and the flowsheet),
`fluidmechanics` (azoth has its own hydraulics) and the support packages — are named at
the end rather than silently dropped. The column-by-column record is
`databank/manifest.toml`, which this file orders; it does not re-inventory it.

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

### Tier 0 — the one data path

The cubic Peng-Robinson core is ported: `ComponentPR`, `AttractiveTermPr`,
`PhasePrEos`, `SystemPrEos`, `TPflash`, `PHflash`, `PSFlash`, the tangent-plane
stability test, `CriticalPointFlash`, the phase envelope, `RachfordRice` and the
classical vdW1f mixing rule, reading 14 of the databank's columns. What remains is the
foundation the rest stands on:

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
solubility, electrolytes, salts and scale, and freezing.

- **The water EOS.** `ComponentWater`, `PhaseWaterIAPWS`, `SystemWaterIF97`,
  `thermo/util/steam/Iapws_if97`.
- **Water content and dehydration.** `WATcalc`, `WaterDewPointTemperatureFlash`,
  `WaterDewPointTemperatureMultiphaseFlash`, `WaterDewPointEquilibriumLine`.
- **Acid gas.** `ComponentSoreideWhitson` and `AttractiveTermSoreideWhitson` with
  `PhaseSoreideWhitson` (`eos.soreide_whitson_phase`) and `SystemSoreideWhitson`;
  `ComponentGEVanLaarAcid` with `PhaseGEVanLaarAcid` (`eos.ge_van_laar_acid_phase`).
  **The brine operations are not ported.** `CO2BrinePhaseEquilibrium`,
  `ReactiveCO2BrinePhaseEquilibrium`, `SaturateWithWater` and `CalcIonicComposition` all
  need the EoS/GE hybrid seam — `SystemEosGE`, `HybridEosGeFlashModel` and
  `TPHybridEosGeFlash`, 1,636 lines between them — which is a flash whose two phases run
  *different models*. azoth's `Mixture` carries one mixing rule for the whole system, so
  closing this is a change in the flash and not in the electrolyte physics.
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
  `HydrateFormationPressureFlash` (`eos.hydrate_formation_pressure`); and `TPHydrateFlash`
  (`eos.hydrate_fraction`), **with the composition taken from both cavities and the material
  balance asserted at the answer**. **Not ported**:
  `PitzerHydrateFlash` with `ComponentHydratePitzer` and `ComponentHydrateGF`, which are
  reachable through the model name and whose coupling to the electrolyte phases is the seam
  this would need; the two inhibitor flashes (`HydrateInhibitorConcentrationFlash`,
  `HydrateInhibitorwtFlash`); `HydrateEquilibriumLine` and `HydrateEquilibriumDiagnostics`;
  and **carried as unreachable upstream** `ComponentHydrateKluda`, which has no construction
  site anywhere in `src/main`, `ComponentHydrateStatoil` and `ComponentHydrateBallard`, whose
  only sites are two commented-out lines in `PhaseHydrate`, and `OLGAhydrateCurveGenerator`,
  which only a test builds. `process/chemistry/hydrate/` is carried with the process tier.
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
  process tier.
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
  closes the defect.
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
- **Reactive and equilibrium.** `thermodynamicoperations/flashops/reactiveflash/`
  (`ReactiveMultiphaseTPflash`, `ReactiveMultiphasePHflash`, `ReactiveStabilityAnalysis`,
  `ModifiedRANDSolver`, `DIISAccelerator`), `ChemicalEquilibrium`, and `chemicalreactions/`.
- **Black oil.** `blackoil/`.
- **Standards.** `standards/` — standard and regulatory calculations.

## Not a port target

The following NeqSim trees are deliberately not ported, so that the mapping above is
complete rather than silent about them:

- **`process/`** (unit operations, equipment, the flowsheet). azoth's specification puts
  the unit-operation tier at tranche P11; the tier was deleted, not deferred piecemeal,
  and it is built on the physics above once that is complete.
- **`fluidmechanics/`** — azoth has its own hydraulics (`hydraulics.*`); this tree is
  NeqSim's parallel one and is not the port source.
- **`statistics/`, `util/`, `mcp/`, `mathlib/`, `integration/`, `datapresentation/`,
  `api/`** — infrastructure and tooling, not physics.

## Beyond the port

The interoperation surface — the middleware that a flowsheet editor, a notebook and an
agent all drive — is built after the port, not beside it. Its shape is
[The middleware](docs/src/architecture/middleware.md), and it gates on the backend closing
in this order:

1. **Tier 0** — the one data path closes: the `not-yet` ideal-gas Cp and reference-state
   columns that `eos.ideal_gas_cp` and `eos.molar_enthalpy_entropy` still take from the
   caller.
2. **P1** — transport properties; the pipe and equipment kernels need density and
   viscosity from a mixture.
3. **P2–P10** — the physics.
4. **P11** — a kernel for every palette entry, not six of 25.
5. **P12** — the flowsheet executor: topological order, a recycle fixed point, a session,
   structured diagnostics, the TOML/JSON round trip and a result codec.

Known blockers that gate P11/P12 kernels:

- `pipe` waits on density and viscosity assembled from a `Mixture` — the P1 transport
  models that read the collision and liquid-viscosity columns the manifest marks
  `not-yet`. The databank fields it also needs (molar mass, critical volume, dipole)
  are carried.
- `compressor`/`expander` wait on the **molar-entropy field `s`** on the stream record;
  the `entropy_at` helper already exists (`ps_flash`).
- `eos.pt_phase_envelope`'s dew branch and critical point need NeqSim's analytic Jacobian
  and `calcCrit` for a tight port; the central-difference version is approximate.
