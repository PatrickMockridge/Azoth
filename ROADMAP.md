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
classical vdW1f mixing rule, reading 14 of the databank's columns. **Both things it was
waiting on are read rather than owed**, and `tools/audit_manifest.py` is what holds the second:

- **the single `databank → keycard → every calculation` path** — demonstrated rather than
  declared: `eos.phase_transport` names its components, resolves them through the databank or a
  card, and reads each one's parameters off the entry. Thirty-three specs still take a component
  property as an *input* — `Tc`, `Pc`, `omega`, `molar_mass` — and those are the low-level terms
  NeqSim also constructs with the constant: an `AttractiveTerm` is built *with* an omega, and it
  is the `Component` that holds one. `eos.ideal_gas_cp` stays a scalar calc by design.
- **the gas diffusivity method** — ported, as `eos.phase_transport`'s gas branch: NeqSim's base
  `Diffusivity`, which `GasPhysicalProperties` assigns to a gas, taking the pair from the two
  components' Lennard-Jones parameters. **The `LJEPS` row said `not-ported` and the tree read
  it**, the same correction `PARACHOR` carries; the row, and the three
  `LIQUIDCONDUCTIVITY` columns beside it, name their reader now. `eos.viscosity` and
  `eos.thermal_conductivity` still do not read any of these — they port
  `PFCTViscosityMethodHeavyOil` and `PFCTConductivityMethod`.

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
  `SystemUNIFACpsrk`, `SystemGEWilson`. The four `eos.ge_*_phase` ids carry
  `ComponentGE.fugcoef`'s two reference states: `gamma_i P0_i / P` for a `solvent`
  component and `(gamma_i / gamma_i^inf) H_i / P` for every other, with the IAPWS table
  preferred for its gases in water. The refusals, and the measurement behind each, are in
  the eos-family bullet below.
- **Component correlations** read from `thermo/component/` (normal boiling point,
  Antoine vapour pressure) and `physicalproperties/` (standard liquid density, heat of
  vaporisation, liquid heat capacity) — the columns the manifest lists under those
  packages.
- **Petroleum-fraction characterisation.** `thermo/characterization/`: `Characterise`,
  `PlusCharacterize`, `TBPCharacterize`, `TBPfractionModel`, `PlusFractionModel`,
  `PedersenPlusModelSolver`, `LumpingModel`, `LumpingConfigBuilder`,
  `PseudoComponentCombiner`, `Recombine`, `OilAssayCharacterisation`, `TbpClosure`,
  `RefineryAssayBlend`, `BioFeedstock`, `BiomassCharacterization`.
- **Transport properties.** `physicalproperties/methods/` and `interfaceproperties/`.
  Ported, each with the id that carries it: `ChungViscosityMethod` and the liquid polynom
  `Viscosity` (`eos.chung_viscosity`, `eos.liquid_viscosity_pure`, `eos.aqueous_viscosity`),
  Wilke's mixing rule (`eos.wilke_viscosity`), `ChungConductivityMethod`,
  `MasonSaxenaConductivityMethod` and the polynom `Conductivity` (`eos.chung_conductivity`,
  `eos.mason_saxena_conductivity`, `eos.liquid_conductivity_polynom`), `Rackett` and
  `Costald` (`eos.rackett_molar_volume`, `eos.costald_molar_volume`), the seven diffusivity
  methods (`eos.chapman_enskog_diffusivity`, `eos.co2_water_diffusivity`,
  `eos.fuller_schettler_giddings_diffusivity`, `eos.hayduk_minhas_diffusivity`,
  `eos.siddiqi_lucas_diffusivity`, `eos.tyn_calus_diffusivity`,
  `eos.wilke_chang_diffusivity`), `ParachorSurfaceTension` (the pure and mixture forms,
  `eos.parachor_surface_tension` / `eos.parachor_mixture_surface_tension`),
  `FirozabadiRamleyInterfaceTension` (`eos.firozabadi_ramley_interface_tension`),
  `FilippovConductivityMethod` (`eos.filippov_conductivity`) and
  `FrictionTheoryViscosityMethod` (`eos.friction_theory_viscosity`), the phase-type
  dispatch that assigns the transport correlations (`eos.phase_transport`), and
  `InterfaceProperties`' tension dispatch (`eos.interface_properties`): model 0, the
  default, installs Parachor for gas-oil and gas-aqueous and Firozabadi-Ramley for
  liquid-liquid, models 3 and 5 are its aliases, and the string selector's two
  calculators are selectable. `LGTSurfaceTension` is `eos.lgt_surface_tension`, model 2's
  calculator, whose 200-step profile runs on `Mixture::helmholtz_gradient`, and
  `GTSurfaceTensionSimple` is `eos.gt_surface_tension_simple`, model 4's gas-oil and
  model 1's gas-aqueous and liquid-liquid slot, which adds the `dmu/dn` inner Newton, and
  `GTSurfaceTension` is `eos.gt_surface_tension`, model 1's gas-oil slot, solved through its
  reachable FullGT branch on a remeshed finite-difference grid (`GTSurfaceTensionODE` is
  unreachable upstream - `useFullGT = 1` has no setter). `CDFTSurfaceTension` is
  `eos.cdft_surface_tension`, the string selector's other calculator: a tanh profile on a
  2048-point grid, varied over its width, taking each phase's molar volume and mole
  fractions because that is all the class reads from the flashed system.
  `eos.interface_properties` now carries all six arms and the string selector's three
  calculators: the gas-liquid, gas-aqueous and liquid-liquid slots install model 0's
  Parachor/Firozabadi pair, model 1's `GTSurfaceTension`/`GTSurfaceTensionSimple` pair,
  model 2's `LGTSurfaceTension` in every slot, and model 4's
  `GTSurfaceTensionSimple`/`LGTSurfaceTension`/Parachor, from the flashed pair and the
  profile parameters. Its Parachor and Firozabadi slots read the databank's `parachor`
  (the `parachor_cpa` column for a CPA mixture, per upstream PR #4004) and `tc`, and both
  phases' translated densities - the SRK Peneloux shift, which the
  databank now carries beside the Peng-Robinson one, is what makes those densities
  NeqSim's rather than the untranslated cubic's. `FrictionTheoryConductivityMethod` is
  `eos.friction_theory_conductivity`, reached the only way the class reaches it - an
  explicit `setConductivityModel("friction theory")` - and taking the phase's own
  viscosity as an input, as the class reads `getViscosity()`.
  **`GTSurfaceTensionSimple` returns a partial integral from a collapsed profile** on
  CO2/water 0.02/0.98 at 313.15 K and 50 bar: densities are zero from step 3 and it answers
  `0.0145737553` over two steps. The port refuses that state; `gt_simple_profile_probe.tsv`
  is the reproduction. The captures `surface_tension_probe.tsv`,
  `surface_tension_solvers_probe.tsv`, `lgt_profile_probe.tsv`, `lgt_local_probe.tsv`,
  `gt_simple_profile_probe.tsv`, `gt_simple_initial_probe.tsv`, `full_gt_probe.tsv`,
  `cdft_probe.tsv` and `friction_conductivity_probe.tsv` record the default pairs, models
  1, 2 and 4, the profiles' negative-pressure segments, the GT-simple profile, the FullGT
  answers, the CDFT states with the parachor beside each so a silent fallback would show,
  and the friction-theory conductivity on both cubics.
- **Flash breadth.** `TVflash`, `PVflash`, `VUflash`, `THflash`, `TSflash`, `TUflash`,
  `PUflash`, `PVFflash`, `TVfractionFlash`, `QfuncFlash`, `ImprovedVUflashQfunc`.
- **Phase envelopes.** `PTphaseEnvelope`, `HPTphaseEnvelope`, `PTPhaseEnvelopeMichelsen`,
  `CricondenBarFlash`, `CricondenThermFlash`, `SysNewtonRhapsonPhaseEnvelope`.
- **Flow assurance.** `pvtsimulation/flowassurance/` is carried as **one named group — 21
  classes and 12,486 lines** — because it is not one family. **Its seven scale calculations are
  ported**, each with the id that carries it: `ScalePredictionCalculator`
  (`eos.scale_prediction`), `ScaleMassCalculator` (`eos.scale_mass`),
  `PitzerScaleActivityModel` (`eos.pitzer_scale_activity`),
  `BariteCelestiteSolidSolution` (`eos.barite_celestite_solid_solution`),
  `MultiMineralScaleEquilibrium` (`eos.multi_mineral_scale_equilibrium`),
  `FlowlineScaleProfile` (`eos.flowline_scale_profile`) and `WaterCompatibilityScreener`
  (`eos.water_compatibility`). **The remaining reachable calculators are ported too**:
  `DeWaardMilliamsCorrosion` (`eos.de_waard_milliams_corrosion`),
  `ErosionPredictionCalculator` (`eos.erosion_prediction`),
  `EmulsionViscosityCalculator` (`eos.emulsion_viscosity`),
  `PipelineCooldownCalculator` (`eos.pipeline_cooldown`),
  `HydrateRiskMapper` (`eos.hydrate_risk_map`) and `SurfCooldownAnalyzer`
  (`eos.surf_cooldown`). **One class is left**: `CO2CorrosionAnalyzer`, whose pH comes from an
  electrolyte-CPA flash this tree does not carry the system route for; it is Tier 4's alone.
  The six asphaltene screens stay carried on NeqSim defect `#3913`, and the one wax
  (`WaxCurveCalculator`) belongs to Tier 3's wax bullet. **The rest of `pvtsimulation/` is not
  here**: the simulation suite,
  the model tuning and the reservoir properties are engineering deliverables, and
  [the specification](docs/src/architecture/specification.md) puts them beyond the port.

### Tier 2 — water and gas-water

Water and the aqueous models: the water EOS, water content and dehydration, acid-gas
solubility, electrolytes, salts and scale, freezing, and the reactions an aqueous phase
reaches equilibrium through.

- **The water EOS.** `ComponentWater`, `PhaseWaterIAPWS`, `SystemWaterIF97`,
  `thermo/util/steam/Iapws_if97` - **and `Iapws_if97` is the tree's already.** Regions 1, 2 and
  4 are `crates/azoth-eos/src/iapws_if97.rs`, with the dimensionless Gibbs functions and their
  derivatives, and `eos.water_phase` is that module as a model. **This is the fourth bullet in
  this register that names as owed what the tree already carries** - after the brine operations,
  the wax bullet's `WATcalc`, and the water-dew-point seeds. What is left of it is small and
  named here rather than left to be re-derived: `w`, the speed of sound, is computed and checked
  against NeqSim in the kernel's own test and is **not** an output of the model; and the Region 4
  saturation curve, which the module carries in both directions, has no id of its own - it is the
  boundary everything else in IF97 is selected by, and it is neither `eos.pure_saturation` (a
  cubic's roots crossing) nor `eos.antoine_vapor_pressure`.
  **The classes the bullet does name are the plumbing, and the same finding is in every bullet of
  this tier.** `ComponentWater`'s `calca`, `calcb` and `alpha` are zero and its `fugcoef` is
  `1.0`; `PhaseWaterIAPWS` reads IF97 and types itself by specific volume, `v_spec < 0.01` being
  `AQUEOUS` and anything above it `GAS`; `SystemWaterIF97` is a pure-water fluid you can flash,
  and it refuses any component but water.
  **The `System*` classes the bullets below name are not one finding, and the distinction is
  which of them a flash reaches.** **One is**: `SystemPitzer`, whose phase is
  `eos.pitzer_phase`, called by `eos.multiphase`, `reactions.reactive_hybrid_eos_ge_flash` and
  `eos.salt_precipitation`. **Four are not**: `eos.kent_eisenberg_phase`,
  `eos.desmukh_mather_phase`, `eos.furst_electrolyte_phase` (with its Mod2004 twin) and
  `eos.soreide_whitson_phase` are each a phase state whose resolver **no flash calls** -
  `soreide_whitson_mixture_of`'s only caller is the phase model itself. So the gap is a route
  from these phases into a flash, which is not the same thing as the phase's arithmetic being
  absent, and the reason it is missing differs by model: the Fürst one's Huron-Vidal rule
  refuses the `d ln phi / d n` Jacobian a flash needs, while the other three have never been
  wired. The one class in the family whose phase model really is not a cubic is
  `SystemWaterIF97`. Owed, with this reason.
- **Water content and dehydration.** Two of the four classes in this bullet are ids now.
  `WaterDewPointTemperatureFlash` is `eos.water_dew_point_temperature`: it **prescribes** its
  incipient phase - `1e-30` on every component but one `water`-or-`MEG` at `1.0` - and walks the
  temperature to the fugacity equality. `WaterDewPointTemperatureMultiphaseFlash` is
  `eos.water_dew_point_multiphase`: it brackets on where a *multiphase flash* reports an aqueous
  phase, and it reproduces the capture to `3e-13` K at all three of its states.
  **The third, `WaterDewPointEquilibriumLine`, is `eos.water_dew_point_line` now, and the
  measurement is why it waited.** It is the multiphase route swept over ten pressures, so its
  only output is the curve - and on the capture's own sweep from `1` to `200` bara **nine of its
  ten points reproduce to `4e-13` K and the first does not**: at `1` bara this library answers
  `224.4312` K where NeqSim answers `224.2818`, a `0.149` K that is a flash-level difference in a
  marginal state - at NeqSim's own answer the flash reports two phases, gas and aqueous, and it
  loses the aqueous phase only `0.149` K warmer. That point rides the model family's `1e-3`
  capture tolerance and is declared in the case's own `source`. The other state that held the
  class back was **not** a divergence but a stall - `111.6` bara, where `eos.tp_multiflash` took
  three hundred steps to a residual of `3.7e-5` against its own `1e-10` - and that is fixed: a
  stalled two-phase flash now warns and seeds from the feed instead of refusing.
  **And this bullet named the wrong blocker for the route, twice over** - the way the entry above
  it named the wrong reason for the brine operations and the wax bullet below named the wrong tier
  for `WATcalc`. `seedHydrocarbonLiquidFromFeed` and `seedAdditionalPhaseFromFeed` seed an
  **OIL** and a **GAS** and each **requires an aqueous phase to exist already**, so neither can
  fire on a lean gas; and the aqueous seed inline in `TPmultiflash.run()` is not what the route
  was missing either. What it was missing is that `TPmultiflash` **overrides** `stabilityAnalysis`
  with a second method of the same name - the same two Wilson trials **plus a pure-component
  fallback** - and this tree had ported the *base* method. That is `eos.multiflash_stability_test`.
  **`WATcalc` was filed in this bullet and is not a water calculation at all** - it is Tier 3's
  wax bullet now, with its own blocker.
- **Acid gas.** `ComponentSoreideWhitson` and `AttractiveTermSoreideWhitson` with
  `PhaseSoreideWhitson` (`eos.soreide_whitson_phase`) and `SystemSoreideWhitson`;
  `ComponentGEVanLaarAcid` with `PhaseGEVanLaarAcid` (`eos.ge_van_laar_acid_phase`).
  **All four brine operations are ported.** `CO2BrinePhaseEquilibrium` is
  `eos.co2_brine_phase_equilibrium`, `CalcIonicComposition` is `eos.calc_ionic_composition`,
  and `SaturateWithWater` is `eos.saturate_with_water`, all on the electrolyte-CPA route
  above. `ReactiveCO2BrinePhaseEquilibrium` is
  `reactions.reactive_co2_brine_phase_equilibrium`: that phase selection and the chemistry
  of `reactions.reactive_phase_equilibrium`, alternated in a bounded loop, with the
  unsymmetric activity coefficients `Phase.initRefPhases` builds. **Both captured states
  reproduce the phase inventory to `1e-9`; the pass count and the accepted residual are
  path-decided**, because this port's chemistry solve fails from the floored post-flash
  composition where NeqSim's converges - the adaptive-derivative refinement the reactions
  bullet below lists as not ported is the reason. The first is not a flash - it holds one aqueous phase and runs its own two trials
  and its own ion-confined split - and the second solves nothing at all: it reads a phase and
  reports molalities and activity coefficients at `Phase`'s two reference states. The third
  floods a feed with water and reads the water dissolved in the non-aqueous phases; **an
  associating fluid is the Statoil CPA**, the generic CPA mixture plus `PhaseSrkCPAs.calc_g`'s
  Chapman distribution and rule 10's temperature-dependent `cpakijT_SRK`. The two source CPA
  states reproduce to 0.2-0.3% once both are in; the multiphase flash's own one-phase answer
  carried a negative flash's extrapolated composition, which is fixed in both kernels. The
  coupling's gate was once recorded here as the EoS/GE hybrid seam, and that
  was wrong: it is `fluid instanceof SystemElectrolyteCPAstatoil`, which extends
  `SystemFurstElectrolyteEos` extends `SystemSrkEos`, so the hybrid route is not in its chain
  at all — that route is ported (`eos.hybrid_eos_ge_flash`,
  `reactions.reactive_hybrid_eos_ge_flash`) and was never the blocker.
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
  `SystemElectrolyteCPA`. `SystemElectrolyteCPAstatoil` and `PhaseElectrolyteCPAstatoil`
  are **ported** as `eos.electrolyte_cpa_statoil_phase` — 145 lines, and the four
  differences from the Fürst fluid are the point of it: the ion covolumes are the
  `furstParamsCPA` set, the attractive term is `AttractiveTermCPAstatoil`, the mixing rule
  is the classical one with the CPA `kij` (rule 10), and the association's cross term is
  the CR-1 table rather than Elliott. The four brine operations that run on it are the
  acid-gas bullet's, above.
- **Salts and scale.** Ported, each with the id that carries it: `CheckScalePotential`
  (`eos.scale_saturation_ratio`, one salt's `IAP/Ksp` against its solubility product) and
  `MultiSaltPrecipitation` over `CalcSaltSatauration` (`eos.salt_precipitation`, one mineral's
  extent, which is the per-mineral half of the complementarity loop). The
  `pvtsimulation/flowassurance/` scale set is the other route and is ported too, its seven ids
  named in Tier 1's flow-assurance bullet: the water-analysis screen, its Pitzer coefficients,
  the solid solution, the mass estimate, the coupled equilibrium and the two composers.
  **Not ported**: `AddIonToScaleSaturation` and `CalciumSulfatePhaseBoundaryQualification`.
  **All eleven `process/chemistry/scale/` classes are ported**, each with the id that carries it:
  `ElectrolyteScaleCalculator` (`eos.electrolyte_scale_calculator`),
  `BrineMixingScaleEvaluator` (`eos.brine_mixing_scale_evaluator`),
  `ScaleInhibitorPerformance` (`eos.scale_inhibitor_performance`), `ScaleKinetics`
  (`eos.scale_kinetics`), `ScaleRemediationAdvisor` (`eos.scale_remediation_advisor`),
  `ScaleControlAssessor` (`eos.scale_control_assessor`), `StreamScaleAnalyzer`
  (`eos.stream_scale_analyzer`), `ValveScaleDrift` (`eos.valve_scale_drift`) and
  `ScaleDepositionAccumulator` (`eos.scale_deposition_accumulator`),
  `ProductionChemicalScaleScenario` (`eos.production_chemical_scale_scenario`) and
  `ClosedLoopDepositionSolver` (`eos.closed_loop_deposition_solver`).
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
  stand on. **The subsystem is ported, as eighteen ids over the tables it works on**:
  `characterization.pedersen_plus_split` and `whitson_gamma_split` divide a plus fraction into
  carbon-number cuts, `lumping` groups a cut table into equal-mass lumps and `delumping` splits
  one back into a finer grid,
  `characterize_to_reference` re-cuts one fluid onto another's slate,
  `tbp_cut_properties` describes a cut by any of the ten models, `tbp_closure` and `tbp_density`
  convert between a cut's boiling point, its molar mass and its gravity, `tbp_cut_boundaries`
  cuts a cumulative TBP curve into an assay, `assay_mass_fractions` resolves one to a mass basis,
  `oil_assay_derivation` derives each cut's molar mass and mole amount, `tbp_grouping` groups a
  phase by boiling point, `set_heavy_tbp_to_plus` picks the row a plus end begins at,
  `characterise_plus_fraction` is the chain from a C7+ end to lumps,
  `combine_reservoir_fluids` blends several fluids onto one slate,
  `characterize_to_common_slate` re-characterizes several onto a shared one,
  `transfer_binary_interaction_parameters` moves an interaction table between two, and
  `tbp_fraction_component` is `addTBPfraction`'s card for one cut - the Rackett compressibility
  among it, from the reference SRK phase's liquid root. `mixture_from_characterization` turns the
  rows any of them answers with into a `Mixture` and its ideal-gas model, and
  `oil_assay_to_mixture` is `OilAssayCharacterisation.apply`'s composition - an assay's cuts
  derived, carded and built into a fluid. What is owed of the class is only its parts that build
  nothing: `Recombine`, whose constructor is empty and whose `runRecombination` returns null.
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
  **The Pitzer route is ported**: `eos.pitzer_hydrate_formation_temperature` and
  `eos.pitzer_hydrate_formation_pressure` carry `PitzerHydrateFlash` with
  `ComponentHydratePitzer`, whose water reference is the pure-water Antoine pressure, over
  `eos.hybrid_eos_ge_flash`. A salty CO2 state runs the manual PHREEQC CO2-chloride dataset
  the upstream test applies, selected by the `pitzer_dataset` input. The tranche also closed
  the two defects it exposed: the hybrid flash's role activation
  (`activateSupportedHybridPhaseRoles`) and the neutral layer's contribution to an ion's
  Pitzer coefficient. `HydrateEquilibriumDiagnostics` is an audit of a state rather than a
  model, and what it asserts belongs in a test. **Carried as unreachable
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
  **`WATcalc` belongs here and not in Tier 2's water bullet, where it was filed.** Its own
  javadoc says *"Wax appearance temperature from a bracket of independent multiphase TP
  flashes"*, its threshold is a **wax** mass fraction of `1e-8` and its predicate is `hasWax`.
  **Not ported, and the blocker is measured rather than assumed**: it answers only on a fluid
  whose wax came from `addTBPWax()` on a *characterized plus fraction*, which is
  `WaxCharacterise` above. On every fluid of the databank's named wax formers it throws
  `Invalid wax mass fraction`, because `getPhaseFraction("wax", "mass")` returns
  `-1.6e-52` for a wax phase whose mole fraction is `0.0707` - **filed upstream as NeqSim
  `#4309`**, with the reproduction in `captures/water_dew_point_probe.tsv`.
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
- **Flow assurance's one remainder.** `CO2CorrosionAnalyzer` couples the de Waard-Williams
  correlation to an **electrolyte-CPA flash for pH**, and this tree carries the
  electrolyte-CPA phase without the system route that solves an aqueous pH - so the coupled
  analyzer has no route here. **Everything else in `pvtsimulation/flowassurance/` is ported**:
  the seven scale calculations (Tier 2's salts bullet), `DeWaardMilliamsCorrosion`,
  `ErosionPredictionCalculator`, `EmulsionViscosityCalculator`, `PipelineCooldownCalculator`,
  `HydrateRiskMapper` and `SurfCooldownAnalyzer` (Tier 1's flow-assurance bullet). The six
  asphaltene screens belong to Tier 3's asphaltene bullet and `WaxCurveCalculator` to its
  wax bullet; neither is here.
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
types, the stream record, the palette loader, the checker and the executor, and `specs/unit_ops/`
declares 31 unit operations of which **30 carry kernels**. The one that does not says why, and
the executor refuses it by name. **`unit_ops.simple_absorber` is refused on measured
evidence**: `SimpleAbsorber` is not the stage-wise absorber its ports describe but a fixed-point
loop over MDEA/CO₂ loading whose `setNumberOfStages` writes a field its `run` never reads, and a
faithful port needs the amine electrolyte chemistry P8 declined — `AmineSystem` and
`AmineKentEisenberg` are the classes that would close it, and P8's tier is where it belongs.
That names what would close it, which is the point: it is not "out of scope".

**The distillation column is no longer parked, and the family it founded is ported.**
`unit_ops.distillation_column` is the one entry that is a solver rather than a composition of
kernels this library already has — 40,058 lines across 46 files in NeqSim's
`process/equipment/distillation/`, `NaphtaliSandholmSolver` alone 5,153 — and its declaration now
carries the feed stage, the two ends and the two product specifications a rigorous column takes.
**Both solves are in**, the sequential-substitution core and Naphtali-Sandholm, and the oracle is
differential because NeqSim carries no absolute reference numbers for a column. Four of the five
entries beside it are in as well: the FUG shortcut column, the absorber and the stripper, which
are one class under two sets of names, and `unit_ops.packed_column`, whose packing turned out to
be a report read after the solve — measured, NeqSim's packed column at a 2.0 m bed and its base
column at four stages are bit-identical on all twenty-two captured quantities. **The fifth is a
measured non-port and a tranche of its own**, as its own bullet below says, and so are the
reactive trays the base class reaches for — and **those are in now too**: their two flashes
take a cubic, the tray flashes the fluid the column carries, and a reactive answer whose two
phases are the same state is refused rather than given an outlet it did not earn. That is the same discipline the rest of
this record follows — each class beside the column is either ported or recorded as a non-port with
the class and the blocker named — and it replaces an earlier sentence here which claimed the five
were "neither owed nor refused", which the fifth entry's own measurement contradicted.

**The trays' side draws and pumparounds are in as well.** `SimpleTray`'s three draw fractions
are a split of the tray's own outlet phase — `getGasOutStream` is the vapour scaled by
`1 - gasSideDrawFraction`, `getLiquidOutStream` the liquid by `1 - liquid - pumparound`, and each
draw that same phase scaled by its own fraction — so a draw carries the tray's composition,
temperature and pressure and a different flow. **Both of the column's closures count them**,
which is what `columnReportsSideDrawAsOutletStream` asserts, and a fraction on an end is refused
because this port's ends are kernels rather than stages. What is owed is not a split but a
recycle: the draw's **flow-rate specification** (`ColumnSideDrawSpecification`,
`DistillationColumn.addSideDrawFlowSpecification`), which moves the fraction until a target flow
is met through the class's own tear loop, and the pumparound's **return** to another tray
(`ColumnPumparound`, a return tray and a temperature drop). Both are iterated by the column's own
solver and belong to P12's executor, and the capture carries the measurement that says so: on the
liquid-and-pumparound row NeqSim's own balance does not close — it answers `RECONCILED_PRODUCTS`
with a `2.4e-3` imbalance and products `0.02` mol/s above this port's — where this port closes
to `1e-8`, because the class is reconciling against a return this port has no tray to return to.

**P11's palette is closed.** The nine kernels the column's unparking deferred have all landed -
`component_splitter`, `ejector`, `flare`, `gas_scrubber`, `plug_flow_reactor`,
`stirred_tank_reactor`, `tank` and `three_phase_separator` - each a registered `process.*` id
with a spec, two implementations, a NeqSim capture and a case set; so have the classes the column
is the source of, beside the column itself, and so has `gibbs_reactor`. **The plug-flow reactor's
axial profile is published rather than summarised** - the summed reaction rate, on a vocabulary
dimension added for it, the per-station molar flows as a matrix, and the residence time - so the
class's `toJson()` has nothing left to serialise that the result does not already carry:

- **`gibbs_reactor` was not the composition the P11 plan expected, and it cost a tier of its
  own.** `GibbsReactor.run` never calls `ChemicalEquilibrium`: the class is 3,163 lines carrying
  its own Lagrange-multiplier Newton solve and its own species database, whose formation
  properties are not the databank's. The tier compiled that database from the class's two CSVs,
  read the solver, and built it - `data/reactors/`, `reactor::gibbs_solver`, and the id. **The
  port reproduces the class to the accuracy the class's conditioning permits**: one iteration
  from the same state agrees to `8e-14` on the major species, and the hundred damped steps the
  class takes separate two runs by the class's own `cond(J) = 6.12e6`. **The solver's own
  machinery is in**: a singular Jacobian falls back to `SimpleMatrix.pseudoInverse`'s
  construction (`solve_least_squares`), the Armijo and adaptive step controls and the Tikhonov
  augmentation are optional inputs at the class's own defaults, and `getConditionNumberHistory()`
  is an output - all held to `process_gibbs_algorithms.tsv`, the capture the class's own
  algorithm test supplies. **`runTransient` is stated as the interface's fallback rather than
  owed**, because neither reactor overrides it.

**And the one that was a second physics is ported.** `unit_ops.rate_based_packed_column` is
4,081 lines of segment model - axial slices through the packing, each with gas and liquid film
coefficients, an interface equilibrium and an interphase heat balance, coupled counter-currently
by a profile solver - and no `src/main` caller. It is `process.rate_based_packed_column` now,
and its **packing is inside the equations** rather than a report read after a solve, which is
why it needed `hydraulics.packing_hydraulics` and the phase-level transport dispatch before it
could be a model at all. On the class's own absorber state re-cased on PR it converges in
fifteen passes at a residual of `6.02e-10` mol/s where NeqSim's own row says fifteen and
`6.009670053264138e-10`, and the products, the segment temperatures, `kGa`, `kLa`, the wetted
area, the flood and the pressure drop agree inside `1e-4`. **Three of the class's own fallbacks
are taken on its own states, and that is measured**: the effective-diffusivity vector is never
populated upstream, so the reference is the constant while the pair matrix is real; and a
gas-and-aqueous interface answers a surface tension of `0.0`, so the constant stands in there.
What stays owed is named with its class: the equation-oriented column solver with its homotopy
continuation, the simultaneous-residual segment solve whose own test NeqSim disables, and the
second mass-transfer correlation, which is a constant multiplier rather than a correlation.
`unit_ops.simple_absorber` stays refused, as above.

**The third was refused for a defect in its own declaration, and the entry was fixed rather than
the executor worked around.** `unit_ops.packed_column` had a kernel and a registered `process.*`
id, and its palette entry declared the *packing* rather than the column - `feed_stage`,
`number_of_stages`, the two pressures and the two ends are what `PackedColumn.run` reads through
`super.run`. The entry now carries them, and the kernel derives the stage count from
`packed_height` in one place instead of the model doing it: **`packed_height` replaces
`number_of_stages`**, which is why the entry declares no stage count at all. It is the palette
failing its own test, caught by the dispatch table's first run - and the test's remedy was the
entry's.

- **`fluidmechanics/`** — azoth has its own hydraulics (`hydraulics.*`); this tree is
  NeqSim's parallel one and is not the port source.
- **`statistics/`, `util/`, `mcp/`, `mathlib/`, `integration/`, `datapresentation/`,
  `api/`** — infrastructure and tooling, not physics.

### The classes the ported entries still owe

**A capability that is missing has no parameter to key a declaration to.** The porting fact now
has a field — a `[[unported]]` row in the model spec, which `tools/check_unported.py` holds to
both implementations — and a row is only available where the port *refuses a value of a declared
input*. Where the class is absent rather than the value refused, the spec's input and output
lists are the statement, and the class that would close the gap is named here instead. A
sentence claiming the same fact in a spec is refused by that same tool, because a sentence is
a claim nothing holds.

- **The mechanical-design tree** — `process/mechanicaldesign/` and the `*MechanicalDesign`
  class beside each equipment class — is where most of them are, and no entry declares a
  parameter any of them read. **`DistillationColumnMechanicalDesign` is ported and wired** onto
  `process.distillation_column`: `crates/azoth-process/src/column/mechanical.rs` is `calcDesign`,
  `python/src/azoth/process/reference/mechanical.py` its twin, and the id now declares the six
  inputs it reads — `tray_efficiency`, `max_flooding_factor`, `tray_type`,
  `contactor_internals_type`, `material_grade` and `max_operation_pressure`, five of them class
  defaults — and publishes its eleven quantities, held to all three binary rows of
  `validation/neqsim/captures/process_column_mechanical_design.tsv` by cases on that id. Two of
  the eleven have no vocabulary unit — the weir loading in m3/(hr·m) and the tray drop in mbar —
  and are declared dimensionless on `total_pressure_drop_mbar`'s precedent. **`packed` is a
  `[[unported]]` row**: the class's packed arm is a second machine, refused by name in both
  languages. **`packed_column` is not one of the ids still owed this, and the reason is the
  refusal above**: its class builds the designer with `internalsType = "packed"`, so driving
  `calcDesign` on it takes exactly the branch this port refuses by name. **The absorber pair is
  landed too**, on `process.absorption_column` and `process.stripping_column`: the same six inputs
  and the same eleven quantities, held to the capture's `absorber_mechanical` row by a case on the
  absorber. **The geometry is not a declared surface on those two** — `AbsorptionColumn` carries
  no tray spacing, weir height or hole diameter of its own, so the two entries declare the six
  *class* inputs alone and the vessel is built at the designer's own defaults, which is what the
  capture's row is. **`stripping_column` has its own row too**: the capture's
  `hydrocarbon_stripper_mechanical` is built from the case's own recipe - it inherits every getter
  and overrides none, and one absorber row was not enough to claim the shape - and the case holds
  the eleven outputs to it.
  **`AbsorberMechanicalDesign` is a measured
  non-port**: its only `src/main` caller is `SimpleAbsorber.getMechanicalDesign()`, and
  `SimpleAbsorber extends Separator`, so driving it on an `AbsorptionColumn` throws — it is the
  separator family's sizing, which this tree already refuses. **`ReactorMechanicalDesign` is
  landed on `process.gibbs_reactor`**: eleven declared inputs and thirteen outputs, held to
  `validation/neqsim/captures/process_reactor_mechanical_design.tsv`'s three rows.
  **`PlugFlowReactor` is not owed it**: the class overrides no `getMechanicalDesign`, so it
  answers the base `MechanicalDesign`, whose `calcDesign` only reads design specifications and
  computes no vessel.
  **`TrayHydraulicsCalculator` is no longer one of them**: its arithmetic is the registered
  `hydraulics.tray_hydraulics` id.
- **`ColumnInternalsDesigner`** is two halves and both are in. Its **packing** side is done
  twice over: the report is `process.packed_column`'s own outputs, and the **sizing** —
  `calculatePacked`, `sizeColumnDiameter`, `roundToStandardDiameter` — is
  `hydraulics.packing_sizing`, which `calcPackingHydraulics`'s branch on a non-positive
  `column_diameter` now takes. Its **trayed** half — `calculateTrayed` and the getters it feeds,
  `getRequiredDiameter`, `getMaxPercentFlood`, `getTotalPressureDrop`, `getAverageTrayEfficiency`,
  `getTrayResults`, `getPackingResult` and `getControllingTrayIndex` — is
  `crates/azoth-process/src/column/designer.rs`, and its quantities are the outputs
  `process.distillation_column` publishes (`required_diameter`, `internals_design_ok`,
  `max_percent_flood`, `min_percent_flood`, `average_tray_efficiency`, `controlling_tray_index`,
  `total_pressure_drop`), held to `validation/neqsim/captures/process_internals_designer.tsv`.
  The pressure-drop coupling `DistillationColumn.updatePressureProfileFromHydraulics` performs
  when `hydraulicPressureDropCouplingEnabled` is set — which changes the tear's shape and not
  only its outputs — is ported too, on the same id.
- **`GravityDrainageMargin` and `PackingFoulingModel` are recorded as non-ports**, and the
  measurement is the record: `git grep` over `f0c7436:src/main` answers **one file for each,
  the class itself**, so neither has a `src/main` caller at the pin, and each has a single test
  of its own in `src/test`. **Neither name appears anywhere in the pinned jar**, which is where
  every capture in this tree came from, so neither is in the artifact the oracle was read from.
  A port would be a class nothing calls.
- **`MechanicalDesign`** — `process.ejector`'s own size, whose areas, diameters, lengths and
  volumes `updateDesign` fills. They are absent rather than guessed.
- **`process.filter`'s transient surfaces** — solids loading, breakthrough, backwash,
  regeneration and the element-integrity check. The four clean-drop models and
  `updateParticleCapturePerformance` are ported; what remains is the loading the transient
  accumulates, which is a dynamic state the five-field record has no slot for.
- **`getCapacityUtilization` is ported**, as a pair of optional inputs — `internal_diameter` and
  `design_gas_load_factor`, **both or neither**, because half a vessel has nothing to check — and
  an optional `capacity_utilization`. **It reads the system `run` entrained rather than the flash
  it started from**: `Separator.run` applies `addPhaseFractionToPhase` before the getter, which
  the capture's pair settles at `0.0033609430039118087` against the dry `0.003341057093792472`.
  Its two rules are its own — `rho_l` is the `oil` phase's, `aqueous`'s where there is no oil and
  `1000.0` where there is neither, which is not the columns' `10.0` floor, and a stream with no
  gas phase answers `0.0` — and NeqSim's `NaN` where the liquid is no denser than the gas is
  refused rather than carried.
- **The energy port** — `process.stirred_tank_reactor` runs `PHflash(inlet + duty)` when an
  `EnergyStream` is connected, and this library's five-field record has no place for one. The
  duty is declared out with the port named rather than guessed at.
- **`process.pipe`'s transient**: `AdiabaticPipe.runTransient` re-runs the steady state, so the
  line participates in a dynamic system without line-pack or wave propagation; this port has
  the steady state and no transient, with the rest of the tree. **`process.mixer`** owes
  the pressure-mismatch flag, which warns that the inlets arrived at different pressures — a
  fact about upstream equipment rather than the record — and, with the rest of the tree,
  `runTransient`.
- **`process.distillation_column`'s cold seed** — `runBostonSullivanRefinement`, declined on a
  measurement: without it the class's first Newton step moves the variables by `1e11` times their
  own size. **`process.shortcut_distillation_column`** owes the tray hydraulics, the tray
  optimization and `getResultsJson`; its own `setNumberOfTrays` override is inert, and says so.
- **`process.throttling_valve`'s laminar branch and controller** — the laminar form is the
  `prod choke` sizing class's and the controller a transient surface. The `ControlValveSizing`
  default, its inverse flow and `isCalcPressure` solve, the choked cap, the fouling fraction and
  the isothermal mode are ported; **the eight other sizing standards are refused by name** in
  the spec's `[[unported]]` rows, each naming its class.
- **The `eos` family's**: `ComponentGE.fugcoef`'s **Henry's-law branch is complete** - the
  four `eos.ge_*_phase` ids resolve both reference states, which the vendored table's
  `REFERENCESTATETYPE` selects on 115 of its 348 rows (113 `solute` and the two literal
  `0.0`). **Three states are refused rather than given a number**: a non-ion row carrying
  the four zeros NeqSim evaluates to `1.802 bar`, a solvent-tagged IAPWS species in water
  (upstream's one-component reference phase throws), and a solvent component with no
  Antoine row. **`SysNewtonRhapsonPHflash`**, which
  NeqSim's `ThrottlingValve` and `Compressor` select and which `eos.ph_flash` and `eos.ps_flash`
  decline on a measured two-phase `(dH/dT)_P` factor of two; **`PHflashSingleComp`** and
  **`PSflashSingleComp`**; **`VHflashQfunc`**'s coupled 2×2 with a cross term, and `VHflash`,
  which has no construction site; **`eos.tp_multiflash`'s seeding is complete** - the
  pure-component fallback, `stabilityAnalysisEnhanced` at its default post-flash call site,
  `seedAdditionalPhaseFromFeed` (a **GAS**), `seedHydrocarbonLiquidFromFeed` (an **OIL**) and
  the inline **aqueous** seed are all ported, and the one piece left is the opt-in pre-solve
  enhanced call, which needs `setEnhancedMultiPhaseCheck(true)` and has no model input; and
  `eos.soreide_whitson_phase`'s two newer parameterisations, **Chabab 2019** and
  **Burgoyne-Nielsen 2026**. The CPA rules and the hydrate coupling are P7's and P9's, each
  gated on a component family or a system type this crate does not carry.
- **`PhaseEos.getSoundSpeed` is `eos.sound_speed`** — the phase value, the beta-weighted system
  value and the phase. The system value doubles on an unflashed system (two phase objects, both
  `beta = 1.0`), which the port records and does not reproduce; `PhaseGE`'s constant `1470.0`
  m/s and `PhaseDefault`'s `NaN` are not models. **`EquilibriumSoundSpeed`**, the re-flash
  finite-difference utility for two-phase screening, is owed.
- **The reactions family's**: the **ionic branch** of `reactions.reactive_tp_flash` and
  `reactions.reactive_ph_flash`, which is why a charged substance is refused at all;
  `reactions.reactive_ph_flash`'s **entropy-specified variant** and its `NR = 0` delegation to
  a standard PH flash; **`useAdaptiveDerivatives`** on `reactions.chemical_equilibrium`, which
  `solveChemEq` sets from the second refinement onward — azoth has that derivative, and a
  procedure over vectors has no phase to reach it through; and **`SimplexSolver`**, whose
  estimated vertex `LinearProgrammingChemicalEquilibrium` seeds through, where this port solves
  the program exactly.
- **The steady-state entries, and the transient and design surface beside each.** Every
  `process.*` class NeqSim carries also carries a `runTransient`, and most carry a mechanical
  design and a performance or control surface the steady `run` never reads. Owed entry by
  entry: `Splitter`'s capacity map; `Flare`'s radiation model and capacity check; `Heater`'s
  electrical and instrument design and its bounded outlet temperatures; `Compressor`'s
  catalogue and readers, the anti-surge recycle calculators and controllers (its surge curve,
  read-outs and run-time recycle block have landed), deposit model and
  `useGERG2008`/`useLeachman`/`useVega` property overrides - **its two alternative map
  lookups, the `CompressorChartKhader2015` map and the speed-solve mode have landed**;
  `Expander`'s `powerSet`/`useOutTemperature` modes, dynamic nozzle and capacity
  constraints; `HeatExchanger`'s
  thermal-hydraulic rating, `useDeltaT` variant and correlation-based UA; `Separator`'s
  capacity check and droplet-performance calculator; `ThreePhaseSeparator`'s levels and
  outlet-valve fractions; `GasScrubber`'s vessel sizing and mesh pad; `Tank`'s wall
  temperature, steel mass and auto-sizing; `ThrottlingValve`'s laminar branch and controller
  path; `Filter`'s solids loading, breakthrough, backwash, regeneration and element-integrity
  check; and `StirredTankReactor`'s vessel volume, residence time and agitator power. **The
  turbo-expander landed**: `TurboExpanderCompressor` and `ExpanderChartKhader` are one coupled
  machine, the chart's only reader. **`MapTurboExpanderCompressor` landed as its own id**,
  `process.map_turbo_expander_compressor`: the second shaft, its speed bisected over the same
  chart's measured lines. What is still owed is the map surface around it:
  `RadialExpanderGeometryMap`, the synthetic-map generator whose consumers
  (`TurboExpanderMapIngestion`, `TurboMachineryChartLibrary`) nothing here calls; both
  classes' `toJson`/`toMap` reports and the map machine's `feasible` flag. **The expander
  chart's two stated sound speeds are closed**: its own surface reads them from
  `eos.sound_speed`, and the class's `830.84` m/s is recorded as its unflashed doubling.
  **Anti-surge landed** - the two boundary splines, the chart's
  generated surge and stonewall curves, the machine's surge read-outs and `run`'s recycle
  block. What it left owed is the steady `AntiSurgeRecycleCalculator`, the legacy
  `Calculator.runAntiSurgeCalc` step, `AntiSurgeController` and the dynamic/application layer
  (`AntiSurgeDynamicBenchmark`, `AdvancedAntiSurgeControlSystem`, `CompressorAntiSurgeApplication`,
  `CompressorOperatingHistory`, `RecycleFlowCoordinator`).
- **Three behaviours of the same shape, in `eos` and `reactions`.** `RachfordRice.calcBeta`
  clamps its root to `[1e-12, 1 - 1e-12]` and declines to iterate where no root exists - the
  port returns the root the equation gives, so the negative flash stays visible;
  `HydrateInhibitorConcentrationFlash` and its `...wtFlash` sibling are the models that add
  methanol or MEG to `eos.hydrate_formation_temperature`; and the `EPS` pinning and
  `getLogInfiniteDiluteFugacity` reference-state correction behind
  `reactions.reactive_tp_flash`'s ionic refusal. `Separator`'s `oilInGas` and `waterInGas` go
  with them: no oil or aqueous phase exists here, and they are the only fields that re-run the
  *vapour*.

- **The characterisation methods that rebuild a fluid rather than a table.** The arithmetic is
  ported and the fluid is built from it: `PseudoComponentCombiner`'s `combineReservoirFluids`,
  `characterizeToCommonSlate` and `transferBinaryInteractionParameters` are
  `characterization.combine_reservoir_fluids`, `characterization.characterize_to_common_slate`
  and `characterization.transfer_binary_interaction_parameters`, each over rows;
  `characterization.tbp_fraction_component` is `addTBPfraction`'s card for one cut, its Rackett
  compressibility a reference SRK liquid root; and `mixture_from_characterization` turns rows
  into a `Mixture` and its ideal-gas model. So the rows are the whole of the fluid's
  description, and three engine-side quirks travel with them: `getCpD`'s `1e-10` offset, the
  acentric factor `setm` rewrites, and the coefficient `init` then re-derives from it.
  `characterization.characterize_to_reference` carries all four of the class's own options -
  `inheritReferenceProperties`, `delumpBeforeRecharacterization` with its resolution, and
  `sharedImaginaryBoundaries`, which is `determineReferenceEqualMassBoundaries`. They are
  composed into a fluid by `mixture_from_characterization`, and an assay reaches it through
  `oil_assay_to_mixture` - the derivation, the cards and the assembler in that order. **What is
  left of the subsystem is what builds nothing**: `Recombine` is unreachable (its constructor is
  empty and `runRecombination` returns null), `normalizeComposition` computes a total and never
  uses it, and `generateValidationReport` reports on data; none is an id.

## Beyond the port

**The interoperation surface is built.** It is `crates/azoth-process`'s `middleware` module and
three bindings of it — `azoth-wasm` for a browser, `azoth.process.Session` for a notebook, and
`azoth edit`/`azoth forms` for a shell. Its shape is
[The middleware](docs/src/architecture/middleware.md), and **what is still not built is named
there** rather than counted here. The editor it feeds is `ui/`.

**Three holes in the checker closed on the way, and each was the same shape.** The checker
verified units, ports, names and connections, and never that a value was *there*, was of the
*right kind*, or *named* something the run could do: `MissingParameter` (closed with P12),
`ParameterKind`, and `Acceleration` — a tear's `acceleration_method` reached `executor::session`
and nothing else, so a misspelt name, or `broyden`, validated clean and then failed the run. The
first two were reachable from a form; the third became reachable when the command model gave a
front-end a way to write one.

**Nothing in this tier is waiting on a prerequisite, and the three that were are now cases.**
The paragraph here used to name them: `pipe` waited on a density and a viscosity assembled from
a `Mixture`, and `compressor`/`expander` waited on a molar-entropy field on the stream record.
All three landed - `pipe` with the two-density formula its own capture pins, and the two
isentropic machines on `Stream::entropy` - so what the paragraph recorded as a blocker is now
three registered ids with oracles.

**P12 is landed.** The executor is `crates/azoth-process`'s `executor` module, and it needed no
calculation this tier owes - its work was ordering, a tear's fixed point, a session and a codec,
each a use of machinery the crate already carried. It runs `specs/flowsheets/demo.toml` and
converges its recycle, held to a NeqSim `ProcessSystem` capture, and every instance's answer is
cross-checked against that unit's own registered model, which is the claim "the executor calls
the kernels" made testable.

**A flowsheet is a self-contained simulation**, and the document's shape is
[the process schema](docs/src/architecture/process-schema.md)'s. What this tier added is that it
runs with no argument but the file: from the library, from `azoth run --flowsheet` and from
`azoth.process.run_flowsheet`, whose `feeds` parameter is an override of a boundary the document
already states. The checker gained the shape rule that holds an input's record to being one, and
the shipped document's own input is held to the NeqSim capture's feed rather than agreeing with
it by memory.

**The low-flow cutoff is ported, and it closed a divergence rather than adding a feature.**
`Recycle.deactivateOnLowFlow` switches a tear off when its inlet carries less than `minimumFlow`
kg/hr: the four residuals are *declared* zero rather than measured, the tear is marked inactive,
and `solved()` answers true for it - so an empty loop stops at one pass instead of running until
the zero-flow floor closes it a pass later. The class's `runUnitProfiled` then skips an inactive
unit, so the extra pass the outer loop still insists on does not advance the tear's own counter.
That is why the capture beside this reads `recycle_iterations=1` for the shipped document while
the loop around it runs twice, and the port now reproduces both numbers.

**The tear's acceleration is carried, and Broyden is refused with the measurement that closes
it.** The two methods and `[[recycles]]`'s own convergence fields are
[the process schema](docs/src/architecture/process-schema.md)'s. What is recorded here is the
evidence. Wegstein is wired where `Recycle.run` applies it, held to
`captures/process_flowsheet_accelerated.tsv` at `1.9e-11` relative, and both methods' arithmetic
is held step by step to `captures/process_acceleration.tsv` - which is where each one's delay
boundary and Wegstein's `q` factors are visible. Broyden's refusal is measured rather than
conservative: the class writes its accelerated composition through `Component.setx` on both
phases, a two-phase system does not read back what it was written, and on the condensing graph at
a tolerance the loop never closes at, NeqSim runs all hundred passes at `8.6243` mol/s where the
transcribed step reports convergence after eighteen at `0.0833`. That capture also records the
class's Broyden step taking the wrong sign against its own comment - `dx = -B^{-1} f` written,
`+B^{-1} f` computed - which walks away from a fixed point one Newton step would reach.

**A `many` outlet is addressable by position, so a splitter can be wired.** The grammar is
[the process schema](docs/src/architecture/process-schema.md)'s. What it buys is a **purge**,
and a purge is what gives a loop a steady state: the test over it wires half a split back to the
mixer and half to a boundary, and at a tolerance that means it, the loop closes at the feed's own
flow where `demo.toml` at any tolerance cannot.

**A parameter's declaration says whether a kernel can run without it, and the checker reads it.**
That was the last hole P12 named, and it is the one that let `specs/flowsheets/demo.toml` ship:
`unit_ops.separator` requires `gas_in_liquid`, the file never supplied it, and `validate` printed
`OK` for as long as the file existed - so the checker's own promise, that a flowsheet printing
`OK` is one an executor can consume, was false. `Param` now carries `required`, as
[the process schema](docs/src/architecture/process-schema.md) states, and the flag's two sources
agree on all twenty-eight runnable entries - which `tests/palette.rs` re-checks on every run.

**What is named as not built, with the class that would close each**: Broyden through
`BroydenAccelerator`, for the measured reason above. The tier's own acceptance is the NeqSim
capture beside it: an executor whose convergence is only checked against its own arithmetic is an
executor nobody has measured.
