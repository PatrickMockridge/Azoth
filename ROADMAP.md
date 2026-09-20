# NeqSim port roadmap

The order in which azoth ports [NeqSim](https://github.com/equinor/neqsim) 3.20.0, from
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
- **PVT and flow assurance.** `pvtsimulation/` (simulation, model tuning, reservoir
  properties).

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
- **Duan-Sun, reachable and unimplementable.** `ComponentGeDuanSun`, `PhaseDuanSun`,
  `SystemDuanSun` and `thermo/util/empiric/DuanSun.java` are **not ported, and no state
  they accept exists**: `SystemDuanSun.addComponent` throws for every name but `CO2`, and
  the phase it would build divides by the moles and molar mass of a component named
  `water`, which a system that admits only `CO2` does not have. There is nothing to
  reproduce (NeqSim issues 3837 and 3839), and this is not the carried-as-unreachable case
  below, because `SystemThermo`'s model-name factory does construct it.
- **Electrolyte-CPA.** `PhaseElectrolyteCPA`, `PhaseElectrolyteCPAAdvanced`,
  `PhaseElectrolyteCPAMM`, `PhaseElectrolyteCPAOld` and their `SystemElectrolyte*` are
  **carried as unreachable upstream**: 5,727 lines with no non-test construction site
  anywhere — each phase is built only by its own `SystemElectrolyte*`, and
  `SystemThermo`'s factory maps `Electrolyte-CPA-EOS` to `CPAstatoil` and never to
  `SystemElectrolyteCPA`. `SystemElectrolyteCPAstatoil` and
  `PhaseElectrolyteCPAstatoil` are **not ported** — 145 lines, eleven `src/main` sites
  that build them, and Fürst plus the Wertheim association this library already has.
- **Salts and scale.** `MultiSaltPrecipitation`, `CalcSaltSatauration`,
  `CheckScalePotential`, `AddIonToScaleSaturation`.
- **Freezing.** `FreezeOut`, `FreezingPointTemperatureFlash`.

### Tier 3 — wax, hydrate, hydrogen, asphaltene

The specialist physics.

- **Hydrate.** `ComponentHydrate`, `ComponentHydrateBallard`, `ComponentHydrateGF`,
  `ComponentHydrateKluda`, `ComponentHydratePitzer`, `ComponentHydratePVTsim`,
  `ComponentHydrateStatoil`, `PhaseHydrate`, `TPHydrateFlash`,
  `HydrateFormationPressureFlash`, `HydrateFormationTemperatureFlash`,
  `HydrateInhibitorConcentrationFlash`, `HydrateInhibitorwtFlash`, `PitzerHydrateFlash`,
  `HydrateEquilibriumLine`, `HydrateEquilibriumDiagnostics`, `OLGAhydrateCurveGenerator`.
- **Wax.** `ComponentWax`, `ComponentWonWax`, `ComponentCoutinhoWax`,
  `ComponentWaxWilson`, `PhaseWax`, `TPmultiflashWAX`, `WaxCharacterise`,
  `WaxModelInterface`, and `pvtsimulation/flowassurance/WaxCurveCalculator`.
- **Asphaltene.** `AsphalteneCharacterization`, `PedersenAsphalteneCharacterization`,
  `AsphalteneOnsetPressureFlash`, `AsphalteneOnsetTemperatureFlash`, and
  `pvtsimulation/flowassurance/FloryHugginsAsphalteneModel`.
- **Hydrogen and cryogenic.** `ComponentGERG2008Eos`, `ComponentGERG2004`,
  `PhaseGERG2008Eos`, `SystemGERG2008Eos`, `PHflashGERG2008`, `PSFlashGERG2008` and
  `thermo/util/gerg/`; `ComponentLeachmanEos`, `PHflashLeachman`, `PSFlashLeachman` and
  `thermo/util/leachman/`; `thermo/util/hydrogen/ParaOrthoH2Correction`;
  `thermo/util/solid/ParaHydrogenSolidHelmholtzEquation`.

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
- **Solids.** `ComponentSolid`, `ComponentSolidHelmholtzEos`, `PhaseSolid`,
  `PhaseSolidComplex`, `PhasePureComponentSolid`, `PhaseSolidHelmholtzEos`,
  `SystemSolidHelmholtzEos`, `SystemArgonSolidHelmholtzEos`, `SolidFlash`, `SolidFlash1`,
  `SolidFlash12`, `PHsolidFlash`, `thermo/util/solid/`.
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
