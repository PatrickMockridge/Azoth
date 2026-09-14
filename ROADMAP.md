# NeqSim port roadmap

The order in which azoth ports [NeqSim](https://github.com/equinor/neqsim) 3.20.0, from
the physics the most callers need to the physics the fewest do. Impact first: a property
or model that ordinary oil, gas and water work uses is ported before one only a
specialist reaches for.

This is the tiered view over `databank/manifest.toml`, which is the column-by-column
record of what is carried and what is not — 35 of NeqSim's resource files, 1,496
columns, 1,459 carried of which 14 are read today. This file orders that backlog; it
does not re-inventory it.

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

The cubic Peng-Robinson core is ported: TP, PH and PS flash, the tangent-plane
stability test, the mixture critical point, the phase envelope, Rachford-Rice and the
classical vdW1f mixing rule, reading 14 of the databank's columns. What remains is the
foundation the rest stands on:

- the single `databank → keycard → every calculation` path, so a calculation names its
  components and every constant it reads comes through one route;
- the `not-yet` columns the manifest names for `eos.ideal_gas_cp` and
  `eos.molar_enthalpy_entropy` — the ideal-gas Cp polynomial and the reference-state
  quantities those two models still take from the caller.

### Tier 1 — oil and gas

The cubic family in breadth, the standard component correlations, petroleum-fraction
characterisation, and the transport properties every process model needs.

- Cubic EOS beyond PR: SRK and RK, Peneloux volume-shift
  (`ComponentSrkPeneloux`, `PhaseSrkPenelouxEos`), and the alpha variants
  (`SystemSrkTwuCoonEos`, `SystemPrMathiasCopeman`, `SystemPrEos1978`,
  `SystemPrLeeKeslerEos`, `SystemSrkSchwartzentruberEos`, `SystemPrGassemEos`), with the
  mixing-rule handler beyond classical vdW (`EosMixingRuleHandler`, Huron-Vidal and
  Wong-Sandler).
- Generic component correlations read from `thermo/component/` (normal boiling point,
  Antoine vapour pressure) and `physicalproperties/` (standard liquid density, heat of
  vaporisation, liquid heat capacity) — the columns the manifest already lists under
  those packages.
- Petroleum-fraction characterisation: `Characterise`, `TBPCharacterize`,
  `PlusFractionModel`, `LumpingModel` in `thermo/characterization/`.
- Transport properties: Chung viscosity and conductivity, Rackett and Costald density,
  and `ParachorSurfaceTension`.

### Tier 2 — water and gas-water

Water and the aqueous models: the water EOS, water content and dehydration, acid-gas
solubility, and the electrolytes and hydrate inhibitors that water work depends on.

- The water EOS: `ComponentWater`, `PhaseWaterIAPWS`, `Iapws_if97`.
- Water content and dehydration: `WATcalc`, `WaterDewPointTemperatureFlash`.
- Acid gas: `ComponentSoreideWhitson`, `CO2BrinePhaseEquilibrium`, `SaturateWithWater`.
- Electrolytes and salts: `PhasePitzer`, `PhaseKentEisenberg`, `PhaseDesmukhMather`,
  `PhaseDuanSun` and the `Pitzer*` machinery, with MEG, methanol and glycol inhibition.

### Tier 3 — wax, hydrate, hydrogen, asphaltene

The specialist physics.

- Hydrate: `PhaseHydrate`, `TPHydrateFlash`, `HydrateFormationPressureFlash`,
  `HydrateInhibitorConcentrationFlash`.
- Wax: `ComponentWax`, `ComponentWonWax`, `ComponentCoutinhoWax`, `TPmultiflashWAX`,
  `WaxCharacterise`.
- Asphaltene: `AsphalteneOnsetPressureFlash`, `FloryHugginsAsphalteneModel`.
- Hydrogen and cryogenic: `ComponentGERG2008Eos` and `GERG2008`, `ComponentLeachmanEos`,
  `ParaOrthoH2Correction`, `ParaHydrogenSolidHelmholtzEquation`.

### Tier 4 — the long tail

Port on demand, as a caller needs them:

- CPA (`ComponentPrCPA`, `CPAMixingRuleHandler`), PC-SAFT and SAFT-VR-Mie
  (`ComponentPCSAFT`, `ComponentSAFTVRMie`), the Span-Wagner reference EOS
  (`ComponentSpanWagnerEos`), solids (`ComponentSolid`, `SolidFlash1`), sulfur
  (`SulfurThermodynamics`), amines (`AmineKentEisenberg`), and the reactive flash.
