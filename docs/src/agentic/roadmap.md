# The skills roadmap

The skill catalog grows in two directions that move at different speeds. The
cross-cutting and `azoth`-basis skills are written as the library implements the
calculation they drive; the rest of the catalog is a set of placeholders that
become real in the same commit that lands their backing tranche. This page records
the mapping, so an agent can tell "validated" from "placeholder" at a glance.

## What is `azoth`-basis, and what is still a placeholder

A skill whose `calculation_basis` is `azoth` drives the validated library. Nineteen
are: the four `library/` skills, four under `eos/`, four under `flow-assurance/`,
two under `hydraulics/`, two under `process/` and one under `thermal/`, driving the
86 calculations and 171 models the library implements.

The remaining 71 are placeholders — 43 `screening`, 18 `advisory` and 10
`data-retrieval` — and each `screening` one names **what would back it**,
which is a tranche where a tranche covers the physics and a named family where none does:

- **What a tranche covers** — 7 skills under P4 (six under `pvt/`,
  `azoth-near-well-and-injectivity` under `subsurface/`); and
  `azoth-water-dewpoint-dehydration-screening` under P1, whose correlation P1's tree
  covers and did not port.
- **What no tranche covers** — the specification puts **mechanical design**, `safety/`,
  `statistics/`, `automation`, `standards/` and `fluidmechanics/` *beyond* P12, so a skill
  whose physics belongs to one of those names the family rather than the nearest
  P-number. 16 are mechanical design (vibration, noise, wall thickness, flexibility,
  seals, turbines, vessel sizing), 9 are `safety/`, 3 `statistics/`, 3
  `automation`, 2 `standards/` and 1 `fluidmechanics/`.
- **What is a unit operation and a flowsheet away** — `azoth-teg-dehydration-modeling`
  keeps the P11 label, because P11's palette is the tranche that built the absorber,
  flash, column, stripper and recycle it composes. What it does not have is the plant:
  no such flowsheet ships, and the unit operations it would stand on are `unverified`.

**A skill is promoted when the library computes its arithmetic, not when a tranche
closes.** Reading the two as the same is what left 39 skills labelled P11 after P11 and
P12 had both closed — a `screening` skill *had* to name a tranche, so a noise screen and
a PSV orifice calculation named the nearest one.

Three physics families that once blocked a set of skills have landed:

- **Activity-coefficient models (NRTL, UNIFAC, UNIQUAC, Wilson, Van Laar)** — P5
  is closed: all five activity models are ported, and all five have a phase
  (`eos.ge_nrtl_phase`, `eos.ge_unifac_phase`, `eos.ge_uniquac_phase`,
  `eos.ge_wilson_phase`, `eos.ge_van_laar_acid_phase`) with the NRTL one also
  carrying a gamma-phi flash. One of the five phases has no differential oracle,
  and the reason is a withdrawal rather than an absence: upstream commit `c5ec5fb`
  gave Wilson a published coefficient, which is now compared against, and withdrew
  `PhaseGEUniquac` outright — its `getGamma` throws and both standalone
  constructors refuse — so UNIQUAC cannot be compared against anything.
- **Associating models (CPA, PC-SAFT, SAFT-VR-Mie)** — P7 has landed the phases:
  `eos.pr_cpa_phase`, `eos.srk_cpa_phase` and `eos.umr_cpa_phase` for the cubic
  association, and `eos.pcsaft_rahmat_phase`, `eos.saft_vr_mie_phase` and
  `eos.tp_flash_saft` for SAFT.
- **Electrolytes** — P8 is closed: `eos.pitzer_phase`, `eos.kent_eisenberg_phase`,
  `eos.desmukh_mather_phase`, `eos.soreide_whitson_phase` and the Fürst pair
  (`eos.furst_electrolyte_phase`, `eos.furst_electrolyte_mod2004_phase`) are
  ported. Both dehydration skills still wait on P11.

P9 has landed hydrate, wax, scale and freezing, and four skills were promoted to
`azoth`-basis on it: `azoth-hydrate-margin-check` and `azoth-hydrate-screening`
(`eos.hydrate_formation_temperature`, `eos.hydrate_formation_pressure`,
`eos.hydrate_fraction`), `azoth-wax-margin-check` (`eos.tp_multiflash_wax`,
`eos.wax_solid_fugacity`) and `azoth-produced-water-scale-screening`
(`eos.scale_saturation_ratio`, `eos.salt_precipitation`).

**The scale family is now two routes, seven water-analysis ids beside the reactive two.**
`eos.scale_prediction`, `eos.scale_mass`, `eos.pitzer_scale_activity`,
`eos.barite_celestite_solid_solution`, `eos.multi_mineral_scale_equilibrium`,
`eos.flowline_scale_profile` and `eos.water_compatibility` take a water analysis - mg/L, a
temperature and a pressure - and form their own Davies, B-dot or binary-Pitzer activities,
where the reactive route takes a speciated brine and the phase's own coefficients. **Their
auto-pH path carries a defect upstream and reproduces it**: NeqSim's CO2 Henry constant puts
the aqueous estimate near `1e48`, so an auto-pH carbonate index lands about `-52` and a
manual pH is the physical route. The specs pin that value, so a future repair upstream is a
visible divergence rather than a silent one.
**The hydrate ids now cover both of NeqSim's routes**: the three fluent-cubic ids above,
and `eos.pitzer_hydrate_formation_temperature` with `eos.pitzer_hydrate_formation_pressure`
for a `SystemPitzer` brine. The latter run `PitzerHydrateFlash` over the hybrid EoS/GE flash,
whose water reference is the pure-water Antoine pressure; a salty CO2 state takes the manual
PHREEQC CO2-chloride dataset through the `pitzer_dataset` input, because the catalogue does
not carry the CO2-cation-Cl zeta rows and the upstream caller supplies them as zero.
**Three of the flow-assurance skills wait on no tranche at all**:
`azoth-two-phase-flow-regime-screening` and `azoth-multiphase-flow-slug-screening`
are `advisory` — their maps are NeqSim's `fluidmechanics/`, which is not a port
source — and `azoth-olga-multiphase-simulator` is `data-retrieval`, because it
drives a commercial simulator rather than computing anything.
**Asphaltene is carried rather than pending**: the reachable pair of onset flashes
scans over a flash that upstream collapses to the feed, so there is nothing to
drive, and [`ROADMAP.md`](../../../ROADMAP.md) records the defect and what would
have to change. No skill is blocked on it.

**Petroleum-fraction characterisation** is the fourth block, and it is not a
tranche's: the `pvt/` pseudocomponent and regression skills wait on Tier 1 of
[`ROADMAP.md`](../../../ROADMAP.md).

## Tranche → domain → skills

The `tranche` field in `skills.toml` is the machine-readable form of this table, and
`tools/validate_skills.py` requires a `screening` skill to name what would back it —
a P-number, or one of the beyond-P12 families the specification names. A domain's
`azoth`-basis skills carry none, because they drive the library rather than wait on it:

| Domain | Basis now | Backed by |
|---|---|---|
| `library` | `azoth` | P0 (done) |
| `eos`, `hydraulics`, `thermal` | `azoth` | P0–P2 (done) |
| `process` | `azoth`, `screening` | P11 (done), and the beyond-P12 families |
| `safety` | `screening` | `safety` — beyond P12, so no tranche |
| `pvt` | `screening` | P4 |
| `subsurface` | `screening` | P4 |
| `flow-assurance` | `azoth`, `screening`, `advisory`, `data-retrieval` | P9 |
| `subsea` | `data-retrieval` | none needed |
| `environment` | `advisory` | none needed |
| `field-development` | `advisory` | none needed |
| `reporting` | `advisory` | none needed |
| `engineering-data` | `data-retrieval` | none needed |

The P-numbers are the P0–P12 order in
[The specification](../architecture/specification.md), mapped one-to-one to NeqSim's
classes in [`ROADMAP.md`](../../../ROADMAP.md); the families are the trees that page puts
beyond P12. `(done)` means the `azoth`-basis skills for that domain are written, not that
every class in the tranche is ported.

## How a placeholder becomes real

- **At `azoth`-basis**: the skills named above.
- **As the calculation that backs it lands**: the `screening` placeholder is promoted to
  `azoth`-basis in that commit, its `tranche` field is removed, and its `SKILL.md` is
  rewritten to drive the ids — a promotion that left the placeholder text in place would
  be a page claiming a number nothing produced.

A `screening` skill is not a promise the library will compute the result — it is an
honest marker that it does not yet. Until its calculation lands, its
`Related Azoth functionality` section names the NeqSim class that does, today.

## The agentic surface

The tool schema and session an agent drives — the reflection surface an editor and a
notebook share with a model — is built: `middleware::tools`, one tool per command,
projected from the command model rather than written beside it. It is stated in
[The middleware](../architecture/middleware.md), not here: a skill is an
instruction an agent reads, a tool schema is the machine-readable form of the same
commands.
