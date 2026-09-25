# The skills roadmap

The skill catalog grows in two directions that move at different speeds. The
cross-cutting and `azoth`-basis skills are written as the library implements the
calculation they drive; the rest of the catalog is a set of placeholders that
become real in the same commit that lands their backing tranche. This page records
the mapping, so an agent can tell "validated" from "placeholder" at a glance.

## What is `azoth`-basis, and what is still a placeholder

A skill whose `calculation_basis` is `azoth` drives the validated library. Fifteen
are: the four `library/` skills, four under `eos/`, four under `flow-assurance/`,
two under `hydraulics/` and one under `thermal/`, driving the 67 calculations and
93 models the library implements.

The remaining 75 are placeholders — 47 `screening`, 18 `advisory` and 10
`data-retrieval` — and each `screening` one names the single tranche that will
back its physics:

- **Unit operations and the flowsheet** — 39 skills wait on P11: 30 under
  `process/` and 9 under `safety/`. The specification places the tier after the
  physics, and the kernels are nearly closed: `specs/unit_ops/` declares 29 palette
  entries and 26 carry a kernel, so the skills are waiting on the P12 executor rather than
  on the tier's existence.
- **Reference equations of state, and the PVT chain** — 7 skills wait on P4: six
  under `pvt/` and `azoth-near-well-and-injectivity` under `subsurface/`.
- **Cooldown** — `azoth-surf-cooldown-screening` is the one skill waiting on P9,
  which it keeps for the hydrate temperature it compares against. The cooldown
  itself is long-tail work under Tier 4 rather than a tranche's.

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
**The hydrate ids carry a caveat worth stating**: they reproduce NeqSim's non-Pitzer
route, and NeqSim answers a `SystemPitzer` brine from `PitzerHydrateFlash` instead
— a path that needs a flash whose two phases run different models, which this
library does not have. So a hydrate margin on a brine is not what these three
compute, and the difference is the electrolyte coupling rather than the hydrate
physics.
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

The `tranche` field in `skills.toml` is the machine-readable form of this table,
and `tools/validate_skills.py` requires a `screening` skill to name exactly one.
A domain's `azoth`-basis skills carry none, because they drive the library rather
than wait on it:

| Domain | Basis now | Backed by tranche |
|---|---|---|
| `library` | `azoth` | P0 (done) |
| `eos`, `hydraulics`, `thermal` | `azoth` | P0–P2 (done) |
| `process` | `screening` | P11 |
| `safety` | `screening` | P11 |
| `pvt` | `screening` | P4 |
| `subsurface` | `screening` | P4 |
| `flow-assurance` | `azoth`, `screening`, `advisory`, `data-retrieval` | P9 |
| `subsea` | `data-retrieval` | none needed |
| `environment` | `advisory` | none needed |
| `field-development` | `advisory` | none needed |
| `reporting` | `advisory` | none needed |
| `engineering-data` | `data-retrieval` | none needed |

The tranches are the P0–P12 order in
[The specification](../architecture/specification.md), mapped one-to-one to NeqSim's
classes in [`ROADMAP.md`](../../../ROADMAP.md). `(done)` means the `azoth`-basis skills
for that domain are written, not that every class in the tranche is ported.

## How a placeholder becomes real

- **At `azoth`-basis**: the fifteen skills named above.
- **As a tranche lands**: the `screening` placeholder for a domain is promoted to
  `azoth`-basis in the same commit, and its `tranche` field is removed.

A `screening` skill is not a promise the library will compute the result — it is an
honest marker that it does not yet. Until its tranche lands, its
`Related Azoth functionality` section names the NeqSim class that does, today.

## The agentic surface

The tool schema and session an agent drives — the reflection surface an editor and a
notebook share with a model — is gated on P12, like the unit-operation tier it describes.
It is stated in [The middleware](../architecture/middleware.md), not here: a skill is an
instruction an agent reads, a tool schema is the machine-readable form of the same
commands.
