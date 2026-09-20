# The skills roadmap

The skill catalog grows in two directions that move at different speeds. The
cross-cutting and `azoth`-basis skills are written as the library implements the
calculation they drive; the rest of NeqSim's catalog is a set of placeholders that
become real in the same commit that lands their backing tranche. This page records
the mapping, so an agent can tell "validated" from "placeholder" at a glance.

## What is blocked, and why

A skill whose `calculation_basis` is `azoth` drives the validated library. Most of
NeqSim's 77 community skills are thermo skills, and the library does not yet compute
the physics behind them:

- **Associating models (CPA, PC-SAFT, SAFT-VR-Mie)** — no flow-assurance hydrate or
  process skill can be `azoth`-basis until tranche P7.
- **Electrolytes** — P8 is closed: `eos.pitzer_phase`, `eos.kent_eisenberg_phase`,
  `eos.desmukh_mather_phase`, `eos.soreide_whitson_phase` and the Fürst pair
  (`eos.furst_electrolyte_phase`, `eos.furst_electrolyte_mod2004_phase`) are ported, and no
  skill drives them. The brine and scale skill waits on P9 and both dehydration skills on
  P11, so this is the one closed tranche that promotes nothing: the physics a scale
  calculation is built on is here, and the calculation is not.
- **Activity-coefficient models (NRTL, UNIFAC, UNIQUAC, Wilson, Van Laar)** — P5 is
  closed: all five activity models are ported, and all five have a phase
  (`eos.ge_nrtl_phase`, `eos.ge_unifac_phase`, `eos.ge_uniquac_phase`,
  `eos.ge_wilson_phase`, `eos.ge_van_laar_acid_phase`) with the NRTL one also carrying a
  gamma-phi flash. One of the five phases has no differential oracle, and the reason is a
  withdrawal rather than an absence: upstream commit `c5ec5fb` gave Wilson a published
  coefficient, which is now compared against, and withdrew `PhaseGEUniquac` outright —
  its `getGamma` throws and both standalone constructors refuse — so UNIQUAC still cannot
  be compared against anything.
- **Unit operations and the flowsheet** — the 30 process skills wait on P11 and P12,
  which the specification places after the physics tiers.
- **Reference equations of state (GERG, IAPWS-95, Span-Wagner, …)** — pvt benchmark
  and reference-data skills wait on P4.
- **Transport properties** beyond heat capacity and heat of vaporisation — viscosity,
  conductivity, diffusivity, surface tension, density — wait on P1.
- **Petroleum-fraction characterisation** — pvt pseudocomponent and regression skills
  wait on Tier 1 of [`ROADMAP.md`](../../../ROADMAP.md).
- **Hydrates, wax, asphaltene and solids** — flow-assurance skills wait on P9.

## Tranche → domain → skills

| Domain | Basis now | Backed by tranche |
|---|---|---|
| `library` | `azoth` | P0 (done) |
| `eos`, `hydraulics`, `thermal` | `azoth` | P0–P2 (done) |
| `pvt` | `screening` | P1, P3 (mixing), P4 (reference EOS), Tier 1 (characterisation) |
| `process` | `screening` | P1, P5, P6, P10, P11 (unit ops), P12 (flowsheet) |
| `flow-assurance` | `screening` | P7, P9 |
| `safety` | `screening` | P1, P6, P11 |
| `subsurface` | `screening` | Tier 1 (PVT) |
| `environment` | `advisory` | none needed |
| `field-development` | `advisory` | none needed |
| `subsea` | `data-retrieval` | none needed |
| `reporting` | `advisory` | none needed |
| `engineering-data` | `data-retrieval` | none needed |

The tranches are the P0–P12 order in
[The specification](../architecture/specification.md), mapped one-to-one to NeqSim's
classes in [`ROADMAP.md`](../../../ROADMAP.md). `(done)` means the `azoth`-basis skills
for that domain are written, not that every class in the tranche is ported — the 25
alpha functions of P2, for example, are still ahead.

## The two speeds

- **Now, at `azoth`-basis**: the four `azoth`-basis `library/` skills, and the `eos/`,
  `hydraulics/` and `thermal/` skills that drive the 43 implemented calculations.
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
