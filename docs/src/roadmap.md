# The roadmap, and the boundary it moves

<!-- Hand-written, and listed in tools/gen_docs.py's STATIC_PAGES. mdBook silently
     drops a page the generated summary does not name, so the registration is what
     makes this page exist rather than a note in a plan file. -->

azoth is porting a **Pareto-efficient core** of [NeqSim](https://github.com/equinor/neqsim)
— the unit operations and user features that make a thermodynamic library usable, not the
specialist physics. This page is where that boundary is written down, so a reader finds
out what is missing on purpose rather than one failed call at a time.

[azoth and NeqSim](./comparison/neqsim.md) is the companion page: what the two libraries
*are*, and how their design choices differ. This one is what is being built and when.

## The boundary this moves

`docs/src/index.md` says *"each calculation is independent, and the composition is done
by the caller"*, and until now that was a **refusal**, not just a description. This page
used to be a paragraph in the comparison explaining why azoth would not grow a flowsheet:
a flowsheet has no single retraceable worked example, which is the guarantee every other
page in this book is built on.

**That guarantee is being traded deliberately, and for something narrower.** What
survives:

- A **unit operation** is still a model with a spec, a source, and a worked example a
  human can retrace. Nothing about the process layer relaxes that.
- A **flowsheet** carries a hand-computable case — a mixer into a flash, no recycle —
  plus conservation checks that hold at every answer rather than at one recorded one.
  That is the honest analogue of a worked example, and it is weaker. Saying so is the
  point of this section.
- **Two implementations, always.** A flowsheet is a document both languages execute and
  compare, with iteration counts required to match. There is no third lane.

What is given up is narrower than it sounds: a plant model still tells you *whether*
its numbers are self-consistent, which a worked example does, and it tells you *whether*
the physics is right only through the unit operations it is built from.

## Where the port is

**One unit operation in NeqSim is a thin wrapper over a flash call.** That is the whole
argument for the scope below, and it is read from NeqSim's source rather than guessed:

| Unit | Its actual physics | NeqSim lines |
|---|---|---|
| `Splitter` | scale component flows, `TPflash` each branch | 827 |
| `Heater` / `Cooler` | `PHflash(oldH + duty)`, or `TPflash` | 1,113 |
| `Mixer` | sum moles, `PHflash(total enthalpy)` | 1,235 |
| `ThrottlingValve` | isenthalpic `PHflash(enthalpy)` | 1,894 |
| `Separator` | one `TPflash`, then split the phases out | 4,511 |
| `Compressor` | `PSflash(entropy)` → efficiency → `PHflash` | 6,593 |

The line counts are performance charts, entrainment models, geometry sizing, drivers and
mechanical design. NeqSim has 669 files under `process/equipment/`; this port needs three
flash models and about ten thin wrappers. The ranking below is NeqSim's own test suite
counting `new X(` — Separator 303, Compressor 202, ThrottlingValve 172, Heater 120,
Cooler 110, Mixer 82, Splitter 74 — not a judgement about which are interesting.

## The programme

Status is stated as of the last commit that touched this page, and a tranche is not
"done" until both implementations agree on it.

| | Tranche | Status |
|---|---|---|
| 0 | Repair the claims this plan falsifies; a roadmap page; an announcement test that covers models | **in progress** |
| 1 | `eos.ph_flash` and `eos.ps_flash` — the critical path, since eight of the ten units above are one of them plus arithmetic | planned |
| 2a | The stream and the flowsheet: a `stream` quantity type, `specs/flowsheets/`, the sequential solver, proved on `process.mixer` and `process.separator` | planned |
| 2b | The rest of the Pareto set: valve, heater, cooler, splitter, compressor, pump, expander, heat exchanger, three-phase separator; plus `eos.viscosity` and `eos.thermal_conductivity` | planned |
| 2c | Recycle convergence (direct substitution and Wegstein) | planned |
| 3 | `azoth report` — Markdown and HTML from Rust, no new dependency | planned |
| 4 | The agent surface: `azoth describe`, generated from the registry, and an MCP server over it | planned |

**Tranche 1 comes first because it is the dependency, not because it is the most
visible.** There is no unit operation without it. Tranches 3 and 4 are the user features
the port was asked for, and they sit on top of something worth reporting about.

## Deliberately not ported

Recorded so the boundary is visible. Each of these is a real capability in NeqSim and a
separate programme here.

| Out | Why |
|---|---|
| **Distillation columns** (15,553 lines in NeqSim) | A full MESH solve with specification homotopy. If ever wanted, start from a Fenske-Underwood-Gilliland shortcut and a simple Newton MESH — not from a port. |
| **Transient pipeline flow** (`TwoFluidPipe`, 11,369 lines) | Finite-volume two-fluid. Steady state is already covered by `hydraulics.*`. |
| **Reactors, networks, energy and power** (63 + 54 + 46 files) | Reaction and power-system equilibria. Each its own domain. |
| **CPA, SAFT, GERG-2008, electrolytes, Helmholtz reference equations** | Different physics from a cubic equation of state. |
| **Mechanical design, cost, field development, safety and risk** (219 files) | Engineering deliverables, not flowsheet physics. |
| **PVT simulation, hydrates, wax, asphaltene, scale, LNG** | Real capabilities, all outside a cubic-EOS library. |
| **Transient and dynamic simulation** | The steady-state path does not need it — only 2 of NeqSim's 669 equipment files touch its dynamics package, and 58 of the rest implement their transient step as `run()` plus a clock advance. |
| **The paperlab paper-writing pipeline** | A separate heavyweight product for driving NeqSim's own development. |
| **Natural-language flowsheet authoring** (`neqsim-studio`) | An agent feature; revisit once tranche 4 exists. |
