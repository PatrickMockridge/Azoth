# azoth and NeqSim

<!-- Hand-written, like copyright.md, and listed in tools/gen_docs.py's STATIC_PAGES
     for the same reason: the summary is generated, and an unregistered page is
     silently dropped by mdBook. -->

[NeqSim](https://github.com/equinor/neqsim) is a Java library that solves the same
problem azoth does: cubic equations of state, phase equilibrium, thermodynamic
properties. It was developed at NTNU and is maintained by Equinor, it is Apache-2.0,
and it is incomparably larger.

This page exists because a side-by-side comparison is the only way to answer a question
this project would otherwise never have to face: **what is azoth for?** Reading another
library's choices makes your own visible, including the ones that are limitations.

*Snapshot: NeqSim 3.20.0, cloned 2026-09-13. NeqSim is under active development and the
counts below will age. The shape is the part worth reading.*

## The scale, first

| | NeqSim | azoth |
|---|---|---|
| Language | Java 8 | Rust + Python |
| Lines in `src/main` | 1,285,392 across 3,371 files | — |
| Thermodynamics alone | 199,598 lines | — |
| Test files | 2,123 | 102 — 57 Python, 45 Rust |
| Registered calculations | 60+ equations of state, 33+ equipment types | 21 calcs, 17 models |
| Component data | 258 rows in `COMP.csv`, 76,705 in `COMP_EXT.csv`, 1,309 kij rows | 173 in `data/components/`, 516 kij rows, vendored from NeqSim |
| Licence | Apache-2.0 | AGPL-3.0 code, CC-BY-4.0 docs and data |

The comparison that matters is not the arithmetic. It is that NeqSim is a **process
simulator** and azoth is a **calculation library**, and every difference below follows
from that.

## Seven axes where the two differ by design

### 1. Where the truth lives

NeqSim's thermodynamics is the Java source, with prose documentation written alongside
it. azoth's thermodynamics is a YAML spec, from which the documentation, the registries,
the range checks and the test cases are all generated — and CI regenerates them and
fails on any difference.

This has a visible consequence in NeqSim's own repository: `docs/BROKEN_API_AUDIT_REPORT.md`
is a hand-run scan finding documentation that references classes which do not exist —
`SimpleWell`, `ChokeValve` — with a table of file and line number for each. It is a
careful and honest document, produced because the same defect kept recurring.

**That defect class cannot occur here.** A page under `docs/src/eos/` is rendered from the
spec that both implementations read. There is no second description of the equation to
drift from the first.

### 2. What data ships

NeqSim ships a component databank: 258 components in `COMP.csv`, 76,705 in
`COMP_EXT.csv`, and 1,309 binary interaction rows, each component carrying around 180
columns — critical properties, ideal-gas heat capacity coefficients, Antoine constants,
association parameters, SAFT parameters, hydrate coefficients, viscosity correlations.

**azoth now ships a databank too, and it is NeqSim's.** When this page was first
written the row above read "none", and the comparison made the point that an azoth user
had to supply every critical constant by hand while a NeqSim user wrote
`addComponent("methane", 1.0)` and got an answer. That was the better half of the
usability trade and it was conceded as such. It is no longer true: `data/components/` is
generated from this repository's `COMP.csv` and `INTER.csv` by
[`tools/gen_databank.py`](https://github.com/PatrickMockridge/Azoth/blob/main/tools/gen_databank.py),
and `azoth.eos.from_names(["methane", "n-butane"])` is the NeqSim call in azoth's
spelling.

Three things about the vendoring are worth stating, because they are where it differs
from what NeqSim ships.

**Not every row came across.** 62 of the 258 components are ions, and a cubic equation
of state has no notion of one — NeqSim fills their critical properties with a shared
default, which shows up as 29 rows carrying the same `Pc`, `omega` and `Vc`. That is
not a coincidence and it is not data, so those rows are excluded, along with `ice`,
`salt`, `seawater`, `asphaltene` and the rows with no type at all. 173 substances came
across. `COMP_EXT.csv` did not: 86 MB of heavy fluids this library cannot characterise.

**The provenance is institutional, and that is now the whole of it.** `COMP.csv` has no
citation column, so a value in it carries a project, a version and a file rather than a
per-value citation. That is real — more than most engineering data has — and it is not a
person having checked anything. The vendored component rows therefore carry that
citation verbatim, naming Equinor and NTNU and the version it came from, rather than an
assertion about how far it can be trusted. The interaction-parameter table is the
exception: its rows carry no citation column at all, so its attribution is the `NOTICE`
file rather than a per-row one.

**The residual risk is inherited, not resolved.** Apache-2.0 permits redistributing
NeqSim's compilation; it does not establish that every value inside was cleanly sourced
upstream. That is not inspectable from here and the mitigation is attribution rather than
inspection. See [Copyright and licensed data](../copyright.md).

### 3. How correctness is claimed

NeqSim has 2,123 test files in the conventional style, plus benchmarks and regression
suites. azoth has two **independent implementations** of every calculation, in Python
and Rust, compared case by case with the **iteration counts required to match exactly** —
the sharpest cheap check that both ran the same procedure — plus a small set of external
validation cases that each carry a `source.verification` status.

The difference is sharpest on the one calculation both libraries have recently built.
NeqSim's `CriticalPointFlash` implements Heidemann & Khalil (1980) correctly as far as
can be told from reading it, and has **no validation of any kind**: no check against a
pure component's known critical point, no check against a mixture critical locus, and a
silent `break` on a `NaN`. Nothing in it is wrong; nothing in it is checked either.

That is where azoth's contribution actually lies. The algorithm is not ours. **What holds
the algorithm is.**

### 4. Units

NeqSim is `double` throughout, with unit conventions that are consistent and documented
by habit rather than by type: `COMP.csv` stores critical temperature in degrees Celsius
and the loader adds `273.15`; pressure is in bara; molar mass is g/mol in the file and
kg/mol in memory; volumes carry an internal factor of `1e5`. None of this is a defect —
it is how a mature library written before dimensional types were common works, and the
conversions are all in one place.

azoth puts `uom` on the Rust boundary and `pint` on the Python one, declares every
input's unit in the spec, and runs a `unit_round_trip` test per calculation. A
Celsius/kelvin confusion is a type error rather than a factor of 273.15.

### 5. Scope

NeqSim is 1.29 million lines: 33+ equipment packages, PVT simulation, pipeline flow,
hydrates, safety and relief, mechanical design, cost estimation, field development
economics — and, more recently, an MCP server and tooling for AI agents.

azoth is 21 calculations and 17 models, eight of which are unit operations. The distance
from NeqSim's 1.29 million lines is
a **scope decision, not a stage of work**: nearly every NeqSim unit operation is a flash
call plus arithmetic, and what surrounds it is performance charts, entrainment models,
geometry sizing and mechanical design. Azoth takes what makes a flowsheet run — which
now includes a process layer, where a calculation becomes a transformation of streams
and the composition is done for you rather than by the caller. What makes a model
acceptable to a detailed-design review is a different product, and
[Roadmap](../roadmap.md) records where the line is drawn.

### 6. Licence, and who each library is for

NeqSim is Apache-2.0 — permissive, embeddable in a proprietary product, no obligation
beyond attribution.

azoth is AGPL-3.0 for code, with documentation and data under CC-BY-4.0. That split is
deliberate and is argued on the copyright page: *"the equations and the coefficients are
the part most people want to reuse or cite, and they should not require adopting a
copyleft obligation to do it."* The two libraries are aimed at different things. NeqSim
can go inside a closed product; azoth's code cannot, and azoth's *numbers* can be cited
by anyone.

### 7. The same values, mechanised differently

This is the most interesting finding, and it is not a difference of principle.

NeqSim cares about exactly what azoth cares about. Its Pitzer documentation says *"No
missing interaction is silently converted to zero"* and describes failing closed on
incomplete parameter coverage. Its reaction-model audit says *"a shared reaction name or
stoichiometry is not evidence that its equilibrium constants, activity convention, or
validity range transfer between thermodynamic models."* Those are the same convictions
as azoth's *"a check that could not run is not a check that passed"* and *"a wrong number
that looks reasonable is the failure this project is organised against."*

The difference is **when** and **by what**. NeqSim enforces its values at **runtime**,
with fail-closed guards and diagnostics, and reviews them with **retrospective audit
tools** — a provenance document here, an audit report there. azoth enforces its values
at **build time**, with schema-required spec fields, generated documentation, and CI
drift checks that fail the build.

Both are honest answers to the same problem. Saying so is more accurate than claiming a
virtue the other library lacks.

## What azoth does not do, and why

The boundary, stated so that it is visible rather than discovered:

| Not implemented | Why |
|---|---|
| Heat-capacity coefficients and the SAFT/CPA association parameters | The databank carries critical constants and `kij` only. The rest arrive when models that read them do. |
| CPA, SAFT, GERG-2008, Helmholtz reference equations, electrolytes | Different physics from a cubic EOS. Each is its own programme. |
| Distillation columns, transient pipeline flow, reactors, networks | Separate programmes, each larger than this one. See `docs/src/roadmap.md`. |
| PVT simulation, hydrates, wax, asphaltene, scale | Real capabilities, all outside a cubic-EOS library. |
| Mechanical design, cost, field development, safety and risk | Engineering deliverables, not flowsheet physics. |
| Relief-valve *sizing* to a standard | The de-rating coefficients are the caller's to compose; `choked_flow_area` is the isentropic basis. |

Two rows left this table, and a third item never was one. Each is worth naming, because
each was a boundary this page argued for and each has since been crossed by the databank
work or by a port:

**A component databank** was "the central question of axis 2" and the answer here was
that this library ships none. It ships one now — 173 substances and 516 `kij` pairs,
generated from NeqSim's `COMP.csv` and `INTER.csv` and carrying NeqSim's provenance
rather than a per-value one. The reasoning that changed is on
[Copyright and licensed data](../copyright.md); the residual risk — a permissive licence
on a compilation does not vouch for every value inside it — is named in `NOTICE` rather
than hidden.

**A mixture critical point** was "not shipped yet", with the mechanical route
`dP/dV = d2P/dv2 = 0` rejected because it returns the same `Z_c` for every mixture. It
ships, on Heidemann & Khalil's conditions, and the test that distinguishes the two
routes is that a mixture's `Z_c` moves with composition — it varies by 0.146 across
methane/n-butane, against a constant.

**Equipment models have arrived; flowsheets have not.** This page used to say they were
"a thesis change, not a feature", and that was the right call when it was written: the
project's guarantee is that every calculation has one phase, one composition and a
worked example a human can retrace, and a flowsheet has none of those. The user has
since asked for the wider port — unit operations, reports and an agent surface — and
eight unit operations now ship, each still a model with a spec, a source and a
retraceable example. A flowsheet is designed and not built, so the boundary moves
deliberately rather than by drift; what replaces the guarantee when it does is stated in
`docs/src/roadmap.md`.

## What we take from it

NeqSim is Apache-2.0, which permits reuse with attribution, and this project ports from
it. Two rules govern what that means, and they are in `CONTRIBUTING.md`:

**A port is never evidence.** Reading someone's Java is not reading the paper it came
from, and a port is a second implementation of a method rather than a second source for
it. A port cites the paper for the *method* and the upstream source for the *port*, and
nothing in the tree records a status claiming otherwise, because a status is a claim
about a person having read a source and no tool can check one. Attribution itself lives
in `NOTICE`, once, rather than restated per spec.

**What a port changed is a limitation report**, and it belongs in
[Required improvements](../required-improvements.md) rather than in the spec — a spec is
a data sheet ([Spec files](../spec-files.md)). The one thing the spec keeps is the
*citation*: the upstream class, and the line ranges taken and not taken. That is what a
licence asks for, and it is a citation rather than an argument.

**A port is accepted on this library's tests, never on its provenance.** That a
well-known library implements something is evidence that it can be implemented. It is
not evidence that it is right, and NeqSim's critical point is the concrete case: correct
by inspection, and validated nowhere.

Every ported algorithm is listed with the NeqSim class it came from, so the boundary
between what is ours and what is borrowed is legible in the specs themselves. The
critical point is the only port that also records a commit; the unit operations record a
class and a line range in the 3.20.0 source tree, and the flash models cite the paper the
method comes from rather than a file.

## Credit

NeqSim is developed at NTNU and maintained by Equinor, and is used in real oil and gas,
carbon capture and hydrogen work. Its documentation is unusually thorough about *why* its
models are built the way they are, and reading it is how this project resolved a problem
it had recorded as unsolvable. The attribution obligations owed when code or data is
reused are listed in the repository's `NOTICE` file.
