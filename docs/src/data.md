# What ships

<!-- Hand-written, unlike the pages under the namespaces, which are generated from
     the specs. It is listed in SUMMARY.md by tools/gen_docs.py, and
     tools/check_links.py fails the build if a page under docs/src is missing
     from the summary. -->

The library carries data, and this page says what and where it came from. That is the
whole of what the library owes on the subject — it does not ask you to record the
provenance of your own values, and it does not tell you which of its own to trust. It
tells you what it has and who produced it, and you judge.

## The component databank

`data/components/components.csv` — **173 substances**, and `kij.csv` — **516 binary
interaction parameter pairs**.

Both are generated from [NeqSim](https://github.com/equinor/neqsim)'s `COMP.csv` and
`INTER.csv` by `tools/gen_databank.py`, and every row cites the source it came from:
NeqSim v3.20.0, Equinor and NTNU, Apache-2.0. The full attribution is in
[`NOTICE`](https://github.com/PatrickMockridge/Azoth/blob/main/NOTICE).

**What is not in it, and why.** NeqSim's `COMP.csv` is larger. The ions are excluded
because a cubic equation of state has no notion of one, and NeqSim fills their
critical properties with a shared default — twenty-nine rows carrying the same `Pc`,
`omega` and `Vc` is not a coincidence, and shipping them would be shipping
plausible-looking wrong numbers. `COMP_EXT.csv` is not vendored at all: 86 MB of heavy
fluids this library cannot characterise.

**The provenance is institutional, and that is the whole of it.** `COMP.csv` has no
per-value citation column, so a value in it carries a project, a version and a file
rather than a source for that number. That is more than most engineering data has, and
it is not a person having checked anything. The vendored rows therefore carry that
citation verbatim rather than an assertion about how far it can be trusted.

**The residual risk is inherited, not resolved.** Apache-2.0 permits redistributing
NeqSim's compilation; it does not establish that every value inside was cleanly sourced
upstream. That is not inspectable from here.

To change any of it, supply your own in a [keycard](./keycard.md).

## The fitting coefficients are placeholders

`data/fittings/crane_k_factors.csv` holds **estimated dummy values**. They are
plausible magnitudes chosen so the software has something to run against. They are not
from Crane TP-410 or any other standard, and **a pressure drop computed from them can
be wrong by a factor of two while looking entirely reasonable.**

The `verify_status` column records that, and a test fails the day someone populates the
file properly. Supply your own with a [keycard](./keycard.md).

## Water and air are real

`data/fluids/water.csv` and `air.csv` are a different case: real published values,
consistent with IAPWS-IF97, NIST and the CRC Handbook. They are marked `unverified`,
which is **not** a statement that the numbers are doubtful — it says that nobody here
has checked them against a primary formulation.

Water's viscosity varies by a factor of six across 0–100 °C, and the provider
interpolates linearly between the points given and refuses to extrapolate: past either
end, a straight-line extension is a confident wrong number rather than a small error.

## Why there is a status column here and not in your keycard

A deliberate asymmetry, not a leftover. The data *this repository ships* records a
derived `verify_status` so that the library can warn you when a result rests on its
placeholders — that is disclosure, and it is the library's own statement about itself.

A value *you* supply is yours. There is no status field in a keycard, because a form
asking you to assert something no tool can check would teach you to fill it in rather
than to know the answer. See
[Specification, S6](./spec.md#s6-provenance-is-the-engineers-job-not-the-librarys).

## Reproducing any of it

```bash
python tools/gen_databank.py --check     # the component databank, from NeqSim
python tools/gen_user_data.py keycard.yaml --check   # your own, from a keycard
mdbook build docs && python tools/check_links.py
```

The data files are written in place, because both languages read them at fixed paths —
Rust embeds them with `include_str!` at compile time. So if you generate a keycard's
data over the shipped files, `git checkout -- data/` puts the shipped placeholders
back.
