# Inherent safety

A relief valve is a protection layer. It contains a hazard the process still has, and it
works until the day it is blocked, undersized or isolated. A chemistry that cannot run
away is a different kind of answer, and process safety has preferred it since Kletz put
the slogan to it — *what you don't have, can't leak*.

This library applies that preference to its own hazard. Naming the hazard is the first
thing, because it is not a crash.

## The failure mode

`SECURITY.md` states it and means it: the worst thing this code can do is not crash,
corrupt memory or leak data. It is **return a wrong number that looks reasonable**.

The asymmetry is the whole of the problem. A crash is loud, arrives at the moment of the
mistake and costs a stack trace. A plausible wrong number is silent, costs nothing to
notice, and goes into a design. Every mechanism on this page is aimed at that one failure,
and none of them is aimed at making the numbers right — that is what the source, the valid
range and the oracle are for.

What the library owes is narrower and harder: a number that cannot be trusted must not be
tellable apart from one that can. "Checked and fine" and "never checked" are different
states and stay different; a value nothing read is not a value in use; an answer that is
two answers is not an answer.

## What is enforced, and by what

| The property | The mechanism | Where |
|---|---|---|
| A failure is loud, and names the input | typed errors, one variant per class, each carrying the offending field | `crates/azoth-core/src/error.rs` |
| An answer that is two answers is refused | a single-phase root, or an error naming the two-phase state | `crates/azoth-process/src/stream.rs` |
| A partial function is total, clamped, or deliberately `NaN` | the rule, and a gate that reads the source for it | [`numerics.md`](./calculus/numerics.md), `tools/check_numerics.py` |
| A value nothing reads is refused, not skipped | the card reader's load-time refusal | `crates/azoth-eos/src/card.rs` |
| A name the data cannot answer is an error, not a default | one choke point for every substance lookup | `crates/azoth-eos/src/databank.rs` |
| The Python path cannot quietly stand in for the Rust one | `AZOTH_REQUIRE_RUST=1` makes a missing core fatal | `python/src/azoth/_dispatch.py` |
| A proof with a gap fails the build | one `#print axioms` line per claimed theorem, against an allow-list of three | `lean/Azoth/Axioms.lean`, `tools/check_lean_axioms.py` |

## What the checks caught

azoth is a port of NeqSim, and porting against an oracle is largely the discovery that the
oracle is wrong in ways nothing inside it can see. These are the defects the port's checks
caught, and each is filed upstream rather than written up here — the issue is the record.

| The property it violates | What it was | Filed |
|---|---|---|
| A value nothing reads is data that looks in use and is not | `PhysicalPropertyMixingRule` never reaches its database lookup, so `Gij` is zero for every fluid | [#3989](https://github.com/equinor/neqsim/issues/3989) |
| | `ComponentSrkCPAMM` reads the plain `PARACHOR` column, because the dispatch tests exact class names | [#3990](https://github.com/equinor/neqsim/issues/3990) |
| A zero that nothing wrote is not a value | `PhaseGEWilson` and `PhaseGEUniquac` never publish an activity coefficient, so their fugacity coefficients are `0` or `NaN` | [#3770](https://github.com/equinor/neqsim/issues/3770) |
| A term absent from both sides of a comparison is invisible to it | `Component.getHID` multiplies the ideal-gas enthalpy of formation by zero | [#3991](https://github.com/equinor/neqsim/issues/3991) |
| A derivative is the derivative | `ComponentSrkCPA.dFCPAdNdN` is not the derivative of `dFCPAdN`; `calc_lngij` carries an extra term equal to `calc_lngi` | [#3992](https://github.com/equinor/neqsim/issues/3992) |
| A check that did not run is not a check that passed | `performGibbsMinimization` performs no minimization, and logs that it completed | [#3993](https://github.com/equinor/neqsim/issues/3993) |
| A default where a measurement was needed is the defect | a zero-filled vapour-pressure row is silently given a curve from its normal boiling point, in the wrong unit | [#3994](https://github.com/equinor/neqsim/issues/3994) |
| | three Antoine defaults cover 63% of `COMP.csv`, and some rows put a vapour pressure on ions | [#3771](https://github.com/equinor/neqsim/issues/3771) |
| A declaration that disagrees with the data under it does not balance | `elementNames` has seven entries for an eight-column table, so the `Ar` balance reads the file's `Na` column, the `Z` balance reads `Ar`, and the charge column is read by nothing; `inlet_mole` is filled per phase, so a two-phase feed solves on phase 0 alone | [#3986](https://github.com/equinor/neqsim/issues/3986) |
| An unreachable branch is a branch nobody has run | the DIPPR-101 branch of `getAntoineVaporPressure` is unreachable, and 20 components get 1e38–1e96 bar | [#3768](https://github.com/equinor/neqsim/issues/3768) |
| A desync is a throw waiting for a reason | `SystemUNIFAC`'s group array and group list are never synchronised | [#3769](https://github.com/equinor/neqsim/issues/3769) |
| A no-op that reports success is the hardest of these to see | `ComponentEos.setAttractiveTerm` logs an error, then logs success, and keeps the previous alpha function | [#3995](https://github.com/equinor/neqsim/issues/3995) |

None of these is a claim about the people who wrote them. The left column is the point:
each defect is a *shape*, and the shape is what the checks are for. A table like this is
evidence that a reader can re-derive — every row names a file, a function and an
observation, and the issue carries the measurement.

## The one that was azoth's own

**A divergence is a finding, never a failure.** Read the other way round — a divergence is
a finding, *so a failure is not one* — the sentence cost a session. A finite difference
that subtracted the wrong quantity reported a NeqSim derivative as wrong; it was
[#3805](https://github.com/equinor/neqsim/issues/3805), the maintainer checked it and
published the comparison, the two agreed to eight figures, and the report was withdrawn on
the issue. That is the whole episode, in public, and it is the rule the exercise runs on: a
divergence starts from [the assumption](./architecture/specification.md) that azoth is the
one that is wrong.

That rule is not a courtesy to the upstream project. It is what makes the table above worth
reading, because a finding that has not survived being wrong about which side is at fault
is not a finding.

The corollary is that a filed finding can be withdrawn too, and one was. #3992 was filed,
and the execution pass run against the pinned build afterwards could not reproduce it in
situ: the flashed liquid state is not the state the filed number came from, and the obvious
finite-difference check returns exactly `0.0` for every derivative — including a volume
derivative whose analytic value is known to be non-zero — because `dFCPAdN` reads the
component's stored site fractions rather than the state it is handed. The control against
the analytic value is what made that visible instead of looking like a confirmation. The
correction is on the issue, which is where it belongs.

## The formal layer, and where it stops

The premise is that a quantity's dimension and the arithmetic over it should be a proved
object rather than a convention everybody remembers. `mm` and `m` differ by a factor; `Pa`
and `Pa*s` are different dimensions; `W` and `J/s` are the same one. Each of those is a
statement about a group, and if the group is proved then a quantity that does not compose
*fails* rather than converting — no test has to be written to catch the case, because there
is no case to catch.

That premise is what the Lean layer is built on, and the honest statement of how far it
currently reaches is the repository's own three-way vocabulary, set out in
[the calculus](./calculus/index.md): **proved**, **specified** and **characterised**. Proved
means a theorem in `lean/Azoth/` that the axiom gate covers. Specified means the layer does
not exist yet, so the claim is what will make the tranche that builds it checkable.
Characterised means the layer exists and the claim describes it rather than guarding it —
nothing this repository could do would violate it, so it is stated and deliberately not
proved.

What is proved today is the dimensional calculus and two load-bearing derivative formulas:
the dimension group (`lean/Azoth/Dim.lean`), the unit vocabulary
(`lean/Azoth/Vocabulary.lean`), the fractional-power conventions and their disagreement
(`lean/Azoth/Pow.lean`), the raw-to-normalised conversion (`Azoth.Normalisation.lean`,
`hasFDerivAt_extensive`), the first-order sensitivity of a solution
(`Azoth/Implicit.lean`, `hasFDerivAt_of_solution`), and the keycard as a capability that
cannot amplify authority (`Azoth/Capability.lean`, `run_derives_in_grant`).

What is deliberately not proved is named rather than glossed. That the balances close is
**specified**: `specs/unit_ops/` declares the channels, `azoth_process::validate` holds a
flowsheet to them in Rust, and the balance is therefore *checked* — what is absent is the
theorem, because a balance is about values crossing channels and a barb records a channel
rather than a magnitude. That a recycle has a fixed point is **characterised**: the layer
exists and the claim describes it, because the conclusion is barbed congruence while the
five constructors of `Process` carry no value — so a loop and a process emitting its fixed
point barb the same channel whatever the loop converges to, and uniqueness, the half that
matters, is invisible to the layer. Proving a claim weakened into something the layer can
carry would be a different claim, and dressing one as a theorem would be a vacuity.
[The calculus](./calculus/index.md) has the argument.

**And the gate is a completeness gate, not a correctness one.** It proves a proof has no
gap; it does not prove the theorem is the one a reader expects. `lean/Azoth/Axioms.lean`
says so about itself, which is the only way that sentence can be trusted.

## Two implementations, mirrored

Every calculation is written twice by hand, one Rust kernel and one Python kernel, from the
same spec — and the two are **mirrors, not a core and a front end**. Neither is a wrapper
over the other, so the two can disagree with each other, and a disagreement is a finding.
See [the architecture](./architecture/index.md).

It is also why a comparison that cannot run is worse than no comparison, because it reports
agreement. If the extension is missing, every cross-language test passes on the Python path
while the guarantee the project promises is exercised by nothing: the suite stays green and
the guarantee evaporates. `AZOTH_REQUIRE_RUST=1` closes that, turning "the extension is
missing" from a quiet fallback into a failure, and CI runs the cross-implementation job with
it set. Were one implementation a binding over the other, a silent fallback would cost
nothing; because the two are mirrors, it would hide the whole comparison.

The same asymmetry shows up in the keycard. A card is read by
`python/src/azoth/keycard.py` and by `azoth_eos::card` independently, so a Rust-native
caller holds a card of their own rather than being handed the values Python resolved —
and one card's text goes to both readers, with the two resolutions compared name by name
and pair by pair. See [the keycard as a capability](./calculus/capability.md).

## Standalone

**The aim is a standalone simulator, and the aim is not the status.** Being honest about
which half is which is the same discipline as everywhere else on this page, so:

What is already true. The Rust implementation has no Python dependency — it compiles to a
library and a CLI. The keycard is a value a caller passes rather than a global in force, so a
Rust-native caller supplies their own authority rather than inheriting a Python session's.
The unit-operation and executor layers compose, and the Rust implementation is built for
large or many concurrent simulations.

What is not. Port coverage is incomplete. Everything ported ships `unverified`, and a port
is accepted on azoth's own tests rather than on its provenance — that NeqSim implements
something is evidence it can be implemented, not evidence it is right. The fitting
coefficients this repository ships are `estimated_dummy` placeholders, which is why
[the front door](./index.md#not-for-design-work-yet) says not for design work yet. A
standalone simulator is a longer-term aim, and `ROADMAP.md` is where the order is
normative.

## What this does not claim

- **The physics is not proved.** The Lean layer covers the dimensional calculus and two
  derivative formulas. No thermodynamic kernel is a theorem, and none is claimed to be.
- **Everything ported ships `unverified`.** Reading someone's Java is not reading the
  paper, so a port never upgrades a verification status.
- **The oracle does not gate the build.** It is a second opinion, not a source of truth,
  and a divergence from it is a finding to explain.
- **Some shipped data is placeholder.** A pressure drop computed from the Crane fitting
  coefficients can be wrong by a factor of two while looking entirely reasonable, and the
  library says so rather than rounding the difference away.
