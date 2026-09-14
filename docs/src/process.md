# Process

How work is done in this repository. Normative, like [the specification](./spec.md): where
a habit or a comment disagrees with this page, this page wins.

It is short because most of the process is already in the design — the specs are the
requirements, the two implementations are independent verification, CI is the gate. What
this page adds is that those are **binding**, and that each rule names the mechanism that
enforces it. A rule with no mechanism is a preference.

Read [S9](./spec.md#s9-what-a-contribution-costs) first if you are adding a calculation.
Read this page second.

---

## P1. Requirements traceability

**The spec is the requirement.** Every registered id traces to a specification, to two
implementations and to verification. Nothing is implemented that no spec asks for, and no
spec asks for anything nothing implements.

**Mechanism.** `tools/spec_lint.py` requires every spec to carry a `worked_example` test —
skippable, but only with a reason — and rejects a case whose inputs or expected outputs are
not declared ones. `test_registry_contract.py` and its model counterpart
`python/tests/models/test_model_contract.py` hold spec ↔ registry ↔ code agreement in both
languages, including `test_registry_matches_the_spec_files`, which fails when a registry
entry and a spec file disagree. Both linters and the link checker run in
`.github/workflows/ci.yml`.

## P2. Verification and validation

Two questions, and they are not the same one.

**Verification — did we build it right?** The two implementations agree case by case; the
spec cases pass; the contracts hold; the gates are green. `docs/src/test-plan.md` is the
plan and CI is the mechanism.

**Validation — did we build the right thing?** A worked example a human can retrace, which
`tools/spec_lint.py` enforces, and agreement with NeqSim where the port claims to be a port.

**Verification is done by something other than the author, at a named commit.** "I ran the
tests" is not a record. A commit hash and the gate results at it is.

## P3. Configuration control

**Main is always a baseline.** Every gate green at every commit on `main`.

**A change set is atomic.** Spec, both implementations, tests and docs move together. A
code fix without its spec change, or a spec change without regenerating, is incomplete.

**Generated files are never hand-edited**, and drift fails the build. The `--check` jobs are
the mechanism.

**Nothing is committed red.** Not "commit then fix". A change set that is not a baseline
does not merge and does not get called finished. If it cannot be finished, it is not merged.

## P4. Defect and change management

**A defect has a lifecycle**, and it is in one of these states at all times: raised,
triaged, root-caused, dispositioned, fixed, verified, closed. [Required
improvements](./required-improvements.md) is the log. An entry carries its status, its
**disposition** (fix now, defer, or won't fix, with the reason), and the test that closes it.

**A change to shipped behaviour carries a reason**, and the reason names the defect it
addresses or the requirement it implements.

**A fix is verified by a test that would have caught it.** A fix with no such test is a fix
that will come back.

---

## P5. The coding standard for prose

Prose in this repository is one of three kinds.

**KEEP** — explains the code in front of the reader: what it does, why this expression, what
a constant means, what an argument must satisfy, what an error means.

**MOVE** — reports a defect, a measurement, a correction, a root cause, a missing
capability, or a rationale for a change. It belongs in
[Required improvements](./required-improvements.md).

**DELETE** — narrates the code's own history, or the author's process. What a thing used to
be called, what it used to do, what did not exist before, what the author was trying to
achieve, what the author thinks is important about it.

Three specifics, because they are where this has gone wrong:

- **Attribution for anything ported lives in `NOTICE`, once.** Not per calc, and not in a
  spec field: a block repeated in twenty specs is a block nobody reads. A spec carries a
  citation (a standard, a DOI, a URL) and nothing beyond it. Rust doc comments, Python
  docstrings and test docstrings are not further places to cite upstream line numbers. This
  is the rule in [S6](./spec.md#s6-the-library-implements-the-engineer-decides), applied to
  prose.
- **"NeqSim does X, this does Y" is a MOVE.** A divergence belongs in the defect log. It is
  a limitation report, and a spec is a data sheet — see [Spec files](./spec-files.md).
- **A test docstring says what the test asserts.** Not why it is valuable, not what it would
  catch, not what it catches that another test does not.

**Mechanism.** `tools/prose_lint.py` fails the build on phrases that are never legitimate
here. It is a **backstop**: a phrase list cannot judge, and the standard has to hold while
writing.
