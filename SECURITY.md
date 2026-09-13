# Security

## What counts as a vulnerability here

This is a calculation library, so the usual list is only half the story. The
worst thing this code can do is not crash, corrupt memory or leak data. It is
**return a wrong number that looks reasonable**.

That means the following are security issues, not merely bugs:

**A result that is wrong without saying so.** The equation is implemented
incorrectly, or a units conversion is applied twice, and the result carries no
warning. This is the highest-severity class: a plausible wrong number goes into
a design.

**A range check that does not run and does not say so.** The library's central
promise is that "checked and fine" and "never checked" are distinguishable.
Anything that breaks that - a check silently skipped, a warning dropped, an
out-of-range value returned clean - is a vulnerability, because it removes the
signal a user relies on.

**Provenance that can be forged or is silently wrong.** A `provenance.json` that
verifies against a tree it does not describe, a signature check that passes when
it should not, a `verify_status` promoted without a verifier. See
[TRUST.md](TRUST.md) for what the trust layer claims.

**A source that is wrong or misattributed.** A spec citing an equation number
that does not exist, or a coefficient recorded against a standard that does not
contain it. This is the one class where the code can be perfect and the result
still unsafe.

Also in scope: memory unsafety, panics reachable from library code, and anything
that makes a badly-formed input produce an undefined result rather than an error.

## What does not count

- **Being outside the validated range.** That is what the range is for, and the
  result says so. If the warning is missing, that *is* in scope.
- **Known placeholder data.** The fitting coefficients are labelled
  `ESTIMATED_DUMMY` and every result using them carries `ESTIMATED_DATA`. They
  are not a vulnerability; they are a documented limitation. Using them for
  design work is a misuse of the library, not a flaw in it.
- **The Rust core being slower than Python for a single call.** Stated in the
  README; not a bug.
- **A dependency's own vulnerability** unless this project's use of it is what
  makes it reachable.

## Reporting

**Do not open a public issue for a security matter.** Use GitHub's private
vulnerability reporting on the repository, or contact the maintainer directly.

Include, as far as you can:

- The calc id and version, or the commit.
- The inputs.
- What you expected, and what you got, with any warnings.
- Where the expectation comes from - a published source, a hand calculation, or
  a cross-check against another tool. For a wrong-number report this is the most
  useful part, and a report without it takes much longer to act on.

## What to expect

This is a small project and there is no security team and no guaranteed response
time. Being honest about that is better than publishing a target that is not met.

- Acknowledgement, best effort, within a week.
- For a wrong-number report, an assessment of whether the spec or the
  implementation is at fault - they are separable here and it changes the fix.
- A fix, a failing test that would have caught it, and credit in the commit
  unless you prefer otherwise.

If a fix requires invalidating published results, that will be said plainly. A
calculation library that quietly changes a number is worse than one that admits
it was wrong.

## Supported versions

Only the latest release. This project is at an early stage and does not maintain
back-ported fixes.
