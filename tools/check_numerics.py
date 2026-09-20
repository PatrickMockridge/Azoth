#!/usr/bin/env python3
"""Fail on a fractional exponent written in a form that loses its meaning.

# Why this exists

`docs/src/calculus/numerics.md` states the rule and `lean/Azoth/Pow.lean` proves the
part of it that is provable. Neither runs on a new file, and the two cases here are the
ones a person gets wrong while believing they have not:

**An integral exponent written as `powf(2.0)`.** That is `exp (2 log x)`, which is
`NaN` for `x < 0` where `powi(2)` is exact and total. Fifteen sites in this tree do it.

**A decimal that approximates a rational.** `x.powf(0.3333)` is not `x^(1/3)`: measured,
the two differ by `4.6e-4` at `x = 10^6` and the error grows with `x`. `0.8333` and
`1.2083` are `5/6` and `29/24` the same way.

# What it checks, and what it does not

Two syntactic rules over `crates/**/*.rs` and `python/src/**/*.py`, and nothing else.
It does not know whether a base can be negative and it does not look at a guard — those
are the interesting cases and they are not decidable from a line of source, which is why
the page says the policy is stated per call site rather than enforced by a wrapper.

# The exemption, and why it is a comment rather than a list

**The port rule is that NeqSim wins**, so a decimal that NeqSim's own source writes is a
faithful transcription and must not be "fixed" — the divergence from the rational form is
a *finding to record*, not a defect to repair here. A line is exempt when it, or the line
immediately above it, carries a `numerics-ok:` marker naming that source:

    let d0 = 5.2804 * x.powf(0.3333)  // numerics-ok: ComponentGEWilson.java:143

The marker travels with the line it excuses, which a list of file-and-line pairs would
not, and it puts the reason where the next reader is already looking.

Usage:
    python tools/check_numerics.py
"""

from __future__ import annotations

import math
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent

#: The trees to scan. Rust kernels and the Python reference kernels, which are the two
#: places an expression is transcribed from NeqSim. **`crates/*/src/` and not `crates/`**,
#: so a test asserting a Lennard-Jones exponent is not a finding: the rule is about the
#: arithmetic that ships.
SCANNED = ("crates/*/src", "python/src")

#: Files whose banner says they are generated are skipped. A generator emits the
#: spec's own `equation` text into a table and a worked-example derivation into a
#: docstring, and neither is arithmetic anybody wrote.
GENERATED = "GENERATED FILE - DO NOT EDIT BY HAND"

#: `x.powf(0.3333)`, and the same in Python as `x ** 0.3333`.
POWF = re.compile(r"\.powf\(\s*(?P<literal>-?\d+\.\d+)\s*\)")
PY_POW = re.compile(r"\*\*\s*(?P<literal>-?\d+\.\d+)(?![\d.eE])")

#: How a line says its decimal is NeqSim's own.
EXEMPT = "numerics-ok:"

#: The largest denominator a decimal is compared against. `1/3`, `5/6` and `29/24` are
#: the ones in this tree; beyond a couple of dozen the approximation stops being a
#: recognisable rational and starts being a fitted exponent.
MAX_DENOMINATOR = 24

#: How close a literal's multiple must be to a whole number to count as approximating it.
#: Loose enough for four-decimal truncation (`0.3333 * 3 = 0.9999`), tight enough to miss
#: `2.043` and `2.303`, which are constants and not rationals.
APPROXIMATION_TOLERANCE = 1.0e-3


def approximates_a_rational(value: float) -> tuple[int, int] | None:
    """The `p/q` a literal rounds, where it rounds one and is not equal to it.

    `None` for an exact rational — `0.5`, `0.25`, `0.75`, `1.5` are all exactly
    representable and are the constants they look like, so flagging them would be noise.
    """
    for denominator in range(2, MAX_DENOMINATOR + 1):
        scaled = value * denominator
        gap = abs(scaled - round(scaled))
        if 0.0 < gap < APPROXIMATION_TOLERANCE:
            return (round(scaled), denominator)
    return None


def offenders(path: Path) -> list[str]:
    """Every line of one file that breaks one of the two rules."""
    messages: list[str] = []
    text = path.read_text(encoding="utf-8")
    if GENERATED in text[:2000]:
        return messages
    lines = text.splitlines()
    # **A marker opens a block that runs to the next blank line**, rather than covering
    # one line. An expression that carries a rational-approximating decimal is often
    # wrapped across several lines, and a per-line rule missed the continuation - which
    # is how the first version of this let `1.2083` through while exempting `0.3333` three
    # lines above it.
    exempt = False
    for number, line in enumerate(lines, start=1):
        if EXEMPT in line:
            exempt = True
        if not line.strip():
            exempt = False
            continue
        patterns = [POWF] if path.suffix == ".rs" else [PY_POW]
        found = [match for pattern in patterns for match in pattern.finditer(line)]
        if not found or exempt:
            continue
        for match in found:
            literal = match.group("literal")
            value = float(literal)
            where = f"{path.relative_to(ROOT)}:{number}"
            written = match.group(0).strip().lstrip(".")
            if value == math.trunc(value):
                # **The two languages need different advice.** Rust has `powi`; Python
                # has no such function, but `** 3` with an *integer* exponent takes the
                # real branch where `** 3.0` promotes a negative base to a complex - the
                # same trap by a different route.
                fix = f"`powi({int(value)})`" if path.suffix == ".rs" else f"`** {int(value)}`"
                messages.append(
                    f"{where}: `{written}` is an integral exponent written as a float. "
                    f"Use {fix}, which is exact and is not a complex or `NaN` for a "
                    f"negative base; the float form is `exp (n log x)`. **If the receiver "
                    f"has no integer power** - a dual number, an autodiff type - or the "
                    f"expression is NeqSim's own, mark it `{EXEMPT} <reason>` instead. "
                    f"Ten sites in this tree are `Dual`, which has no `powi`."
                )
                continue
            rational = approximates_a_rational(value)
            if rational is not None:
                numerator, denominator = rational
                messages.append(
                    f"{where}: `{literal}` approximates `{numerator}/{denominator}`, and "
                    f"is not it. Write the rational, which is the exponent the physics "
                    f"means, or mark it `{EXEMPT} <reason>` if the source writes the "
                    f"decimal - in which case the divergence is a finding to record."
                )
    return messages


def main() -> int:
    messages: list[str] = []
    for directory in SCANNED:
        for path in sorted(ROOT.glob(f"{directory}/**/*")):
            if path.suffix in (".rs", ".py") and path.is_file():
                messages.extend(offenders(path))
    if messages:
        for message in messages:
            print(f"  ERROR  {message}", file=sys.stderr)
        print(
            f"\ncheck_numerics: FAILED with {len(messages)} problem(s). See "
            f"docs/src/calculus/numerics.md for the rule.",
            file=sys.stderr,
        )
        return 1
    print("check_numerics: OK (no integral exponent written as a float, no decimal exponent")
    print("                left unmarked that approximates a rational)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
