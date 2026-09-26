#!/usr/bin/env python3
"""Fail on the ways an absence gets written as a number.

# Why this exists

`docs/src/calculus/numerics.md` states the rule and `lean/Azoth/Pow.lean` proves the
part of it that is provable. Neither runs on a new file, and the cases here are the ones
a person gets wrong while believing they have not.

**An integral exponent written as `powf(2.0)`.** That is `exp (2 log x)`, which is
`NaN` for `x < 0` where `powi(2)` is exact and total. Fifteen sites in this tree do it.

**A decimal that approximates a rational.** `x.powf(0.3333)` is not `x^(1/3)`: measured,
the two differ by `4.6e-4` at `x = 10^6` and the error grows with `x`. `0.8333` and
`1.2083` are `5/6` and `29/24` the same way.

**A physical field defaulted by `unwrap_or_default()`.** `over.tc.unwrap_or_default()`
gives a missing critical temperature the value zero, which is a number every model
computes with - so a card that omits one is a state the arithmetic runs on rather than
refuses.

**A `f64::MAX` outside a `let mut` seed.** The largest finite float is three things a
reader cannot tell apart from the line: a bound somebody meant, a ceiling that saturates
an infinity, and a `min` fold's seed.

**An index lookup defaulted to zero.** `names.iter().position(|n| n == want).unwrap_or(0)`
answers zero for a name the list does not carry, and zero is a *different* entry - so a
caller asking for a field the table does not declare reads the first one. This is the one
rule that reads a whole statement rather than a line, because the lookup and its default
are usually two lines apart.

# What it checks, and what it does not

Syntactic rules over `crates/**/*.rs` and `python/src/**/*.py`, and nothing else. A rule
does not know whether a base can be negative and it does not look at a guard — those are
the interesting cases and they are not decidable from a line of source, which is why the
page says the policy is stated per call site rather than enforced by a wrapper.

**A comment is not arithmetic**, so a line whose first non-space is `//` or `#` is not
matched. The pages and the doc comments name these forms on purpose — this file's own
docstring does — and flagging prose would make the rule about the documentation.

# The exemption, and why it is a comment rather than a list

**The port rule is that NeqSim wins**, so a decimal that NeqSim's own source writes is a
faithful transcription and must not be "fixed" — the divergence from the rational form is
a *finding to record*, not a defect to repair here. A line is exempt when it, or the line
immediately above it, carries a `numerics-ok:` marker naming that source:

    let d0 = 5.2804 * x.powf(0.3333)  // numerics-ok: ComponentGEWilson.java:143

The marker travels with the line it excuses, which a list of file-and-line pairs would
not, and it puts the reason where the next reader is already looking. **The count of
excused lines is printed per rule**, so an exemption that widens is visible in the output
rather than being a silence.

Usage:
    python tools/check_numerics.py
"""

from __future__ import annotations

import json
import math
import re
import sys
from collections.abc import Callable
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

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

#: `f64::MAX`, and the one shape where it is unambiguously a fold's seed.
F64_MAX = re.compile(r"\bf64::MAX\b")
SEED = re.compile(r"^\s*let mut \w+ = f64::MAX;\s*$")

#: An index lookup that defaults to zero: an index *value* defaulted, which the statement
#: has to also contain an index idiom for. `map_or(0.0, |index| ...)` is deliberately not
#: matched: that computes a value the absent case contributes, and this is about a number
#: that is later used to address an entry.
ZERO_INDEX_DEFAULT = re.compile(r"\.unwrap_or\(0\)")
INDEX_IDIOM = re.compile(r"\.(?:position|rposition|binary_search)\(|\bindex_of\(")

#: How a line says its number is deliberate.
EXEMPT = "numerics-ok:"

#: A comment line, which is prose about arithmetic rather than arithmetic.
COMMENT = re.compile(r"^\s*(//|#)")

#: The component parameters a card may state, from the schema rather than from a list
#: here. `specs/schema/component.schema.json` is what a card's component block is checked
#: against, so a parameter added there is one this rule starts looking for in the same
#: edit - which a hand-written list would not do.
COMPONENT_SCHEMA = ROOT / "specs" / "schema" / "component.schema.json"

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


def physical_fields() -> tuple[str, ...]:
    """The component parameter names the schema declares, sorted.

    Read rather than written here, so the list cannot rot. A missing schema answers the
    empty tuple rather than raising: the synthetic trees the tests compile have no
    `specs/`, and a rule that cannot be built is one that silently stops running.
    """
    if not COMPONENT_SCHEMA.exists():
        return ()
    schema: dict[str, Any] = json.loads(COMPONENT_SCHEMA.read_text(encoding="utf-8"))
    return tuple(sorted(schema.get("properties", {})))


@dataclass(frozen=True)
class Rule:
    """One syntactic rule: a name, the patterns that recognise it, and what a match means.

    `explain` returns the message for one match, or `None` where the pattern matched
    something the rule does not mean — which is how the exponent rules separate an exact
    rational from a truncated one without a second regex. `suffix` is the file's
    extension, because the advice differs between the two languages.
    """

    name: str
    #: Keyed by the file's suffix, because the two languages spell the same defect
    #: differently and a `.py` file's docstring quotes the Rust it mirrors - so a pattern
    #: applied to both would report prose.
    patterns: dict[str, re.Pattern[str]]
    explain: Callable[[re.Match[str], str, str, str], str | None]

    #: Whether the defect spans a statement rather than fitting on a line. A chain that
    #: names a field on one line and defaults it on the next is one expression, so a
    #: per-line rule would miss half the sites.
    multiline: bool = False

    def messages(self, path: Path, line: str, number: int) -> list[str]:
        """Every finding this rule makes on one line of one file."""
        pattern = self.patterns.get(path.suffix)
        if pattern is None:
            return []
        where = f"{path.relative_to(ROOT)}:{number}"
        return [
            message
            for match in pattern.finditer(line)
            if (message := self.explain(match, line, where, path.suffix)) is not None
        ]


def exponent_is_integral(match: re.Match[str], line: str, where: str, suffix: str) -> str | None:
    """`powf(2.0)`, which is `exp (2 log x)`."""
    _ = line
    literal = match.group("literal")
    value = float(literal)
    if value != math.trunc(value):
        return None
    # **The two languages need different advice.** Rust has `powi`; Python has no such
    # function, but `** 3` with an *integer* exponent takes the real branch where `** 3.0`
    # promotes a negative base to a complex - the same trap by a different route.
    fix = f"`powi({int(value)})`" if suffix == ".rs" else f"`** {int(value)}`"
    return (
        f"{where}: `{match.group(0).strip().lstrip('.')}` is an integral exponent written "
        f"as a float. Use {fix}, which is exact and is not a complex or `NaN` for a "
        f"negative base; the float form is `exp (n log x)`. **If the receiver has no "
        f"integer power** - a dual number, an autodiff type - or the expression is NeqSim's "
        f"own, mark it `{EXEMPT} <reason>` instead. Ten sites in this tree are `Dual`, "
        f"which has no `powi`."
    )


def exponent_approximates_a_rational(
    match: re.Match[str], line: str, where: str, suffix: str
) -> str | None:
    """`powf(0.3333)`, which is not `x^(1/3)`."""
    _ = (line, suffix)
    literal = match.group("literal")
    rational = approximates_a_rational(float(literal))
    if rational is None:
        return None
    numerator, denominator = rational
    return (
        f"{where}: `{literal}` approximates `{numerator}/{denominator}`, and is not it. "
        f"Write the rational, which is the exponent the physics means, or mark it "
        f"`{EXEMPT} <reason>` if the source writes the decimal - in which case the "
        f"divergence is a finding to record."
    )


def physical_field_defaulted(
    match: re.Match[str], line: str, where: str, suffix: str
) -> str | None:
    """`over.tc.unwrap_or_default()`, which is a missing critical temperature as zero."""
    _ = (line, suffix)
    parameter = match.group("field")
    return (
        f"{where}: `{parameter}.unwrap_or_default()` gives an absent {parameter} the value "
        f"zero, and zero is a number every model computes with rather than an absence. "
        f"Refuse it, or mark it `{EXEMPT} <reason>` where zero is the deliberate statement - "
        f"an ion's missing critical constants are the case in this tree."
    )


def maximum_float_outside_a_seed(
    match: re.Match[str], line: str, where: str, suffix: str
) -> str | None:
    """`f64::MAX` where it is not a `min` fold's seed."""
    _ = (match, suffix)
    if SEED.match(line):
        return None
    return (
        f"{where}: `f64::MAX` is the largest finite float, and a reader cannot tell a bound "
        f'somebody meant from "no bound" - nor a *ceiling*, which it also is: an infinite '
        f"norm against it scales to zero. Say which it is: a named constant where the class "
        f"writes it, a refusal where the region is unreachable, or mark it `{EXEMPT} <reason>` "
        f"where it seeds a paired `min`/`max`."
    )


def index_defaulted_to_zero(match: re.Match[str], line: str, where: str, suffix: str) -> str | None:
    """A lookup whose index defaults to zero, and the statement it sits in."""
    _ = (suffix,)
    if not INDEX_IDIOM.search(line):
        return None
    return (
        f"{where}: `{match.group(0)}` reads a name the lookup did not find as **index zero**. "
        f"Zero is a different component, field or phase rather than an absence, and nothing "
        f"about the number says which was meant. Refuse it, or mark it `{EXEMPT} <reason>` "
        f"where zero is the deliberate answer - a substance a phase does not carry "
        f"contributing nothing is the case in this tree."
    )


def build_rules() -> tuple[Rule, ...]:
    """The rules, in the order their findings are reported."""
    fields = "|".join(re.escape(name) for name in physical_fields())
    # A pattern that matches nothing where the schema is missing, rather than a rule that
    # raises: the tests compile synthetic trees with no `specs/`.
    # Case-insensitive: the schema names the parameters as a *card* writes them (`Tc`) and
    # the Rust field is the same name lowercased (`tc`), so a pattern that matched one
    # spelling would miss half the sites.
    defaults = re.compile(
        rf"\b(?P<field>{fields})\.unwrap_or_default\(\)" if fields else r"(?!)", re.IGNORECASE
    )
    return (
        Rule(
            "integral-exponent-as-float",
            {".rs": POWF, ".py": PY_POW},
            exponent_is_integral,
        ),
        Rule(
            "decimal-exponent-approximating-a-rational",
            {".rs": POWF, ".py": PY_POW},
            exponent_approximates_a_rational,
        ),
        Rule("physical-field-defaulted-to-zero", {".rs": defaults}, physical_field_defaulted),
        Rule("maximum-float-outside-a-seed", {".rs": F64_MAX}, maximum_float_outside_a_seed),
        Rule(
            "index-defaulted-to-zero",
            {".rs": ZERO_INDEX_DEFAULT},
            index_defaulted_to_zero,
            multiline=True,
        ),
    )


@dataclass
class Findings:
    """What one file produced, and what was excused."""

    messages: list[str] = field(default_factory=list)
    excused: dict[str, int] = field(default_factory=dict)

    def excuse(self, rule: str, count: int) -> None:
        self.excused[rule] = self.excused.get(rule, 0) + count


def exempt(lines: list[str], number: int) -> bool:
    """Whether a `numerics-ok:` marker covers one line.

    **A marker opens a block that runs to the next blank line**, rather than covering one
    line. An expression that carries a rational-approximating decimal is often wrapped
    across several lines, and a per-line rule missed the continuation - which is how the
    first version of this let `1.2083` through while exempting `0.3333` three lines above
    it. So this walks back to the previous blank line, which is the block's own start.
    """
    for index in range(number - 1, -1, -1):
        if not lines[index].strip():
            return False
        if EXEMPT in lines[index]:
            return True
    return False


def statements(lines: list[str]) -> list[tuple[int, str]]:
    """A file's statements as `(first line number, comment-stripped text)`.

    A statement ends at a line whose code ends in `;` or `}`, or at a blank line - which is
    also where a marker's block ends, so the two agree on what a unit is.
    """
    out: list[tuple[int, str]] = []
    start = 0
    buffer: list[str] = []
    for number, line in enumerate(lines, start=1):
        # A `//` inside a string truncates the statement rather than the comment. That can
        # only shorten what a rule sees, and the alternative is a lexer this tool is not.
        code = line.split("//", 1)[0]
        if not line.strip():
            if buffer:
                out.append((start, "\n".join(buffer)))
                buffer = []
            continue
        if not buffer:
            start = number
        if code.strip():
            buffer.append(code)
        if code.rstrip().endswith((";", "}")):
            out.append((start, "\n".join(buffer)))
            buffer = []
    if buffer:
        out.append((start, "\n".join(buffer)))
    return out


def scan(path: Path) -> Findings:
    """Every line of one file that breaks a rule, and how many a marker excused."""
    findings = Findings()
    text = path.read_text(encoding="utf-8")
    if GENERATED in text[:2000]:
        return findings
    lines = text.splitlines()
    for number, line in enumerate(lines, start=1):
        if not line.strip() or COMMENT.match(line):
            continue
        for rule in RULES:
            if rule.multiline:
                continue
            found = rule.messages(path, line, number)
            if not found:
                continue
            if exempt(lines, number):
                findings.excuse(rule.name, len(found))
            else:
                findings.messages.extend(found)
    for start, statement in statements(lines):
        for rule in RULES:
            if not rule.multiline:
                continue
            found = rule.messages(path, statement, start)
            if not found:
                continue
            if exempt(lines, start):
                findings.excuse(rule.name, len(found))
            else:
                findings.messages.extend(found)
    return findings


#: Built once, at import: one of them reads the component schema, and the tests redirect
#: `ROOT` to a synthetic tree that has no `specs/`.
RULES: tuple[Rule, ...] = build_rules()


def offenders(path: Path) -> list[str]:
    """Every line of one file that breaks one of the rules."""
    return scan(path).messages


def main() -> int:
    messages: list[str] = []
    excused: dict[str, int] = {}
    for directory in SCANNED:
        for path in sorted(ROOT.glob(f"{directory}/**/*")):
            if path.suffix in (".rs", ".py") and path.is_file():
                findings = scan(path)
                messages.extend(findings.messages)
                for rule, count in findings.excused.items():
                    excused[rule] = excused.get(rule, 0) + count
    names = sorted({rule.name for rule in RULES})
    for name in names:
        print(f"check_numerics: {name}: {excused.get(name, 0)} exemption(s) marked")
    if messages:
        for message in messages:
            print(f"  ERROR  {message}", file=sys.stderr)
        print(
            f"\ncheck_numerics: FAILED with {len(messages)} problem(s). See "
            f"docs/src/calculus/numerics.md for the rule.",
            file=sys.stderr,
        )
        return 1
    print(
        f"check_numerics: OK ({len(names)} rule(s), {sum(excused.values())} exemption(s) "
        f"marked, no unmarked finding)"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
