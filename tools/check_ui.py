#!/usr/bin/env python3
"""Fail on the shapes that let a front end hold a second answer, or invent a number.

# Why this exists

`docs/src/architecture/middleware.md` says the editor "holds no copy of the flowsheet" and that
"a unit set is a reading of a run, not a second run". Both are true, and both were true only
because a person kept them so: `ui/src` is not read by `check_numerics.py` (which scans the two
kernel trees), by `prose_lint.py` (which reads the UI for *history*, not for arithmetic) or by
anything else. The vocabulary page's rule — "no conversion factor appears in the table, or in any
generated file, or **anywhere else in this repository**" — had three trees behind it and a fourth
nobody looked at.

# The three rules

Each is decidable from a line of source, which is the standard the numeric-safety rules are held
to, and each was verified by breaking the thing it guards.

**A cast on the boundary.** `JSON.parse(text) as Envelope` asserts a shape where the wire gave an
untyped document. The symptom of a wrong one is `undefined` in a widget three files from the field
that moved, and the fix is a parse that says which path was wrong — `ui/src/wire/decode.ts`, which
is the one file allowed to cast because a cast inside a checked function is the whole technique.
A second exempt file is a second parser, and the count of them is printed.

**A measurement written in the front end.** `6894.757`, `273.15`, `0.0254` — a decimal with three
or more significant digits, or an integer with four or more, that no display computed. A factor is
a number `uom` and `pint` each already know and the catalogue carries; a number written here is a
fourth copy that can disagree with all of them, which is the defect the vocabulary's gate exists to
refuse. `1`, `0`, `4`, `0.5` and `0.75` are structure and fractions rather than units' worths, and
the threshold is where the two stop being distinguishable from the source.

**A component holding a part of the document.** `useState<Graph>`, `useState<GraphNode[]>`,
`useState<Diagnostic[]>` — a *part* of what the envelope carries, kept beside it. The envelope
itself in state is the design (`App.tsx` is one document and no reducer over the graph); a graph
in state is the second copy the claim says does not exist, and the two then disagree after any
edit that only one of them saw.

# The exemption, and why it is a comment rather than a list

A line is exempt when it, or the line immediately above it in the same block, carries a `ui-ok:`
marker naming why:

    const value = thing as Role;  // ui-ok: the projection guarantees the prefix

A whole file is exempt when its banner carries `ui-ok-file:` and a reason. The marker travels with
what it excuses, which a list of file-and-line pairs does not, and **the count is printed per
rule**, so an exemption that widens is visible in the output rather than being a silence.

# What it does not do

- **Not the stylesheets or the drawing.** `ui/src/styles/**` and the glyph geometry are a
  presentation of a machine, not a quantity in it, and a rule that flagged a radius would be a
  rule about SVG.
- **Not a type checker.** `tsc --noEmit` is the type gate and runs in the same job; this reads for
  the three shapes a type cannot see, which is why it is a lint and not a fourth type.
- **Not the tests.** `ui/test/**` is where a wrong number is *supposed* to appear.

Usage:
    python tools/check_ui.py
"""

from __future__ import annotations

import re
import sys
from collections.abc import Callable
from dataclasses import dataclass, field
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent

#: The tree the editor's own code lives in. **`ui/src` and not `ui/`**, so the tests are out of
#: it: a test asserting a conversion is the place a literal belongs.
SCANNED = "ui/src"

#: The module `wasm-pack` writes, which is generated, is not TypeScript anybody wrote, and is
#: gitignored — but it exists in a built working tree and would be read otherwise.
GENERATED = ("/wasm/",)

#: A cast to a type: `x as Role`, `x as unknown as Catalogue`, `x as Record<string, unknown>`.
#: A re-export alias (`export type { EditorCommand as Command }`) is not one, which is what the
#: import/export guard below is for.
CAST = re.compile(r"\bas\s+(?:unknown\s+as\s+)?[A-Z]\w*")

#: Import and export lines, where `as` renames a binding rather than asserting a type.
REEXPORT = re.compile(r"^\s*(?:import|export)\b")

#: A literal that looks measured: a decimal with three or more significant digits, or an integer
#: with four or more. `6894.757`, `0.0254` and `273.15` are the shapes; `1`, `0`, `4`, `0.5` and
#: `0.75` are structure and a fraction rather than a unit's worth. **Significant digits and not
#: decimal places**: `273.15` has two of the latter and five of the former, which is the whole
#: reason the earlier form of this pattern missed it.
MEASUREMENT = re.compile(r"(?<![\w.])(?P<literal>\d+\.\d+|[1-9]\d{3,})(?![\w.])")

#: The count a literal needs to be one of those.
SIGNIFICANT = 3

#: Local state whose type is a *part* of the envelope rather than the envelope.
#:
#: `Envelope` is the whole document and holding it is the design; `Session` is a door's handle and
#: `Catalogue` is the library's own answer about itself, neither of which an edit changes. These
#: are the things an edit *does* change, so a copy of one is a copy that goes stale.
PARTS = (
    "Graph",
    "GraphNode",
    "GraphEdge",
    "NodeData",
    "Ports",
    "Handle",
    "InputRecord",
    "Diagnostic",
    "StreamRecord",
    "TearRecord",
    "SessionReport",
    "Quantity",
    "Form",
    "FormPort",
    "FormParameter",
)
COPY_IN_STATE = re.compile(
    rf"\buse(?:State|Ref)<\s*(?P<part>{'|'.join(PARTS)})\b(?:\[\])?\s*\|?\s*(?:null)?\s*>"
)

#: How a line or a file says its shape is deliberate.
EXEMPT = "ui-ok:"
FILE_EXEMPT = "ui-ok-file:"

#: A comment line, in either of the two languages' styles plus the `*` continuation of a block.
COMMENT = re.compile(r"^\s*(?://|/\*|\*)")


@dataclass(frozen=True)
class Rule:
    """One syntactic rule: a name, what recognises it, and what a match means."""

    name: str
    pattern: re.Pattern[str]
    explain: Callable[[re.Match[str], str, str], str | None]


def cast_on_the_boundary(match: re.Match[str], line: str, where: str) -> str | None:
    """`JSON.parse(...) as Envelope`, which asserts the shape the wire did not promise."""
    _ = line
    return (
        f"{where}: `{match.group(0)}` asserts a shape. What crossed the wire is an untyped "
        f"document, and a field that moved produces `undefined` here rather than a refusal — "
        f"decode it instead (`ui/src/wire/decode.ts`), or mark the line `{EXEMPT} <reason>` "
        f"where the value was already decoded and this narrows it back."
    )


def measurement_written_down(match: re.Match[str], line: str, where: str) -> str | None:
    """`6894.757`, which is a unit's worth written in the fourth place."""
    _ = line
    literal = match.group("literal")
    digits = len(literal.replace(".", "").lstrip("0"))
    # An integer is a measurement at four digits and a count below it; a decimal at three
    # significant digits, which is where a fraction stops and a fitted constant starts.
    if digits < SIGNIFICANT + (0 if "." in literal else 1):
        return None
    return (
        f"{where}: `{literal}` is a magnitude with {digits} significant digits. A "
        f"factor belongs to `uom` and `pint` and reaches the editor through the catalogue's own "
        f"`factor`/`offset`; a number written here is a copy that can disagree with both, which "
        f"is the defect `docs/src/calculus/vocabulary.md` refuses. Name what it is - a count, a "
        f"limit, a formatting width - or mark it `{EXEMPT} <reason>`."
    )


def a_part_of_the_document_in_state(match: re.Match[str], line: str, where: str) -> str | None:
    """`useState<Graph>(...)`, which is a second copy of what the envelope carries."""
    _ = line
    return (
        f"{where}: `{match.group(0)}` keeps a `{match.group('part')}` in local state, and the "
        f"envelope already carries it. **The document is the state** - a part held beside it goes "
        f"stale the first time an edit lands on one and not the other. Read it from the envelope, "
        f"or mark it `{EXEMPT} <reason>` where the value is a display's own and not the "
        f"document's."
    )


RULES: tuple[Rule, ...] = (
    Rule("cast-on-the-wire-boundary", CAST, cast_on_the_boundary),
    Rule("measurement-written-in-the-front-end", MEASUREMENT, measurement_written_down),
    Rule("part-of-the-document-in-local-state", COPY_IN_STATE, a_part_of_the_document_in_state),
)


@dataclass
class Findings:
    """What one file produced, and what was excused."""

    messages: list[str] = field(default_factory=list)
    excused: dict[str, int] = field(default_factory=dict)

    def excuse(self, rule: str, count: int) -> None:
        self.excused[rule] = self.excused.get(rule, 0) + count


def file_exempt(text: str) -> bool:
    """Whether a file's banner declares it a decoder or otherwise out of the rules' reach."""
    return FILE_EXEMPT in "\n".join(text.splitlines()[:60])


def exempt(lines: list[str], number: int) -> bool:
    """Whether a `ui-ok:` marker covers one line.

    **A marker opens a block that runs back to the previous blank line**, the way the
    numeric-safety rules' does: a cast is often the tail of an expression wrapped across several
    lines, and a per-line rule would miss the continuation it deliberately explains.
    """
    for index in range(number - 1, -1, -1):
        if not lines[index].strip():
            return False
        if EXEMPT in lines[index]:
            return True
    return False


def scan(path: Path) -> Findings:
    """Every line of one file that breaks a rule, and how many a marker excused."""
    findings = Findings()
    text = path.read_text(encoding="utf-8")
    if file_exempt(text):
        return findings
    lines = text.splitlines()
    for number, line in enumerate(lines, start=1):
        if not line.strip() or COMMENT.match(line) or REEXPORT.match(line):
            continue
        for rule in RULES:
            found = [
                message
                for match in rule.pattern.finditer(line)
                if (message := rule.explain(match, line, f"{path.relative_to(ROOT)}:{number}"))
                is not None
            ]
            if not found:
                continue
            if exempt(lines, number):
                findings.excuse(rule.name, len(found))
            else:
                findings.messages.extend(found)
    return findings


def scanned_files() -> list[Path]:
    """Every file the rules read: the editor's own TypeScript, and nothing generated."""
    root = ROOT / SCANNED
    if not root.exists():
        return []
    return [
        path
        for path in sorted(root.rglob("*"))
        if path.is_file()
        and path.suffix in (".ts", ".tsx")
        and not any(part in f"/{path.relative_to(root)}" for part in GENERATED)
    ]


def offenders(path: Path) -> list[str]:
    """Every line of one file that breaks one of the rules."""
    return scan(path).messages


def main() -> int:
    files = scanned_files()
    if not files:
        print(
            f"check_ui: FAILED - no TypeScript under {SCANNED}/. The tree is not where this "
            f"looks, and a rule that reads nothing passes.",
            file=sys.stderr,
        )
        return 1
    messages: list[str] = []
    excused: dict[str, int] = {}
    decoders = 0
    for path in files:
        text = path.read_text(encoding="utf-8")
        if file_exempt(text):
            decoders += 1
            continue
        findings = scan(path)
        messages.extend(findings.messages)
        for rule, count in findings.excused.items():
            excused[rule] = excused.get(rule, 0) + count
    for rule in RULES:
        print(f"check_ui: {rule.name}: {excused.get(rule.name, 0)} exemption(s) marked")
    if messages:
        for message in messages:
            print(f"  ERROR  {message}", file=sys.stderr)
        print(
            f"\ncheck_ui: FAILED with {len(messages)} problem(s). See "
            f"docs/src/architecture/middleware.md for what the editor claims.",
            file=sys.stderr,
        )
        return 1
    print(
        f"check_ui: OK ({len(RULES)} rule(s), {len(files)} file(s), "
        f"{sum(excused.values())} exemption(s) marked, {decoders} file(s) declared a decoder)"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
