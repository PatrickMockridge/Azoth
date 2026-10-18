#!/usr/bin/env python3
"""Hold every executable claim to a test that would falsify it.

# Why this exists

A `panic!`, an `assert!` or a guarantee-worded `expect` is a statement that a state
cannot happen. Those statements were prose in the tree: the claim was in the code, the
reason was in a comment, and nothing connected either to a test. So a state the code
says is unreachable is an assertion no run can falsify, and a comment that says a branch
is dead is a claim nothing checks.

This tool makes the claim a marker. Every executable claim carries
``claim-ok: <test-name>`` on its line or in the block above it, and the named test must
exist. The marker is the forcing function and the named test is the behavioural half;
**the checker cannot prove the test drives the branch** - it says so rather than implying
otherwise - but it can prove a name exists and that the marker was not left behind.

# What it checks

Rust, in ``crates/*/src/**/*.rs`` up to the first ``#[cfg(test)]`` and never in
``azoth-test-support`` - whose panics are its test verdict rather than a claim about the
library:

* ``panic!``, ``unreachable!``, ``assert!``, ``debug_assert!`` and ``debug_assert_eq!``;
* ``.expect("...")`` whose message carries a guarantee word, because the other expects
  are error plumbing rather than claims.

Python, in ``python/src/azoth/**/*.py`` except the generated ``*_gen.py``:

* an ``assert`` statement, and ``raise AssertionError``. Read through ``ast``, so a
  docstring or a comment that quotes one is prose rather than a claim.

**A marker names a test that exists.** Rust test functions are read from
``crates/*/{src,tests}/**/*.rs`` and Python ones from ``python/tests/**/test_*.py``. A
marker naming nothing is an error, and so is a marker with no claim in its block - the
check runs in both directions, so neither a claim nor an exemption can be left behind.

**A marker opens a block that runs to the next blank line**, the same rule
``numerics-ok:`` uses: a claim wrapped across several lines is one statement, and its
marker sits above it.

# Carve-outs, and why each is a distinction

**A comment is not a claim.** A line whose first non-space is ``//`` is prose - the tree
quotes the forms on purpose - and flagging it would make the rule about the comment.

**A generated file is not swept.** Its text is a spec's text, and a claim emitted from a
spec is held by the spec.

Usage:
    python tools/check_claims.py
"""

from __future__ import annotations

import ast
import re
import sys
from collections.abc import Iterable, Iterator, Sequence
from dataclasses import dataclass
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent

#: The trees a claim can live in.
RUST_TREES = ("crates/*/src",)
PYTHON_TREE = "python/src/azoth"

#: Where a marker's named test can be found.
RUST_TEST_TREES = ("crates/*/src", "crates/*/tests")
PYTHON_TEST_TREE = "python/tests"

#: **The one crate whose panics are not claims.** `azoth-test-support` is the test
#: harness: a panic there is the verdict it exists to deliver, not a statement about the
#: library's reachable states.
EXEMPT_CRATES = ("azoth-test-support",)

#: The banner a generator puts on its output.
GENERATED = "GENERATED FILE - DO NOT EDIT BY HAND"

#: The marker, and the test name it must carry.
MARKER = re.compile(r"claim-ok:\s*(?P<test>[A-Za-z_][A-Za-z0-9_:]*)")

#: The Rust forms that are executable claims.
RUST_CLAIMS = (
    re.compile(r"\bpanic!\s*\("),
    re.compile(r"\bunreachable!\s*\("),
    re.compile(r"\bassert!\s*\("),
    re.compile(r"\bdebug_assert!\s*\("),
    re.compile(r"\bdebug_assert_eq!\s*\("),
)

#: An `.expect` whose message says the state cannot happen. The literal is captured so
#: the guarantee word can be read from the message rather than from the call.
EXPECT = re.compile(r'\.expect\(\s*"(?P<message>(?:[^"\\]|\\.)*)"')

#: The phrases that make an `expect` a claim. A message using one of these says the branch
#: is unreachable, impossible or guaranteed - which is a statement to falsify. A message
#: that states the caller's requirement ("the field must be positive") is error plumbing
#: rather than a claim about this code's reachable states, so the bare words `must`,
#: `always` and `never` are not enough.
GUARANTEE_WORDS = (
    "unreachable",
    "impossible",
    "cannot happen",
    "can't happen",
    "cannot fail",
    "cannot be empty",
    "must hold",
    "must never",
    "must exist",
    "always holds",
    "always true",
    "never happens",
    "can never",
    "is never empty",
    "is never zero",
    "checked above",
    "checked present",
    "checked before",
    "the caller checked",
    "asked for above",
    "was asked for",
    "every row exists",
    "guaranteed",
    "invariant",
    "by construction",
    "by definition",
)

#: A Rust comment line, which is prose about a claim rather than a claim.
RUST_COMMENT = re.compile(r"^\s*//")

#: Where a source file's tests begin, and therefore where the sweep stops.
CFG_TEST = "#[cfg(test)]"

#: A Rust function name, which is what a marker names.
RUST_FUNCTION = re.compile(r"\bfn\s+(?P<name>[A-Za-z_][A-Za-z0-9_]*)")


def display(path: Path) -> str:
    """A path for a message: relative to the checkout where it is under it."""
    try:
        return str(path.relative_to(ROOT))
    except ValueError:
        return str(path)


@dataclass(frozen=True)
class Claim:
    """One executable claim: where it is, and what it says."""

    path: Path
    line: int
    kind: str
    text: str


def rust_files(root: Path) -> Iterator[Path]:
    """Every non-test Rust source file under a tree, in path order."""
    for directory in RUST_TREES:
        for path in sorted(root.glob(f"{directory}/**/*.rs")):
            if any(part in path.parts for part in EXEMPT_CRATES):
                continue
            if path.is_file():
                yield path


def python_files(root: Path) -> Iterator[Path]:
    """Every hand-written Python file under the package tree, in path order."""
    for path in sorted(root.glob(f"{PYTHON_TREE}/**/*.py")):
        if path.name.endswith("_gen.py"):
            continue
        if path.is_file():
            yield path


def generated(path: Path) -> bool:
    """Whether a file's own banner says a generator wrote it."""
    head = path.read_text(encoding="utf-8")[:2000]
    return GENERATED in head


def marker(lines: Sequence[str], number: int) -> str | None:
    """The test a line's block names, walking back to the previous blank line.

    `number` is 1-based and the line itself is included: a marker may sit on the claim's
    line or anywhere in the block above it.
    """
    for index in range(number - 1, -1, -1):
        if not lines[index].strip():
            return None
        found = MARKER.search(lines[index])
        if found is not None:
            return found.group("test")
    return None


def statement_cut(lines: Sequence[str]) -> int:
    """The line where a Rust file's test module begins, or the file's length.

    The sweep stops there because everything after is a test's own assertion, where the
    test *is* the falsifier and a marker would name itself.
    """
    for index, line in enumerate(lines):
        if line.strip().startswith(CFG_TEST):
            return index
    return len(lines)


def rust_claims(path: Path) -> list[Claim]:
    """Every executable claim in one Rust source file."""
    if generated(path):
        return []
    lines = path.read_text(encoding="utf-8").splitlines()
    cut = statement_cut(lines)
    out: list[Claim] = []
    for index in range(cut):
        line = lines[index]
        if not line.strip() or RUST_COMMENT.match(line):
            continue
        claim = next((pattern for pattern in RUST_CLAIMS if pattern.search(line)), None)
        if claim is not None:
            out.append(Claim(path, index + 1, "rust", line.strip()))
            continue
        expectation = EXPECT.search(line)
        if expectation is not None:
            message = expectation.group("message").lower()
            if any(word in message for word in GUARANTEE_WORDS):
                out.append(Claim(path, index + 1, "rust", line.strip()))
    return out


def python_claims(path: Path) -> list[Claim]:
    """Every executable claim in one Python source file, read through the AST."""
    if generated(path):
        return []
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    lines = path.read_text(encoding="utf-8").splitlines()
    out: list[Claim] = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Assert):
            out.append(Claim(path, node.lineno, "python", lines[node.lineno - 1].strip()))
        elif isinstance(node, ast.Raise) and node.exc is not None:
            raised = node.exc.func if isinstance(node.exc, ast.Call) else node.exc
            name = (
                raised.id
                if isinstance(raised, ast.Name)
                else raised.attr
                if isinstance(raised, ast.Attribute)
                else ""
            )
            if name == "AssertionError":
                out.append(Claim(path, node.lineno, "python", lines[node.lineno - 1].strip()))
    return sorted(out, key=lambda claim: claim.line)


def claims(rust_root: Path, python_root: Path) -> list[Claim]:
    """Every claim under the two source trees."""
    out: list[Claim] = []
    for path in rust_files(rust_root):
        out.extend(rust_claims(path))
    for path in python_files(python_root):
        out.extend(python_claims(path))
    return out


def rust_tests(root: Path) -> set[str]:
    """Every Rust function name a marker may name, from the source and the test files."""
    names: set[str] = set()
    for tree in RUST_TEST_TREES:
        for path in root.glob(f"{tree}/**/*.rs"):
            if not path.is_file():
                continue
            names.update(m.group("name") for m in RUST_FUNCTION.finditer(path.read_text("utf-8")))
    return names


def python_tests(root: Path) -> set[str]:
    """Every Python test function name a marker may name, from `python/tests`."""
    names: set[str] = set()
    for path in root.glob(f"{PYTHON_TEST_TREE}/**/test_*.py"):
        if not path.is_file():
            continue
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        names.update(
            node.name
            for node in ast.walk(tree)
            if isinstance(node, ast.FunctionDef | ast.AsyncFunctionDef)
        )
    return names


def markers_without_claims(path: Path, found: Sequence[Claim]) -> list[str]:
    """Markers whose block contains no claim, in both directions.

    A claim looks upward for its marker; this looks downward from every marker, so an
    exemption left behind by a deleted claim or a moved line is an error rather than a
    silence.
    """
    lines = path.read_text(encoding="utf-8").splitlines()
    claim_lines = {claim.line for claim in found}
    messages: list[str] = []
    for number, line in enumerate(lines, start=1):
        matches = list(MARKER.finditer(line))
        if not matches:
            continue
        covered = False
        for index in range(number - 1, len(lines)):
            if index > number - 1 and not lines[index].strip():
                break
            if index + 1 in claim_lines:
                covered = True
                break
        if not covered:
            messages.append(
                f"{display(path)}:{number}: `claim-ok:` names no claim in its "
                f"block; the claim it excused was removed or moved"
            )
    return messages


def check(
    rust_root: Path,
    python_root: Path,
    test_root: Path,
    *,
    rust_names: Iterable[str] | None = None,
    python_names: Iterable[str] | None = None,
) -> list[str]:
    """Every unmarked claim and every marker that names nothing.

    `rust_names` and `python_names` are passed by the synthetic-tree tests, whose test
    files live under their own root rather than this checkout's.
    """
    rust_known = set(rust_names) if rust_names is not None else rust_tests(test_root)
    python_known = set(python_names) if python_names is not None else python_tests(test_root)
    messages: list[str] = []
    by_file: dict[Path, list[Claim]] = {}
    for claim in claims(rust_root, python_root):
        by_file.setdefault(claim.path, []).append(claim)

    for path, found in by_file.items():
        lines = path.read_text(encoding="utf-8").splitlines()
        for claim in found:
            named = marker(lines, claim.line)
            where = f"{display(path)}:{claim.line}"
            if named is None:
                messages.append(
                    f"{where}: executable claim `{claim.text[:90]}` carries no "
                    f"`claim-ok: <test>` marker naming the test that would falsify it"
                )
            elif named not in rust_known and named not in python_known:
                messages.append(
                    f"{where}: `claim-ok: {named}` names no test. A Rust test must be a "
                    f"function under crates/*/{{src,tests}} and a Python one under "
                    f"python/tests"
                )
        messages.extend(markers_without_claims(path, found))

    # A marker in a file with no claim at all is still a marker with nothing to excuse.
    swept = set(by_file)
    for path in [*rust_files(rust_root), *python_files(python_root)]:
        if path in swept or generated(path):
            continue
        lines = path.read_text(encoding="utf-8").splitlines()
        for number, line in enumerate(lines, start=1):
            if MARKER.search(line):
                messages.append(
                    f"{display(path)}:{number}: `claim-ok:` names no claim in its block"
                )
    return messages


def main() -> int:
    messages = check(ROOT, ROOT, ROOT)
    found = claims(ROOT, ROOT)
    marked = 0
    for claim in found:
        lines = claim.path.read_text(encoding="utf-8").splitlines()
        if marker(lines, claim.line) is not None:
            marked += 1
    print(
        f"check_claims: {len(found)} executable claim(s), {marked} marked, "
        f"{len(found) - marked} unmarked"
    )
    if messages:
        for message in messages:
            print(f"  ERROR  {message}", file=sys.stderr)
        print(
            f"\ncheck_claims: FAILED with {len(messages)} problem(s). A claim names the "
            f"test that would falsify it, or it is converted to a typed refusal.",
            file=sys.stderr,
        )
        return 1
    print("check_claims: OK (every claim names an existing test)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
