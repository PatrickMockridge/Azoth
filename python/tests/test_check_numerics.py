"""The numerical-safety lint: the two rules, and the exemption that keeps them honest.

`tools/check_numerics.py` runs in CI and over a tree that is currently clean, which is
the shape of a check that can stop working without anyone noticing - an exemption that
widened, or a rule that stopped matching. So the rules are tested against synthetic
source rather than only against the tree.

**The exemption is the part most likely to rot.** It is a `numerics-ok:` marker that
opens a block running to the next blank line, and the first version of it covered one
line - which exempted `0.3333` while letting `1.2083` through on the statement's
continuation three lines below. That is why the block case is tested.
"""

from __future__ import annotations

import importlib
import sys
from pathlib import Path
from types import ModuleType
from typing import Any

REPO_ROOT = Path(__file__).resolve().parents[2]


def lint() -> ModuleType:
    """`tools/check_numerics.py`, imported by name. The tools are not a package."""
    sys.path.insert(0, str(REPO_ROOT / "tools"))
    try:
        return importlib.import_module("check_numerics")
    finally:
        sys.path.pop(0)


def offenders(tmp_path: Path, name: str, text: str) -> list[str]:
    """The lint's verdict on one synthetic file.

    The root is redirected because the message is relative to it, and a synthetic file
    under `/tmp` against the real root would raise from `relative_to` rather than report.
    `Any` at the point the tool is reached rather than a widened signature: a `ModuleType`
    declares neither `ROOT` nor `offenders`, and the tools are outside mypy's `files`.
    """
    path = tmp_path / name
    path.write_text(text, encoding="utf-8")
    tool: Any = lint()
    original = tool.ROOT
    tool.ROOT = tmp_path
    try:
        found: list[str] = tool.offenders(path)
        return found
    finally:
        tool.ROOT = original


def test_an_integral_exponent_written_as_a_float_is_reported(tmp_path: Path) -> None:
    """The Rust advice is `powi`, which exists and is total."""
    messages = offenders(tmp_path, "a.rs", "let y = x.powf(2.0);\n")
    assert len(messages) == 1, messages
    assert "powi(2)" in messages[0]

    # And `powi` itself is not reported.
    assert offenders(tmp_path, "b.rs", "let y = x.powi(2);\n") == []


def test_python_is_told_the_python_fix(tmp_path: Path) -> None:
    """**`powi` does not exist in Python**, and the message must not say it does.

    `x ** 3` with an *integer* exponent takes the real branch where `x ** 3.0` promotes a
    negative base to a complex, so the fix is the integer literal rather than a function.
    """
    messages = offenders(tmp_path, "a.py", "y = x ** 3.0\n")
    assert len(messages) == 1, messages
    assert "`** 3`" in messages[0]
    assert "powi" not in messages[0].split("If the receiver")[0]


def test_a_decimal_that_approximates_a_rational_is_reported(tmp_path: Path) -> None:
    """`0.3333` is not `1/3`, and the message names the rational it rounds."""
    messages = offenders(tmp_path, "a.rs", "let y = x.powf(0.3333);\n")
    assert len(messages) == 1, messages
    assert "1/3" in messages[0]


def test_an_exact_rational_is_left_alone(tmp_path: Path) -> None:
    """`0.5`, `0.25` and `1.5` are the constants they look like, not approximations."""
    for literal in ("0.5", "0.25", "0.75", "1.5", "0.1"):
        assert offenders(tmp_path, "a.rs", f"let y = x.powf({literal});\n") == [], literal
    # And a fitted constant that is not a rational at all.
    assert offenders(tmp_path, "a.rs", "let y = x.powf(2.303);\n") == []


def test_a_marker_exempts_its_whole_statement(tmp_path: Path) -> None:
    """**The block case, which the first version got wrong.**

    The marker covers every following line up to a blank one, so a wrapped expression is
    exempt once rather than once per physical line.
    """
    text = (
        "// numerics-ok: NeqSim ComponentGEWilson.java:146 writes 0.3333\n"
        "let d0 = 5.2804 * x.powf(0.3333)\n"
        "    + 12.865 * x.powf(0.8333)\n"
        "    + 1.171 * x.powf(1.2083);\n"
    )
    assert offenders(tmp_path, "a.rs", text) == []

    # The block ends at a blank line, so the next statement is checked again.
    text += "\nlet unmarked = x.powf(0.3333);\n"
    messages = offenders(tmp_path, "a.rs", text)
    assert len(messages) == 1, messages
    assert ":6:" in messages[0], messages[0]


def test_a_generated_file_is_not_scanned(tmp_path: Path) -> None:
    """A generator emits an `equation` string and a worked example, not arithmetic."""
    text = "//! GENERATED FILE - DO NOT EDIT BY HAND.\nlet y = x.powf(2.0);\n"
    assert offenders(tmp_path, "a.rs", text) == []
