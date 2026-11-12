"""The sweep's stale-literal classifier, on the cases its report depends on.

A moved key is `(name, was, is_now)`. `is_stale` answers whether a literal in a Rust test
is the *old* value, which is only decidable when the old value is nonzero: `relative` is
absolute at a zero reference, so without that guard every near-zero constant in the tree
matches every key that moved off zero.
"""

from __future__ import annotations

import importlib
import sys
from pathlib import Path
from types import ModuleType

REPO_ROOT = Path(__file__).resolve().parents[2]
TOOLS = REPO_ROOT / "tools"


def sweep() -> ModuleType:
    """`tools/oracle_sweep.py`, imported by name; the tools are not a package."""
    sys.path.insert(0, str(TOOLS))
    try:
        return importlib.import_module("oracle_sweep")
    finally:
        sys.path.pop(0)


def test_a_zero_old_value_never_matches() -> None:
    """The false positive that made the 26,299-line report unusable."""
    assert not sweep().is_stale(0.0, ("k", 0.0, 5.0))
    assert not sweep().is_stale(1e-50, ("k", 0.0, 5.0))


def test_the_old_literal_is_stale() -> None:
    assert sweep().is_stale(2.5, ("k", 2.5, 2.6))


def test_the_new_literal_is_not_stale() -> None:
    assert not sweep().is_stale(2.6, ("k", 2.5, 2.6))


def test_a_literal_equal_to_neither_is_not_stale() -> None:
    assert not sweep().is_stale(9.0, ("k", 2.5, 2.6))
