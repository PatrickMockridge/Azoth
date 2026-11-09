"""The shipped comparison policy is the tests' policy, checked rather than assumed."""

from __future__ import annotations

import importlib
from typing import Any

from azoth.rgen import policy

# Imported through the import system so mypy does not demand a stub for the tests' own
# helper module, which is a top-level module rather than part of the package.
_helpers: Any = importlib.import_module("_helpers")


def test_the_uncompared_field_sets_are_equal() -> None:
    assert set(policy.UNCOMPARED_FIELDS) == set(_helpers._UNCOMPARED_FIELDS)


def test_the_diagnostic_field_sets_are_equal() -> None:
    assert set(policy.DIAGNOSTIC_FIELDS) == set(_helpers._DIAGNOSTIC_FIELDS)
