"""The unit vocabulary is one set, agreed on by three separate places.

A calc spec names its units as strings. Three artifacts have to agree about which
strings are legal and what each one means:

* ``specs/schema/calc.schema.json``'s ``$defs.unit.enum``, which is what a spec is
  validated against;
* :data:`azoth.core.units.CANONICAL_UNITS`, which turns a name into a ``pint`` unit
  and is what the Python side actually converts with;
* ``crates/azoth-core/src/units.rs``'s ``UNIT_NAMES``, the Rust vocabulary, reachable
  from here through ``azoth._core.unit_names()``.

Nothing checked the Rust leg before, and it showed: ``K`` was a permitted unit the
schema allowed and ``CANONICAL_UNITS`` knew about, while ``azoth-core`` had no
temperature type at all. No calc used it, so no test touched it. A vocabulary is
only a contract if something compares the copies, which is what this file does.

The other thing checked here is that no entry is *dead* - each one converts in both
directions rather than merely being listed. Whether the numbers those conversions
produce are correct is a separate question and lives in
``test_units_conversion.py``; this file is about the vocabulary being one set and
every member of it working.
"""

from __future__ import annotations

import importlib
import json
from pathlib import Path
from types import ModuleType
from typing import Any

import pytest

from azoth._registry_gen import CALCS
from azoth.core.units import CANONICAL_UNITS, from_si, quantity, to_si

REPO_ROOT = Path(__file__).resolve().parents[2]
SCHEMA_PATH = REPO_ROOT / "specs" / "schema" / "calc.schema.json"


def _extension() -> ModuleType:
    """The compiled extension, imported by name.

    The same accessor `test_cross_impl.py` uses, and for the same reason: the
    module does not exist until the bindings are built, and an attribute mypy
    cannot resolve is a worse trade than a lookup that fails clearly.
    """
    try:
        return importlib.import_module("azoth._core")
    except ImportError as exc:  # pragma: no cover - the marker skips these
        raise AssertionError("azoth._core is not built; run `maturin develop`") from exc


def schema_units() -> set[str]:
    """The unit strings the spec schema permits."""
    schema: dict[str, Any] = json.loads(SCHEMA_PATH.read_text(encoding="utf-8"))
    return set(schema["$defs"]["unit"]["enum"])


def test_the_schema_and_the_python_vocabulary_agree() -> None:
    """The schema's enum and ``CANONICAL_UNITS`` are the same set.

    A name the schema allows but this dict lacks fails at runtime, in
    ``to_si`` -> ``unit_for``, and only for whichever calc happens to use it -
    ``spec_lint`` does not check units at all. This is the cheap version of that
    check, and it runs without the extension built.
    """
    schema, python_side = schema_units(), set(CANONICAL_UNITS)
    assert schema == python_side, (
        f"unit vocabularies differ\n"
        f"  only in the schema: {sorted(schema - python_side)}\n"
        f"  only in CANONICAL_UNITS: {sorted(python_side - schema)}"
    )


@pytest.mark.requires_rust
def test_the_rust_vocabulary_agrees_with_both() -> None:
    """All three lists are one set.

    The Rust leg is the one that was unchecked, and the one that was wrong: `K`
    was permitted by the schema and known to `CANONICAL_UNITS`, while azoth-core
    had no temperature type at all. `unit_names()` exists so this can be asserted
    from Python instead of by reading Rust source.
    """
    rust_side = set(_extension().unit_names())
    assert schema_units() == rust_side, (
        f"the Rust vocabulary differs from the schema\n"
        f"  only in the schema: {sorted(schema_units() - rust_side)}\n"
        f"  only in Rust: {sorted(rust_side - schema_units())}"
    )
    assert rust_side == set(CANONICAL_UNITS)


def test_every_vocabulary_entry_has_a_working_conversion() -> None:
    """Each name converts, in both directions, without raising.

    The vocabulary is a promise that a spec may declare this unit. A name that
    ``pint`` cannot parse, or that ``CANONICAL_UNITS`` maps to something that is not
    a unit, would break that promise at runtime and only for whoever declared it.

    This says nothing about whether the *numbers* are right - that is
    ``test_units_conversion.py``, which is a different question and a different
    file. Here the concern is only that the vocabulary has no dead entries.
    """
    for unit in CANONICAL_UNITS:
        # A round trip through both directions, which is the path a real value
        # takes: a caller's quantity in, a result quantity out.
        si = to_si(quantity(1.0, unit), unit, "x")
        assert isinstance(si, float), f"{unit}: to_si did not return a float"
        rebuilt = from_si(si, unit)
        assert rebuilt.dimensionality == quantity(1.0, unit).dimensionality, (
            f"{unit}: from_si(to_si(...)) changed the dimensionality"
        )


def test_every_unit_a_spec_declares_is_in_the_vocabulary() -> None:
    """No spec may name a unit the vocabulary does not carry.

    The schema already enforces this, so this test is a guard on the guard: it
    proves the schema is actually being applied to the specs in the registry, and
    that the generated registry is current.
    """
    assert CALCS, "no calcs in the registry - this test would pass vacuously"

    declared: set[str] = set()
    for calc in CALCS:
        for section in ("inputs", "outputs"):
            for declaration in calc[section].values():
                unit = declaration.get("unit")
                if unit is not None:
                    declared.add(unit)

    unknown = declared - schema_units()
    assert not unknown, f"spec(s) declare unit(s) {sorted(unknown)} which the vocabulary lacks"
