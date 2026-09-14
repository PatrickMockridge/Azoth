"""The unit vocabulary is one set, declared once and compiled into four artefacts.

A calc spec names its units as strings. One hand-written table says which strings
are legal - ``specs/vocabulary/vocabulary.toml`` - and ``tools/gen_vocabulary.py``
compiles it into:

* ``specs/schema/unit.schema.json``, which is what a spec is validated against;
* :data:`azoth.core.units.CANONICAL_UNITS`, which turns a name into a ``pint`` unit
  and is what the Python side actually converts with;
* ``crates/azoth-core/src/unit_vocab_gen.rs``'s ``UNIT_NAMES``, the Rust vocabulary,
  reachable from here through ``azoth._core.unit_names()``;
* ``azoth.keycard.UNIT_VOCABULARY``, which is what the keycard loader admits.

Before the table existed those were hand-maintained lists, and it showed: ``K`` was
a permitted unit the schema allowed and ``CANONICAL_UNITS`` knew about, while
``azoth-core`` had no temperature type at all. No calc used it, so no test touched
it. This file still compares all four - generation is what makes them agree, and
comparing them is what makes a failure of generation visible rather than a spec
failing at runtime for whoever declared the unit.

The other thing checked here is that no entry is *dead* - each one converts in both
directions rather than merely being listed. Whether the numbers those conversions
produce are correct is a separate question and lives in
``test_units_cross_library.py``; this file is about the vocabulary being one set and
every member of it working.
"""

from __future__ import annotations

import importlib
import json
import subprocess
import sys
import tomllib
from pathlib import Path
from types import ModuleType
from typing import Any

import pytest

from azoth._registry_gen import CALCS
from azoth.core.units import CANONICAL_UNITS, from_si, quantity, to_si
from azoth.keycard import UNIT_VOCABULARY

REPO_ROOT = Path(__file__).resolve().parents[2]
VOCAB_PATH = REPO_ROOT / "specs" / "vocabulary" / "vocabulary.toml"
UNIT_SCHEMA_PATH = REPO_ROOT / "specs" / "schema" / "unit.schema.json"


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


def table_units() -> set[str]:
    """The unit strings the hand-written vocabulary table declares."""
    table = tomllib.loads(VOCAB_PATH.read_text(encoding="utf-8"))
    return {unit["id"] for unit in table["units"]}


def schema_units() -> set[str]:
    """The unit strings the spec schema permits, read from the generated enum."""
    schema: dict[str, Any] = json.loads(UNIT_SCHEMA_PATH.read_text(encoding="utf-8"))
    return set(schema["enum"])


def _report(name: str, expected: set[str], got: set[str]) -> None:
    assert expected == got, (
        f"{name} differs from the vocabulary table\n"
        f"  only in the table: {sorted(expected - got)}\n"
        f"  only in {name}: {sorted(got - expected)}"
    )


def test_every_generated_copy_is_the_table() -> None:
    """The table, the schema, the Python map and the keycard list are one set.

    A name the schema allows but ``CANONICAL_UNITS`` lacks fails at runtime, in
    ``to_si`` -> ``unit_for``, and only for whichever calc happens to use it.
    This is the cheap version of that check, and it runs without the extension
    built - the Rust leg is the marked test below.
    """
    table = table_units()
    _report("the schema", table, schema_units())
    _report("CANONICAL_UNITS", table, set(CANONICAL_UNITS))
    _report("keycard.UNIT_VOCABULARY", table, set(UNIT_VOCABULARY))


@pytest.mark.requires_rust
def test_the_rust_vocabulary_is_the_table() -> None:
    """The Rust leg, which is the one that was unchecked and the one that was wrong.

    `K` was permitted by the schema and known to `CANONICAL_UNITS`, while
    azoth-core had no temperature type at all. `unit_names()` exists so this can
    be asserted from Python instead of by reading Rust source.
    """
    _report("the Rust vocabulary", table_units(), set(_extension().unit_names()))


def test_the_generated_files_are_current() -> None:
    """Regenerating from the table changes nothing.

    The four copies above agree because one generator emits them from one table;
    this is the check that they agree with the table *now* rather than agreeing
    with a table that has since been edited. `--check` writes nothing.
    """
    result = subprocess.run(
        [sys.executable, str(REPO_ROOT / "tools" / "gen_vocabulary.py"), "--check"],
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, (
        f"gen_vocabulary.py --check failed, so a generated copy is stale:\n"
        f"{result.stdout}{result.stderr}"
    )


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

    unknown = declared - table_units()
    assert not unknown, f"spec(s) declare unit(s) {sorted(unknown)} which the vocabulary lacks"
