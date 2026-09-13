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

The other thing checked here is the *convention*, which is subtler than the set. The
Python side works in the magnitude of the unit's canonical form; the Rust side's
``uom`` quantities always report ``.value`` as the SI **base** magnitude. Those are
the same number only when a unit's canonical form is its SI base unit. Every entry
here satisfies that except ``mm``, which is recorded rather than hidden - see
``NON_SI_BASE_UNITS``.
"""

from __future__ import annotations

import importlib
import json
from pathlib import Path
from types import ModuleType
from typing import Any, Final

import pytest

from azoth._registry_gen import CALCS
from azoth.core.units import CANONICAL_UNITS, quantity

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


#: Units whose canonical form is not the SI base unit, mapped to the SI base
#: magnitude of one canonical unit.
#:
#: ``mm`` is the only entry and it is a **known defect, not a decision**. A spec
#: declaring ``mm`` would have the Python side working in millimetres and the Rust
#: side in metres, and the cross-language agreement test would fail by 1000x rather
#: than by a rounding error. It has never mattered because no spec uses ``mm``.
#:
#: It is pinned rather than fixed or deleted for two reasons. It cannot be deleted
#: on the evidence - pipe diameters are conventionally quoted in millimetres, so
#: this is the unit the next hydraulics calc is most likely to want. And it cannot
#: be fixed here: making the two sides agree means separating "the unit the spec
#: declares" from "the unit the calculation works in", which is a change to
#: ``to_si``, ``quantity``, both test helpers and both languages.
#:
#: Recording it as an exception means a *second* unit with the same problem fails
#: this test instead of joining the set. When the convention is made coherent this
#: dict should become empty, and the test will say so.
NON_SI_BASE_UNITS: Final[dict[str, float]] = {"mm": 1.0e-3}


def schema_units() -> set[str]:
    """The unit strings the spec schema permits."""
    schema: dict[str, Any] = json.loads(SCHEMA_PATH.read_text(encoding="utf-8"))
    return set(schema["$defs"]["unit"]["enum"])


def si_base_magnitude(spec_unit: str) -> float:
    """The SI base magnitude of one unit of ``spec_unit``.

    ``1.0`` for a unit that *is* its own SI base unit, which is every entry the
    vocabulary is allowed to contain.
    """
    base = quantity(1.0, spec_unit).to_base_units()
    return float(base.magnitude)


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


def test_the_canonical_unit_is_the_si_base_unit() -> None:
    """Every unit is its own SI base unit, except the recorded exceptions.

    This is the convention that makes the two implementations comparable. The
    Python side hands a calculation the magnitude of ``spec_unit``; the Rust side
    hands it ``.value``, which is always the SI base magnitude. The two agree
    exactly when ``spec_unit`` is the SI base unit for that dimension.

    Stated as an equality against a named exception set, so a new unit that breaks
    the rule fails here rather than by 1000x somewhere downstream.
    """
    factors = {unit: si_base_magnitude(unit) for unit in CANONICAL_UNITS}
    not_base = {unit: factor for unit, factor in factors.items() if factor != 1.0}
    assert not_base == NON_SI_BASE_UNITS, (
        f"the set of units that are not their own SI base unit has changed:\n"
        f"  observed: {not_base}\n"
        f"  recorded: {NON_SI_BASE_UNITS}\n"
        f"A new entry means a unit whose two implementations disagree. An entry "
        f"that has gone means the convention was fixed - update this list."
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
