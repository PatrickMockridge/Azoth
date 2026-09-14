"""The vocabulary table, held to `pint` and to `uom` at once.

`specs/vocabulary/vocabulary.toml` declares what each canonical unit is: its
dimension, its name in `pint`, and its conversion path in `uom`. Neither of those
libraries exposes a dimension as data - `uom` carries it in the type system and
`pint` in its own registry - so the table is the only place the dimension is
written down, and these tests are what keep it honest about both.

There are two directions, and they are different questions:

* **The dimension.** `pint` *can* be asked what a unit's dimensionality is, so
  :func:`test_every_unit_has_the_dimension_the_table_declares` asks it and compares
  the answer with the table's exponents. A table entry whose exponents were edited
  to something wrong fails here.
* **The factor.** `pint` and `uom` each know what one of their units is worth, and
  the two have to agree - but neither number is written down anywhere in this
  repository. :func:`test_the_two_units_libraries_agree_on_every_factor` compares
  them: `pint`'s answer for a unit name against the result of running the very
  conversion a calculation runs, which is reached through
  `azoth._core.unit_si_factor`.

That second one is the check this repository needed and did not have. A prior
version of the Rust vocabulary paired each unit with a hand-typed SI magnitude, and
`mm`'s was wrong: it said the two implementations agreed while putting a factor of
1000 between them, squared by the bore diameter of an orifice. Nothing here is
typed; every number in this file is one a units library reported.
"""

from __future__ import annotations

import importlib
import tomllib
from pathlib import Path
from types import ModuleType
from typing import Any, cast

import pytest

from azoth.core._units_gen import SLOTS
from azoth.core.units import ureg

REPO_ROOT = Path(__file__).resolve().parents[2]
VOCAB_PATH = REPO_ROOT / "specs" / "vocabulary" / "vocabulary.toml"

#: How `pint` spells each of the table's slots. A closed map rather than a lookup
#: on `pint`'s internal names, so a `pint` release that renames a base dimension
#: fails the bijection test below instead of silently matching nothing.
PINT_BASE = {
    "L": "[length]",
    "M": "[mass]",
    "T": "[time]",
    "I": "[current]",
    "Th": "[temperature]",
    "N": "[substance]",
    "J": "[luminosity]",
}

#: One `pint` unit per slot, used to check the slot really is a base dimension.
SLOT_UNITS = {
    "L": "meter",
    "M": "kilogram",
    "T": "second",
    "I": "ampere",
    "Th": "kelvin",
    "N": "mole",
    "J": "candela",
}


def _extension() -> ModuleType:
    """The compiled extension, imported by name. See `test_units_contract.py`."""
    try:
        return importlib.import_module("azoth._core")
    except ImportError as exc:  # pragma: no cover - the marker skips these
        raise AssertionError("azoth._core is not built; run `maturin develop`") from exc


def table() -> dict[str, Any]:
    """The hand-written vocabulary table."""
    return tomllib.loads(VOCAB_PATH.read_text(encoding="utf-8"))


def units() -> list[dict[str, Any]]:
    """Every unit in the table, in table order."""
    return cast("list[dict[str, Any]]", table()["units"])


def _exponents(unit: dict[str, Any]) -> tuple[int, ...]:
    """The exponents the table gives one unit, through its named dimension."""
    dimensions = {d["id"]: tuple(d["exponents"]) for d in table()["dimensions"]}
    return dimensions[unit["dimension"]]


def exponents_of(dimensionality: Any) -> tuple[int, ...]:
    """`pint`'s dimensionality as exponents in the table's slot order.

    A dimension `pint` reports that the table has no slot for is an error rather
    than a zero: a unit whose dimensionality includes a base dimension this table
    does not carry is a unit this vocabulary cannot express, and reading its absent
    slot as zero would make it look like one it can.
    """
    known = {PINT_BASE[slot] for slot in SLOTS}
    reported = {str(part) for part in dimensionality}
    unknown = reported - known
    assert not unknown, (
        f"pint reports base dimension(s) {sorted(unknown)}, which the vocabulary table "
        f"has no slot for. Slots: {SLOTS}."
    )
    return tuple(int(dimensionality[PINT_BASE[slot]]) for slot in SLOTS)


def test_every_slot_is_a_pint_base_dimension() -> None:
    """Each slot is a `pint` base dimension, and not a combination of the others.

    Two things are asserted per slot, and they are different. The representative
    unit's *dimensionality* is the slot's own name - so if `pint` renamed
    `[substance]`, this fails rather than the map quietly matching nothing. And
    that unit is its own base unit, so `to_base_units` leaves it alone - a derived
    dimension would reduce, `[energy]` to `[length]**2*[mass]/[time]**2`, and the
    seven exponents in the table would stop being independent.

    The companion guard, against `pint` *adding* a base dimension this vocabulary
    cannot express, is in :func:`exponents_of`: it refuses a dimensionality
    containing a base dimension the table has no slot for, rather than reading the
    absent slot as a zero exponent.
    """
    assert set(SLOTS) == set(PINT_BASE), (
        f"the table's slots are {SLOTS} and this file maps {sorted(PINT_BASE)}; the "
        f"generated table and this map have drifted"
    )
    assert len(set(PINT_BASE.values())) == len(PINT_BASE), "two slots map to one pint dimension"
    for slot, representative in SLOT_UNITS.items():
        quantity = ureg.Quantity(1.0, representative)
        assert str(quantity.dimensionality) == PINT_BASE[slot], (
            f"slot {slot}: {representative!r} has dimensionality "
            f"{quantity.dimensionality}, not {PINT_BASE[slot]}"
        )
        assert quantity.to_base_units().units == quantity.units, (
            f"slot {slot}: {representative!r} is not its own base unit - it reduces to "
            f"{quantity.to_base_units().units}, so it is a derived dimension"
        )


def test_every_unit_has_the_dimension_the_table_declares() -> None:
    """`pint`'s dimensionality for a unit equals the exponents the table gives it.

    This is the only place a dimension is checked against a units library, because
    `pint` is the only one of the two that can report one.
    """
    for unit in units():
        declared = _exponents(unit)
        reported = exponents_of(ureg.Quantity(1.0, unit["pint"]).dimensionality)
        assert reported == declared, (
            f"{unit['id']}: the table declares {declared} but pint reports {reported} "
            f"for {unit['pint']!r}"
        )


@pytest.mark.requires_rust
def test_the_two_units_libraries_agree_on_every_factor() -> None:
    """`pint`'s factor and the conversion a calculation runs are the same number.

    Neither side is a literal. The left is `pint` reading its own definition; the
    right is `azoth._core.unit_si_factor`, which runs the generated conversion
    table - the same function `to_si`'s Rust counterpart reaches. `mm` is the entry
    that would break first, by 1000x, which is what happened when the factor was
    written down instead of asked for.
    """
    factors = _extension().unit_si_factor
    for unit in units():
        expected = float(ureg.Quantity(1.0, unit["pint"]).to_base_units().magnitude)
        got = factors(unit["id"])
        assert got is not None, f"{unit['id']}: the Rust vocabulary has no conversion for it"
        assert got == pytest.approx(expected, rel=1.0e-15), (
            f"{unit['id']}: pint says one of it is {expected} in SI base, and the Rust "
            f"conversion this library actually runs gives {got}"
        )


@pytest.mark.requires_rust
def test_the_slot_order_is_the_same_in_both_languages() -> None:
    """The Rust and Python sides write exponent tuples in one order.

    The table states the order once, and both generated files are built from it - so
    this is a check on the generator, not on a person. It matters because the order
    is unobservable from either side alone: a permuted tuple in Rust and a matching
    permutation in Python would compare equal to each other and to nothing the table
    says, and every unit whose exponents are not symmetric would be silently wrong.
    """
    assert tuple(_extension().unit_slots()) == SLOTS


@pytest.mark.requires_rust
def test_the_rust_dimension_agrees_with_the_table() -> None:
    """The generated Rust tuple for each unit is the table's, in the table's order.

    The Rust assertion that catches a wrong exponent is a compile error, and it
    needs the `uom` quantity to differ - which it cannot when two dimensions share
    a quantity type. This compares all seven exponents for every unit, so an
    ordering mistake that the compiler cannot see fails here.
    """
    declared = {unit["id"]: _exponents(unit) for unit in units()}
    reported = _extension().unit_dimensions
    for name, exponents in declared.items():
        rust = reported(name)
        assert rust is not None, f"{name}: the Rust vocabulary has no dimension for it"
        assert tuple(rust) == exponents, (
            f"{name}: the table says {exponents} and the generated Rust says {tuple(rust)}"
        )
