"""Unit conversion: does the number a calculation receives mean what it says.

``test_units_contract.py`` checks that the three vocabularies agree on *which*
units exist. This file checks what the conversion does with them, which is a
different question and the one where a silent factor hides.

The convention is that a calculation works in SI base magnitudes. That is what the
Rust side's ``uom`` ``.value`` always is, so it is what makes the two
implementations run the same arithmetic on the same numbers. Everything here is
either an assertion that the convention holds or a demonstration of what goes wrong
when it is confused with a neighbouring function.

The expectations are written out rather than derived by calling the code under
test. ``to_si`` could be compared against ``pint``'s own ``to_base_units`` and
would pass by construction, which proves nothing; the table below states what the
answer physically is, and ``mm`` is in it precisely because it is the entry where
"the number in the spec" and "the number the calculation uses" differ.
"""

from __future__ import annotations

import pytest

import _helpers as h
from azoth.core.errors import UnitMismatchError
from azoth.core.units import CANONICAL_UNITS, from_si, quantity, to_si, unit_for, ureg

#: The SI base magnitude of one of each unit, where it is not simply 1.0.
#:
#: Every vocabulary entry is its own SI base unit except ``mm``, so the default is
#: 1.0 and this table carries only the exceptions. Having them spelled out means
#: the assertion below is a statement about the world rather than a restatement of
#: the implementation - and it is what would catch a change of convention, or a new
#: entry whose base magnitude nobody thought about.
NOT_ITS_OWN_SI_BASE_UNIT: dict[str, float] = {"mm": 1.0e-3}


@pytest.mark.parametrize("spec_unit", sorted(CANONICAL_UNITS))
def test_to_si_yields_the_si_base_magnitude(spec_unit: str) -> None:
    """One declared unit converts to the advertised SI base magnitude.

    This is the invariant the whole boundary rule rests on. For ``mm`` it means a
    calculation handed 500 mm receives 0.5, not 500 - the same number the Rust side
    would hold, whose ``.value`` is always metres.
    """
    expected = NOT_ITS_OWN_SI_BASE_UNIT.get(spec_unit, 1.0)
    assert to_si(quantity(1.0, spec_unit), spec_unit, "x") == pytest.approx(expected)


def test_mm_is_where_the_spec_number_and_the_internal_number_differ() -> None:
    """The discriminating case, kept as executable documentation.

    ``quantity`` and ``to_si`` both take a number and a unit string, and for every
    unit but one they agree. Here they differ by 1000x, which is exactly the
    confusion that makes a unit bug silent: 500 mm handed to a calculation that
    expected millimetres, in a codebase where everything else is metres, is a
    plausible number that is wrong by a factor of a thousand.

    Nothing in the library is allowed to use the wrong one, and no existing calc
    takes an ``mm`` input, so this is the only thing standing between the
    distinction and a future caller getting it backwards.
    """
    five_hundred_mm = ureg.Quantity(500.0, "millimeter")
    assert to_si(five_hundred_mm, "mm", "D") == pytest.approx(0.5)
    assert quantity(500.0, "mm").magnitude == pytest.approx(500.0)
    assert to_si(five_hundred_mm, "mm", "D") != pytest.approx(quantity(500.0, "mm").magnitude)


@pytest.mark.parametrize("spec_unit", sorted(CANONICAL_UNITS))
def test_the_round_trip_is_the_identity(spec_unit: str) -> None:
    """``from_si`` undoes ``to_si``, for every unit in the vocabulary.

    A result travels out the way an input came in: the calc produces an SI base
    magnitude, and the caller is handed a quantity in the unit the spec declares.
    If that path is not the exact inverse of the inbound one, a calc returns a
    number that does not mean what its own documentation says.
    """
    original = quantity(7.5, spec_unit)
    round_tripped = from_si(to_si(original, spec_unit, "x"), spec_unit)
    # Compared as floats in base units: a value comparison between quantities is
    # the physically right thing, but `pytest.approx` iterates its argument and
    # pint's `Quantity.__iter__` only accepts array-backed magnitudes.
    assert float(round_tripped.to_base_units().magnitude) == pytest.approx(
        float(original.to_base_units().magnitude), rel=1e-12
    )


def test_from_si_expresses_the_result_in_the_declared_unit() -> None:
    """The rebuilt quantity is in the spec's unit, not in SI base.

    Round-tripping values is not enough on its own: 0.5 metres and 500 millimetres
    are the same physical quantity, so a value comparison passes even when the unit
    label is wrong. The label is what a caller reads, so it is asserted too.
    """
    rebuilt = from_si(0.5, "mm")
    assert rebuilt.units == ureg.Unit("millimeter")
    assert float(rebuilt.to("millimeter").magnitude) == pytest.approx(500.0)
    assert float(rebuilt.to_base_units().magnitude) == pytest.approx(0.5)


def test_a_bare_number_is_rejected_rather_than_assumed_to_be_si() -> None:
    """The failure this library exists to make impossible.

    A float where a length is expected is the mistake a units library is for. It is
    rejected loudly instead of being read as metres, because a silent
    thousand-fold error looks entirely reasonable in the output.
    """
    with pytest.raises(UnitMismatchError, match="must be a quantity"):
        to_si(100.0, "m", "L")


def test_wrong_dimensions_are_rejected() -> None:
    """A length where a pressure belongs is an error, not a conversion."""
    with pytest.raises(UnitMismatchError):
        to_si(ureg.Quantity(100.0, "meter"), "Pa", "dp")


def test_an_unknown_unit_is_rejected() -> None:
    """A unit outside the vocabulary cannot be converted, and says so."""
    with pytest.raises(UnitMismatchError, match="unknown spec unit"):
        unit_for("furlong")


def test_incoming_values_are_converted_not_merely_accepted() -> None:
    """A caller may work in any unit convertible to the spec's.

    The schema says so explicitly - "API callers may pass any convertible unit" -
    and it is what keeps the spec's own units a documentation choice rather than a
    constraint on whoever calls it.
    """
    # 1 bar in pascals.
    assert to_si(ureg.Quantity(1.0, "bar"), "Pa", "p") == pytest.approx(1.0e5)
    # 1 foot in metres, exactly.
    assert to_si(ureg.Quantity(1.0, "foot"), "m", "L") == pytest.approx(0.3048)
    # 1 US gallon per minute in cubic metres per second.
    assert to_si(ureg.Quantity(1.0, "gallon/minute"), "m**3/s", "Q") == pytest.approx(
        6.30901964e-5, rel=1e-6
    )


def test_dimensionless_values_pass_through_unchanged() -> None:
    """A dimensionless quantity is its own base unit, and stays a float."""
    stated = to_si(ureg.Quantity(0.023, "dimensionless"), "dimensionless", "f")
    bare = to_si(ureg.Quantity(0.023), "dimensionless", "f")
    assert stated == pytest.approx(0.023)
    assert bare == pytest.approx(0.023)


def test_an_offset_unit_is_refused_where_a_difference_is_meant() -> None:
    """`interval=True` refuses `degC`/`degF`, and only those.

    A temperature difference and an absolute temperature are the same dimension, so
    nothing about `pint`'s conversion can tell them apart: it will convert an
    absolute ``Q(30, "degC")`` to 303.15 K for an input that means "a 30 kelvin
    difference", and the caller gets a plausible number wrong by 273.15.

    The unit is not guessed at - it is read the way `has_offset` reads it, by asking
    whether zero of it is zero of its base. That is what makes `delta_degF` pass
    while `degF` is refused: both scale, only one offsets.
    """
    from azoth.core.units import has_offset

    assert has_offset("degC") and has_offset("degF")
    assert not has_offset("delta_degC")
    assert not has_offset("delta_degF")
    # A scaling unit is not an offset unit, which is why "one of this unit is not
    # one kelvin" would have been the wrong test: it would reject delta_degF too.
    assert not has_offset("mm")
    assert not has_offset("K")

    for absolute in ("degC", "degF"):
        with pytest.raises(UnitMismatchError) as excinfo:
            to_si(ureg.Quantity(30.0, absolute), "K", "dT", interval=True)
        assert excinfo.value.field() == "dT"
        assert "delta_degC" in str(excinfo.value)

    # A difference in any of its spellings converts normally.
    assert to_si(ureg.Quantity(30.0, "delta_degC"), "K", "dT", interval=True) == 30.0
    assert to_si(ureg.Quantity(30.0, "K"), "K", "dT", interval=True) == 30.0
    assert to_si(ureg.Quantity(54.0, "delta_degF"), "K", "dT", interval=True) == pytest.approx(30.0)

    # And without the flag, an absolute temperature still converts - the fluid
    # property tables need exactly that.
    assert to_si(ureg.Quantity(20.0, "degC"), "K", "temperature") == pytest.approx(293.15)


def test_every_interval_input_in_the_registry_is_honoured_on_both_backends() -> None:
    """The flag is only worth having if the calcs carrying it act on it.

    Walks the registry rather than naming `conduction_plane_wall`, so a calc that
    declares ``interval: true`` and never passes the flag through is caught here
    instead of at a caller's expense. The failure this guards is a spec that says
    one thing and an implementation that does another, which no amount of reading a
    single calc's test would reveal.

    Both backends, because the refusal happens at the Python boundary: the Rust
    binding converts inputs too, and a fix applied only to the reference would leave
    the default backend silently offsetting.
    """
    from azoth._dispatch import available, resolve, use_backend
    from azoth._registry_gen import CALCS

    marked = [
        (calc, name)
        for calc in CALCS
        for name, declaration in calc["inputs"].items()
        if declaration.get("interval", False)
    ]
    assert marked, "no input declares interval: true; this test would pass vacuously"

    for calc, name in marked:
        case = next(
            (c for c in h.all_tests(calc) if c["status"] == "active" and "inputs" in c),
            None,
        )
        assert case is not None, f"{calc['id']} has no runnable case to build from"
        kwargs = h.kwargs_for(calc, case["inputs"])
        absolute = ureg.Quantity(float(case["inputs"][name]), "degC")
        for backend in sorted(available()):
            with use_backend(backend), pytest.raises(UnitMismatchError) as excinfo:
                resolve(calc["id"])(**{**kwargs, name: absolute})
            assert excinfo.value.field() == name, f"{calc['id']} on {backend}"
