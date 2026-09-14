"""Units, and the rule that governs how they are used.

# The boundary rule

Dimensioned values cross the public API as `pint` quantities. Inside a
calculation everything is a plain `float`, and that float is always the **SI base
magnitude**: millimetres become metres, grams become kilograms, and degrees
Celsius become kelvin, before the number reaches any arithmetic.

That mirrors the Rust side, where `uom` quantities are used at the boundary and
`.value` is used inside. `uom`'s `.value` is always the SI base magnitude, so
converting to SI base here is what makes the two implementations do the same
arithmetic on the same numbers, rather than each trusting its own units library
to arrive there by a different route.

The cost is that dimensional correctness is enforced on the way in and out but
not step by step through an equation. The benefit is that the bodies of the
calculations match the published equations line by line, which is the property a
reviewer can actually check.

# The three functions, which are not interchangeable

Confusing these is easy, and the failure is silent - a wrong number that looks
entirely reasonable - so they have distinct names and the distinction is pinned
by `test_units_conversion.py`:

* :func:`to_si` - **in**: a caller's quantity. **out**: the SI base magnitude as a
  `float`. This is the public-to-internal direction, and every dimensioned input
  of every calc goes through it.

* :func:`from_si` - **in**: an SI base magnitude. **out**: a quantity expressed in
  the spec's own unit. This is the internal-to-public direction, used wherever an
  internal number is handed back to a caller and by the Rust bridge, whose
  transport objects carry an SI magnitude alongside a unit name.

* :func:`quantity` - builds a quantity whose magnitude is *already stated in* the
  named unit. This is what spec files and test data hold: a worked example that
  says `f_t: 0.018` or `D: 100` means 100 of the declared unit, not 100 SI base.

For a unit that is its own SI base unit the three agree, and the distinction looks
academic. For `mm` they do not: 500 mm is `to_si` 0.5 and `quantity` 500.0. That
difference is the whole reason this module exists as its own thing.

# Why the conversion does not depend on the vocabulary

`to_si` converts to base units rather than to whatever `CANONICAL_UNITS` names, so
it is correct for any unit the schema permits - including one whose canonical form
is not its SI base unit, which `mm` is. Conversion is to the SI base unit rather than
to the named one, so correctness does not rest on every entry being SI-coherent -
`mm` is in the vocabulary and would otherwise make the two implementations disagree
by 1000x. A test asserts the round trip for every entry.

# Dimensionless quantities

Reynolds number, relative roughness, friction factor and resistance coefficient
are plain `float`, here and in Rust. They carry no unit to be safe about, and a
newtype around a ratio would be an abstraction with exactly one implementation.

# On pint and mypy --strict

`pint`'s `Quantity` is generic, so a bare `Quantity` annotation is a
``[type-arg]`` error under strict mode. `Q` below is the single alias used
throughout the package, which confines the problem to this module.
"""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any, Final

import pint

from azoth.core._units_gen import CANONICAL_UNITS as CANONICAL_UNITS
from azoth.core.errors import UnitMismatchError

#: The one registry. A second one would have its own unit cache and its own
#: definition overrides, which is a subtle way for two parts of a program to
#: disagree about what a unit means.
ureg: Final[pint.UnitRegistry] = pint.UnitRegistry()

#: Quantity alias. Always annotate `Q`, never a bare `Quantity`.
type Q = pint.Quantity[float]

#: Canonical unit strings, keyed by the strings the spec schema allows, mapped to
#: the unit's name in `pint`'s registry. **Generated** from
#: `specs/vocabulary/vocabulary.toml` by `tools/gen_vocabulary.py`, which is why
#: the table rather than this module is where a unit is added.
#:
#: This is the unit the spec and the generated docs *speak in* - what a worked
#: example's numbers mean, and what a caller gets back. It is deliberately not the
#: unit anything is calculated in: calculations work in SI base magnitudes, and
#: `to_si`/`from_si` are the only places the two representations meet.
#:
#: Entries are SI-first where the unit is a free choice, so that the vocabulary and
#: the internal representation agree for most quantities. But that is a convention
#: about how these are written down, not a load-bearing property: `mm` is not SI
#: base, and the conversion handles it correctly because it converts to base units
#: rather than to this name. A spec may declare any entry here without risking a
#: silent factor.


def unit_for(spec_unit: str) -> str:
    """Map a spec unit string to the pint unit name used to express it."""
    try:
        return CANONICAL_UNITS[spec_unit]
    except KeyError:
        raise UnitMismatchError(
            "unit", spec_unit, f"unknown spec unit; known: {sorted(CANONICAL_UNITS)}"
        ) from None


def quantity(magnitude: float, spec_unit: str) -> Q:
    """A quantity whose magnitude is stated in the spec's own unit.

    This is the constructor for values as they appear in a spec, a worked example
    or a test case. It is *not* the inverse of :func:`to_si`.

    ``quantity(100.0, "mm")`` is 100 millimetres, while ``to_si`` of that same
    quantity is ``0.1``. Use :func:`from_si` to rebuild a quantity from a number a
    calculation produced.
    """
    return ureg.Quantity(magnitude, unit_for(spec_unit))


def has_offset(unit: Any) -> bool:
    """Whether a unit carries an additive offset, like ``degC`` or ``degF``.

    Read out of `pint` rather than listed, so a unit this library has never heard of
    is classified correctly. The test is "is zero of this unit zero of its base
    unit": a `degC` is 273.15 K at zero, a `delta_degC` is 0 K, and a `degF` is
    255.37 K. Scaling alone does not make a unit offset - `mm` and `delta_degF` both
    map zero to zero - which is why `one != 1` would be the wrong test: it would
    call `delta_degF` offset as well, and that is the unit a caller is *supposed* to
    reach for.
    """
    return float(ureg.Quantity(0.0, unit).to_base_units().magnitude) != 0.0


def to_si(value: Any, spec_unit: str, field: str, *, interval: bool = False) -> float:
    """Convert a caller's quantity to the SI base magnitude the calcs work in.

    The input is first converted to the spec's declared unit, which is what checks
    the dimensions, and then to base units, which is what the calculation needs.
    Both steps matter: the first is why a bare number where a length is expected is
    rejected rather than silently read as metres, and the second is why a spec may
    declare `mm` without the two implementations disagreeing by 1000x.

    Args:
        value: the caller's quantity.
        spec_unit: the unit the spec declares, which fixes the dimension.
        field: the input's name, for the error message.
        interval: set for an input the spec marks ``interval: true`` - a temperature
            *difference* rather than an absolute temperature. Both are
            kelvin-dimensioned, so nothing about the dimensions can tell them
            apart, and `pint` will happily convert an absolute ``Q(30, "degC")`` to
            303.15 K for an input that means "a 30 kelvin difference". That is a
            plausible number wrong by 273.15, which is the failure this library
            exists to make impossible, so an offset unit is **refused** instead.

    Raises:
        UnitMismatchError: if `value` is not a quantity, carries the wrong
            dimensions, or is an offset-unit quantity where `interval` says a
            difference is meant. All three are caller errors: a bare number where a
            length is expected is exactly the mistake this library exists to make
            impossible, so it is rejected loudly rather than assumed to be SI.
    """
    if not isinstance(value, pint.Quantity):
        raise UnitMismatchError(field, spec_unit, f"{type(value).__name__} {value!r}")
    target = unit_for(spec_unit)
    if interval and has_offset(value.units):
        # Deliberately not converted. A caller with an absolute temperature who
        # wants a difference should write the conversion, `q.to("delta_degC")` or
        # `q.to("K")`, so the 273.15 lands in their code where it is visible rather
        # than inside this function where it is not.
        raise UnitMismatchError(
            field,
            spec_unit,
            f"{value.units} (an absolute temperature; this input is a temperature "
            f"difference, so pass `delta_degC`/`delta_degF`, or convert explicitly "
            f"with `.to('K')` first)",
        )
    try:
        converted: Q = value.to(target)
    except pint.errors.DimensionalityError as exc:
        raise UnitMismatchError(field, spec_unit, str(value.units)) from exc
    base: Q = converted.to_base_units()
    return float(base.magnitude)


def input_to_si(spec: Mapping[str, Any], field: str, value: Any) -> float:
    """Convert one declared input, reading its unit and flags from the spec.

    The point is the `interval` flag and not the unit string. A call site that
    writes `to_si(dT, "K", "dT")` is restating what the spec already declares, and a
    flag the spec holds but the call site ignores is a check that exists on paper
    and not in the code - exactly the failure this project is organised against. The
    unit is read for the same reason: one fewer literal that can drift.

    Reads the spec's own `inputs` block rather than a generated copy of it, so a
    spec edited without regenerating cannot leave this function agreeing with a
    stale table.

    Raises:
        UnitMismatchError: as :func:`to_si`.
        KeyError: if the spec does not declare `field`. A programming error in the
            implementation, not a caller condition - the whole point of the spec is
            that the input names are the contract.
    """
    declaration = spec["inputs"][field]
    return to_si(
        value,
        declaration["unit"],
        field,
        interval=bool(declaration.get("interval", False)),
    )


def from_si(si_magnitude: float, spec_unit: str) -> Q:
    """Rebuild a quantity in the spec's own unit from an SI base magnitude.

    The inverse of :func:`to_si`, and the direction a result travels. A calc
    computes in SI base magnitudes, and the Rust bridge receives them as
    ``magnitude_si`` alongside the unit name to present them in; both go through
    here so that the number a caller sees is expressed the way the spec describes
    it.

    Round-trips with :func:`to_si` for every unit in the vocabulary, including the
    ones that are not their own SI base unit.
    """
    declared = unit_for(spec_unit)
    # One unit of the declared unit, as an SI base quantity. Its units are the ones
    # the internal magnitude is in, so that is the unit to read the number as.
    reference: Q = ureg.Quantity(1.0, declared).to_base_units()
    rebuilt: Q = ureg.Quantity(si_magnitude, reference.units).to(declared)
    return rebuilt
