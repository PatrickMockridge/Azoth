"""Units, and the rule that governs how they are used.

# The boundary rule

Dimensioned values cross the public API as `pint` quantities and are converted
to SI on the way in. Inside a calculation everything is a plain `float`.

That mirrors the Rust side, where `uom` quantities are used at the boundary and
`.value` is used inside, and it is what makes the two implementations comparable:
both do the same arithmetic on the same numbers rather than each trusting its own
units library to reach the same result by a different route.

The cost is that dimensional correctness is enforced on the way in and out but
not step by step through an equation. The benefit is that the bodies of the
calculations match the published equations line by line, which is the property a
reviewer can actually check.

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

from typing import Any, Final

import pint

from azoth.core.errors import UnitMismatchError

#: The one registry. A second one would have its own unit cache and its own
#: definition overrides, which is a subtle way for two parts of a program to
#: disagree about what a unit means.
ureg: Final[pint.UnitRegistry] = pint.UnitRegistry()

#: Quantity alias. Always annotate `Q`, never a bare `Quantity`.
type Q = pint.Quantity[float]

#: Canonical unit strings, keyed by the strings the spec schema allows. The
#: schema restricts `unit` to this set, so a spec cannot name a unit the code has
#: no conversion path for.
#:
#: # The SI-base rule
#:
#: `to_si` returns the magnitude expressed in the unit named here, and that number
#: is what crosses into a calculation. The Rust side's `uom` quantities always
#: report `.value` as the SI *base* magnitude. For those two to be the same number
#: - which is what the cross-language agreement test rests on - every entry below
#: has to be the SI base unit, so that "the unit named here" and "the SI base unit"
#: are the same thing and no conversion factor can hide between the two languages.
#:
#: `mm` is the one entry that breaks that rule, and it is a known defect rather
#: than a decision. A spec declaring `mm` would have this side work in millimetres
#: while the Rust side works in metres, and the two would disagree by 1000x. It is
#: recorded rather than quietly fixed or quietly deleted because pipe diameters are
#: conventionally quoted in millimetres, so this unit is likely to be wanted - and
#: it must not be reached for before the convention is made coherent.
#: `test_units_contract.py` pins the exception so that a second unit added with the
#: same problem fails rather than joining it.
CANONICAL_UNITS: Final[dict[str, str]] = {
    "dimensionless": "dimensionless",
    "m": "meter",
    # NOT SI base - see the note above. The only such entry, and unused by every
    # spec in the registry, which is the only reason it has never mattered.
    "mm": "millimeter",
    "m**2": "meter**2",
    "m**3/s": "meter**3/second",
    "kg/s": "kilogram/second",
    "kg/m**3": "kilogram/meter**3",
    "m/s": "meter/second",
    "Pa": "pascal",
    "Pa*s": "pascal*second",
    "K": "kelvin",
    "W": "watt",
    "J/(kg*K)": "joule/(kilogram*kelvin)",
    "W/(m*K)": "watt/(meter*kelvin)",
    "W/(m**2*K)": "watt/(meter**2*kelvin)",
    "kg/mol": "kilogram/mole",
}


def unit_for(spec_unit: str) -> str:
    """Map a spec unit string to the pint unit name used to express it."""
    try:
        return CANONICAL_UNITS[spec_unit]
    except KeyError:
        raise UnitMismatchError(
            "unit", spec_unit, f"unknown spec unit; known: {sorted(CANONICAL_UNITS)}"
        ) from None


def quantity(magnitude: float, spec_unit: str) -> Q:
    """Build a quantity in a spec's canonical unit."""
    return ureg.Quantity(magnitude, unit_for(spec_unit))


def to_si(value: Any, spec_unit: str, field: str) -> float:
    """Convert a dimensioned input to its SI magnitude in the spec's unit.

    Returns a plain float, which is what the calculations work in.

    Raises:
        UnitMismatchError: if `value` is not a quantity, or carries the wrong
            dimensions. Both are caller errors: a bare number where a length is
            expected is exactly the mistake this library exists to make
            impossible, so it is rejected loudly rather than assumed to be SI.
    """
    if not isinstance(value, pint.Quantity):
        raise UnitMismatchError(field, spec_unit, f"{type(value).__name__} {value!r}")
    target = unit_for(spec_unit)
    try:
        converted: Q = value.to(target)
    except pint.errors.DimensionalityError as exc:
        raise UnitMismatchError(field, spec_unit, str(value.units)) from exc
    return float(converted.magnitude)


def same_physical_value(a: Any, b: Any) -> bool:
    """Whether two quantities describe the same physical value.

    Converts both to base units before comparing. Comparing pint quantities
    directly mostly works, but does not raise on incompatible dimensions - it
    returns False - so a test asserting equality would pass for the wrong reason
    when both sides are wrong in the same way.
    """
    if not isinstance(a, pint.Quantity) or not isinstance(b, pint.Quantity):
        return bool(a == b)
    if a.dimensionality != b.dimensionality:
        return False
    return bool(abs(a.to_base_units().magnitude - b.to_base_units().magnitude) <= 0.0)
