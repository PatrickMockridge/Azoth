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

from chemeng.core.errors import UnitMismatchError

#: The one registry. A second one would have its own unit cache and its own
#: definition overrides, which is a subtle way for two parts of a program to
#: disagree about what a unit means.
ureg: Final[pint.UnitRegistry] = pint.UnitRegistry()

#: Quantity alias. Always annotate `Q`, never a bare `Quantity`.
type Q = pint.Quantity[float]

#: Canonical unit strings, keyed by the strings the spec schema allows. The
#: schema restricts `unit` to this set, so a spec cannot name a unit the code has
#: no conversion path for.
CANONICAL_UNITS: Final[dict[str, str]] = {
    "dimensionless": "dimensionless",
    "m": "meter",
    "mm": "millimeter",
    "kg/m**3": "kilogram/meter**3",
    "m/s": "meter/second",
    "Pa": "pascal",
    "Pa*s": "pascal*second",
    "K": "kelvin",
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
