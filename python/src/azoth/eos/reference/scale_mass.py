"""``eos.scale_mass`` - NeqSim's ``ScaleMassCalculator``.

The class's approximation: the limiting ion's concentration times
``min(1 - 10^-SI, 0.99)``, one arm per mineral and differing only in molar mass.

Spec: ``specs/calcs/eos/scale_mass.toml``.
"""

from __future__ import annotations

import math
from typing import Any

from azoth.core.errors import InvalidInputError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import ScaleMassResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning

CALC_ID = "eos.scale_mass"

#: The five minerals' molar masses, in g/mol, as the class writes them.
MOLAR_MASSES = {
    "caco3": 100.09,
    "baso4": 233.39,
    "srso4": 183.68,
    "caso4": 136.14,
    "feco3": 115.85,
}


def scale_mass(
    mineral: str,
    cation_concentration: Q,
    anion_concentration: Q,
    saturation_index: float,
) -> ScaleMassResult:
    """The mass a supersaturated water precipitates, by the class's approximation.

    Args:
        mineral: one of ``caco3``, ``baso4``, ``srso4``, ``caso4``, ``feco3``.
        cation_concentration: the cation's concentration.
        anion_concentration: the anion's concentration.
        saturation_index: the mineral's SI at the state.

    Returns:
        The mass, the moles and the fraction the class uses.

    Raises:
        InvalidInputError: for a mineral outside the five.
        OutOfRangeError: if a concentration is negative.

    Example:
        >>> import azoth
        >>> q = azoth.ureg.Quantity
        >>> r = scale_mass("caco3", q(10.0, "mol/m**3"), q(20.0, "mol/m**3"), 2.5)
        >>> abs(r.precipitated_mass.magnitude - 0.990891) < 1e-12
        True
    """
    from azoth._registry_gen import spec as _spec_for

    spec: dict[str, Any] = _spec_for(CALC_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    cation = input_to_si(spec, "cation_concentration", cation_concentration)
    anion = input_to_si(spec, "anion_concentration", anion_concentration)
    apply_checks(
        checks.on_input,
        {"cation_concentration": cation, "anion_concentration": anion}.get,
        warnings,
    )

    if mineral not in MOLAR_MASSES:
        raise InvalidInputError("mineral", f"`{mineral}` is not one of {sorted(MOLAR_MASSES)}")
    molar_mass = MOLAR_MASSES[mineral]

    # The class's own units: mol/L and g/mol.
    c_cation = cation / 1000.0
    c_anion = anion / 1000.0
    limiting = min(c_cation, c_anion)
    if saturation_index <= 0.0:
        fraction, moles = 0.0, 0.0
    else:
        fraction = min(1.0 - math.pow(10.0, -saturation_index), 0.99)
        moles = limiting * fraction
    mass_mg_per_l = moles * molar_mass * 1000.0

    return ScaleMassResult(
        precipitated_mass=from_si(mass_mg_per_l * 1.0e-3, "kg/m**3"),
        precipitated_moles=from_si(moles * 1000.0, "mol/m**3"),
        fraction_excess=fraction,
        limiting_concentration=from_si(limiting * 1000.0, "mol/m**3"),
        warnings=tuple(warnings),
    )
