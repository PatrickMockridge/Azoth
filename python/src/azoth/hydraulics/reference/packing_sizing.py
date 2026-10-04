"""``hydraulics.packing_sizing`` - the internal diameter a packed bed needs at a chosen fraction
of flood.

```text
u_design = u_flood * f
q_v      = m_v / max(rho_v, 0.01)
A        = q_v / u_design
d_raw    = sqrt(4 A / pi)
d        = the first standard vessel size >= d_raw
```

Spec: ``specs/calcs/hydraulics/packing_sizing.toml``

# The trial diameter is inert, and that is what makes this a closed form

``PackingHydraulicsCalculator.sizeColumnDiameter`` sets ``columnDiameter = 1.0`` and then calls
``calculateFloodingVelocity()`` - which reads **no diameter and no area**, only the packing's
factor, the flows, the densities and the liquid's viscosity. So the trial is set and never used,
and sizing is not the iteration the name ``trial`` suggests: one fit evaluation, one area, one
rounding.

# The rounding is up, and its floor is reached

``roundToStandardDiameter`` returns the first size in a fixed table of thirty-one that is
``>= d_raw``, and only above ``8.0`` m falls back to half-metre steps. Nothing smaller than
``0.3`` m is ever returned, so a requirement of ``0.1195`` m and one of ``0.2973`` m both answer
``0.3`` m - which is why the raw diameter is reported beside the rounded one rather than instead
of it.
"""

from __future__ import annotations

import math
from typing import Final

from azoth._registry_gen import spec as _spec_for
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import PackingSizingResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning
from azoth.hydraulics.packing import packing_or_default
from azoth.hydraulics.reference.packing_hydraulics import eckert_flooding_velocity

CALC_ID = "hydraulics.packing_sizing"

#: The standard vessel sizes ``roundToStandardDiameter`` chooses from, in metres.
#:
#: Thirty-one of them, in hundredths to ``1.2``, then in fifths and tenths to ``4.0``, then in
#: halves to ``8.0``. Above the table the class falls back to ``ceil(2 d)/2``, which is the same
#: half-metre step continuing - so the table is a list of what is *bought* rather than a
#: discretisation of what is possible.
STANDARD_DIAMETERS_M: Final[tuple[float, ...]] = (
    0.3,
    0.4,
    0.5,
    0.6,
    0.7,
    0.8,
    0.9,
    1.0,
    1.1,
    1.2,
    1.4,
    1.5,
    1.6,
    1.8,
    2.0,
    2.2,
    2.4,
    2.6,
    2.8,
    3.0,
    3.2,
    3.4,
    3.6,
    3.8,
    4.0,
    4.5,
    5.0,
    5.5,
    6.0,
    7.0,
    8.0,
)

#: The diameter ``sizeColumnDiameter`` starts from, and the one it answers with when the design
#: velocity is not positive.
TRIAL_DIAMETER_M: Final[float] = 1.0

#: The floor the vapour density is held at in the volumetric flow, in kg/m**3.
MIN_VAPOR_DENSITY: Final[float] = 0.01


def round_to_standard_diameter(diameter: float) -> float:
    """``roundToStandardDiameter``: the first standard size at or above ``diameter``."""
    for size in STANDARD_DIAMETERS_M:
        if size >= diameter:
            return size
    return math.ceil(diameter * 2.0) / 2.0


def packing_sizing(
    packing: str,
    design_flood_fraction: float,
    vapor_mass_flow: Q,
    liquid_mass_flow: Q,
    vapor_density: Q,
    liquid_density: Q,
    liquid_viscosity: Q,
    hydraulic_capacity_factor: float,
) -> PackingSizingResult:
    """The internal diameter a packed bed needs at a chosen fraction of flood.

    Args:
        packing: the packing's name or an alias; a name nothing matches is ``Pall-Ring-50``.
            Only its factor is read here.
        design_flood_fraction: the fraction of the flooding velocity to design at.
        vapor_mass_flow: the vapour's mass flow at the sizing state.
        liquid_mass_flow: the liquid's mass flow, which enters only through the flow parameter.
        vapor_density: the vapour's mass density.
        liquid_density: the liquid's mass density.
        liquid_viscosity: the liquid's dynamic viscosity.
        hydraulic_capacity_factor: the packing's relative hydraulic capacity.

    Returns:
        The flooding and design velocities, the vapour's volumetric flow, the area and raw
        diameter the vapour needs, and the standard vessel size that meets it.

    Raises:
        OutOfRangeError: if the design flood fraction, either density or the liquid viscosity is
            not positive.

    Example:
        >>> import azoth
        >>> q = azoth.ureg.Quantity
        >>> r = packing_sizing(
        ...     "Pall-Ring-50", 0.7, q(0.35, "kg/s"), q(3.5, "kg/s"),
        ...     q(45.0, "kg/m**3"), q(990.0, "kg/m**3"), q(6.5e-4, "Pa*s"), 1.0,
        ... )
        >>> round(r.required_diameter.magnitude, 6)
        0.371718
        >>> round(r.column_diameter.magnitude, 6)
        0.4
    """
    spec = _spec_for(CALC_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    values = {
        "vapor_mass_flow": input_to_si(spec, "vapor_mass_flow", vapor_mass_flow),
        "liquid_mass_flow": input_to_si(spec, "liquid_mass_flow", liquid_mass_flow),
        "vapor_density": input_to_si(spec, "vapor_density", vapor_density),
        "liquid_density": input_to_si(spec, "liquid_density", liquid_density),
        "liquid_viscosity": input_to_si(spec, "liquid_viscosity", liquid_viscosity),
    }
    apply_checks(
        checks.on_input,
        {
            **values,
            "design_flood_fraction": design_flood_fraction,
            "hydraulic_capacity_factor": hydraulic_capacity_factor,
        }.get,
        warnings,
    )

    resolved = packing_or_default(packing)

    flooding_velocity = eckert_flooding_velocity(
        resolved.packing_factor,
        hydraulic_capacity_factor,
        values["vapor_mass_flow"],
        values["liquid_mass_flow"],
        values["vapor_density"],
        values["liquid_density"],
        values["liquid_viscosity"],
    )

    design_velocity = flooding_velocity * design_flood_fraction
    # **The class answers its trial diameter here rather than dividing by zero.** ``if (uDesign
    # <= 0) return 1.0;`` is a guard, not a physical range, so the four intermediates are the
    # trial's own state and the answer is the trial rather than an error.
    if design_velocity <= 0.0:
        vapor_volumetric_flow = 0.0
        required_area = 0.0
        required_diameter = TRIAL_DIAMETER_M
        column_diameter = TRIAL_DIAMETER_M
    else:
        vapor_volumetric_flow = values["vapor_mass_flow"] / max(
            values["vapor_density"], MIN_VAPOR_DENSITY
        )
        required_area = vapor_volumetric_flow / design_velocity
        required_diameter = math.sqrt(4.0 * required_area / math.pi)
        column_diameter = round_to_standard_diameter(required_diameter)

    return PackingSizingResult(
        packing_name=resolved.name,
        packing_factor=resolved.packing_factor,
        flooding_velocity=from_si(flooding_velocity, "m/s"),
        design_velocity=from_si(design_velocity, "m/s"),
        vapor_volumetric_flow=from_si(vapor_volumetric_flow, "m**3/s"),
        required_area=from_si(required_area, "m**2"),
        required_diameter=from_si(required_diameter, "m"),
        column_diameter=from_si(column_diameter, "m"),
        warnings=tuple(warnings),
    )
