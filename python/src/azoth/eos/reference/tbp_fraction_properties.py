"""``eos.tbp_fraction_properties`` - a TBP pseudo-component's properties from two numbers.

```text
tc, pc, tb, m = Pedersen's TBP correlations(molar_mass, density)
omega         = 3/7 log10(pc/1.01325) / (tc/tb - 1) - 1
```

Spec: ``specs/calcs/eos/tbp_fraction_properties.toml``, which carries the coefficient sets,
where the two branches switch, and the measurement of what the heavy set does above
``mw = 1120`` g/mol.

NeqSim's ``addTBPfraction`` is this and nothing else: a plus fraction is split into cuts
described by a molar mass and a normal liquid density, and these correlations are what turn a
cut into a component a cubic can take. The units are the correlation's own inside - **g/mol
and g/cm3** - and the conversion is at the boundary, because that is where NeqSim puts it.
"""

from __future__ import annotations

import math

from azoth._registry_gen import spec as _spec_for
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import TbpFractionPropertiesResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning

CALC_ID = "eos.tbp_fraction_properties"

#: The molar mass at and above which the critical-property coefficients change set, g/mol.
HEAVY_CUT_MOLAR_MASS = 1120.0

#: The molar mass at and above which the boiling point takes the power law, g/mol.
POWER_LAW_MOLAR_MASS = 540.0

#: The reference pressure the acentric factor is reduced against, bar.
REFERENCE_PRESSURE_BAR = 1.01325

#: ``PedersenTBPModelSRK.TBPfractionCoefOil``, as ``[tc, pc, m]`` of ``[c0, c1, c2, c3, c4]``.
OIL: tuple[tuple[float, float, float, float, float], ...] = (
    (163.12, 86.052, 0.43475, -1877.4, 0.0),
    (-0.13408, 2.5019, 208.46, -3987.2, 1.0),
    (0.7431, 0.0048122, 0.0096707, -3.7184e-6, 0.0),
)

#: ``TBPfractionCoefsHeavyOil``, in force at and above ``mw = 1120`` g/mol.
HEAVY: tuple[tuple[float, float, float, float, float], ...] = (
    (8.3063e2, 1.75228e1, 4.55911e-2, -1.13484e4, 0.0),
    (8.02988e-1, 1.78396, 1.56740e2, -6.96559e3, 0.25),
    (-4.7268e-2, 6.02931e-2, 1.21051, -5.76676e-3, 0.0),
)


def tbp_fraction_properties(molar_mass: Q, density: Q) -> TbpFractionPropertiesResult:
    """A TBP cut's critical properties, boiling point, acentric factor and alpha exponent.

    Args:
        molar_mass: the cut's molar mass. A pseudo-component has no databank row, so this is
            what stands in for its identity.
        density: the cut's normal liquid density, at 15 C and 1 atm.

    Returns:
        The five properties, and any caveats. An acentric factor outside ``[-1, 2]`` -
        which the heavy coefficient set produces - comes back carrying
        ``OUT_OF_VALID_RANGE`` rather than as an error.

    Raises:
        OutOfRangeError: if the molar mass or the density is not positive, each of which the
            correlation divides by or raises to a fractional power.
    """
    spec = _spec_for(CALC_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []
    mass_si = input_to_si(spec, "molar_mass", molar_mass)
    density_si = input_to_si(spec, "density", density)
    apply_checks(checks.on_input, {"molar_mass": mass_si, "density": density_si}.get, warnings)

    # NeqSim's own units, which is what the coefficients are fitted in.
    mw = mass_si * 1000.0
    d = density_si / 1000.0
    coefficients = OIL if mw < HEAVY_CUT_MOLAR_MASS else HEAVY

    tc = (
        coefficients[0][0] * d
        + coefficients[0][1] * math.log(mw)
        + coefficients[0][2] * mw
        + coefficients[0][3] / mw
    )
    pc_bar = math.exp(
        0.01325
        + coefficients[1][0]
        + coefficients[1][1] * d ** coefficients[1][4]
        + coefficients[1][2] / mw
        + coefficients[1][3] / mw**2
    )
    if mw < POWER_LAW_MOLAR_MASS:
        tb = 2.0e-6 * mw**3 - 0.0035 * mw**2 + 2.4003 * mw + 171.74
    else:
        tb = 97.58 * mw**0.3323 * d**0.04609
    exponent = (
        coefficients[2][0]
        + coefficients[2][1] * mw
        + coefficients[2][2] * d
        + coefficients[2][3] * mw**2
    )
    acentric = 3.0 / 7.0 * math.log10(pc_bar / REFERENCE_PRESSURE_BAR) / (tc / tb - 1.0) - 1.0

    apply_checks(checks.derived, {"acentric_factor": acentric}.get, warnings)

    return TbpFractionPropertiesResult(
        tc=from_si(tc, "K"),
        pc=from_si(pc_bar * 1.0e5, "Pa"),
        boiling_temperature=from_si(tb, "K"),
        acentric_factor=acentric,
        attraction_exponent=exponent,
        warnings=tuple(warnings),
    )


__all__ = ["tbp_fraction_properties"]
