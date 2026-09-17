"""``eos.nitric_sulfuric_acid_vapor_pressure`` - the three pure-component vapour
pressures of the water-nitric-sulfuric acid system.

Spec: ``specs/calcs/eos/nitric_sulfuric_acid_vapor_pressure.toml``

This is the pure-Python reference: a second, independent expression of the same physics
as the Rust kernel, held to it by ``test_cross_impl.py``.
"""

from __future__ import annotations

import math

from azoth.core.errors import OutOfRangeError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import NitricSulfuricAcidVaporPressureResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning

CALC_ID = "eos.nitric_sulfuric_acid_vapor_pressure"

#: Pascal per millibar, the unit the water correlation is written in.
MBAR_TO_PA = 100.0

#: Pascal per torr, the unit the nitric-acid Antoine is written in.
TORR_TO_PA = 133.322368421

#: Pascal per atmosphere, the unit the sulfuric-acid correlation is written in.
ATM_TO_PA = 101325.0

#: The nitric-acid Antoine coefficients, NeqSim's refit of Pennington's pair.
HNO3_ANTOINE_A = 7.57628
HNO3_ANTOINE_B = 1470.385
HNO3_ANTOINE_C = 43.0


def nitric_sulfuric_acid_vapor_pressure(
    T: Q,
) -> NitricSulfuricAcidVaporPressureResult:
    """The pure-component vapour pressures of water, nitric acid and sulfuric acid.

    The three come back together because a phase over this system needs all three at one
    temperature. Each is its own correlation: water a ``log10 P/mbar`` polynomial in
    ``1/T``, nitric acid an Antoine in ``log10 P/torr``, sulfuric acid a ``ln P/atm``
    straight line.

    Args:
        T: absolute temperature in kelvin.

    Returns:
        The three pure-component saturation pressures, in pascals, and any caveats.

    Raises:
        OutOfRangeError: if ``T`` is not positive, or at or below the nitric-acid form's
            43 K pole.

    An out-of-range ``T`` carries ``OUT_OF_VALID_RANGE`` rather than failing: the
    arithmetic is well defined outside 190-298 K.

    Example:
        >>> import azoth
        >>> q = azoth.ureg.Quantity
        >>> r = nitric_sulfuric_acid_vapor_pressure(q(273.15, "K"))
        >>> round(r.p_water.to("Pa").magnitude, 6)
        610.359225
    """
    from azoth import _registry_gen

    spec = _registry_gen.BY_ID[CALC_ID]
    checks = checks_for(spec)
    warnings: list[Warning] = []

    t_si = input_to_si(spec, "T", T)
    apply_checks(checks.on_input, {"T": t_si}.get, warnings)

    if t_si <= HNO3_ANTOINE_C:
        raise OutOfRangeError(
            "T",
            t_si,
            "the nitric-acid Antoine form is `B/(T - 43.0)`, so 43 K and below are a "
            "pole rather than a pressure",
        )

    apply_checks(checks.derived, {"T": t_si}.get, warnings)

    # Water: `log10(P/mbar) = 8.42926609 - 1827.17843/T - 71208.271/T^2`.
    log10_mbar = 8.42926609 - 1827.17843 / t_si - 71208.271 / (t_si * t_si)
    p_water = 10.0**log10_mbar * MBAR_TO_PA

    # Nitric acid: `log10(P/torr) = A - B/(T - C)`.
    log10_torr = HNO3_ANTOINE_A - HNO3_ANTOINE_B / (t_si - HNO3_ANTOINE_C)
    p_nitric_acid = 10.0**log10_torr * TORR_TO_PA

    # Sulfuric acid: `ln(P/atm) = -10156.0/T + 16.259`.
    p_sulfuric_acid = math.exp(-10156.0 / t_si + 16.259) * ATM_TO_PA

    return NitricSulfuricAcidVaporPressureResult(
        p_water=from_si(p_water, "Pa"),
        p_nitric_acid=from_si(p_nitric_acid, "Pa"),
        p_sulfuric_acid=from_si(p_sulfuric_acid, "Pa"),
        warnings=tuple(warnings),
    )
