"""``eos.co2_water_diffusivity`` - the CO2-in-water binary diffusivity from the
Tammi correlation.

Spec: ``specs/calcs/eos/co2_water_diffusivity.toml``, which records the exponential
correlation.
"""

from __future__ import annotations

import math

from azoth._registry_gen import spec as _spec_for
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import Co2WaterDiffusivityResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning

CALC_ID = "eos.co2_water_diffusivity"


def co2_water_diffusivity(T: Q) -> Co2WaterDiffusivityResult:
    """The CO2-in-water binary diffusion coefficient, from NeqSim's `CO2water`.

    The correlation is temperature-only: it carries no solute or solvent argument.

    Args:
        T: absolute temperature.

    Returns:
        The binary diffusion coefficient, in ``m**2/s``.

    Raises:
        OutOfRangeError: if ``T`` is not positive.

    Example:
        >>> import azoth
        >>> q = azoth.ureg.Quantity
        >>> r = co2_water_diffusivity(q(298.15, "K"))
        >>> round(r.d.magnitude, 16)
        2.02082125395796e-09
    """
    spec = _spec_for(CALC_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    values = {
        "T": input_to_si(spec, "T", T),
    }
    apply_checks(checks.on_input, values.get, warnings)

    d = 0.03389 * math.exp(-2213.7 / values["T"]) * 1.0e-4

    apply_checks(checks.derived, lambda name: d if name == "d" else None, warnings)

    return Co2WaterDiffusivityResult(d=from_si(d, "m**2/s"), warnings=tuple(warnings))
