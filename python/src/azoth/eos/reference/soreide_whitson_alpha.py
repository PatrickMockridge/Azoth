"""``eos.soreide_whitson_alpha`` - the Soreide-Whitson alpha function for water.

Spec: ``specs/calcs/eos/soreide_whitson_alpha.toml``, which carries the provenance and
why the salinity is molality and the water-only scope.
"""

from __future__ import annotations

import math

from azoth._registry_gen import spec as _spec_for
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import SoreideWhitsonAlphaResult
from azoth.core.warnings import Warning

CALC_ID = "eos.soreide_whitson_alpha"


def soreide_whitson_alpha(salinity: float, Tr: float) -> SoreideWhitsonAlphaResult:
    """The Soreide-Whitson alpha function for water.

    Args:
        salinity: the molality (mol NaCl / kg H2O).
        Tr: the reduced temperature.

    Returns:
        The temperature-dependent alpha function.

    Raises:
        OutOfRangeError: if ``Tr`` is not positive or ``salinity`` is negative.

    Example:
        >>> soreide_whitson_alpha(1.0, 0.7).alpha
        1.3125796067429514
    """
    spec = _spec_for(CALC_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    values = {"salinity": salinity, "Tr": Tr}
    apply_checks(checks.on_input, values.get, warnings)

    a = (
        1.0
        + 0.453 * (1.0 - Tr * (1.0 - 0.0103 * math.pow(salinity, 1.1)))
        + 0.0034 * ((1.0 / Tr) ** 3 - 1.0)
    )
    alpha = a * a

    apply_checks(checks.derived, lambda name: alpha if name == "alpha" else None, warnings)

    return SoreideWhitsonAlphaResult(alpha=alpha, warnings=tuple(warnings))


def bracket(salinity: float, reduced_temperature: float) -> float:
    """``A(Tr)``, the bracket the alpha squares.

    Factored out because both derivatives need it and the three must agree by construction
    rather than by three copies of the same expression.
    """
    return (
        1.0
        + 0.453 * (1.0 - reduced_temperature * (1.0 - 0.0103 * math.pow(salinity, 1.1)))
        + 0.0034 * ((1.0 / reduced_temperature) ** 3 - 1.0)
    )


def diff_alpha_t(salinity: float, reduced_temperature: float, critical_temperature: float) -> float:
    """``d alpha / dT`` in **kelvin**, from ``AttractiveTermSoreideWhitson.diffalphaT``.

    ``alpha = A(Tr)^2`` and ``Tr = T / Tc``, so the chain rule carries one ``1/Tc``. `Tc` is
    an argument and not a constant, which is why this is not a registered calculation: the
    one beside it takes the reduced temperature and has no use for `Tc`, while a derivative
    in kelvin cannot drop it.
    """
    slope = (
        -0.453 * (1.0 - 0.0103 * math.pow(salinity, 1.1))
        - 3.0 * 0.0034 * (1.0 / reduced_temperature) ** 4
    )
    return 2.0 * bracket(salinity, reduced_temperature) * slope / critical_temperature


def diff2_alpha_t(
    salinity: float, reduced_temperature: float, critical_temperature: float
) -> float:
    """``d^2 alpha / dT^2`` in kelvin squared, from ``diffdiffalphaT``."""
    slope = (
        -0.453 * (1.0 - 0.0103 * math.pow(salinity, 1.1))
        - 3.0 * 0.0034 * (1.0 / reduced_temperature) ** 4
    )
    curvature = 12.0 * 0.0034 * (1.0 / reduced_temperature) ** 5
    scale = critical_temperature * critical_temperature
    return (2.0 * slope * slope + 2.0 * bracket(salinity, reduced_temperature) * curvature) / scale
