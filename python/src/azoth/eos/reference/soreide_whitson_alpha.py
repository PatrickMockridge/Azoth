"""``eos.soreide_whitson_alpha`` - the Soreide-Whitson alpha function for water.

Spec: ``specs/calcs/eos/soreide_whitson_alpha.toml``, which carries the provenance and
why the salinity is molality and the water-only scope.
"""

from __future__ import annotations

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

    a = 1.0 + 0.453 * (1.0 - Tr * (1.0 - 0.0103 * salinity**1.1)) + 0.0034 * ((1.0 / Tr) ** 3 - 1.0)
    alpha = a * a

    apply_checks(checks.derived, lambda name: alpha if name == "alpha" else None, warnings)

    return SoreideWhitsonAlphaResult(alpha=alpha, warnings=tuple(warnings))
