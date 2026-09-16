"""``eos.pr_lee_kesler_alpha`` - the Peng-Robinson alpha with a Soave-form m-factor.

Spec: ``specs/calcs/eos/pr_lee_kesler_alpha.toml``, which carries the provenance and why
the m-factor is Soave's rather than the Lee-Kesler method the class name suggests.
"""

from __future__ import annotations

import math

from azoth._registry_gen import spec as _spec_for
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import PrLeeKeslerAlphaResult
from azoth.core.warnings import Warning

CALC_ID = "eos.pr_lee_kesler_alpha"


def pr_lee_kesler_alpha(omega: float, Tr: float) -> PrLeeKeslerAlphaResult:
    """The Peng-Robinson alpha function with a Soave-form m-factor.

    Args:
        omega: the acentric factor.
        Tr: the reduced temperature.

    Returns:
        The temperature-dependent alpha function.

    Raises:
        OutOfRangeError: if ``Tr`` is not positive.

    Example:
        >>> pr_lee_kesler_alpha(0.1, 0.7).alpha
        1.2184305594583278
    """
    spec = _spec_for(CALC_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    values = {"omega": omega, "Tr": Tr}
    apply_checks(checks.on_input, values.get, warnings)

    m = 0.480 + 1.574 * omega - 0.176 * omega * omega
    t = 1.0 + m * (1.0 - math.sqrt(Tr))
    alpha = t * t

    apply_checks(checks.derived, lambda name: alpha if name == "alpha" else None, warnings)

    return PrLeeKeslerAlphaResult(alpha=alpha, warnings=tuple(warnings))
