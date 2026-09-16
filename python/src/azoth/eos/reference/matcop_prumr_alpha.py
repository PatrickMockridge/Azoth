"""``eos.matcop_prumr_alpha`` - the Mathias-Copeman alpha with the UMR-PR fallback.

Spec: ``specs/calcs/eos/matcop_prumr_alpha.toml``, which carries the provenance and why
the three coefficients are the caller's rather than derived.
"""

from __future__ import annotations

import math

from azoth._registry_gen import spec as _spec_for
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import MatcopPrumrAlphaResult
from azoth.core.warnings import Warning

CALC_ID = "eos.matcop_prumr_alpha"


def matcop_prumr_alpha(
    omega: float, mc1: float, mc2: float, mc3: float, Tr: float
) -> MatcopPrumrAlphaResult:
    """The Mathias-Copeman alpha function for a pure component, with the standard
    Peng-Robinson alpha when the first coefficient is essentially zero.

    Args:
        omega: the acentric factor, used only for the fallback.
        mc1: the first Mathias-Copeman coefficient.
        mc2: the second Mathias-Copeman coefficient.
        mc3: the third Mathias-Copeman coefficient.
        Tr: the reduced temperature.

    Returns:
        The temperature-dependent alpha function.

    Raises:
        OutOfRangeError: if ``Tr`` is not positive.

    Example:
        >>> matcop_prumr_alpha(0.1, 0.5, 0.2, -0.1, 0.7).alpha
        1.180634768967036
    """
    spec = _spec_for(CALC_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    values = {"omega": omega, "mc1": mc1, "mc2": mc2, "mc3": mc3, "Tr": Tr}
    apply_checks(checks.on_input, values.get, warnings)

    if mc1 < 1e-20:
        m = (
            0.384401
            + 1.52276 * omega
            - 0.213808 * omega * omega
            + 0.034616 * omega * omega * omega
            - 0.001976 * omega * omega * omega * omega
        )
        t = 1.0 + m * (1.0 - math.sqrt(Tr))
        alpha = t * t
    else:
        root = 1.0 - math.sqrt(Tr)
        t = 1.0 + mc1 * root + mc2 * root * root + mc3 * root * root * root
        alpha = t * t

    apply_checks(checks.derived, lambda name: alpha if name == "alpha" else None, warnings)

    return MatcopPrumrAlphaResult(alpha=alpha, warnings=tuple(warnings))
