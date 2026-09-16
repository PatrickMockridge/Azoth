"""``eos.matcop_prumr_new_alpha`` - the five-parameter Mathias-Copeman alpha, UMR-PR new.

Spec: ``specs/calcs/eos/matcop_prumr_new_alpha.toml``, which carries the provenance and
why the five coefficients are the caller's rather than derived.
"""

from __future__ import annotations

import math

from azoth._registry_gen import spec as _spec_for
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import MatcopPrumrNewAlphaResult
from azoth.core.warnings import Warning

CALC_ID = "eos.matcop_prumr_new_alpha"


def matcop_prumr_new_alpha(
    omega: float,
    mc1: float,
    mc2: float,
    mc3: float,
    mc4: float,
    mc5: float,
    Tr: float,
) -> MatcopPrumrNewAlphaResult:
    """The five-parameter Mathias-Copeman alpha function for a pure component, UMR-PR
    new variant.

    Args:
        omega: the acentric factor, used only for the fallback.
        mc1: the first Mathias-Copeman coefficient.
        mc2: the second Mathias-Copeman coefficient.
        mc3: the third Mathias-Copeman coefficient.
        mc4: the fourth Mathias-Copeman coefficient.
        mc5: the fifth Mathias-Copeman coefficient.
        Tr: the reduced temperature.

    Returns:
        The temperature-dependent alpha function.

    Raises:
        OutOfRangeError: if ``Tr`` is not positive.

    Example:
        >>> matcop_prumr_new_alpha(0.1, 0.5, 0.2, -0.1, 0.05, 0.01, 0.7).alpha
        1.1807146411895904
    """
    spec = _spec_for(CALC_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    values = {
        "omega": omega,
        "mc1": mc1,
        "mc2": mc2,
        "mc3": mc3,
        "mc4": mc4,
        "mc5": mc5,
        "Tr": Tr,
    }
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
        u = 1.0 - math.sqrt(Tr)
        s = (
            1.0
            + mc1 * u
            + mc2 * u * u
            + mc3 * u * u * u
            + mc4 * u * u * u * u
            + mc5 * u * u * u * u * u
        )
        alpha = s * s

    apply_checks(checks.derived, lambda name: alpha if name == "alpha" else None, warnings)

    return MatcopPrumrNewAlphaResult(alpha=alpha, warnings=tuple(warnings))
