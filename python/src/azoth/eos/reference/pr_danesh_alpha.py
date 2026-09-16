"""``eos.pr_danesh_alpha`` - the Danesh alpha function.

Spec: ``specs/calcs/eos/pr_danesh_alpha.toml``, which carries the provenance and why the
1.21 factor damps the attraction past the critical temperature.
"""

from __future__ import annotations

import math

from azoth._registry_gen import spec as _spec_for
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import PrDaneshAlphaResult
from azoth.core.warnings import Warning

CALC_ID = "eos.pr_danesh_alpha"


def pr_danesh_alpha(omega: float, Tr: float) -> PrDaneshAlphaResult:
    """The Danesh alpha function for a pure component.

    Args:
        omega: the acentric factor.
        Tr: the reduced temperature.

    Returns:
        The temperature-dependent alpha function.

    Raises:
        OutOfRangeError: if ``Tr`` is not positive.

    Example:
        >>> pr_danesh_alpha(0.152, 0.7).alpha
        1.2066270990034567
    """
    spec = _spec_for(CALC_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    values = {"omega": omega, "Tr": Tr}
    apply_checks(checks.on_input, values.get, warnings)

    m = 0.37464 + 1.54226 * omega - 0.26992 * omega * omega
    m_mod = 1.21 * m if Tr > 1.0 else m
    one_minus = 1.0 - math.sqrt(Tr)
    alpha = (1.0 + m_mod * one_minus) * (1.0 + m_mod * one_minus)

    apply_checks(checks.derived, lambda name: alpha if name == "alpha" else None, warnings)

    return PrDaneshAlphaResult(alpha=alpha, warnings=tuple(warnings))
