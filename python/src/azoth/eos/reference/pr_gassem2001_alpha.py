"""``eos.pr_gassem2001_alpha`` - the Gassem (2001) alpha function.

Spec: ``specs/calcs/eos/pr_gassem2001_alpha.toml``, which carries the five fitted
constants and the provenance.
"""

from __future__ import annotations

import math

from azoth._registry_gen import spec as _spec_for
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import PrGassem2001AlphaResult
from azoth.core.warnings import Warning

CALC_ID = "eos.pr_gassem2001_alpha"

_A = 2.0
_B = 0.836
_C = 0.134
_D = 0.508
_E = -0.0467


def pr_gassem2001_alpha(omega: float, Tr: float) -> PrGassem2001AlphaResult:
    """The Gassem (2001) alpha function for a pure component.

    Args:
        omega: the acentric factor.
        Tr: the reduced temperature.

    Returns:
        The temperature-dependent alpha function.

    Raises:
        OutOfRangeError: if ``Tr`` is not positive.

    Example:
        >>> pr_gassem2001_alpha(0.152, 0.7).alpha
        1.2052404595821262
    """
    spec = _spec_for(CALC_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    values = {"omega": omega, "Tr": Tr}
    apply_checks(checks.on_input, values.get, warnings)

    exponent = _C + _D * omega + _E * omega * omega
    alpha = math.exp((_A + _B * Tr) * (1.0 - Tr**exponent))

    apply_checks(checks.derived, lambda name: alpha if name == "alpha" else None, warnings)

    return PrGassem2001AlphaResult(alpha=alpha, warnings=tuple(warnings))
