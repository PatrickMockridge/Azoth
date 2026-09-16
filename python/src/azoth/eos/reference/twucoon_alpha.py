"""``eos.twucoon_alpha`` - the Twu-Coon alpha function.

Spec: ``specs/calcs/eos/twucoon_alpha.toml``, which carries the six fitted constants
and the provenance.
"""

from __future__ import annotations

import math

from azoth._registry_gen import spec as _spec_for
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import TwucoonAlphaResult
from azoth.core.warnings import Warning

CALC_ID = "eos.twucoon_alpha"

_A = -0.201158
_B = 0.141599
_C = 2.29528
_D = -0.660145
_E = 0.500315
_F = 2.63165


def twucoon_alpha(omega: float, Tr: float) -> TwucoonAlphaResult:
    """The Twu-Coon alpha function for a pure component.

    Args:
        omega: the acentric factor.
        Tr: the reduced temperature.

    Returns:
        The temperature-dependent alpha function.

    Raises:
        OutOfRangeError: if ``Tr`` is not positive.

    Example:
        >>> twucoon_alpha(0.152, 0.7).alpha
        1.2469731129671788
    """
    spec = _spec_for(CALC_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    values = {"omega": omega, "Tr": Tr}
    apply_checks(checks.on_input, values.get, warnings)

    low = Tr**_A * math.exp(_B * (1.0 - Tr**_C))
    high = Tr**_D * math.exp(_E * (1.0 - Tr**_F))
    alpha = low + omega * (high - low)

    apply_checks(checks.derived, lambda name: alpha if name == "alpha" else None, warnings)

    return TwucoonAlphaResult(alpha=alpha, warnings=tuple(warnings))
