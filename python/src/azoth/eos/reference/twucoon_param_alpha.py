"""``eos.twucoon_param_alpha`` - the Twu-Coon parameter alpha function.

Spec: ``specs/calcs/eos/twucoon_param_alpha.toml``, which carries the provenance and why
the three parameters are the caller's rather than derived.
"""

from __future__ import annotations

import math

from azoth._registry_gen import spec as _spec_for
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import TwucoonParamAlphaResult
from azoth.core.warnings import Warning

CALC_ID = "eos.twucoon_param_alpha"


def twucoon_param_alpha(a: float, b: float, c: float, Tr: float) -> TwucoonParamAlphaResult:
    """The Twu-Coon parameter alpha function for a pure component.

    Args:
        a: the first Twu-Coon parameter.
        b: the second Twu-Coon parameter.
        c: the third Twu-Coon parameter.
        Tr: the reduced temperature.

    Returns:
        The temperature-dependent alpha function.

    Raises:
        OutOfRangeError: if ``Tr`` is not positive.

    Example:
        >>> twucoon_param_alpha(0.1, 0.5, 2.0, 0.7).alpha
        1.4720779056478814
    """
    spec = _spec_for(CALC_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    values = {"a": a, "b": b, "c": c, "Tr": Tr}
    apply_checks(checks.on_input, values.get, warnings)

    alpha = Tr ** (c * (b - 1.0)) * math.exp(a * (1.0 - Tr ** (b * c)))

    apply_checks(checks.derived, lambda name: alpha if name == "alpha" else None, warnings)

    return TwucoonParamAlphaResult(alpha=alpha, warnings=tuple(warnings))
