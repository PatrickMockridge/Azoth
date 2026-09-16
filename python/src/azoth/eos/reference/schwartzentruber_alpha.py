"""``eos.schwartzentruber_alpha`` - the Schwartzentruber-Renon alpha function.

Spec: ``specs/calcs/eos/schwartzentruber_alpha.toml``, which carries the provenance and
why the three parameters are the caller's rather than derived.
"""

from __future__ import annotations

import math

from azoth._registry_gen import spec as _spec_for
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import SchwartzentruberAlphaResult
from azoth.core.warnings import Warning

CALC_ID = "eos.schwartzentruber_alpha"


def schwartzentruber_alpha(
    omega: float, p1: float, p2: float, p3: float, Tr: float
) -> SchwartzentruberAlphaResult:
    """The Schwartzentruber-Renon alpha function for a pure component.

    Args:
        omega: the acentric factor.
        p1: the first Schwartzentruber-Renon parameter.
        p2: the second Schwartzentruber-Renon parameter.
        p3: the third Schwartzentruber-Renon parameter.
        Tr: the reduced temperature.

    Returns:
        The temperature-dependent alpha function.

    Raises:
        OutOfRangeError: if ``Tr`` is not positive.

    Example:
        >>> schwartzentruber_alpha(0.1, 0.3, 0.2, 0.1, 0.7).alpha
        0.9946408503265206
    """
    spec = _spec_for(CALC_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    values = {"omega": omega, "p1": p1, "p2": p2, "p3": p3, "Tr": Tr}
    apply_checks(checks.on_input, values.get, warnings)

    m = 0.48508 + 1.55191 * omega - 0.15613 * omega * omega
    inner = 1.0 + m * (1.0 - math.sqrt(Tr)) - p1 * (1.0 - Tr) * (1.0 + p2 * Tr + p3 * Tr * Tr)
    alpha = inner * inner

    apply_checks(checks.derived, lambda name: alpha if name == "alpha" else None, warnings)

    return SchwartzentruberAlphaResult(alpha=alpha, warnings=tuple(warnings))
