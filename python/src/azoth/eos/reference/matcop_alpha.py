"""``eos.matcop_alpha`` - the Mathias-Copeman alpha function.

Spec: ``specs/calcs/eos/matcop_alpha.toml``, which carries the provenance and why the
three coefficients are the caller's rather than derived.
"""

from __future__ import annotations

import math

from azoth._registry_gen import spec as _spec_for
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import MatcopAlphaResult
from azoth.core.warnings import Warning

CALC_ID = "eos.matcop_alpha"


def matcop_alpha(mc1: float, mc2: float, mc3: float, Tr: float) -> MatcopAlphaResult:
    """The Mathias-Copeman alpha function for a pure component.

    Args:
        mc1: the first Mathias-Copeman coefficient.
        mc2: the second Mathias-Copeman coefficient.
        mc3: the third Mathias-Copeman coefficient.
        Tr: the reduced temperature.

    Returns:
        The temperature-dependent alpha function.

    Raises:
        OutOfRangeError: if ``Tr`` is not positive.

    Example:
        >>> matcop_alpha(0.1, 0.05, 0.02, 0.7).alpha
        1.0358255509077803
    """
    spec = _spec_for(CALC_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    values = {"mc1": mc1, "mc2": mc2, "mc3": mc3, "Tr": Tr}
    apply_checks(checks.on_input, values.get, warnings)

    one_minus = 1.0 - math.sqrt(Tr)
    alpha = (1.0 + mc1 * one_minus + mc2 * one_minus**2 + mc3 * one_minus**3) ** 2

    apply_checks(checks.derived, lambda name: alpha if name == "alpha" else None, warnings)

    return MatcopAlphaResult(alpha=alpha, warnings=tuple(warnings))
