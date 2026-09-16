"""``eos.umrpr_alpha`` - the UMR-PR alpha function.

Spec: ``specs/calcs/eos/umrpr_alpha.toml``, which carries the provenance and why the
m-factor is the UMR-PR kappa.
"""

from __future__ import annotations

import math

from azoth._registry_gen import spec as _spec_for
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import UmrprAlphaResult
from azoth.core.warnings import Warning

CALC_ID = "eos.umrpr_alpha"


def umrpr_alpha(omega: float, Tr: float) -> UmrprAlphaResult:
    """The UMR-PR alpha function for a pure component.

    Args:
        omega: the acentric factor.
        Tr: the reduced temperature.

    Returns:
        The temperature-dependent alpha function.

    Raises:
        OutOfRangeError: if ``Tr`` is not positive.

    Example:
        >>> umrpr_alpha(0.1, 0.7).alpha
        1.1822586823466172
    """
    spec = _spec_for(CALC_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    values = {"omega": omega, "Tr": Tr}
    apply_checks(checks.on_input, values.get, warnings)

    m = (
        0.384401
        + 1.52276 * omega
        - 0.213808 * omega * omega
        + 0.034616 * omega * omega * omega
        - 0.001976 * omega * omega * omega * omega
    )
    t = 1.0 + m * (1.0 - math.sqrt(Tr))
    alpha = t * t

    apply_checks(checks.derived, lambda name: alpha if name == "alpha" else None, warnings)

    return UmrprAlphaResult(alpha=alpha, warnings=tuple(warnings))
