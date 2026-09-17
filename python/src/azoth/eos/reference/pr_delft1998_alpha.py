"""``eos.pr_delft1998_alpha`` - the Peng-Robinson alpha function, Delft (1998).

Spec: ``specs/calcs/eos/pr_delft1998_alpha.toml``, which carries the provenance and the
methane-specific branch this general form does not express.
"""

from __future__ import annotations

import math

from azoth._registry_gen import spec as _spec_for
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import PrDelft1998AlphaResult
from azoth.core.warnings import Warning

CALC_ID = "eos.pr_delft1998_alpha"


def pr_delft1998_alpha(omega: float, Tr: float) -> PrDelft1998AlphaResult:
    """The Peng-Robinson alpha function, Delft (1998), for a pure component.

    Args:
        omega: the acentric factor.
        Tr: the reduced temperature.

    Returns:
        The temperature-dependent alpha function.

    Raises:
        OutOfRangeError: if ``Tr`` is not positive.

    Example:
        >>> pr_delft1998_alpha(0.1, 0.7).alpha
        1.1792745256672486
    """
    spec = _spec_for(CALC_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    values = {"omega": omega, "Tr": Tr}
    apply_checks(checks.on_input, values.get, warnings)

    if omega > 0.49:
        m = 0.379642 + 1.48503 * omega - 0.164423 * omega * omega + 0.01666 * omega * omega * omega
    else:
        m = 0.37464 + 1.54226 * omega - 0.26992 * omega * omega
    t = 1.0 + m * (1.0 - math.sqrt(Tr))
    alpha = t * t

    apply_checks(checks.derived, lambda name: alpha if name == "alpha" else None, warnings)

    return PrDelft1998AlphaResult(alpha=alpha, warnings=tuple(warnings))
