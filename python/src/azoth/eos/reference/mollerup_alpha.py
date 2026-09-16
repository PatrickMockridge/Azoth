"""``eos.mollerup_alpha`` - the Mollerup alpha function.

Spec: ``specs/calcs/eos/mollerup_alpha.toml``, which carries the provenance and why the
three parameters are the caller's rather than derived.
"""

from __future__ import annotations

import math

from azoth._registry_gen import spec as _spec_for
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import MollerupAlphaResult
from azoth.core.warnings import Warning

CALC_ID = "eos.mollerup_alpha"


def mollerup_alpha(p1: float, p2: float, p3: float, Tr: float) -> MollerupAlphaResult:
    """The Mollerup alpha function for a pure component.

    Args:
        p1: the first Mollerup parameter.
        p2: the second Mollerup parameter.
        p3: the third Mollerup parameter.
        Tr: the reduced temperature.

    Returns:
        The temperature-dependent alpha function.

    Raises:
        OutOfRangeError: if ``Tr`` is not positive.

    Example:
        >>> mollerup_alpha(0.1, 0.05, 0.02, 0.7).alpha
        1.0243735198192874
    """
    spec = _spec_for(CALC_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    values = {"p1": p1, "p2": p2, "p3": p3, "Tr": Tr}
    apply_checks(checks.on_input, values.get, warnings)

    alpha = 1.0 + p1 * (1.0 / Tr - 1.0) + p2 * Tr * math.log(Tr) + p3 * (Tr - 1.0)

    apply_checks(checks.derived, lambda name: alpha if name == "alpha" else None, warnings)

    return MollerupAlphaResult(alpha=alpha, warnings=tuple(warnings))
