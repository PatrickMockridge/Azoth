"""``hydraulics.friction_factor_swamee_jain`` - the explicit Swamee-Jain friction
factor.

```text
f = 0.25 / (log10(relative_roughness / 3.7 + 5.74 / Re**0.9))**2
```

Swamee, P. K.; Jain, A. K. (1976). "Explicit equations for pipe-flow problems."
Journal of the Hydraulics Division, ASCE, 102(5), 657-664.

Spec: ``specs/calcs/hydraulics/friction_factor_swamee_jain.toml``

This is an *approximation to* Colebrook, not a more correct alternative to it. It
is offered so a caller can trade about 1% accuracy for the absence of an
iteration, and the wording matters: presenting the two as equivalent methods
would be wrong.
"""

from __future__ import annotations

import math

from azoth._registry_gen import spec as _spec_for
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import SwameeJainResult
from azoth.core.warnings import Warning

CALC_ID = "hydraulics.friction_factor_swamee_jain"


def friction_factor_swamee_jain(re: float, relative_roughness: float) -> SwameeJainResult:
    """Swamee-Jain explicit friction factor.

    Returns immediately: no iteration, and so no convergence report. Where the
    Colebrook result carries ``iterations`` and ``converged``, this carries only
    ``f``.

    Args:
        re: Reynolds number. Dimensionless.
        relative_roughness: ``epsilon / D``. Dimensionless.

    Raises:
        OutOfRangeError: if ``re <= 0`` (``Re**0.9`` is zero at the origin,
            making the ``5.74/Re**0.9`` term singular) or if
            ``relative_roughness < 0``.

    Outside the paper's stated range the value is still returned, carrying an
    ``OUT_OF_VALID_RANGE`` warning: the equation is well defined there, it simply
    is not covered by the accuracy claim.

    Example:
        >>> r = friction_factor_swamee_jain(100_000.0, 4.6e-4)
        >>> round(r.f, 8)
        0.02024003
    """
    spec = _spec_for(CALC_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    values = {"re": re, "relative_roughness": relative_roughness}
    apply_checks(checks.on_input, values.get, warnings)

    inner = math.log10(relative_roughness / 3.7 + 5.74 / re**0.9)
    f = 0.25 / (inner * inner)

    apply_checks(checks.derived, lambda name: f if name == "f" else None, warnings)

    return SwameeJainResult(f=f, warnings=tuple(warnings))
