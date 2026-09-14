"""``hydraulics.friction_factor_haaland`` - the explicit Haaland friction factor.

```text
f = (1.0 / (-1.8 * log10((relative_roughness / 3.7)**1.11 + 6.9 / Re)))**2
```

Haaland, S. E. (1983). "Simple and explicit formulas for the friction factor in
turbulent pipe flow." Journal of Fluids Engineering, 105(1), 89-90.

Spec: ``specs/calcs/hydraulics/friction_factor_haaland.toml``

This is an *approximation to* Colebrook, not a more correct alternative to it, and
it is the second such approximation in this registry alongside Swamee-Jain. Having
two is not redundancy: they were fitted differently, they are accurate to about 1%
and 2% respectively, and a caller comparing them at the same inputs gets a cheap
sense of how much the choice of explicit form matters. Where they disagree by more
than either claims against Colebrook, that is worth knowing.

The citation is `unverified` - see the spec's ``verification`` block. The equation
is not in doubt; nobody has opened the paper and checked this against it.
"""

from __future__ import annotations

import math

from azoth._registry_gen import spec as _spec_for
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import HaalandResult
from azoth.core.warnings import Warning

CALC_ID = "hydraulics.friction_factor_haaland"


def friction_factor_haaland(re: float, relative_roughness: float) -> HaalandResult:
    """Haaland explicit friction factor.

    Returns immediately: no iteration, and so no convergence report.

    Args:
        re: Reynolds number. Dimensionless.
        relative_roughness: ``epsilon / D``. Dimensionless.

    Raises:
        OutOfRangeError: if ``re <= 0`` (the ``6.9/Re`` term is singular at the
            origin) or if ``relative_roughness < 0``.

    Outside the range in which the correlation describes turbulent pipe flow the
    value is still returned, carrying an ``OUT_OF_VALID_RANGE`` warning: the
    equation is well defined there, it simply is not covered by the accuracy
    claim. That range is the Moody/Colebrook framework's, not a figure read from
    Haaland's paper - the spec says so, and the reason is in its ``verification``.

    Example:
        >>> r = friction_factor_haaland(100_000.0, 4.6e-4)
        >>> round(r.f, 9)
        0.019898058
    """
    spec = _spec_for(CALC_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    values = {"re": re, "relative_roughness": relative_roughness}
    apply_checks(checks.on_input, values.get, warnings)

    # `1/sqrt(f) = -1.8 * log10(...)`, so `f` is the reciprocal of that, squared.
    # Written in that order rather than expanded, so the line matches the published
    # form a reader will be checking it against.
    #
    # The square is `x * x` rather than `x ** 2`, matching Swamee-Jain and the Rust
    # side: `pow` is not required to be correctly rounded for an integer exponent on
    # every libm, and the two implementations should not differ in the last bits for
    # a reason that is only about how the exponent was spelled.
    inner = (relative_roughness / 3.7) ** 1.11 + 6.9 / re
    inv_sqrt_f = -1.8 * math.log10(inner)
    f = 1.0 / (inv_sqrt_f * inv_sqrt_f)

    apply_checks(checks.derived, lambda name: f if name == "f" else None, warnings)

    return HaalandResult(f=f, warnings=tuple(warnings))
