"""``hydraulics.friction_factor_colebrook`` - the implicit Colebrook-White
friction factor.

```text
1 / sqrt(f) = -2 * log10(relative_roughness / 3.7 + 2.51 / (Re * sqrt(f)))
```

Colebrook, C. F. (1939). "Turbulent flow in pipes, with particular reference to
the transition region between the smooth and rough pipe laws." J. Inst. Civil
Engineers 11(4), 133-156. DOI 10.1680/ijoti.1939.13150

Spec: ``specs/calcs/hydraulics/friction_factor_colebrook.toml``

The equation is implicit in ``f``, so it is solved rather than evaluated.
Substituting ``x = 1/sqrt(f)`` makes it a fixed point in ``x``, which converges
quickly and - crucially - converges the *same way* in both languages, because the
scheme is fixed by the spec rather than chosen per implementation.
"""

from __future__ import annotations

import math

from azoth._registry_gen import spec as _spec_for
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import ColebrookResult
from azoth.core.solver import (
    Convergence,
    fixed_point,
    require_converged,
)
from azoth.core.warnings import Warning

CALC_ID = "hydraulics.friction_factor_colebrook"


def friction_factor_colebrook(re: float, relative_roughness: float) -> ColebrookResult:
    """Solve the Colebrook-White equation for the Darcy friction factor.

    Args:
        re: Reynolds number. Dimensionless.
        relative_roughness: ``epsilon / D``. Dimensionless.

    Returns:
        The friction factor and the solver's own report. That report is not
        decoration: when ``converged`` is false, ``f`` is the last iterate of an
        iteration that did not meet tolerance and does not solve the equation.
        This function only returns such a value when it converged - a
        non-converged run raises.

    Raises:
        OutOfRangeError: if ``re <= 0`` (the ``2.51/(re*sqrt(f))`` term is
            singular) or if ``relative_roughness < 0`` (unphysical).
        SolverNotConvergedError: if the iteration hits its cap. With the spec's
            tolerance and cap this does not happen for physical inputs, and a
            test asserts so.

    Below ``Re = 4000`` the result is out of the range the correlation was fitted
    to; the value is still returned, carrying an ``OUT_OF_VALID_RANGE`` warning.
    Use ``64/Re`` for laminar flow - a different equation, deliberately not part
    of this calc.

    Example:
        >>> # Commercial steel: epsilon = 0.046 mm in a 100 mm pipe.
        >>> r = friction_factor_colebrook(100_000.0, 4.6e-4)
        >>> round(r.f, 8)
        0.02016203
        >>> r.converged
        True
    """
    spec = _spec_for(CALC_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    values = {"re": re, "relative_roughness": relative_roughness}
    apply_checks(checks.on_input, values.get, warnings)

    solver = spec["solver"]
    convergence = Convergence.parse(solver["convergence"])

    # x = 1/sqrt(f). The initial guess is a friction factor, so its reciprocal
    # square root is the starting x.
    x0 = 1.0 / math.sqrt(solver["initial_guess"])
    outcome = fixed_point(
        x0=x0,
        tolerance=solver["tolerance"],
        max_iterations=solver["max_iterations"],
        convergence=convergence,
        f=lambda x: -2.0 * math.log10(relative_roughness / 3.7 + 2.51 * x / re),
    )
    outcome = require_converged(outcome, solver["tolerance"])

    f = 1.0 / (outcome.x * outcome.x)
    apply_checks(checks.derived, lambda name: f if name == "f" else None, warnings)

    return ColebrookResult(
        f=f,
        iterations=outcome.iterations,
        converged=outcome.converged,
        residual=outcome.residual,
        warnings=tuple(warnings),
    )


def fully_rough_limit(relative_roughness: float) -> float:
    """The fully-rough asymptote of the Colebrook equation, which is explicit.

    As ``Re -> infinity`` the roughness term dominates and
    ``1/sqrt(f) -> -2*log10(relative_roughness/3.7)``. Exposed because it is a
    useful sanity bound, and because a test uses it to confirm the solver
    converges to the correct asymptote rather than to a nearby fixed point.
    """
    x = -2.0 * math.log10(relative_roughness / 3.7)
    return 1.0 / (x * x)
