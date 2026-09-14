"""``eos.pr_z_factor`` - the Peng-Robinson compressibility factor.

```text
z**3 - (1 - B)*z**2 + (A - 3*B**2 - 2*B)*z - (A*B - B**2 - B**3) = 0
```

Spec: ``specs/calcs/eos/pr_z_factor.toml``, which carries the provenance, the proof
that the admissible root count is one or three, and what the answer does near the
critical point.

This calc returns the outermost two roots that are admissible - ``z > B`` - and
discards the middle one, which lies on the unstable branch between the spinodals.
"""

from __future__ import annotations

from azoth._registry_gen import spec as _spec_for
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import PrZFactorResult, RootStructure
from azoth.core.solver import (
    Convergence,
    cubic_roots,
    require_cubic_converged,
)
from azoth.core.warnings import Warning

CALC_ID = "eos.pr_z_factor"


def pr_z_factor(a_reduced: float, b_reduced: float) -> PrZFactorResult:
    """The Peng-Robinson compressibility factor, for one state.

    Args:
        a_reduced: the cubic's attraction parameter ``A``, from
            :func:`azoth.eos.pr_alpha_ab`. Dimensionless, and already carrying the
            composition - this calc does not apply a mixing rule.
        b_reduced: the cubic's repulsion parameter ``B``, likewise.

    Returns:
        The smallest and largest admissible roots, and how many there were. The
        middle root is deliberately not returned; see the module documentation.

    Raises:
        OutOfRangeError: if ``b_reduced <= 0`` (the cubic degenerates to a trivial
            double root at zero pressure) or ``a_reduced < 0`` (which no
            Peng-Robinson state produces).
        SolverNotConvergedError: if the polish hits its cap.

    Example:
        >>> r = pr_z_factor(0.20206500174625697, 0.02431127309496514)
        >>> round(r.z_max, 12)
        0.790778966297
        >>> str(r.root_structure)
        'three_roots'
    """
    spec = _spec_for(CALC_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    values = {"a_reduced": a_reduced, "b_reduced": b_reduced}

    apply_checks(checks.on_input, values.get, warnings)

    solver = spec["solver"]
    convergence = Convergence.parse(solver["convergence"])

    # The monic cubic z**3 + c2*z**2 + c1*z + c0, coefficients straight from the
    # published form so a reader can check them against the equation above.
    c2 = -(1.0 - b_reduced)
    c1 = a_reduced - 3.0 * b_reduced * b_reduced - 2.0 * b_reduced
    c0 = -(a_reduced * b_reduced - b_reduced * b_reduced - b_reduced * b_reduced * b_reduced)

    outcome = cubic_roots(
        c2,
        c1,
        c0,
        solver["tolerance"],
        solver["max_iterations"],
        convergence,
    )
    outcome = require_cubic_converged(outcome, solver["tolerance"])

    # Admissible roots only. For B > 0 there is always at least one: the polynomial
    # equals -2*B**2 at z = B and tends to positive infinity, so it crosses zero
    # above B. `min`/`max` rather than the first and last element, so that the
    # filter's correctness does not depend on an ordering the solver happens to
    # guarantee.
    admissible = [z for z in outcome.roots if z > b_reduced]
    if not admissible:  # pragma: no cover - guarded by the input checks
        raise AssertionError(
            f"no admissible root for b_reduced = {b_reduced}; the bounds should have "
            f"refused this before the solver ran"
        )
    z_min = min(admissible)
    z_max = max(admissible)

    root_structure = RootStructure.ONE_ROOT if len(admissible) == 1 else RootStructure.THREE_ROOTS

    computed = {"z_min": z_min, "z_max": z_max}
    apply_checks(checks.derived, computed.get, warnings)

    return PrZFactorResult(
        z_min=z_min,
        z_max=z_max,
        root_structure=root_structure,
        iterations=outcome.iterations,
        converged=outcome.converged,
        residual=outcome.residual,
        warnings=tuple(warnings),
    )
