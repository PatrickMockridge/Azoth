"""``eos.pr_z_factor`` - the Peng-Robinson compressibility factor.

```text
z**3 - (1 - B)*z**2 + (A - 3*B**2 - 2*B)*z - (A*B - B**2 - B**3) = 0
```

Peng, D. Y.; Robinson, D. B. (1976). "A New Two-Constant Equation of State."
Ind. Eng. Chem. Fundam. 15(1), 59-64. DOI 10.1021/i160057a011

Spec: ``specs/calcs/eos/pr_z_factor.yaml``

# Which roots come back, and which do not

The cubic has up to three real roots. This calc returns the **outermost two of
those that are admissible**, where admissible means ``z > B``: ``z = B`` is the
zero-volume limit, and a root below it makes ``ln(z - B)`` the logarithm of a
negative number.

The middle root is discarded. It is a genuine root of the polynomial and it is not
a state the equation describes - it lies on the unstable branch between the
spinodals - and the spec asserts it is *absent* rather than merely unmentioned.

# Why only one or three, never two

``f(B)`` is exactly ``-2*B**2``, by cancellation of every other term:

```text
f(B) = B**3 - (1 - B)B**2 + (A - 3B**2 - 2B)B - (AB - B**2 - B**3)
     = -2B**2
```

which is negative for every permitted ``B``. A cubic with a positive leading
coefficient is negative below its smallest root and between its middle and largest
ones, so ``B`` lies either below all three roots or between the middle and largest.
One admissible root or three, and no third case.

# Near the critical point

The constants exist to place a triple root at the critical point, and a triple root
is cubically ill-conditioned. At ``A = OMEGA_A, B = OMEGA_B`` the polynomial stays
within ``1e-12`` of zero across a window about ``2e-4`` wide in ``z``, changing
sign repeatedly inside it on rounding noise alone - so every value in that window
is a root as far as double precision can tell. This calc cannot pin the last four
digits there, and cannot detect that it is in that region, because the reduction
happened two calcs upstream and it never sees ``Tr`` or ``Pr``.
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
    # published form so a reader can check them against the equation above. Written
    # exactly as the Rust arm writes them, statement for statement.
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
