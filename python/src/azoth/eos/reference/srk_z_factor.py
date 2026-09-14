"""``eos.srk_z_factor`` - the Soave-Redlich-Kwong compressibility factor.

```text
z**3 - z**2 + (a_reduced - b_reduced - b_reduced**2)*z - a_reduced*b_reduced = 0
```

Spec: ``specs/calcs/eos/srk_z_factor.toml``, which carries the provenance and the same
admissible-root theorem as the Peng-Robinson form.

This calc returns the outermost two roots that are admissible - ``z > B`` - and
discards the middle one, exactly as ``eos.pr_z_factor`` does.
"""

from __future__ import annotations

from azoth._registry_gen import spec as _spec_for
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import RootStructure, SrkZFactorResult
from azoth.core.solver import (
    Convergence,
    cubic_roots,
    require_cubic_converged,
)
from azoth.core.warnings import Warning
from azoth.eos.cubic import SRK

CALC_ID = "eos.srk_z_factor"


def srk_z_factor(a_reduced: float, b_reduced: float) -> SrkZFactorResult:
    """The Soave-Redlich-Kwong compressibility factor, for one state.

    Args:
        a_reduced: the cubic's attraction parameter ``A``, from
            :func:`azoth.eos.srk_alpha_ab`.
        b_reduced: the cubic's repulsion parameter ``B``, likewise.

    Raises:
        OutOfRangeError: if ``b_reduced <= 0`` or ``a_reduced < 0``.
        SolverNotConvergedError: if the polish hits its cap.

    Example:
        >>> r = srk_z_factor(0.19315231739218255, 0.02707510936404929)
        >>> round(r.z_max, 12)
        0.801955197256
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

    c2, c1, c0 = SRK.z_coefficients(a_reduced, b_reduced)

    outcome = cubic_roots(
        c2,
        c1,
        c0,
        solver["tolerance"],
        solver["max_iterations"],
        convergence,
    )
    outcome = require_cubic_converged(outcome, solver["tolerance"])

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

    return SrkZFactorResult(
        z_min=z_min,
        z_max=z_max,
        root_structure=root_structure,
        iterations=outcome.iterations,
        converged=outcome.converged,
        residual=outcome.residual,
        warnings=tuple(warnings),
    )
