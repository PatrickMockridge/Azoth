"""``eos.pr_alpha_ab`` - the Peng-Robinson alpha function and reduced parameters.

```text
alpha     = (1 + kappa*(1 - Tr**0.5))**2
a_reduced = Omega_a * alpha * Pr / Tr**2
b_reduced = Omega_b * Pr / Tr
```

Spec: ``specs/calcs/eos/pr_alpha_ab.yaml``, which carries the provenance, the
derivation of the two Omega constants from the triple-root condition, and what is
not claimed about the source.
"""

from __future__ import annotations

from azoth._registry_gen import spec as _spec_for
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import PrAlphaAbResult
from azoth.core.warnings import Warning

CALC_ID = "eos.pr_alpha_ab"

#: ``Omega_a``, the attraction constant of the Peng-Robinson cubic.
#:
#: Full precision, and load-bearing: the paper prints ``0.45724``, which is this
#: rounded, and which puts the cubic's critical point 4.55% out. See the spec's
#: ``notes``.
OMEGA_A = 0.4572355289213822

#: ``Omega_b``, the repulsion constant of the Peng-Robinson cubic.
#:
#: Full precision. The paper prints ``0.07780``.
OMEGA_B = 0.07779607390388846


def pr_alpha_ab(kappa: float, Tr: float, Pr: float) -> PrAlphaAbResult:
    """The Peng-Robinson alpha function and the reduced attraction parameters.

    Args:
        kappa: the alpha-function coefficient, from :func:`azoth.eos.pr_kappa`.
            Dimensionless.
        Tr: reduced temperature, ``T / Tc``. Dimensionless.
        Pr: reduced pressure, ``P / Pc``. Dimensionless.

    Returns:
        ``alpha``, ``a_reduced`` (the cubic's ``A``) and ``b_reduced`` (its ``B``).
        All dimensionless.

    Raises:
        OutOfRangeError: if ``Tr <= 0`` (a square root and a squared divisor) or
            ``Pr <= 0`` (``B`` is a divisor in the fugacity expression downstream).

    Example:
        >>> r = pr_alpha_ab(0.60282728832, 0.8, 0.25)
        >>> round(r.a_reduced, 17)
        0.20206500174625697
        >>> round(r.b_reduced, 17)
        0.02431127309496514
    """
    spec = _spec_for(CALC_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    values = {"kappa": kappa, "Tr": Tr, "Pr": Pr}

    apply_checks(checks.on_input, values.get, warnings)

    # A square, so never negative; exactly 1 at Tr = 1 whatever kappa is, which is
    # what makes the critical point special.
    #
    # Written as `t ** 2` to mirror the Rust side's `t * t`, following the same
    # convention `darcy_weisbach` uses for `v**2` / `v * v`.
    attraction = 1.0 + kappa * (1.0 - Tr**0.5)
    alpha = attraction**2
    # Guarded by the checks above, so Tr is positive here.
    a_reduced = OMEGA_A * alpha * Pr / Tr**2
    b_reduced = OMEGA_B * Pr / Tr

    computed = {"alpha": alpha, "a_reduced": a_reduced, "b_reduced": b_reduced}
    apply_checks(checks.derived, computed.get, warnings)

    return PrAlphaAbResult(
        alpha=alpha,
        a_reduced=a_reduced,
        b_reduced=b_reduced,
        warnings=tuple(warnings),
    )
