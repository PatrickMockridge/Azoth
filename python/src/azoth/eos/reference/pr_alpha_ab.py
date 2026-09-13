"""``eos.pr_alpha_ab`` - the Peng-Robinson alpha function and reduced parameters.

```text
alpha     = (1 + kappa*(1 - Tr**0.5))**2
a_reduced = Omega_a * alpha * Pr / Tr**2
b_reduced = Omega_b * Pr / Tr
```

Peng, D. Y.; Robinson, D. B. (1976). "A New Two-Constant Equation of State."
Ind. Eng. Chem. Fundam. 15(1), 59-64. DOI 10.1021/i160057a011

Spec: ``specs/calcs/eos/pr_alpha_ab.yaml``

# Why the Omega constants are not the printed ones

The paper prints ``0.45724`` and ``0.07780``. Those are roundings, and using them
makes the cubic fail its own critical point: at ``Tr = Pr = 1`` the real root comes
out 0.321379025174 instead of Peng-Robinson's critical compressibility
0.307401308699 - 4.55% wrong exactly where the equation is anchored.

The cause is conditioning. The constants exist to place a *triple* root at the
critical point, and a triple root is cubically ill-conditioned: perturbing the
coefficients by ``epsilon`` moves the roots by about ``epsilon ** (1/3)``.

So this module carries the full-precision pair, which is the unique solution of the
triple-root condition rather than a value transcribed from anywhere. The Rust half
asserts both - that these satisfy it, and that the printed pair does not. See the
spec's ``verification`` notes for what remains unconfirmed.
"""

from __future__ import annotations

from azoth._registry_gen import spec as _spec_for
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import PrAlphaAbResult
from azoth.core.warnings import Warning

CALC_ID = "eos.pr_alpha_ab"

#: ``Omega_a``, the attraction constant of the Peng-Robinson cubic.
#:
#: Full precision, and load-bearing. The paper prints ``0.45724``, which is this
#: rounded, and which puts the cubic's critical point 4.55% out. Mirrored exactly
#: in ``crates/azoth-eos/src/pr_alpha_ab.rs``.
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

    # Guarded by the checks above, so Tr is positive here.
    #
    # Written as `attraction ** 2` to mirror the Rust side's `attraction * attraction`,
    # following the same convention `darcy_weisbach` uses for `v**2` / `v * v`.
    #
    # Measured: the two backends give bit-identical results for all three outputs on
    # the worked example. Not *guaranteed* - IEEE-754 pins `+ - * /` and `sqrt`, and
    # says nothing about how a language lowers `x ** 2` - so the claim this calc makes
    # is agreement within the spec's tolerance, not bit-equality.
    attraction = 1.0 + kappa * (1.0 - Tr**0.5)
    alpha = attraction**2
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
