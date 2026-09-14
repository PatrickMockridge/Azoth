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
#: **NeqSim 3.20.0's value, not the paper's, and carrying it is the point.** This
#: library is a port, and NeqSim's ``ComponentPR`` constructor sets
#: ``a = .45724333333 * R**2 * Tc**2 / Pc``. The Peng-Robinson paper prints
#: ``0.45724`` and the cubic's triple-root condition gives ``0.4572355289213822``
#: exactly; NeqSim's is neither. It is the paper's printed value plus 3.3333e-6 -
#: and so is its ``Omega_b``, by the same offset - which is what makes the pair read
#: as a transcription artefact carried forward rather than as a refit.
#:
#: Substituting NeqSim's pair for the exact one reproduces NeqSim's own ``TPflash``
#: to twelve significant figures, where the exact pair leaves a 1.4e-4 residue - see
#: ``validation/eos/methane_butane_flash_against_neqsim.json``, which records both
#: the measurement and how to repeat it. A port that corrected its upstream would
#: disagree with it by 1.4e-4 forever. The spec's ``notes`` carries the argument in
#: full.
OMEGA_A = 0.45724333333

#: ``Omega_b``, the repulsion constant of the Peng-Robinson cubic.
#:
#: NeqSim's value, for the reason :data:`OMEGA_A` gives: ``ComponentPR`` sets
#: ``b = .077803333 * R * Tc / Pc``, against the cubic's exact
#: ``0.07779607390388846`` and the paper's printed ``0.07780``. It is 7.26e-6 above
#: the exact value and 3.3333e-6 above the printed one - the same offset
#: ``Omega_a`` carries.
OMEGA_B = 0.077803333


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
