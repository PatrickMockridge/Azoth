"""``eos.srk_alpha_ab`` - the Soave-Redlich-Kwong alpha function and reduced parameters.

```text
alpha     = (1 + kappa*(1 - Tr**0.5))**2
a_reduced = Omega_a * alpha * Pr / Tr**2
b_reduced = Omega_b * Pr / Tr
```

Spec: ``specs/calcs/eos/srk_alpha_ab.toml``, which carries the provenance and the two
Omega constants' derivation from the cubic's triple-root condition.
"""

from __future__ import annotations

from azoth._registry_gen import spec as _spec_for
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import SrkAlphaAbResult
from azoth.core.warnings import Warning
from azoth.eos.cubic import SRK

CALC_ID = "eos.srk_alpha_ab"


def srk_alpha_ab(kappa: float, Tr: float, Pr: float) -> SrkAlphaAbResult:
    """The Soave-Redlich-Kwong alpha function and the reduced attraction parameters.

    Args:
        kappa: the alpha-function coefficient, from :func:`azoth.eos.srk_kappa`.
        Tr: reduced temperature, ``T / Tc``.
        Pr: reduced pressure, ``P / Pc``.

    Raises:
        OutOfRangeError: if ``Tr <= 0`` or ``Pr <= 0``.

    Example:
        >>> r = srk_alpha_ab(0.715181696, 0.8, 0.25)
        >>> round(r.a_reduced, 17)
        0.19315231739218255
        >>> round(r.b_reduced, 17)
        0.02707510936404929
    """
    spec = _spec_for(CALC_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    values = {"kappa": kappa, "Tr": Tr, "Pr": Pr}

    apply_checks(checks.on_input, values.get, warnings)

    attraction = 1.0 + kappa * (1.0 - Tr**0.5)
    alpha = attraction**2
    a_reduced = SRK.omega_a * alpha * Pr / Tr**2
    b_reduced = SRK.omega_b * Pr / Tr

    computed = {"alpha": alpha, "a_reduced": a_reduced, "b_reduced": b_reduced}
    apply_checks(checks.derived, computed.get, warnings)

    return SrkAlphaAbResult(
        alpha=alpha,
        a_reduced=a_reduced,
        b_reduced=b_reduced,
        warnings=tuple(warnings),
    )
