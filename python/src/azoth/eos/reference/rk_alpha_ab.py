"""``eos.rk_alpha_ab`` - the Redlich-Kwong alpha function and reduced parameters.

```text
alpha     = 1/Tr**0.5
a_reduced = Omega_a * alpha * Pr / Tr**2
b_reduced = Omega_b * Pr / Tr
```

Spec: ``specs/calcs/eos/rk_alpha_ab.toml``. The original Redlich-Kwong temperature
dependence, kappa-free.
"""

from __future__ import annotations

from azoth._registry_gen import spec as _spec_for
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import RkAlphaAbResult
from azoth.core.warnings import Warning
from azoth.eos.cubic import RK

CALC_ID = "eos.rk_alpha_ab"


def rk_alpha_ab(Tr: float, Pr: float) -> RkAlphaAbResult:
    """The Redlich-Kwong alpha function and the reduced attraction parameters.

    Args:
        Tr: reduced temperature, ``T / Tc``.
        Pr: reduced pressure, ``P / Pc``.

    Raises:
        OutOfRangeError: if ``Tr <= 0`` or ``Pr <= 0``.

    Example:
        >>> r = rk_alpha_ab(0.8, 0.25)
        >>> round(r.a_reduced, 17)
        0.1866943088347048
        >>> round(r.b_reduced, 17)
        0.02707510936404929
    """
    spec = _spec_for(CALC_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    values = {"Tr": Tr, "Pr": Pr}

    apply_checks(checks.on_input, values.get, warnings)

    alpha = 1.0 / Tr**0.5
    a_reduced = RK.omega_a * alpha * Pr / Tr**2
    b_reduced = RK.omega_b * Pr / Tr

    computed = {"alpha": alpha, "a_reduced": a_reduced, "b_reduced": b_reduced}
    apply_checks(checks.derived, computed.get, warnings)

    return RkAlphaAbResult(
        alpha=alpha,
        a_reduced=a_reduced,
        b_reduced=b_reduced,
        warnings=tuple(warnings),
    )
