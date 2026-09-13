"""``eos.prsv_kappa`` - the Peng-Robinson-Stryjek-Vera alpha-function coefficient.

```text
kappa = 0.378893 + 1.4897153*omega - 0.17131848*omega**2 + 0.0196554*omega**3
      + kappa1*(1 + Tr**0.5)*(0.7 - Tr)
```

Stryjek, R.; Vera, J. H. (1986). "PRSV: An improved Peng-Robinson equation of
state for pure compounds and mixtures." Can. J. Chem. Eng. 64(2), 323-333.
DOI 10.1002/cjce.5450640224

Spec: ``specs/calcs/eos/prsv_kappa.yaml``

# What PRSV changes, and what it does not

It changes the temperature dependence of the attraction coefficient and nothing
else. The alpha function it feeds is Peng-Robinson's, unchanged, which is why
:func:`azoth.eos.pr_alpha_ab` serves both and why this is a coefficient rather than
a second equation of state.

Peng-Robinson's coefficient is a constant per substance. This one is not: the
``(1 + sqrt(Tr))*(0.7 - Tr)`` term varies with temperature, and ``kappa1`` is fitted
per component so the variation matches that component's vapour pressure.

# The anchor at Tr = 0.7

The temperature factor is exactly zero at ``Tr = 0.7``, for any ``kappa1``, because
``0.7 - Tr`` is. So the coefficient collapses to its acentric-only part there - the
temperature at which the acentric factor is defined. PRSV is pinned to the
acentric-only fit at Tr = 0.7 and departs from it elsewhere by whatever ``kappa1``
says.

# The parameter this library does not ship

``kappa1`` is fitted per substance, and a table of fitted parameters is the databank
this library deliberately has none of. It is an input. A caller who expected PRSV's
published accuracy for free will not get it.
"""

from __future__ import annotations

from azoth._registry_gen import spec as _spec_for
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import PrsvKappaResult
from azoth.core.warnings import Warning

CALC_ID = "eos.prsv_kappa"


def prsv_kappa(omega: float, Tr: float, kappa1: float) -> PrsvKappaResult:
    """The PRSV alpha-function coefficient for a pure component.

    Args:
        omega: the acentric factor. Dimensionless.
        Tr: reduced temperature, ``T / Tc``. Dimensionless.
        kappa1: the pure-compound parameter PRSV adds, fitted per substance.
            Supplied by the caller, and this library ships no values for it - see
            the module documentation. Its sign is unconstrained.

    Returns:
        The coefficient. Unlike Peng-Robinson's it varies with temperature, so it
        must be recomputed at each ``Tr`` rather than passed on as a constant.

    Raises:
        OutOfRangeError: if ``Tr <= 0`` - the temperature factor takes a square
            root and a reduced temperature at or below zero is not a state.

    Example:
        >>> at_anchor = prsv_kappa(0.152, 0.7, 0.05)
        >>> round(at_anchor.kappa, 13)
        0.601440609429
        >>> r = prsv_kappa(0.152, 0.8, 0.05)
        >>> round(r.kappa, 16)
        0.5919684734740435
    """
    spec = _spec_for(CALC_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    values = {"omega": omega, "Tr": Tr, "kappa1": kappa1}

    apply_checks(checks.on_input, values.get, warnings)

    # Written in the equation's own order, with `omega * omega` for the square
    # rather than `omega ** 2`. The two agree bit-for-bit for these inputs -
    # measured - but the explicit form is the one the Rust side can match without
    # relying on how each language lowers `**`.
    acentric_only = (
        0.378893
        + 1.4897153 * omega
        - 0.17131848 * (omega * omega)
        + 0.0196554 * (omega * omega * omega)
    )
    kappa = acentric_only + kappa1 * (1.0 + Tr**0.5) * (0.7 - Tr)

    apply_checks(checks.derived, lambda name: kappa if name == "kappa" else None, warnings)

    return PrsvKappaResult(kappa=kappa, warnings=tuple(warnings))
