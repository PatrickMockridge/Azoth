"""``eos.prsv_kappa`` - the Peng-Robinson-Stryjek-Vera alpha-function coefficient.

```text
kappa = 0.378893 + 1.4897153*omega - 0.17131848*omega**2 + 0.0196554*omega**3
      + kappa1*(1 + Tr**0.5)*(0.7 - Tr)
```

Spec: ``specs/calcs/eos/prsv_kappa.toml``, which carries the provenance, the
parameter ``kappa1`` this library does not ship, and why no bound asserts the fitted
range.

PRSV changes the temperature dependence of the attraction coefficient and nothing
else, so :func:`azoth.eos.pr_alpha_ab` serves it unchanged. Unlike Peng-Robinson's
coefficient this one varies with ``Tr``, so it has to be recomputed at each
temperature rather than passed as a fixed number.
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
    # `0.7 - Tr` is zero at Tr = 0.7 whatever `kappa1` is, which is the reduced
    # temperature the acentric factor is defined at - the correlation's anchor.
    kappa = acentric_only + kappa1 * (1.0 + Tr**0.5) * (0.7 - Tr)

    apply_checks(checks.derived, lambda name: kappa if name == "kappa" else None, warnings)

    return PrsvKappaResult(kappa=kappa, warnings=tuple(warnings))
