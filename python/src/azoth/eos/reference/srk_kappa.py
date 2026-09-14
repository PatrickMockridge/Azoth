"""``eos.srk_kappa`` - the Soave-Redlich-Kwong attraction-parameter coefficient.

```text
kappa = 0.48 + 1.574*omega - 0.176*omega**2
```

Spec: ``specs/calcs/eos/srk_kappa.toml``, which carries the citation and why a negative
``kappa`` is returned with a warning rather than refused.

Soave's coefficient is the same *form* as Peng-Robinson's with three different
constants, which is why the alpha functions share ``azoth.eos.alpha_term``.
"""

from __future__ import annotations

from azoth._registry_gen import spec as _spec_for
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import SrkKappaResult
from azoth.core.warnings import Warning

CALC_ID = "eos.srk_kappa"


def srk_kappa(omega: float) -> SrkKappaResult:
    """The Soave-Redlich-Kwong alpha-function coefficient for a pure component.

    Args:
        omega: the Pitzer acentric factor. Dimensionless, and supplied by the caller.

    Returns:
        The coefficient, and any caveats. A negative ``kappa`` comes back carrying
        ``OUT_OF_VALID_RANGE`` rather than as an error.

    Example:
        >>> r = srk_kappa(0.152)
        >>> round(r.kappa, 9)
        0.715181696
        >>> r.is_clean
        True
    """
    spec = _spec_for(CALC_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    values = {"omega": omega}

    apply_checks(checks.on_input, values.get, warnings)

    kappa = 0.48 + 1.574 * omega - 0.176 * omega * omega

    apply_checks(checks.derived, lambda name: kappa if name == "kappa" else None, warnings)

    return SrkKappaResult(kappa=kappa, warnings=tuple(warnings))
