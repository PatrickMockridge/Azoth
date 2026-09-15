"""``eos.pr78_kappa`` - the Peng-Robinson (1978) attraction-parameter coefficient.

```text
kappa = 0.379642 + 1.48503*omega - 0.164423*omega**2 + 0.01666*omega**3   if omega > 0.49
kappa = 0.37464 + 1.54226*omega - 0.26992*omega**2                        otherwise
```

Spec: ``specs/calcs/eos/pr78_kappa.toml``. The 1978 revision of the Peng-Robinson
coefficient: the 1976 form for light components, a heavier-acentric branch for
``omega`` above 0.49.
"""

from __future__ import annotations

from azoth._registry_gen import spec as _spec_for
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import Pr78KappaResult
from azoth.core.warnings import Warning

CALC_ID = "eos.pr78_kappa"


def pr78_kappa(omega: float) -> Pr78KappaResult:
    """The 1978 Peng-Robinson alpha-function coefficient for a pure component.

    Args:
        omega: the Pitzer acentric factor.

    Example:
        >>> r = pr78_kappa(0.6)
        >>> round(r.kappa, 11)
        1.21506628
        >>> pr78_kappa(0.152).kappa == 0.60282728832
        True
    """
    spec = _spec_for(CALC_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    values = {"omega": omega}

    apply_checks(checks.on_input, values.get, warnings)

    if omega > 0.49:
        kappa = (
            0.379642 + 1.48503 * omega - 0.164423 * omega * omega + 0.01666 * omega * omega * omega
        )
    else:
        kappa = 0.37464 + 1.54226 * omega - 0.26992 * omega * omega

    apply_checks(checks.derived, lambda name: kappa if name == "kappa" else None, warnings)

    return Pr78KappaResult(kappa=kappa, warnings=tuple(warnings))
