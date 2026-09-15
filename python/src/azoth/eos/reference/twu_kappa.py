"""``eos.twu_kappa`` - the Twu attraction-parameter coefficient.

```text
kappa = 0.48 + 1.574*omega - 0.175*omega**2
```

Spec: ``specs/calcs/eos/twu_kappa.toml``. Twu's coefficient, which differs from
Soave's only in the last constant - 0.175 against 0.176.
"""

from __future__ import annotations

from azoth._registry_gen import spec as _spec_for
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import TwuKappaResult
from azoth.core.warnings import Warning

CALC_ID = "eos.twu_kappa"


def twu_kappa(omega: float) -> TwuKappaResult:
    """Twu's alpha-function coefficient for a pure component.

    Args:
        omega: the Pitzer acentric factor.

    Example:
        >>> r = twu_kappa(0.152)
        >>> round(r.kappa, 7)
        0.7152048
        >>> r.is_clean
        True
    """
    spec = _spec_for(CALC_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    values = {"omega": omega}

    apply_checks(checks.on_input, values.get, warnings)

    kappa = 0.48 + 1.574 * omega - 0.175 * omega * omega

    apply_checks(checks.derived, lambda name: kappa if name == "kappa" else None, warnings)

    return TwuKappaResult(kappa=kappa, warnings=tuple(warnings))
