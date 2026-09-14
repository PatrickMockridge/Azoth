"""``eos.pr_kappa`` - the Peng-Robinson attraction-parameter coefficient.

```text
kappa = 0.37464 + 1.54226*omega - 0.26992*omega**2
```

Spec: ``specs/calcs/eos/pr_kappa.yaml``, which carries the citation, what a
transposed digit in the polynomial costs, and why a negative ``kappa`` is returned
with a warning rather than refused.

``kappa`` is the whole temperature dependence of the Peng-Robinson attraction term, in
one number:

```text
alpha(T) = (1 + kappa*(1 - sqrt(Tr)))**2
a(T)     = 0.45724 * R**2 * Tc**2 / Pc * alpha(T)
```

It is a property of the substance alone - no temperature, no pressure - which is what
makes it a calculation of its own rather than a line inside one.
"""

from __future__ import annotations

from azoth._registry_gen import spec as _spec_for
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import PrKappaResult
from azoth.core.warnings import Warning

CALC_ID = "eos.pr_kappa"


def pr_kappa(omega: float) -> PrKappaResult:
    """The Peng-Robinson alpha-function coefficient for a pure component.

    Args:
        omega: the Pitzer acentric factor. Dimensionless, and supplied by the
            caller - this library ships no component databank.

    Returns:
        The coefficient, and any caveats. A negative ``kappa`` is a real answer to
        a real question and comes back with an ``OUT_OF_VALID_RANGE`` warning
        rather than as an error.

    Example:
        >>> r = pr_kappa(0.152)
        >>> round(r.kappa, 11)
        0.60282728832
        >>> r.is_clean
        True
    """
    spec = _spec_for(CALC_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    values = {"omega": omega}

    apply_checks(checks.on_input, values.get, warnings)

    kappa = 0.37464 + 1.54226 * omega - 0.26992 * omega * omega

    apply_checks(checks.derived, lambda name: kappa if name == "kappa" else None, warnings)

    return PrKappaResult(kappa=kappa, warnings=tuple(warnings))
