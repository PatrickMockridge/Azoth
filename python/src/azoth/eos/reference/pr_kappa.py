"""``eos.pr_kappa`` - the Peng-Robinson attraction-parameter coefficient.

```text
kappa = 0.37464 + 1.54226*omega - 0.26992*omega**2
```

Peng, D. Y.; Robinson, D. B. (1976). "A New Two-Constant Equation of State."
Ind. Eng. Chem. Fundam. 15(1), 59-64. DOI 10.1021/i160057a011

Spec: ``specs/calcs/eos/pr_kappa.yaml``

# Where this sits in the equation of state

``kappa`` is the whole temperature dependence of the Peng-Robinson attraction term,
in one number:

```text
alpha(T) = (1 + kappa*(1 - sqrt(Tr)))**2
a(T)     = 0.45724 * R**2 * Tc**2 / Pc * alpha(T)
```

It is a property of the substance alone - no temperature, no pressure - which is
what makes it worth a calculation of its own rather than a line inside one. Nothing
downstream can catch a transposed digit in 0.26992: the Z factor, the fugacity
coefficient and the phase split it feeds all still converge and all still pass
their own consistency checks, and all are slightly wrong.

# The sign of kappa

``kappa`` is negative for ``omega`` below -0.23338349942403008, which the quadratic
term makes reachable - helium is at -0.385. A negative coefficient makes ``alpha``
*grow* with temperature, which is not a statement about any fluid. The value is
returned rather than refused, carrying ``OUT_OF_VALID_RANGE``: the arithmetic is
well defined and a caller inspecting the limit deliberately should not be stopped.
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
        omega: the Pitzer acentric factor. A plain float, not a ``pint`` quantity,
            because it is genuinely dimensionless - the same rule that makes
            ``Re``, ``f`` and ``epsilon/D`` plain floats in the hydraulics calcs.

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

    # Nothing is converted here because nothing carries a unit. Both sides of this
    # calc are dimensionless, so there is no `to_si` on the way in and no `from_si`
    # on the way out - which is what writing an equation of state in reduced
    # variables buys at the boundary, and the reason this namespace has no unit
    # handling to get wrong.
    return PrKappaResult(kappa=kappa, warnings=tuple(warnings))
