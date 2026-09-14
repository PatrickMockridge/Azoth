"""``eos.rachford_rice_binary`` - the Rachford-Rice vapour fraction, for a binary.

```text
beta = -(z1*(K1 - 1) + z2*(K2 - 1)) / ((K1 - 1)*(K2 - 1))     z2 = 1 - z1
```

Spec: ``specs/calcs/eos/rachford_rice_binary.yaml``, which carries the provenance, the
derivation of the closed form, and why a ``beta`` outside ``[0, 1]`` is a warning
rather than an error.
"""

from __future__ import annotations

from azoth._registry_gen import spec as _spec_for
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import RachfordRiceBinaryResult
from azoth.core.warnings import Warning

CALC_ID = "eos.rachford_rice_binary"


def rachford_rice_binary(z1: float, K1: float, K2: float) -> RachfordRiceBinaryResult:
    """The vapour fraction that solves the Rachford-Rice equation for two components.

    Args:
        z1: overall mole fraction of component 1. Component 2 is ``1 - z1``.
        K1: the K-value of component 1, ``y1/x1``. Normally the caller's, computed
            from fugacity coefficients.
        K2: the K-value of component 2.

    Returns:
        ``beta``, the vapour fraction. Outside ``[0, 1]`` the feed is single phase
        and the result carries an ``OUT_OF_VALID_RANGE`` warning - the value is
        still meaningful, but it is not a vapour fraction.

    Raises:
        OutOfRangeError: if ``z1`` is outside ``[0, 1]``, if either K-value is not
            positive, or if either equals 1 - which makes ``K - 1`` a divisor and
            drops that component out of the sum entirely.

    Example:
        >>> r = rachford_rice_binary(0.6, 4.0, 0.25)
        >>> round(r.beta, 16)
        0.6666666666666665
        >>> r.is_clean
        True
    """
    spec = _spec_for(CALC_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    values = {"z1": z1, "K1": K1, "K2": K2}

    apply_checks(checks.on_input, values.get, warnings)

    # Hoisted but one operation each, so this is the same double as writing
    # `K1 - 1.0` inline four times - which is what the spec's equation shows.
    a = K1 - 1.0
    b = K2 - 1.0
    # Guarded by the `equals: 1` bounds above, so neither divisor is zero.
    beta = -(z1 * a + (1.0 - z1) * b) / (a * b)

    apply_checks(checks.derived, lambda name: beta if name == "beta" else None, warnings)

    return RachfordRiceBinaryResult(beta=beta, warnings=tuple(warnings))
