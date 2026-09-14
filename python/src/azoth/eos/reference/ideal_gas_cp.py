"""``eos.ideal_gas_cp`` - ideal-gas heat capacity from a four-term polynomial.

```text
Cp/R = a + b*theta + c*theta**2 + d*theta**3,   theta = T / (1000 K)
```

Spec: ``specs/calcs/eos/ideal_gas_cp.yaml``, which carries why the four coefficients
are the caller's, the derivation of the ``T/(1000 K)`` substitution, and the warning
a non-positive ``Cp`` carries.
"""

from __future__ import annotations

from typing import Any

from azoth.core.range import apply_checks, checks_for
from azoth.core.result import IdealGasCpResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning
from azoth.eos.reference.pr_molar_volume import MOLAR_GAS_CONSTANT

MODEL_ID = "eos.ideal_gas_cp"

#: The reference temperature the polynomial is written against, in kelvin.
#:
#: A stated constant of the correlation's form rather than a fitted quantity: it is
#: what makes a table's four printed numbers dimensionless, and a caller whose
#: coefficients are quoted against another reference rescales them once.
REFERENCE_TEMPERATURE = 1000.0


def ideal_gas_cp(a: float, b: float, c: float, d: float, T: Q) -> IdealGasCpResult:
    """The ideal-gas heat capacity at a temperature.

    The coefficients are dimensionless and describe ``Cp/R``; the result is
    dimensioned.

    Args:
        a: the constant term of ``Cp/R``.
        b: the coefficient of ``theta``.
        c: the coefficient of ``theta**2``.
        d: the coefficient of ``theta**3``.
        T: absolute temperature. Must lie inside the range the coefficients were
            fitted over, which is **not checked**.

    Returns:
        ``cp_over_r`` and ``cp``. A non-positive ``cp`` comes back carrying
        ``OUT_OF_VALID_RANGE`` rather than raising: it means the polynomial has been
        evaluated outside its fitted range, and the arithmetic is well defined.

    Raises:
        OutOfRangeError: if ``T`` is not positive.

    Example:
        >>> import azoth
        >>> q = azoth.ureg.Quantity
        >>> r = ideal_gas_cp(4.0, 1.0, -0.5, 0.1, q(500.0, "K"))
        >>> round(r.cp_over_r, 4)
        4.3875
    """
    spec = _spec()
    checks = checks_for(spec)
    warnings: list[Warning] = []

    values = {"a": a, "b": b, "c": c, "d": d, "T": input_to_si(spec, "T", T)}
    apply_checks(checks.on_input, values.get, warnings)

    theta = values["T"] / REFERENCE_TEMPERATURE
    cp_over_r = a + b * theta + c * theta * theta + d * theta * theta * theta
    cp = cp_over_r * MOLAR_GAS_CONSTANT

    apply_checks(
        checks.derived,
        lambda name: cp if name == "cp" else None,
        warnings,
    )

    return IdealGasCpResult(
        cp_over_r=cp_over_r,
        cp=from_si(cp, "J/(mol*K)"),
        warnings=tuple(warnings),
    )


def _spec() -> dict[str, Any]:
    from azoth._registry_gen import spec

    return spec(MODEL_ID)
