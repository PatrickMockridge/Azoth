"""``eos.liquid_heat_capacity`` - NeqSim's pure-component liquid heat-capacity
polynomial.

```text
cp = 1e-3*(c0 + c1*T + c2*T**2 + c3*T**3 + c4*T**4)
```

Spec: ``specs/calcs/eos/liquid_heat_capacity.toml``.
"""

from __future__ import annotations

from azoth._registry_gen import spec as _spec_for
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import LiquidHeatCapacityResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning

CALC_ID = "eos.liquid_heat_capacity"


def liquid_heat_capacity(
    c0: float, c1: float, c2: float, c3: float, c4: float, T: Q
) -> LiquidHeatCapacityResult:
    """The pure-component liquid heat capacity at a temperature, from NeqSim's polynomial.

    ``c0..c4`` are the raw ``cpliquid1..5``, caller-supplied.

    Args:
        c0..c4: the five coefficients, NeqSim's internal scale.
        T: absolute temperature.

    Returns:
        The liquid heat capacity, in ``J/(mol*K)``.

    Raises:
        OutOfRangeError: if ``T`` is not positive.

    Example:
        >>> import azoth
        >>> q = azoth.ureg.Quantity
        >>> r = liquid_heat_capacity(276370.0, -2090.1, 8.125, -0.014116, 9.37e-6, q(300.0, "K"))
        >>> round(r.cp.magnitude, 9)
        75.355
    """
    spec = _spec_for(CALC_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    values = {
        "c0": c0,
        "c1": c1,
        "c2": c2,
        "c3": c3,
        "c4": c4,
        "T": input_to_si(spec, "T", T),
    }

    apply_checks(checks.on_input, values.get, warnings)

    t = values["T"]
    cp = 1e-3 * (c0 + c1 * t + c2 * t * t + c3 * t * t * t + c4 * t * t * t * t)

    apply_checks(checks.derived, lambda name: cp if name == "cp" else None, warnings)

    return LiquidHeatCapacityResult(cp=from_si(cp, "J/(mol*K)"), warnings=tuple(warnings))
