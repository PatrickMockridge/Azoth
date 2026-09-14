"""``eos.ideal_gas_cp`` - ideal-gas heat capacity from a polynomial.

```text
cp = cp_a + cp_b*T + cp_c*T**2 + cp_d*T**3 + cp_e*T**4
```

A port of ``neqsim.thermo.component.Component.getCp0(double)``, and dimensional
throughout: the coefficients carry the powers of temperature in their units, so
``cp_a`` is a heat capacity and ``cp_b`` is one per kelvin. That is what lets the five
``CPA``-``CPE`` columns of NeqSim's ``COMP.csv`` be read as they are stored.

Spec: ``specs/calcs/eos/ideal_gas_cp.toml``.
"""

from __future__ import annotations

from typing import Any

from azoth.core.range import apply_checks, checks_for
from azoth.core.result import IdealGasCpResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning

MODEL_ID = "eos.ideal_gas_cp"


def ideal_gas_cp(
    cp_a: float, cp_b: float, cp_c: float, cp_d: float, cp_e: float, T: Q
) -> IdealGasCpResult:
    """The ideal-gas heat capacity at a temperature.

    Args:
        cp_a: the constant term, in J/(mol*K).
        cp_b: the coefficient of ``T``, in J/(mol*K**2).
        cp_c: the coefficient of ``T**2``, in J/(mol*K**3).
        cp_d: the coefficient of ``T**3``, in J/(mol*K**4).
        cp_e: the coefficient of ``T**4``, in J/(mol*K**5). Zero for a component whose
            heat capacity is described by four terms.
        T: absolute temperature.

    Returns:
        ``cp``, the ideal-gas heat capacity in J/(mol*K), carrying a warning when it
        comes out non-positive.

    Raises:
        OutOfRangeError: if ``T`` is not positive.

    Example:
        >>> import azoth
        >>> q = azoth.ureg.Quantity
        >>> r = ideal_gas_cp(
        ...     37.978352, -0.07461815, 0.000301881, -2.83e-07, 9.070574e-11, q(300.0, "K")
        ... )
        >>> round(r.cp.to("J/(mol*K)").magnitude, 6)
        35.855913
    """
    spec = _spec()
    checks = checks_for(spec)
    warnings: list[Warning] = []

    # Every argument becomes an SI magnitude before any arithmetic. The coefficients
    # carry units of their own - a heat capacity and one per kelvin per degree - so
    # leaving them as quantities and mixing them with a plain temperature leaves each
    # term with a different power of kelvin and the sum does not add up.
    values = {
        "cp_a": input_to_si(spec, "cp_a", cp_a),
        "cp_b": input_to_si(spec, "cp_b", cp_b),
        "cp_c": input_to_si(spec, "cp_c", cp_c),
        "cp_d": input_to_si(spec, "cp_d", cp_d),
        "cp_e": input_to_si(spec, "cp_e", cp_e),
        "T": input_to_si(spec, "T", T),
    }
    apply_checks(checks.on_input, values.get, warnings)

    # Abbreviated multiplication rather than `**`, because these are the terms as the
    # source writes them and a reader checking one against the other should not have to
    # expand an operator.
    (a, b, c, d, e) = (
        values["cp_a"],
        values["cp_b"],
        values["cp_c"],
        values["cp_d"],
        values["cp_e"],
    )
    temperature = values["T"]
    cp = (
        a
        + b * temperature
        + c * temperature * temperature
        + d * temperature * temperature * temperature
        + e * temperature * temperature * temperature * temperature
    )

    apply_checks(
        checks.derived,
        lambda name: cp if name == "cp" else None,
        warnings,
    )

    return IdealGasCpResult(
        cp=from_si(cp, "J/(mol*K)"),
        warnings=tuple(warnings),
    )


def _spec() -> dict[str, Any]:
    from azoth._registry_gen import spec

    return spec(MODEL_ID)
