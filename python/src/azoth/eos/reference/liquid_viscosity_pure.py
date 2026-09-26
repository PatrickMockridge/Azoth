"""``eos.liquid_viscosity_pure`` - one component's pure-liquid viscosity, from the LIQVISC
correlation.

Spec: ``specs/calcs/eos/liquid_viscosity_pure.toml``, which records the four polynomial forms, the
pressure correction and **the two NeqSim classes whose ladders disagree**.
"""

from __future__ import annotations

import math
from typing import Any

from azoth.core.errors import InvalidInputError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import LiquidViscosityPureResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning

CALC_ID = "eos.liquid_viscosity_pure"

#: The value NeqSim answers above a component's critical temperature, in cP.
ABOVE_CRITICAL_CP = 0.5
#: The value it answers for a component with no LIQVISC model, in cP.
NO_MODEL_CP = 0.7
#: NeqSim's value is in cP; the port reports Pa*s.
CP_TO_PA_S = 1.0e-3

LADDERS = ("common_phase", "liquid")


def pressure_correction(
    temperature: float, pressure: float, tc: float, pc: float, omega: float
) -> float:
    """``getViscosityPressureCorrection``: NeqSim's own four-coefficient form."""
    reduced_t = temperature / tc
    if reduced_t > 1.0:
        return 1.0
    delta_pr = pressure / pc
    a = 0.9991 - (4.674e-4 / (1.0523 * math.pow(reduced_t, -0.03877) - 1.0513))
    d = (0.3257 / math.pow(1.0039 - reduced_t**2.573, 0.2906)) - 0.2086
    c = (
        -0.07921
        + 2.1616 * reduced_t
        - 13.4040 * reduced_t**2
        + 44.1706 * reduced_t**3
        - 84.8291 * reduced_t**4
        + 96.1209 * reduced_t**5
        - 59.8127 * reduced_t**6
        + 15.6719 * reduced_t**7
    )
    return (1.0 + d * math.pow(delta_pr / 2.118, a)) / (1.0 + c * omega * delta_pr)


def liquid_viscosity_pure(
    form: str,
    model: int,
    l1: float,
    l2: float,
    l3: float,
    l4: float,
    Tc: Q,
    Pc: Q,
    omega: float,
    T: Q,
    P: Q,
) -> LiquidViscosityPureResult:
    """One component's pure-liquid viscosity, from its LIQVISC correlation.

    ``form`` picks which of NeqSim's two ladders runs: ``common_phase``, whose LIQVISC model 2
    branch is empty so the component comes out at zero, or ``liquid``, which implements it.
    A phase gets whichever class its own viscosity model extends, and both are on the rate-based
    column's path.

    Args:
        form: ``common_phase`` or ``liquid``.
        model: the LIQVISC model index, 1 to 4; anything else answers the 0.7 cP sentinel.
        l1: ``LIQVISC1``.
        l2: ``LIQVISC2``.
        l3: ``LIQVISC3``, read by model 3 only.
        l4: ``LIQVISC4``, read by model 3 only.
        Tc: the component's critical temperature.
        Pc: its critical pressure.
        omega: its acentric factor.
        T: absolute temperature.
        P: absolute pressure.

    Returns:
        The component's pure-liquid viscosity, in ``Pa*s``.

    Raises:
        InvalidInputError: for a ``form`` that names neither ladder.
        OutOfRangeError: if ``T``, ``Tc`` or ``Pc`` is not positive.

    Example:
        >>> import azoth
        >>> q = azoth.ureg.Quantity
        >>> r = liquid_viscosity_pure(
        ...     "liquid", 3, -27.952757828, 4665.22592993, 0.052323342, -3.8356e-5,
        ...     q(647.3, "K"), q(22089000.0, "Pa"), 0.344,
        ...     q(313.15, "K"), q(50.0e5, "Pa"))
        >>> round(r.mu.magnitude, 18)
        0.0006528922494381046
    """
    if form not in LADDERS:
        raise InvalidInputError(
            "form",
            f"`{form}` is not one of the two ladders: {', '.join(LADDERS)}",
        )

    from azoth._registry_gen import spec as _spec_for

    spec: dict[str, Any] = _spec_for(CALC_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    values = {
        "Tc": input_to_si(spec, "Tc", Tc),
        "Pc": input_to_si(spec, "Pc", Pc),
        "T": input_to_si(spec, "T", T),
        "P": input_to_si(spec, "P", P),
    }
    apply_checks(checks.on_input, values.get, warnings)

    temperature = values["T"]
    if temperature > values["Tc"]:
        uncorrected = ABOVE_CRITICAL_CP
    elif model == 1:
        uncorrected = l1 * temperature**l2
    elif model == 2:
        # **The one branch the two ladders disagree on.**
        uncorrected = math.exp(l1 + l2 / temperature) if form == "liquid" else 0.0
    elif model == 3:
        uncorrected = math.exp(l1 + l2 / temperature + l3 * temperature + l4 * temperature**2)
    elif model == 4:
        uncorrected = 10 ** (l1 * (1.0 / temperature - 1.0 / l2))
    else:
        uncorrected = NO_MODEL_CP

    mu_cp = (
        uncorrected
        * (pressure_correction(temperature, values["P"], values["Tc"], values["Pc"], omega) + 1.0)
        / 2.0
    )
    mu = mu_cp * CP_TO_PA_S

    apply_checks(checks.derived, lambda name: mu if name == "mu" else None, warnings)

    return LiquidViscosityPureResult(mu=from_si(mu, "Pa*s"), warnings=tuple(warnings))
