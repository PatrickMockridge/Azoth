"""``eos.wilson_activity_coefficients`` - the activity coefficients from the
paraffin-wax Wilson model.

Spec: ``specs/models/eos/wilson_activity_coefficients.toml``. A *direct* model: no
iteration, so no algorithm block.
"""

from __future__ import annotations

import math
from collections.abc import Sequence

from azoth import _models_gen
from azoth.core.errors import InvalidInputError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import WilsonActivityCoefficientsResult
from azoth.core.units import Q, input_to_si
from azoth.core.warnings import Warning

MODEL_ID = "eos.wilson_activity_coefficients"

R = 8.3144621


def _interaction_energy(m: float, tc: float, t: float) -> float:
    coordination = 6.0
    carbon = m / 0.014
    x = 1.0 - t / tc
    d0 = (
        5.2804 * math.pow(x, 0.3333)
        + 12.865 * math.pow(x, 0.8333)
        + 1.171 * math.pow(x, 1.2083)
        - 13.166 * x
        + 0.4858 * x**2
        - 1.088 * x**3
    )
    d1 = (
        0.80022 * math.pow(x, 0.3333)
        + 273.23 * math.pow(x, 0.8333)
        + 465.08 * math.pow(x, 1.2083)
        - 638.51 * x
        - 145.12 * x**2
        - 74.049 * x**3
    )
    d2 = (
        7.2543 * math.pow(x, 0.3333)
        - 346.45 * math.pow(x, 0.8333)
        - 610.48 * math.pow(x, 1.2083)
        + 839.89 * x
        + 160.05 * x**2
        - 50.711 * x**3
    )
    omega = 0.0520750 + 0.0448946 * carbon - 0.000185397 * carbon * carbon
    dh_vap = R * tc * (d0 + omega * d1 + omega * omega * d2) * 4.1868
    dh_tot = (3.7791 * carbon - 12.654) * 1000.0
    tf = 374.5 + 0.2617 * m - 20.172 / m
    dh_f = 0.1426 * m * tf * 4.1868
    dh_sub = dh_vap + dh_f + (dh_tot - dh_f)
    return -2.0 / coordination * (dh_sub - R * t)


def _char_energy(m: list[float], tc: list[float], t: float, i: int, j: int) -> float:
    if i == j or m[i] > m[j]:
        return 1.0
    li = _interaction_energy(m[i], tc[i], t)
    lj = _interaction_energy(m[j], tc[j], t)
    return math.exp(-(lj - li) / (R * t))


def wilson_activity_coefficients(
    T: Q,
    x: Sequence[float],
    M: Sequence[Q],
    Tc: Sequence[Q],
) -> WilsonActivityCoefficientsResult:
    """The activity coefficients of a mixture, from the paraffin-wax Wilson model.

    ``M`` is molar mass in kg/mol and ``Tc`` the critical temperature in K; both feed
    the ``lambda_i`` correlation, and ``Lambda_ij`` is ``1.0`` for ``i == j`` or the
    heavier first component, else ``exp(-(lambda_j - lambda_i) / (R T))``. The activity
    coefficient is ``ln gamma_c = 1 - ln(sum_i x_i Lambda_c_i) - sum_i x_i Lambda_i_c /
    sum_j x_j Lambda_i_j``. ``x`` is checked rather than renormalised; a supercritical
    component yields NaN, as in NeqSim.

    Args:
        T: absolute temperature.
        x: mole fractions; non-negative and summing to one.
        M: molar mass of each component, in kg/mol.
        Tc: critical temperature of each component.

    Returns:
        The natural logarithm and the value of each activity coefficient.

    Raises:
        InvalidInputError: if the vectors disagree in length, or ``x`` is not a
            composition.
        OutOfRangeError: if ``T`` is not positive.

    Example:
        >>> import azoth
        >>> q = azoth.ureg.Quantity
        >>> m = [q(0.058123, "kg/mol"), q(0.1703, "kg/mol")]
        >>> tc = [q(425.12, "K"), q(658.0, "K")]
        >>> r = wilson_activity_coefficients(q(298.15, "K"), [0.5, 0.5], m, tc)
        >>> round(r.gamma[0], 14)
        1.21306131942502
    """
    spec = _models_gen.model(MODEL_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    values = {"T": input_to_si(spec, "T", T)}
    apply_checks(checks.on_input, values.get, warnings)

    t = input_to_si(spec, "T", T)
    x = list(x)
    m = [input_to_si(spec, "M", value) for value in M]
    tc = [input_to_si(spec, "Tc", value) for value in Tc]

    n = len(x)
    if n == 0:
        raise InvalidInputError("x", "a mixture of zero components has no activity coefficient")
    if len(m) != n or len(tc) != n:
        raise InvalidInputError(
            "M",
            f"the per-component vectors disagree in length: `x` has {n} entries, `M` "
            f"{len(m)} and `Tc` {len(tc)}",
        )
    if any(value < 0.0 for value in x):
        bad = next(i for i, value in enumerate(x) if value < 0.0)
        raise InvalidInputError("x", f"x[{bad}] is {x[bad]} but a mole fraction cannot be negative")
    if abs(sum(x) - 1.0) > 1.0e-09:
        raise InvalidInputError(
            "x",
            f"the mole fractions sum to {sum(x)}, not to one. Renormalising them here "
            "would make a composition error invisible in every number downstream, so "
            "it is refused instead",
        )

    ln_gamma: list[float] = []
    gamma: list[float] = []
    for c in range(n):
        s1 = sum(x[i] * _char_energy(m, tc, t, c, i) for i in range(n))
        s2 = 0.0
        for i in range(n):
            temp = sum(x[j] * _char_energy(m, tc, t, i, j) for j in range(n))
            s2 += x[i] * _char_energy(m, tc, t, i, c) / temp
        lng = 1.0 - math.log(s1) - s2
        ln_gamma.append(lng)
        gamma.append(math.exp(lng))

    return WilsonActivityCoefficientsResult(
        ln_gamma=tuple(ln_gamma),
        gamma=tuple(gamma),
        warnings=tuple(warnings),
    )
