"""``eos.wilson_activity_coefficients`` - the activity coefficients from the
paraffin-wax Wilson model.

Spec: ``specs/models/eos/wilson_activity_coefficients.toml``. A *direct* model: no
iteration, so no algorithm block.
"""

from __future__ import annotations

import math
from collections.abc import Sequence

from azoth import _models_gen
from azoth.core.errors import InvalidInputError, PropertyUnavailableError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import WilsonActivityCoefficientsResult
from azoth.core.units import Q, input_to_si
from azoth.core.warnings import Warning
from azoth.eos.mixture import Mixture

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
    mixture: Mixture,
    T: Q,
    x: Sequence[float],
) -> WilsonActivityCoefficientsResult:
    """The activity coefficients of a mixture, from the paraffin-wax Wilson model.

    ``mixture`` carries each component's molar mass (kg/mol) and critical temperature
    (K), both of which feed the ``lambda_i`` correlation; ``Lambda_ij`` is ``1.0`` for
    ``i == j`` or the heavier first component, else
    ``exp(-(lambda_j - lambda_i) / (R T))``. The activity coefficient is
    ``ln gamma_c = 1 - ln(sum_i x_i Lambda_c_i) - sum_i x_i Lambda_i_c /
    sum_j x_j Lambda_i_j``. ``x`` is checked rather than renormalised; a supercritical
    component yields NaN, as in NeqSim.

    Args:
        mixture: the components and their interaction parameters.
        T: absolute temperature.
        x: mole fractions; non-negative and summing to one.

    Returns:
        The natural logarithm and the value of each activity coefficient.

    Raises:
        PropertyUnavailableError: if a component carries no molar mass.
        InvalidInputError: if ``x`` is not a composition of ``mixture``.
        OutOfRangeError: if ``T`` is not positive.

    Example:
        >>> import azoth
        >>> from azoth.eos.components import mixture_of
        >>> q = azoth.ureg.Quantity
        >>> mixture = mixture_of(["n-butane", "nc12"])[0]
        >>> r = wilson_activity_coefficients(mixture, q(298.15, "K"), [0.5, 0.5])
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

    n = len(x)
    if n == 0:
        raise InvalidInputError("x", "a mixture of zero components has no activity coefficient")
    if len(mixture) != n:
        raise InvalidInputError(
            "components",
            f"the mixture has {len(mixture)} components but `x` has {n} entries",
        )
    # The correlation is over carbon number and the fusion temperature, both of which
    # are masses; a component without one is refused rather than given a zero, which
    # would put a zero on the wrong side of the `M_i > M_j` branch as well as in the
    # energy.
    m: list[float] = []
    for component in mixture.components:
        if component.molar_mass is None:
            raise PropertyUnavailableError(
                "component",
                "molar mass",
                "the paraffin-wax Wilson correlation needs a molar mass for the carbon "
                "number, and a card-added component carries none",
            )
        m.append(component.molar_mass.to_base_units().magnitude)
    tc = [component.Tc.to_base_units().magnitude for component in mixture.components]
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
