"""``eos.unifac_activity_coefficients`` - the activity coefficients from the UNIFAC
group-contribution model.

Spec: ``specs/models/eos/unifac_activity_coefficients.toml``. A *direct* model: no
iteration, so no algorithm block.
"""

from __future__ import annotations

import math
from collections.abc import Sequence

from azoth import _models_gen
from azoth.core.errors import InvalidInputError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import UnifacActivityCoefficientsResult
from azoth.core.units import Q, input_to_si
from azoth.core.warnings import Warning

MODEL_ID = "eos.unifac_activity_coefficients"


def _ln_gamma_group(
    k: int, theta: list[float], group_q: list[float], aij: list[list[float]], t: float
) -> float:
    g = len(group_q)
    s1 = sum(theta[m] * math.exp(-aij[m][k] / t) for m in range(g))
    s3 = 0.0
    for m in range(g):
        s2 = sum(theta[n] * math.exp(-aij[n][m] / t) for n in range(g))
        s3 += theta[m] * math.exp(-aij[k][m] / t) / s2
    return group_q[k] * (1.0 - math.log(s1) - s3)


def unifac_activity_coefficients(
    T: Q,
    x: Sequence[float],
    groups: Sequence[Sequence[float]],
    group_r: Sequence[float],
    group_q: Sequence[float],
    aij: Sequence[Sequence[Q]],
) -> UnifacActivityCoefficientsResult:
    """The activity coefficients of a mixture, from UNIFAC.

    ``groups[i][k]`` is the count of group ``k`` in component ``i`` over the union of
    the named components' groups, ``group_r``/``group_q`` are the per-group volume and
    surface area, and ``aij[m][n] = a_{main(m), main(n)}`` in Kelvin. The component
    volume and area follow from the group sums, the combinatorial term uses ``Z = 10``,
    and the residual is the standard ``ln gamma^R_i = sum_k nu_k^i (ln Gamma_k^mix -
    ln Gamma_k^pure)``. ``x`` is checked rather than renormalised.

    Args:
        T: absolute temperature.
        x: mole fractions; non-negative and summing to one.
        groups: the group counts, one row per component.
        group_r: the volume ``R`` of each group.
        group_q: the surface area ``Q`` of each group.
        aij: the main-group interaction matrix, in Kelvin.

    Returns:
        The natural logarithm and the value of each activity coefficient.

    Raises:
        InvalidInputError: if the shapes disagree, a group count is negative, or ``x``
            is not a composition.
        OutOfRangeError: if ``T`` is not positive.

    Example:
        >>> import azoth
        >>> q = azoth.ureg.Quantity
        >>> aij = [[q(0.0, "K"), q(-181.0, "K")], [q(289.6, "K"), q(0.0, "K")]]
        >>> r = unifac_activity_coefficients(
        ...     q(298.15, "K"), [0.5, 0.5], [[1.0, 0.0], [0.0, 1.0]],
        ...     [1.4311, 0.92], [1.432, 1.4], aij,
        ... )
        >>> round(r.gamma[0], 14)
        1.1156815062468
    """
    spec = _models_gen.model(MODEL_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    values = {"T": input_to_si(spec, "T", T)}
    apply_checks(checks.on_input, values.get, warnings)

    t = input_to_si(spec, "T", T)
    x = list(x)
    groups = [list(row) for row in groups]
    group_r = list(group_r)
    group_q = list(group_q)
    aij_si = [[input_to_si(spec, "aij", value) for value in row] for row in aij]

    n = len(x)
    g = len(group_r)
    if n == 0:
        raise InvalidInputError("x", "a mixture of zero components has no activity coefficient")
    if g == 0 or len(group_q) != g:
        raise InvalidInputError(
            "group_r",
            f"`group_r` has {g} entries and `group_q` has {len(group_q)}; both must be "
            "the same non-empty length",
        )
    if len(groups) != n:
        raise InvalidInputError(
            "groups", f"`groups` has {len(groups)} rows but `x` has {n} entries; it must be N x G"
        )
    for i in range(n):
        if len(groups[i]) != g:
            raise InvalidInputError(
                "groups", f"row {i} of `groups` has {len(groups[i])} entries but must be {g}"
            )
        if any(count < 0.0 for count in groups[i]):
            raise InvalidInputError("groups", f"row {i} of `groups` has a negative group count")
    if len(aij_si) != g:
        raise InvalidInputError("aij", f"`aij` has {len(aij_si)} rows but must be {g} x {g}")
    for i in range(g):
        if len(aij_si[i]) != g:
            raise InvalidInputError(
                "aij", f"row {i} of `aij` has {len(aij_si[i])} entries but must be {g}"
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

    ri = [sum(groups[i][k] * group_r[k] for k in range(g)) for i in range(n)]
    qi = [sum(groups[i][k] * group_q[k] for k in range(g)) for i in range(n)]

    ln_gamma: list[float] = []
    gamma: list[float] = []
    for i in range(n):
        t1 = sum(x[j] * ri[j] for j in range(n))
        t2 = sum(x[j] * qi[j] for j in range(n))
        suml = sum(x[j] * (5.0 * (ri[j] - qi[j]) - (ri[j] - 1.0)) for j in range(n))
        v = x[i] * ri[i] / t1
        f = x[i] * qi[i] / t2
        li = 5.0 * (ri[i] - qi[i]) - (ri[i] - 1.0)
        lng_c = math.log(v / x[i]) + 5.0 * qi[i] * math.log(f / v) + li - (v / x[i]) * suml

        denom = sum(x[j] * qi[j] for j in range(n))
        qmix = [group_q[l] * sum(x[j] * groups[j][l] for j in range(n)) / denom for l in range(g)]
        qcomp = [group_q[l] * groups[i][l] / qi[i] for l in range(g)]
        lng_r = sum(
            groups[i][k]
            * (
                _ln_gamma_group(k, qmix, group_q, aij_si, t)
                - _ln_gamma_group(k, qcomp, group_q, aij_si, t)
            )
            for k in range(g)
        )

        lng = lng_c + lng_r
        ln_gamma.append(lng)
        gamma.append(math.exp(lng))

    return UnifacActivityCoefficientsResult(
        ln_gamma=tuple(ln_gamma),
        gamma=tuple(gamma),
        warnings=tuple(warnings),
    )
