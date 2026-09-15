"""``eos.nrtl_activity_coefficients`` - the activity coefficients from the NRTL
local-composition model.

Spec: ``specs/models/eos/nrtl_activity_coefficients.toml``. A *direct* model: no
iteration, so no algorithm block.
"""

from __future__ import annotations

import math
from collections.abc import Sequence

from azoth import _models_gen
from azoth.core.errors import InvalidInputError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import NrtlActivityCoefficientsResult
from azoth.core.units import Q, input_to_si
from azoth.core.warnings import Warning

MODEL_ID = "eos.nrtl_activity_coefficients"


def nrtl_activity_coefficients(
    T: Q,
    x: Sequence[float],
    Dij: Sequence[Sequence[Q]],
    alpha: Sequence[Sequence[float]],
) -> NrtlActivityCoefficientsResult:
    """The activity coefficients of a mixture, from NRTL (Renon-Prausnitz).

    ``Dij[i][j] = g_ij`` is in Kelvin so that ``tau_ij = Dij[i][j] / T`` is
    dimensionless, ``G_ij = exp(-alpha_ij tau_ij)``, and ``ln gamma_i`` is the
    standard symmetric sum: ``(sum_j tau_ji G_ji x_j) / (sum_j G_ji x_j)`` plus
    ``sum_j (x_j G_ij / C_j) (tau_ij - D_j / C_j)``. ``alpha`` is symmetric and both
    matrices have a zero diagonal; ``x`` is checked rather than renormalised.

    Args:
        T: absolute temperature.
        x: mole fractions; non-negative and summing to one.
        Dij: the NRTL energy matrix, in Kelvin.
        alpha: the NRTL non-randomness matrix, symmetric with a zero diagonal.

    Returns:
        The natural logarithm and the value of each activity coefficient.

    Raises:
        InvalidInputError: if the matrices are not ``N x N`` for the ``N`` of ``x``,
            ``alpha`` is not symmetric, a diagonal is not zero, or ``x`` is not a
            composition.
        OutOfRangeError: if ``T`` is not positive.

    Example:
        >>> import azoth
        >>> q = azoth.ureg.Quantity
        >>> t = q(350.0, "K")
        >>> dij = [[q(0.0, "K"), q(-48.68, "K")], [q(610.6, "K"), q(0.0, "K")]]
        >>> alpha = [[0.0, 0.303], [0.303, 0.0]]
        >>> r = nrtl_activity_coefficients(t, [0.5, 0.5], dij, alpha)
        >>> round(r.gamma[0], 14)
        1.22772692900497
    """
    spec = _models_gen.model(MODEL_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    values = {"T": input_to_si(spec, "T", T)}
    apply_checks(checks.on_input, values.get, warnings)

    t = input_to_si(spec, "T", T)
    x = list(x)
    dij = [[input_to_si(spec, "Dij", value) for value in row] for row in Dij]
    alpha = [list(row) for row in alpha]

    n = len(x)
    if n == 0:
        raise InvalidInputError("x", "a mixture of zero components has no activity coefficient")
    if len(dij) != n or len(alpha) != n:
        raise InvalidInputError(
            "Dij",
            f"the matrices are {len(dij)} x {len(alpha)} rows but `x` has {n} entries; "
            "a matrix input must be N x N for the same N as `x`",
        )
    for i in range(n):
        if len(dij[i]) != n or len(alpha[i]) != n:
            raise InvalidInputError(
                "Dij",
                f"row {i} is {len(dij[i])} (Dij) / {len(alpha[i])} (alpha) entries but must be {n}",
            )
    for i in range(n):
        for j in range(n):
            if i == j:
                if dij[i][j] != 0.0 or alpha[i][j] != 0.0:
                    raise InvalidInputError(
                        "Dij",
                        f"the diagonal must be zero, but Dij[{i}][{i}] = {dij[i][j]} and "
                        f"alpha[{i}][{i}] = {alpha[i][j]}",
                    )
            elif alpha[i][j] != alpha[j][i]:
                raise InvalidInputError(
                    "alpha",
                    f"`alpha` must be symmetric, but alpha[{i}][{j}] = {alpha[i][j]} and "
                    f"alpha[{j}][{i}] = {alpha[j][i]}",
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

    tau = [[dij[i][j] / t for j in range(n)] for i in range(n)]
    g = [[math.exp(-alpha[i][j] * tau[i][j]) for j in range(n)] for i in range(n)]

    ln_gamma: list[float] = []
    gamma: list[float] = []
    for i in range(n):
        a = sum(tau[j][i] * g[j][i] * x[j] for j in range(n))
        b = sum(g[j][i] * x[j] for j in range(n))
        f = 0.0
        for j in range(n):
            c = sum(g[l][j] * x[l] for l in range(n))
            d = sum(tau[l][j] * g[l][j] * x[l] for l in range(n))
            f += (x[j] * g[i][j] / c) * (tau[i][j] - d / c)
        lng = a / b + f
        ln_gamma.append(lng)
        gamma.append(math.exp(lng))

    return NrtlActivityCoefficientsResult(
        ln_gamma=tuple(ln_gamma),
        gamma=tuple(gamma),
        warnings=tuple(warnings),
    )
