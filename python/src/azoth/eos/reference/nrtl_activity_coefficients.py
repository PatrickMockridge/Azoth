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
from azoth.eos.components import NrtlParameters

MODEL_ID = "eos.nrtl_activity_coefficients"


def nrtl_activity_coefficients(
    params: NrtlParameters,
    T: Q,
    x: Sequence[float],
) -> NrtlActivityCoefficientsResult:
    """The activity coefficients of a mixture, from NRTL (Renon-Prausnitz).

    ``params.dij[i][j] = g_ij`` is in Kelvin so that ``tau_ij = params.dij[i][j] / T``
    is dimensionless, ``G_ij = exp(-params.alpha[i][j] tau_ij)``, and ``ln gamma_i`` is
    the standard symmetric sum: ``(sum_j tau_ji G_ji x_j) / (sum_j G_ji x_j)`` plus
    ``sum_j (x_j G_ij / C_j) (tau_ij - D_j / C_j)``. ``params`` is resolved by name
    through :func:`azoth.eos.components.nrtl_parameters`, which is symmetric with a
    zero diagonal in ``alpha`` by construction; ``x`` is checked rather than
    renormalised.

    Args:
        params: the resolved NRTL parameters, by component.
        T: absolute temperature.
        x: mole fractions; non-negative and summing to one.

    Returns:
        The natural logarithm and the value of each activity coefficient.

    Raises:
        InvalidInputError: if either matrix is not ``N x N`` for the ``N`` of ``x``, or
            ``x`` is not a composition.
        OutOfRangeError: if ``T`` is not positive.

    Example:
        >>> import azoth
        >>> from azoth.eos.components import nrtl_parameters
        >>> q = azoth.ureg.Quantity
        >>> r = nrtl_activity_coefficients(
        ...     nrtl_parameters(["methanol", "water"]), q(298.15, "K"), [0.5, 0.5]
        ... )
        >>> round(r.gamma[0], 12)
        1.233178856168
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
    for field, matrix in (("alpha", params.alpha), ("dij", params.dij)):
        if len(matrix) != n * n:
            raise InvalidInputError(
                "components",
                f"the {field} matrix has {len(matrix)} entries but `x` has {n} "
                f"components, which needs {n * n} - an N x N matrix row-major",
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

    tau = [[params.dij[i * n + j] / t for j in range(n)] for i in range(n)]
    g = [[math.exp(-params.alpha[i * n + j] * tau[i][j]) for j in range(n)] for i in range(n)]

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
