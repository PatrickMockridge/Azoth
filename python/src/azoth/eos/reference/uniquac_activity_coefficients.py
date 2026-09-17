"""``eos.uniquac_activity_coefficients`` - the activity coefficients from the UNIQUAC
model.

Spec: ``specs/models/eos/uniquac_activity_coefficients.toml``. A *direct* model: no
iteration, so no algorithm block.
"""

from __future__ import annotations

import math
from collections.abc import Sequence

from azoth import _models_gen
from azoth.core.errors import InvalidInputError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import UniquacActivityCoefficientsResult
from azoth.core.units import Q, input_to_si
from azoth.core.warnings import Warning
from azoth.eos.components import UniquacParameters

MODEL_ID = "eos.uniquac_activity_coefficients"


def uniquac_activity_coefficients(
    params: UniquacParameters,
    T: Q,
    x: Sequence[float],
    aij: Sequence[Sequence[Q]],
) -> UniquacActivityCoefficientsResult:
    """The activity coefficients of a mixture, from UNIQUAC (Abrams-Prausnitz).

    ``phi_i = r_i x_i / sum r_j x_j`` and ``theta_i = q_i x_i / sum q_j x_j`` are the
    volume and surface fractions, ``l_i = 5 (r_i - q_i) - (r_i - 1)``, and
    ``tau_ij = exp(-aij[i][j] / T)``. The combinatorial term is
    ``ln gamma^C_i = ln(phi_i/x_i) + 5 q_i ln(theta_i/phi_i) + l_i - (phi_i/x_i) sum
    x_j l_j`` and the residual is ``q_i (1 - ln(sum_j theta_j tau_ji) - sum_j theta_j
    tau_ij / sum_k theta_k tau_kj)``. ``params`` is resolved by name through
    :func:`azoth.eos.components.uniquac_parameters`; ``aij`` stays the caller's,
    because no upstream table carries a UNIQUAC interaction matrix. It is directional
    with a zero diagonal, and ``x`` is checked rather than renormalised.

    Args:
        params: the resolved volume and surface parameters, by component.
        T: absolute temperature.
        x: mole fractions; non-negative and summing to one.
        aij: the interaction energy matrix, in Kelvin.

    Returns:
        The natural logarithm and the value of each activity coefficient.

    Raises:
        InvalidInputError: if the vectors disagree in length, ``aij`` is not ``N x N``
            with a zero diagonal, or ``x`` is not a composition.
        OutOfRangeError: if ``T`` is not positive.

    Example:
        >>> import azoth
        >>> from azoth.eos.components import uniquac_parameters
        >>> q = azoth.ureg.Quantity
        >>> aij = [[q(0.0, "K"), q(-71.0, "K")], [q(209.0, "K"), q(0.0, "K")]]
        >>> r = uniquac_activity_coefficients(
        ...     uniquac_parameters(["methanol", "water"]), q(298.15, "K"),
        ...     [0.5, 0.5], aij,
        ... )
        >>> round(r.gamma[0], 14)
        1.21854418487282
    """
    spec = _models_gen.model(MODEL_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    values = {"T": input_to_si(spec, "T", T)}
    apply_checks(checks.on_input, values.get, warnings)

    t = input_to_si(spec, "T", T)
    x = list(x)
    r = list(params.r)
    q = list(params.q)
    aij_si = [[input_to_si(spec, "aij", value) for value in row] for row in aij]

    n = len(x)
    if n == 0:
        raise InvalidInputError("x", "a mixture of zero components has no activity coefficient")
    if len(r) != n or len(q) != n:
        raise InvalidInputError(
            "components",
            f"the per-component vectors disagree in length: `x` has {n} entries, `r` "
            f"{len(r)} and `q` {len(q)}",
        )
    if len(aij_si) != n:
        raise InvalidInputError("aij", f"`aij` has {len(aij_si)} rows but must be {n} x {n}")
    for i in range(n):
        if len(aij_si[i]) != n:
            raise InvalidInputError(
                "aij", f"row {i} of `aij` has {len(aij_si[i])} entries but must be {n}"
            )
        if aij_si[i][i] != 0.0:
            raise InvalidInputError(
                "aij", f"the diagonal must be zero, but aij[{i}][{i}] = {aij_si[i][i]}"
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

    sum_r = sum(x[i] * r[i] for i in range(n))
    sum_q = sum(x[i] * q[i] for i in range(n))
    l = [5.0 * (r[i] - q[i]) - (r[i] - 1.0) for i in range(n)]
    sum_l = sum(x[i] * l[i] for i in range(n))
    theta = [q[i] * x[i] / sum_q for i in range(n)]
    tau = [[math.exp(-aij_si[m][k] / t) for k in range(n)] for m in range(n)]

    ln_gamma: list[float] = []
    gamma: list[float] = []
    for i in range(n):
        phi = r[i] * x[i] / sum_r
        lng_c = (
            math.log(phi / x[i])
            + 5.0 * q[i] * math.log(theta[i] / phi)
            + l[i]
            - (phi / x[i]) * sum_l
        )

        s1 = sum(theta[j] * tau[j][i] for j in range(n))
        s3 = 0.0
        for j in range(n):
            denom = sum(theta[k] * tau[k][j] for k in range(n))
            s3 += theta[j] * tau[i][j] / denom
        lng_r = q[i] * (1.0 - math.log(s1) - s3)

        lng = lng_c + lng_r
        ln_gamma.append(lng)
        gamma.append(math.exp(lng))

    return UniquacActivityCoefficientsResult(
        ln_gamma=tuple(ln_gamma),
        gamma=tuple(gamma),
        warnings=tuple(warnings),
    )
