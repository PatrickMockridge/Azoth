"""``eos.unifac_psrk_activity_coefficients`` - the activity coefficients from UNIFAC
with PSRK's temperature-dependent interaction parameters.

Spec: ``specs/models/eos/unifac_psrk_activity_coefficients.toml``. A *direct* model:
no iteration, so no algorithm block.

This is the pure-Python reference: a second, independent expression of the same
physics as the Rust kernel, held to it by ``test_cross_impl.py``.
"""

from __future__ import annotations

import math
from collections.abc import Sequence

from azoth import _models_gen
from azoth.core.errors import InvalidInputError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import UnifacPsrkActivityCoefficientsResult
from azoth.core.units import Q, input_to_si
from azoth.core.warnings import Warning
from azoth.eos.components import UnifacPsrkParameters

MODEL_ID = "eos.unifac_psrk_activity_coefficients"


def _ln_gamma_group(
    k: int, theta: list[float], group_q: Sequence[float], aij: Sequence[float], t: float
) -> float:
    """The residual ``ln Gamma_k`` for one group at a surface-fraction distribution.

    ``aij`` is ``G x G`` row-major, so ``aij[m][n]`` is ``aij[m * g + n]``.
    """
    g = len(group_q)
    s1 = sum(theta[m] * math.exp(-aij[m * g + k] / t) for m in range(g))
    s3 = 0.0
    for m in range(g):
        s2 = sum(theta[n] * math.exp(-aij[n * g + m] / t) for n in range(g))
        s3 += theta[m] * math.exp(-aij[k * g + m] / t) / s2
    return group_q[k] * (1.0 - math.log(s1) - s3)


def unifac_psrk_activity_coefficients(
    params: UnifacPsrkParameters,
    T: Q,
    x: Sequence[float],
) -> UnifacPsrkActivityCoefficientsResult:
    """The activity coefficients of a mixture, from UNIFAC with PSRK's interaction
    parameters.

    Identical to :func:`azoth.eos.reference.unifac_activity_coefficients` except that
    the main-group interaction is a function of temperature:
    ``a_mn(T) = a_mn + b_mn T + c_mn T**2``, NeqSim's
    ``ComponentGEUnifacPSRK.calcaij``. The group basis, the combinatorial term and the
    residual are the same, so this evaluates the interaction matrix at ``T`` and runs
    the shared body. ``params`` is resolved by name through
    :func:`azoth.eos.components.unifac_psrk_parameters`; ``x`` is checked rather than
    renormalised.

    Args:
        params: the resolved UNIFAC-PSRK parameters, by component.
        T: absolute temperature.
        x: mole fractions; non-negative and summing to one.

    Returns:
        The natural logarithm and the value of each activity coefficient.

    Raises:
        InvalidInputError: if the shapes disagree, a group count is negative, or ``x``
            is not a composition.
        OutOfRangeError: if ``T`` is not positive.

    Example:
        >>> import azoth
        >>> from azoth.eos.components import unifac_psrk_parameters
        >>> q = azoth.ureg.Quantity
        >>> cold = unifac_psrk_activity_coefficients(
        ...     unifac_psrk_parameters(["water", "methane"]), q(250.0, "K"), [0.5, 0.5]
        ... )
        >>> round(cold.gamma[0], 14)
        2.31218279323519
    """
    spec = _models_gen.model(MODEL_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    values = {"T": input_to_si(spec, "T", T)}
    apply_checks(checks.on_input, values.get, warnings)

    t = input_to_si(spec, "T", T)
    x = list(x)
    groups = params.groups
    group_r = params.group_r
    group_q = params.group_q

    n = len(x)
    g = len(group_r)
    if n == 0:
        raise InvalidInputError("x", "a mixture of zero components has no activity coefficient")
    if g == 0 or len(group_q) != g:
        raise InvalidInputError(
            "components",
            f"`group_r` has {g} entries and `group_q` has {len(group_q)}; both must be "
            "the same non-empty length",
        )
    if len(groups) != n * g:
        raise InvalidInputError(
            "components",
            f"`groups` has {len(groups)} entries but `x` has {n} components and there "
            f"are {g} groups, which needs {n * g}",
        )
    if any(count < 0.0 for count in groups):
        raise InvalidInputError("components", "a group count cannot be negative")
    for field, matrix in (("aij", params.aij), ("bij", params.bij), ("cij", params.cij)):
        if len(matrix) != g * g:
            raise InvalidInputError(
                "components", f"`{field}` has {len(matrix)} entries but must be {g} x {g}"
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

    aij = [
        a + b * t + c * t * t for a, b, c in zip(params.aij, params.bij, params.cij, strict=True)
    ]

    ri = [sum(groups[i * g + k] * group_r[k] for k in range(g)) for i in range(n)]
    qi = [sum(groups[i * g + k] * group_q[k] for k in range(g)) for i in range(n)]

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
        qmix = [
            group_q[l] * sum(x[j] * groups[j * g + l] for j in range(n)) / denom for l in range(g)
        ]
        qcomp = [group_q[l] * groups[i * g + l] / qi[i] for l in range(g)]
        lng_r = sum(
            groups[i * g + k]
            * (
                _ln_gamma_group(k, qmix, group_q, aij, t)
                - _ln_gamma_group(k, qcomp, group_q, aij, t)
            )
            for k in range(g)
        )

        lng = lng_c + lng_r
        ln_gamma.append(lng)
        gamma.append(math.exp(lng))

    return UnifacPsrkActivityCoefficientsResult(
        ln_gamma=tuple(ln_gamma),
        gamma=tuple(gamma),
        warnings=tuple(warnings),
    )
