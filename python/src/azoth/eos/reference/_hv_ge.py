"""The Huron-Vidal excess-Gibbs model: a co-volume-weighted NRTL.

The Python reference, the mirror of ``crates/azoth-eos/src/hv_ge.rs``. NeqSim's
``PhaseGENRTLmodifiedHV`` feeds the Huron-Vidal mixing rule, and its activity coefficient is
the NRTL formula with two twists: ``G_ij = B_i exp(-alpha_ij tau_ij)`` - the weight is the
row's co-volume, which does not cancel - and the ``tau`` matrix is mixed, with pairs the
database marks ``HV`` carrying the fitted NRTL parameters and every other pair the cubic's
own excess energy, so the rule stays consistent with the pure-component equation of state.
"""

from __future__ import annotations

import math

from azoth.core.errors import InvalidInputError


def _matrices(
    x: list[float],
    t_kelvin: float,
    a: list[float],
    b: list[float],
    kij: list[float],
    hv_gij: list[float],
    hv_gij_t: list[float],
    hv_alpha: list[float],
    hv_pairs: list[bool],
    lambda_: float,
) -> tuple[list[float], list[float]]:
    """The ``tau`` and ``G`` matrices the NRTL sums read, flattened row-major."""
    n = len(x)
    tau = [0.0] * (n * n)
    g = [0.0] * (n * n)
    for i in range(n):
        for j in range(n):
            at = i * n + j
            if hv_pairs[at]:
                tau_ij = hv_gij[at] / t_kelvin + hv_gij_t[at]
                alpha = hv_alpha[at]
            else:
                tau_ij = lambda_ * (
                    a[j] / b[j] - 2.0 * math.sqrt(a[i] * a[j]) / (b[i] + b[j]) * (1.0 - kij[at])
                )
                alpha = 0.0
            tau[at] = tau_ij
            g[at] = b[i] * math.exp(-alpha * tau_ij)
    return tau, g


def hv_ln_gamma(
    x: list[float],
    t_kelvin: float,
    a: list[float],
    b: list[float],
    kij: list[float],
    hv_gij: list[float],
    hv_gij_t: list[float],
    hv_alpha: list[float],
    hv_pairs: list[bool],
    lambda_: float,
) -> list[float]:
    """The natural logarithms of the activity coefficients for the Huron-Vidal model.

    ``a`` and ``b`` are the *reduced* attraction and repulsion, ``kij`` the interaction
    matrix, ``hv_gij`` the fitted ``Dij`` in kelvin, ``hv_gij_t`` its temperature coefficient
    ``DijT``, ``hv_alpha`` the fitted non-randomness, and ``lambda_`` the cubic's Huron-Vidal
    constant.

    Raises:
        InvalidInputError: if the matrices are not ``N x N`` for ``N = len(x)``.
    """
    n = len(x)
    for name, matrix in (
        ("kij", kij),
        ("hv_gij", hv_gij),
        ("hv_gij_t", hv_gij_t),
        ("hv_alpha", hv_alpha),
        ("hv_pairs", hv_pairs),
    ):
        if len(matrix) != n * n:
            raise InvalidInputError(
                "huron_vidal",
                f"{name} is {len(matrix)} long and the mixture has {n} components, so the "
                f"matrix is {n * n}",
            )
    tau, g = _matrices(x, t_kelvin, a, b, kij, hv_gij, hv_gij_t, hv_alpha, hv_pairs, lambda_)
    out: list[float] = []
    for i in range(n):
        first_num = 0.0
        first_den = 0.0
        for j in range(n):
            first_num += tau[j * n + i] * g[j * n + i] * x[j]
            first_den += g[j * n + i] * x[j]
        second = 0.0
        for j in range(n):
            den = 0.0
            num = 0.0
            for l in range(n):
                den += g[l * n + j] * x[l]
                num += g[l * n + j] * tau[l * n + j] * x[l]
            second += g[i * n + j] * x[j] / den * (tau[i * n + j] - num / den)
        out.append(first_num / first_den + second)
    return out
