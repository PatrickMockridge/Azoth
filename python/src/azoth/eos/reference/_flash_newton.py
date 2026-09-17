"""The second-order step the slowly-converging flashes hand over to.

The mirror of ``crates/azoth-eos/src/flash_newton.rs``, and NeqSim's
``SysNewtonRhapsonTPflash``: a Newton step on the isofugacity residuals in the *vapour
mole numbers* ``u_i = beta y_i``, which is the variable set that makes the residuals'
Jacobian the Hessian of Michelsen's reduced Gibbs energy

```text
Q(u) = sum_i u_i ln(y_i phi_i^V) + (z_i - u_i) ln(x_i phi_i^L)
```

and the step a descent direction for it. That is why the line search below is on ``Q``
rather than on the residual: the residual's gradient *is* ``Q``'s, so ``Q`` is the
function whose minimisation the equations describe, and a step that does not decrease
it is a step away from the answer.

Every derivative is analytic, from
:func:`azoth.eos.reference._mixture_state.phase_derivatives`. Nothing is differenced.
"""

from __future__ import annotations

import math
from typing import Any, NamedTuple

from azoth.core.errors import SolverNotConvergedError
from azoth.eos.mixture import Mixture
from azoth.eos.reference._mixture_state import (
    ReducedParameters,
    phase_derivatives,
    phase_state,
)


class Step(NamedTuple):
    """One Newton step's outcome."""

    #: The vapour mole numbers the step lands on.
    u: list[float]
    #: ``max_i |g_i|`` at the iterate the step landed on.
    #:
    #: The **stopping measure**, and it is the residual rather than NeqSim's
    #: ``deviation``. NeqSim stops on the relative step, which is a different quantity
    #: from the one the scheme it replaced stops on - an rms change in ``ln K`` - so the
    #: two answers stop at different distances from the same root. The residual is what
    #: the equations actually say, so stopping on it makes the two stop in the same place.
    residual: float


def split_at(feed: list[float], u: list[float]) -> tuple[list[float], list[float]]:
    """The two phases at a set of vapour mole numbers."""
    beta = sum(u)
    y = [value / beta for value in u]
    x = [(feed[i] - u[i]) / (1.0 - beta) for i in range(len(feed))]
    return x, y


def _gibbs_energy(
    mixture: Mixture,
    reduced: ReducedParameters,
    feed: list[float],
    u: list[float],
) -> float:
    """Michelsen's reduced Gibbs energy at a set of vapour mole numbers."""
    x, y = split_at(feed, u)
    liquid = phase_state(reduced, mixture.kij, x, liquid=True)
    vapour = phase_state(reduced, mixture.kij, y, liquid=False)
    total = 0.0
    for i in range(len(u)):
        total += u[i] * (vapour.ln_phi[i] + math.log(y[i]))
        total += (feed[i] - u[i]) * (liquid.ln_phi[i] + math.log(x[i]))
    return total


def _is_feasible(feed: list[float], u: list[float]) -> bool:
    """Whether a trial ``u`` describes two phases at all.

    NeqSim's ``isFeasible``: every ``u_i`` strictly between zero and ``z_i``, and
    ``beta`` strictly inside ``(0, 1)``. A step that leaves this region has asked for a
    negative mole number, and the compositions would take its logarithm.
    """
    edge = 1.0e-15
    if any(u[i] < edge or u[i] > feed[i] - edge for i in range(len(u))):
        return False
    beta = sum(u)
    return 1.0e-12 < beta < 1.0 - 1.0e-12


def _solve(j: list[list[float]], g: list[float], tolerance: float) -> list[float]:
    """``J dx = g`` by Gaussian elimination with partial pivoting."""
    n = len(g)
    a = [[*j[i], g[i]] for i in range(n)]
    for column in range(n):
        pivot = max(range(column, n), key=lambda r: abs(a[r][column]))
        a[column], a[pivot] = a[pivot], a[column]
        if a[column][column] == 0.0:
            raise SolverNotConvergedError(0, float("nan"), tolerance)
        for row in range(column + 1, n):
            factor = a[row][column] / a[column][column]
            for k in range(column, n + 1):
                a[row][k] -= factor * a[column][k]
    dx = [0.0] * n
    for row in range(n - 1, -1, -1):
        total = sum(a[row][k] * dx[k] for k in range(row + 1, n))
        dx[row] = (a[row][n] - total) / a[row][row]
    return dx


def step(
    mixture: Mixture,
    reduced: ReducedParameters,
    feed: list[float],
    u: list[float],
    algorithm: dict[str, Any],
    temperature: float,
    pressure: float,
) -> Step:
    """One Newton step, from an iterate the caller's own scheme produced."""
    n = len(feed)
    x, y = split_at(feed, u)
    liquid = phase_state(reduced, mixture.kij, x, liquid=True)
    vapour = phase_state(reduced, mixture.kij, y, liquid=False)
    d_vapour = phase_derivatives(
        reduced, mixture.kij, y, vapour.z, temperature=temperature, pressure=pressure
    )
    d_liquid = phase_derivatives(
        reduced, mixture.kij, x, liquid.z, temperature=temperature, pressure=pressure
    )
    beta = sum(u)

    residual = [
        vapour.ln_phi[i] + math.log(y[i]) - liquid.ln_phi[i] - math.log(x[i]) for i in range(n)
    ]
    current = _gibbs_energy(mixture, reduced, feed, u)

    # NeqSim's Jacobian, and it is the same matrix the derivation gives. The two `- 1`
    # terms are the normalisation of each phase's fractions, and the two composition
    # derivatives are `d ln phi_i / d n_j` at unit total, which is why each is divided
    # by its own phase's mole number.
    jacobian = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            dij = 1.0 if i == j else 0.0
            jacobian[i][j] = (d_vapour.d_ln_phi_dn[i][j] + dij / y[i] - 1.0) / beta + (
                d_liquid.d_ln_phi_dn[i][j] + dij / x[i] - 1.0
            ) / (1.0 - beta)
    # Levenberg-Marquardt regularisation, NeqSim's and at NeqSim's magnitude.
    trace = sum(abs(jacobian[i][i]) for i in range(n))
    lam = 1.0e-8 * trace / n
    for i in range(n):
        jacobian[i][i] += lam

    delta = _solve(jacobian, residual, algorithm["tolerance"])
    slope = sum(residual[i] * delta[i] for i in range(n))

    # Armijo backtracking on `Q`, NeqSim's `c1 = 1e-4` and eight halvings. The halved
    # step is taken even when none of the eight satisfies Armijo, which is NeqSim's
    # behaviour and not an oversight in it: `alpha` has been divided down to `1/256` by
    # then, so its direction matters more than its descent, and refusing it outright
    # stops the solve over a line search that has already done its job.
    alpha = 1.0
    for _ in range(8):
        trial = [u[i] - alpha * delta[i] for i in range(n)]
        if _is_feasible(feed, trial):
            try:
                if _gibbs_energy(mixture, reduced, feed, trial) <= current - 1.0e-4 * alpha * slope:
                    break
            except Exception:
                pass
        alpha *= 0.5

    stepped = [u[i] - alpha * delta[i] for i in range(n)]

    # The residual at the iterate the step landed on, which is what the solve stops on.
    landed_x, landed_y = split_at(feed, stepped)
    landed_liquid = phase_state(reduced, mixture.kij, landed_x, liquid=True)
    landed_vapour = phase_state(reduced, mixture.kij, landed_y, liquid=False)
    landed = max(
        abs(
            landed_vapour.ln_phi[i]
            + math.log(landed_y[i])
            - landed_liquid.ln_phi[i]
            - math.log(landed_x[i])
        )
        for i in range(n)
    )

    return Step(u=stepped, residual=landed)
