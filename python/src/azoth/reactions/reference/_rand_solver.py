"""The modified-RAND reactive solve, the Python twin of
``crates/azoth-reactions/src/rand_solver.rs``.

Simultaneous chemical and phase equilibrium by minimising the Gibbs energy subject to
the element balances: the unknowns are the moles ``n[j][i]`` and one Lagrange multiplier
``lambda_k`` per element row, and each pass corrects them together.

.. code-block:: text

    e[j][i] = g0_i + ln x_j,i + ln phi_j,i - sum_k lambda_k A_ki     the potential error
    C_kl    = sum_j sum_i A_ki A_li n_j,i                            the RAND matrix
    rhs_k   = (b_k - sum A n) + sum A n e                            the inventory plus the error
    n_j,i  *= exp(alpha (-e_j,i + sum_k A_ki dlambda_k))             the damped step
    lambda += alpha dlambda

**The rules that only apply above one phase** are the class's own and they matter: a
multiphase step starts at a tenth rather than the whole Newton direction, a ten-iteration
*sliding window* replaces the iteration-to-iteration damping rule, a residual of ``1e-4``
counts as converged rather than ``TOL``, and DIIS may extrapolate. The single-phase branch
keeps the strict rule, which is why the same fluid can report ``1.4e-14`` on one path and
``2.0e-5`` on the other.

**Where two phases converge to the same composition the split is not determined**, and
nothing here pretends otherwise: every split on that line satisfies the equilibrium
conditions, and which one a run reports is decided by its path.

The element residual is the class's **scaled root-mean-square** and not a bare difference:
each row's deviation is divided by ``max(|b_k|, max(totalMoles 1e-6, 1e-10))``, so a charge
row whose inventory is near zero reports a bounded number rather than amplifying its own
round-off.
"""

from __future__ import annotations

import math
from collections.abc import Callable
from typing import NamedTuple

from azoth.core.errors import InvalidInputError
from azoth.reactions.reference._diis import DiisAccelerator
from azoth.reactions.reference._linalg import rank_of_integer_matrix

#: The residual both the potential error and the element balance must come in under, from
#: ``ModifiedRANDSolver.TOL``.
TOL = 1.0e-9

#: The pass cap, from ``MAX_ITER``.
MAX_ITERATIONS = 500

#: The floor a mole number is kept above, from ``EPS``.
EPS = 1.0e-30

#: The gas constant the solver carries, from ``R_GAS``. **Not a units library's**: this is a
#: literal in the class, and the standard potentials are computed with it.
R_GAS = 8.314462

#: The reference pressure the ``ln(P/P_ref)`` term is taken against, from ``P_REF``.
P_REF = 1.0

#: The temperature the Cp polynomials are integrated from, from ``computeG0``'s ``T0``.
T0 = 298.15

#: The singular-pivot floor ``solveLinear`` refuses below.
LINEAR_PIVOT_FLOOR = 1.0e-30

#: The DIIS history's length and the pass before it may extrapolate.
DIIS_DEPTH = 6
DIIS_START = 5

#: The sliding window the multiphase damping rule averages over, from the class's ``WINDOW``.
DAMPING_WINDOW = 10

#: The potential error a multiphase solve accepts.
MULTIPHASE_ERROR_TOLERANCE = 1.0e-4

#: The element residual a *neutral* multiphase solve accepts - ``1e-8`` when a charge row is
#: present and ``1e-4`` when it is not.
MULTIPHASE_ELEMENT_TOLERANCE = 1.0e-4


class ThermoData(NamedTuple):
    """A component's ideal-gas heat-capacity polynomial and its two formation properties."""

    enthalpy_of_formation: float
    absolute_entropy: float
    gibbs_energy_of_formation: float
    cp: tuple[float, float, float, float, float]


class PhaseFeed(NamedTuple):
    """One phase's starting state, which ``initialize`` reads off the phase object."""

    fractions: list[float]
    beta: float


class RandSolution(NamedTuple):
    """What the solve answers with."""

    #: The overall moles, summed over the phases.
    moles: list[float]
    #: The moles in each phase, ``n[j][i]``.
    phase_moles: list[list[float]]
    #: Each phase's share, ``nPhase[j] / totalMoles`` - **the solver's own split**, which is
    #: not the phase objects' beta after the driver writes them back.
    phase_amounts: list[float]
    #: ``totalMoles``: ``nPhase`` summed.
    total_moles: float
    #: The element Lagrange multipliers at the answer.
    lambda_: list[float]
    #: Passes taken.
    iterations: int
    #: The largest absolute potential error over the species that matter.
    max_error: float
    #: The scaled root-mean-square element deviation.
    element_residual: float
    #: ``max(max_error, element_residual)``, lowered by any extrapolated step DIIS kept.
    final_residual: float
    #: Whether the residual passed the test that applied.
    converged: bool
    #: ``getDiisStepsAccepted``.
    diis_steps: int


def standard_potentials(data: list[ThermoData], temperature: float, pressure: float) -> list[float]:
    """``computeG0``'s neutral branch: ``g0_i = (hT - T sT) / (R T) + ln(P/P_ref)``.

    The three fallbacks are the class's own, in its order: with ``dHf`` and ``S0`` but no Cp
    data the polynomial is dropped, with only ``dGf298`` the Gibbs energy is used directly,
    and with nothing at all the potential is ``ln(P/P_ref)`` alone.
    """
    rt = R_GAS * temperature
    ln_p = math.log(pressure / P_REF)
    out: list[float] = []
    for entry in data:
        cp_a, cp_b, cp_c, cp_d, cp_e = entry.cp
        has_cp_data = abs(cp_a) > 1.0e-10 or abs(cp_b) > 1.0e-10
        has_thermo = (
            abs(entry.enthalpy_of_formation) > 1.0e-10 or abs(entry.absolute_entropy) > 1.0e-10
        )
        if has_thermo and has_cp_data:
            dt = temperature - T0
            dt2 = temperature**2 - T0**2
            dt3 = temperature**3 - T0**3
            dt4 = temperature**4 - T0**4
            dt5 = temperature**5 - T0**5
            ln_ratio = math.log(temperature / T0)
            delta_h = (
                cp_a * dt
                + cp_b / 2.0 * dt2
                + cp_c / 3.0 * dt3
                + cp_d / 4.0 * dt4
                + cp_e / 5.0 * dt5
            )
            delta_s = (
                cp_a * ln_ratio + cp_b * dt + cp_c / 2.0 * dt2 + cp_d / 3.0 * dt3 + cp_e / 4.0 * dt4
            )
            h_t = entry.enthalpy_of_formation + delta_h
            s_t = entry.absolute_entropy + delta_s
            out.append((h_t - temperature * s_t) / rt + ln_p)
        elif has_thermo:
            out.append(
                (entry.enthalpy_of_formation - temperature * entry.absolute_entropy) / rt + ln_p
            )
        elif abs(entry.gibbs_energy_of_formation) > 1.0e-10:
            out.append(entry.gibbs_energy_of_formation / rt + ln_p)
        else:
            out.append(ln_p)
    return out


def element_residual(
    a_matrix: list[list[float]], b: list[float], n: list[list[float]], total: float
) -> float:
    """``computeElementResidual``: the scaled root-mean-square deviation of ``A n`` from ``b``."""
    scale_floor = max(total * 1.0e-6, 1.0e-10)
    total_sq = 0.0
    for k, b_k in enumerate(b):
        element_sum = 0.0
        for phase in n:
            for i, moles in enumerate(phase):
                element_sum += a_matrix[k][i] * moles
        scale = max(abs(b_k), scale_floor)
        deviation = (element_sum - b_k) / scale
        total_sq += deviation * deviation
    return math.sqrt(total_sq)


def element_residual_vector(
    a_matrix: list[list[float]], b: list[float], n: list[list[float]], total: float
) -> list[float]:
    """``computeElementResidualVector``: the same deviations, one per element row."""
    scale_floor = max(total * 1.0e-6, 1.0e-10)
    out: list[float] = []
    for k, b_k in enumerate(b):
        element_sum = 0.0
        for phase in n:
            for i, moles in enumerate(phase):
                element_sum += a_matrix[k][i] * moles
        out.append((element_sum - b_k) / max(abs(b_k), scale_floor))
    return out


def solve_linear(matrix: list[list[float]], rhs: list[float]) -> list[float] | None:
    """``solveLinear``: Gaussian elimination with partial pivoting, refusing a small pivot."""
    dim = len(rhs)
    augment = [[*matrix[i][:dim], rhs[i]] for i in range(dim)]
    for column in range(dim):
        pivot = column
        largest = abs(augment[column][column])
        for row in range(column + 1, dim):
            if abs(augment[row][column]) > largest:
                largest = abs(augment[row][column])
                pivot = row
        if largest < LINEAR_PIVOT_FLOOR:
            return None
        augment[column], augment[pivot] = augment[pivot], augment[column]
        for row in range(column + 1, dim):
            factor = augment[row][column] / augment[column][column]
            for k in range(column, dim + 1):
                augment[row][k] -= factor * augment[column][k]
    solution = [0.0] * dim
    for i in range(dim - 1, -1, -1):
        total = augment[i][dim]
        for k in range(i + 1, dim):
            total -= augment[i][k] * solution[k]
        solution[i] = total / augment[i][i]
    return solution


def _initial_lambda(
    a_matrix: list[list[float]], g0: list[float], fractions: list[float], ln_phi: list[float]
) -> list[float]:
    """``initializeLambda``: ``lambda = (A A^T)^-1 A h`` with ``h_i = g0_i + ln x_i + ln phi_i``."""
    ne = len(a_matrix)
    nc = len(g0)
    h = [g0[i] + math.log(max(fractions[i], EPS)) + ln_phi[i] for i in range(nc)]
    ata = [[0.0] * ne for _ in range(ne)]
    ath = [0.0] * ne
    for k in range(ne):
        for l in range(ne):
            for i in range(nc):
                ata[k][l] += a_matrix[k][i] * a_matrix[l][i]
        for i in range(nc):
            ath[k] += a_matrix[k][i] * h[i]
    solved = solve_linear(ata, ath)
    return solved if solved is not None else [0.0] * ne


def _recalc_totals(n: list[list[float]]) -> tuple[list[float], float]:
    """``recalcTotals``: the phase amounts from the moles, floored, and the total."""
    n_phase = [max(sum(phase), EPS) for phase in n]
    total = max(sum(n_phase), EPS)
    return n_phase, total


def solve(
    a_matrix: list[list[float]],
    g0: list[float],
    b: list[float],
    total_moles: float,
    phases: list[PhaseFeed],
    ln_phi: Callable[[int, list[float]], list[float]],
) -> RandSolution:
    """The neutral RAND solve over a phase list, ``ModifiedRANDSolver.solve``.

    ``ln_phi`` is called as ``ln_phi(phase_index, fractions)`` and answers each component's
    ``ln(phi_i)`` there; the index is what a caller needs to pick the phase's own root.

    Raises:
        InvalidInputError: where the shapes disagree.
    """
    ne = len(b)
    nc = len(g0)
    np_ = len(phases)
    if len(a_matrix) != ne or np_ == 0:
        raise InvalidInputError(
            "a_matrix", f"{ne} element row(s), {len(g0)} potential(s) and {np_} phase(s)"
        )
    for row in a_matrix:
        if len(row) != nc:
            raise InvalidInputError(
                "a_matrix", f"a row has {len(row)} entries against {nc} component(s)"
            )
    for phase in phases:
        if len(phase.fractions) != nc:
            raise InvalidInputError(
                "phases", f"a phase has {len(phase.fractions)} entries against {nc} component(s)"
            )

    # `initialize`: the moles and fractions the phase objects hold, floored.
    n: list[list[float]] = [[0.0] * nc for _ in range(np_)]
    fractions: list[list[float]] = [[0.0] * nc for _ in range(np_)]
    for j, phase in enumerate(phases):
        beta = max(phase.beta, EPS)
        for i in range(nc):
            fractions[j][i] = phase.fractions[i]
            n[j][i] = max(fractions[j][i] * beta * total_moles, EPS)
    n_phase = [max(phase.beta, EPS) * total_moles for phase in phases]
    total = sum(n_phase)

    ln_phi_here = [list(ln_phi(j, fractions[j])) for j in range(np_)]

    # `NR = 0` is not an iteration: the element balance alone fixes the composition.
    independent_reactions = nc - rank_of_integer_matrix(a_matrix)
    if independent_reactions == 0:
        element = element_residual(a_matrix, b, n, total)
        return RandSolution(
            moles=[sum(n[j][i] for j in range(np_)) for i in range(nc)],
            phase_moles=n,
            phase_amounts=[n_phase[j] / total for j in range(np_)],
            total_moles=total,
            lambda_=[0.0] * ne,
            iterations=0,
            max_error=element,
            element_residual=element,
            final_residual=element,
            converged=True,
            diis_steps=0,
        )

    lambda_ = _initial_lambda(a_matrix, g0, fractions[0], ln_phi_here[0])
    diis = DiisAccelerator(ne, DIIS_DEPTH)
    diis_steps = 0

    # The class starts a multiphase solve at a tenth of the direction.
    damping = 0.1 if np_ > 1 else 1.0
    previous_residual = float("inf")
    stagnation = 0
    window_residual = float("inf")
    iterations_since_window = 0

    final_error = float("inf")
    final_element = float("inf")
    final_residual = float("inf")
    converged = False
    iterations = 0

    for iteration in range(MAX_ITERATIONS):
        iterations = iteration + 1

        error = [[0.0] * nc for _ in range(np_)]
        for j in range(np_):
            for i in range(nc):
                x_i = max(fractions[j][i], EPS)
                row_sum = sum(lambda_[k] * a_matrix[k][i] for k in range(ne))
                value = g0[i] + math.log(x_i) + ln_phi_here[j][i] - row_sum
                error[j][i] = value if math.isfinite(value) else 0.0

        c = [[0.0] * ne for _ in range(ne)]
        rhs = [0.0] * ne
        for k in range(ne):
            element_sum = 0.0
            element_error = 0.0
            for j in range(np_):
                for i in range(nc):
                    element_sum += a_matrix[k][i] * n[j][i]
                    element_error += a_matrix[k][i] * n[j][i] * error[j][i]
            rhs[k] = (b[k] - element_sum) + element_error
            for l in range(ne):
                value = 0.0
                for j in range(np_):
                    for i in range(nc):
                        value += a_matrix[k][i] * a_matrix[l][i] * n[j][i]
                c[k][l] = value
            c[k][k] += max(1.0e-10 * abs(c[k][k]), 1.0e-14)

        scale = [
            1.0 / math.sqrt(abs(c[k][k])) if abs(c[k][k]) > 1.0e-30 else 1.0 for k in range(ne)
        ]
        for k in range(ne):
            for l in range(ne):
                c[k][l] *= scale[k] * scale[l]
            rhs[k] *= scale[k]

        delta = solve_linear(c, rhs)
        if delta is None:
            break
        delta = [delta[k] * scale[k] for k in range(ne)]

        n_old = [list(phase) for phase in n]
        lambda_old = list(lambda_)

        alpha = damping
        accepted = False
        for _ in range(5):
            n = [list(phase) for phase in n_old]
            lambda_ = list(lambda_old)
            for j in range(np_):
                for i in range(nc):
                    correction = alpha * -error[j][i]
                    for k in range(ne):
                        correction += alpha * a_matrix[k][i] * delta[k]
                    if not math.isfinite(correction):
                        correction = 0.0
                    correction = max(-3.0, min(3.0, correction))
                    stepped = n[j][i] * math.exp(correction)
                    n[j][i] = stepped if math.isfinite(stepped) and stepped >= EPS else EPS
            for k in range(ne):
                lambda_[k] += alpha * delta[k]
            n_phase, total = _recalc_totals(n)
            if element_residual(a_matrix, b, n, total) < previous_residual * 1.5 or alpha < 0.05:
                accepted = True
                break
            alpha *= 0.5

        if not accepted:
            n = [list(phase) for phase in n_old]
            lambda_ = list(lambda_old)
            alpha = 0.1
            for j in range(np_):
                for i in range(nc):
                    correction = alpha * -error[j][i]
                    for k in range(ne):
                        correction += alpha * a_matrix[k][i] * delta[k]
                    if not math.isfinite(correction):
                        correction = 0.0
                    correction = max(-3.0, min(3.0, correction))
                    stepped = n[j][i] * math.exp(correction)
                    n[j][i] = stepped if math.isfinite(stepped) and stepped >= EPS else EPS
            for k in range(ne):
                lambda_[k] += alpha * delta[k]
            n_phase, total = _recalc_totals(n)

        for j in range(np_):
            fractions[j] = [n[j][i] / n_phase[j] for i in range(nc)]
        ln_phi_here = [list(ln_phi(j, fractions[j])) for j in range(np_)]

        max_error = 0.0
        for j in range(np_):
            for i in range(nc):
                if n[j][i] > 1.0e-10 * total:
                    max_error = max(max_error, abs(error[j][i]))
        element = element_residual(a_matrix, b, n, total)
        final_error = max_error
        final_element = element
        residual = max(max_error, element)
        final_residual = residual

        if max_error < TOL and element < TOL:
            converged = True
            break
        if (
            np_ > 1
            and max_error < MULTIPHASE_ERROR_TOLERANCE
            and element < MULTIPHASE_ELEMENT_TOLERANCE
        ):
            converged = True
            break

        if np_ == 1:
            if residual < previous_residual * 0.9:
                damping = min(1.0, damping * 1.5)
                stagnation = 0
            elif residual > previous_residual * 1.1:
                damping = max(0.01, damping * 0.5)
                stagnation += 1
            else:
                stagnation += 1
        else:
            if residual > previous_residual * 1.5:
                damping = max(0.01, damping * 0.5)
                stagnation += 1
            elif residual > previous_residual * 1.05:
                stagnation += 1
            iterations_since_window += 1
            if iterations_since_window >= DAMPING_WINDOW:
                if residual < window_residual * 0.5:
                    ceiling = 0.3 if residual > 10.0 else (0.5 if residual > 1.0 else 1.0)
                    damping = min(ceiling, damping * 2.0)
                    stagnation = 0
                elif residual > window_residual * 2.0:
                    damping = max(0.01, damping * 0.25)
                    stagnation += DAMPING_WINDOW
                window_residual = residual
                iterations_since_window = 0

        if stagnation > 50 and max_error < 1.0e-3 and element < MULTIPHASE_ELEMENT_TOLERANCE:
            converged = True
            break

        # DIIS, on the multipliers, with the element deviation as its error. **Every** pass is
        # recorded and the gate only decides whether an extrapolation is tried.
        diis.add_entry(lambda_, element_residual_vector(a_matrix, b, n, total))
        if iteration >= DIIS_START and diis.can_extrapolate():
            extrapolated = diis.extrapolate()
            if extrapolated is not None:
                saved_moles = [list(phase) for phase in n]
                saved_lambda = list(lambda_)
                for j in range(np_):
                    for i in range(nc):
                        correction = sum(
                            (extrapolated[k] - lambda_[k]) * a_matrix[k][i] for k in range(ne)
                        )
                        correction = max(-3.0, min(3.0, correction))
                        stepped = n[j][i] * math.exp(correction)
                        n[j][i] = stepped if math.isfinite(stepped) and stepped >= EPS else EPS
                lambda_ = list(extrapolated)
                n_phase, total = _recalc_totals(n)
                for j in range(np_):
                    fractions[j] = [n[j][i] / n_phase[j] for i in range(nc)]
                ln_phi_here = [list(ln_phi(j, fractions[j])) for j in range(np_)]

                diis_error = 0.0
                for j in range(np_):
                    for i in range(nc):
                        if n[j][i] > 1.0e-10 * total:
                            x_i = max(fractions[j][i], EPS)
                            row_sum = sum(lambda_[k] * a_matrix[k][i] for k in range(ne))
                            value = g0[i] + math.log(x_i) + ln_phi_here[j][i] - row_sum
                            diis_error = max(diis_error, abs(value))
                diis_element = element_residual(a_matrix, b, n, total)
                diis_residual = max(diis_error, diis_element)
                if diis_residual < residual * 1.1:
                    final_residual = min(final_residual, diis_residual)
                    diis_steps += 1
                else:
                    n = saved_moles
                    lambda_ = saved_lambda
                    n_phase, total = _recalc_totals(n)
                    for j in range(np_):
                        fractions[j] = [n[j][i] / n_phase[j] for i in range(nc)]
                    ln_phi_here = [list(ln_phi(j, fractions[j])) for j in range(np_)]

        previous_residual = final_residual

    return RandSolution(
        moles=[sum(n[j][i] for j in range(np_)) for i in range(nc)],
        phase_moles=n,
        phase_amounts=[n_phase[j] / total for j in range(np_)],
        total_moles=total,
        lambda_=lambda_,
        iterations=iterations,
        max_error=final_error,
        element_residual=final_element,
        final_residual=final_residual,
        converged=converged,
        diis_steps=diis_steps,
    )


def _one_phase(
    ln_phi: Callable[[list[float]], list[float]],
) -> Callable[[int, list[float]], list[float]]:
    """The one-phase closure: the phase index is not asked of a model that has one phase."""
    return lambda _phase, x: ln_phi(x)


def solve_single_phase(
    a_matrix: list[list[float]],
    g0: list[float],
    b: list[float],
    feed_moles: list[float],
    ln_phi: Callable[[list[float]], list[float]],
) -> RandSolution:
    """The single-phase case of :func:`solve`, named so a caller need not build a list."""
    total = sum(feed_moles)
    if total <= 0.0:
        raise InvalidInputError("feed_moles", "the feed holds no moles")
    fractions = [moles / total for moles in feed_moles]
    return solve(
        a_matrix,
        g0,
        b,
        total,
        [PhaseFeed(fractions=fractions, beta=1.0)],
        _one_phase(ln_phi),
    )
