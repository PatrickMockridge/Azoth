"""``eos.tp_multiflash`` - how many phases a feed splits into, and how much of each.

Spec: ``specs/models/eos/tp_multiflash.toml``. NeqSim's ``TPmultiflash``, reached through
``TPflash.runInternal()`` when ``system.doMultiPhaseCheck()`` is true - so
:func:`azoth.eos.reference.pt_flash.pt_flash` is this model with the stability seeding left
out, and the two agree everywhere the seeding has nothing to add.

Three steps, in the order upstream runs them: the two-phase flash; the tangent-plane trial
from :func:`azoth.eos.reference.stability_test.stability_test`, whose stationary point at a
composition that is not already a phase is a phase the flash has not found; and the fraction
solve below, which is ``multiphase.rs``'s ``solve_phase_fractions`` in Python.

The Rust kernel carries the reasoning; this is the same arithmetic. The two are compared case
by case, and the case files record NeqSim's own numbers beside them.

**One known divergence.** NeqSim also runs a pure-component fallback - one trial per
component, started from that component - which this does not, so a phase only such a start
reaches is not found. Measured at N2/CO2/n-octane at 180 K and 10 bar, where NeqSim reaches
three phases and this stops at two.
"""

from __future__ import annotations

import math
from typing import Any, NamedTuple

from azoth import _models_gen
from azoth.core.errors import InvalidInputError, SolverNotConvergedError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import Phase, TpMultiflashResult, TpMultiflashSeed
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning, WarningCode
from azoth.eos.mixture import Mixture
from azoth.eos.reference._mixture_state import (
    ReducedParameters,
    helmholtz_energy,
    mixture_parameters,
    phase_state,
    reduced_parameters,
)
from azoth.eos.reference.pr_z_factor import pr_z_factor
from azoth.eos.reference.pt_flash import pt_flash
from azoth.eos.reference.srk_z_factor import srk_z_factor
from azoth.eos.reference.stability_test import TM_LIMIT, stability_test

MODEL_ID = "eos.tp_multiflash"

#: A trial within this of an existing phase, in sum-absolute composition, *is* that phase.
#: A trial that reproduces an existing phase has a negative tangent-plane distance by
#: construction, because the plane passes through it.
TRIVIAL_TOLERANCE = 1.0e-4

#: The floor under a phase fraction, upstream's ``phaseFractionMinimumLimit``.
FRACTION_FLOOR = 1.0e-12

#: A phase below this is dropped: above the floor by a hair, so a phase the Newton *pinned*
#: at the floor is dropped while one it settled just above is kept.
FRACTION_DROP = 1.1e-12

#: Two phases at the same composition are one phase, however their roots were labelled.
DUPLICATE_TOLERANCE = 1.0e-6

#: The diagonal regulariser. The Hessian is a sum of rank-one terms, so it is singular
#: whenever two phases sit at the same composition - the state the iteration passes through
#: on its way to merging them.
REGULARISER = 1.0e-3

#: The gradient norm the iteration must also reach, upstream's second convergence test.
GRADIENT_TOLERANCE = 1.0e-10

#: The fewest steps taken, upstream's ``|| iter < 3``.
MINIMUM_STEPS = 2


class _Phase(NamedTuple):
    """One phase of the set being solved: an amount, a composition and a root."""

    fraction: float
    composition: list[float]
    liquid: bool


def _distance(a: list[float], b: list[float]) -> float:
    """``sum_i |a_i - b_i|``, the trivial-solution distance upstream tests against ``1e-4``."""
    return sum(abs(x - y) for x, y in zip(a, b, strict=True))


def _is_liquid(reduced: ReducedParameters, kij: Any, x: list[float]) -> bool:
    """Whether a composition's lower-Gibbs root is the liquid one.

    The comparison is ``A^R / RT - ln Z + Z``, which is :func:`_feed_state`'s rule and the same
    one, so the side a seeded phase is placed on cannot disagree with the reference the trial
    was measured from.
    """
    a_mix, b_mix = mixture_parameters(reduced.a, reduced.b, kij, x)
    roots = (
        srk_z_factor(a_mix, b_mix)
        if reduced.cubic.name in ("srk", "rk")
        else pr_z_factor(a_mix, b_mix)
    )
    candidates = [roots.z_min] if roots.z_min == roots.z_max else [roots.z_min, roots.z_max]
    best = None
    best_gibbs = math.inf
    for compressibility in candidates:
        gibbs = (
            helmholtz_energy(reduced, kij, x, compressibility)
            - math.log(compressibility)
            + compressibility
        )
        if gibbs < best_gibbs:
            best_gibbs = gibbs
            best = compressibility
    assert best is not None  # `candidates` is never empty
    return best == roots.z_min


def _phases_of(flash: Any, reduced: ReducedParameters, kij: Any, z: list[float]) -> list[_Phase]:
    """The phases a two-phase flash reported, as a starting set for the fraction solve."""
    if flash.beta is not None:
        return [
            _Phase(1.0 - flash.beta, list(flash.x), liquid=True),
            _Phase(flash.beta, list(flash.y), liquid=False),
        ]
    # No fraction: either every K-value was on one side, or the iteration reached `x = y = z`.
    # The first names the phase; the second does not, so the feed goes on whichever root has
    # the lower Gibbs energy.
    if flash.phase is Phase.ALL_LIQUID:
        liquid = True
    elif flash.phase is Phase.ALL_VAPOUR:
        liquid = False
    else:
        liquid = _is_liquid(reduced, kij, list(z))
    return [_Phase(1.0, list(z), liquid=liquid)]


def _solve(hessian: list[list[float]], gradient: list[float]) -> list[float] | None:
    """``H x = g`` by Gaussian elimination with partial pivoting, or ``None`` if singular.

    The pivot takes the **last** largest magnitude, not the first. That is what Rust's
    ``max_by`` does - it returns the last element among equal maxima - and the two kernels
    have to break a tie the same way or they take different elimination paths. Measured, a
    first-maximum rule here against Rust's last made `eos.tp_multiflash`'s N2/CO2/n-octane
    case take 32 steps on one side and 34 on the other while converging to the same answer.
    """
    n = len(gradient)
    a = [[*hessian[i][:], gradient[i]] for i in range(n)]
    for column in range(n):
        pivot = column
        for row in range(column + 1, n):
            if abs(a[row][column]) >= abs(a[pivot][column]):
                pivot = row
        a[column], a[pivot] = a[pivot], a[column]
        if a[column][column] == 0.0:
            return None
        pivot_value = a[column][column]
        for row in range(column + 1, n):
            factor = a[row][column] / pivot_value
            for entry in range(column, n + 1):
                a[row][entry] -= factor * a[column][entry]
    x = [0.0] * n
    for row in range(n - 1, -1, -1):
        total = sum(a[row][k] * x[k] for k in range(row + 1, n))
        x[row] = (a[row][n] - total) / a[row][row]
    return x


def _solve_phase_fractions(
    reduced: ReducedParameters,
    kij: Any,
    feed: list[float],
    phases: list[_Phase],
    tolerance: float,
    cap: int,
) -> tuple[list[_Phase], int, float, bool]:
    """The fractions of a set of phases at a state, and their compositions.

    ``Q(beta) = sum_k beta_k - sum_i z_i ln E_i`` with ``E_i = sum_k beta_k / phi_ik``,
    whichever way the coefficients are held: the gradient and the Hessian are the two sums
    below, so one step is a dense solve of size *n* - the phase count - rather than a bracketed
    root find. The composition the material balance gives at those fractions is
    ``x_ik = z_i / (E_i phi_ik)``, rebuilt each step and normalised.

    Raises:
        SolverNotConvergedError: if the Hessian is singular at an iterate, which means two
            phases have met. Reaching the cap is not an error: upstream returns whatever
            residual it left with, and so does this, with ``converged`` False.
    """
    n = len(feed)
    count = len(phases)
    residual = math.nan
    gradient_norm = math.inf
    converged = False
    iterations = 0

    for step in range(1, cap + 1):
        iterations = step
        phi = [
            [
                math.exp(value)
                for value in phase_state(
                    reduced, kij, phase.composition, liquid=phase.liquid
                ).ln_phi
            ]
            for phase in phases
        ]

        e = [0.0] * n
        for k, phase in enumerate(phases):
            for i in range(n):
                e[i] += phase.fraction / phi[k][i]
        for i in range(n):
            if e[i] < 1.0e-100:
                e[i] = 1.0e-100

        gradient = [1.0] * count
        for k in range(count):
            for i in range(n):
                gradient[k] -= feed[i] / (e[i] * phi[k][i])
        hessian = [[0.0] * count for _ in range(count)]
        for j in range(count):
            for k in range(count):
                for i in range(n):
                    hessian[j][k] += feed[i] / (e[i] * e[i] * phi[j][i] * phi[k][i])
                if j == k:
                    hessian[j][k] += REGULARISER

        gradient_norm = math.sqrt(sum(value * value for value in gradient))
        correction = _solve(hessian, gradient)
        if correction is None:
            raise SolverNotConvergedError(iterations, math.nan, tolerance)
        residual = math.sqrt(sum(value * value for value in correction))

        # Damped by `n/(n+3)` - a quarter at the first step, approaching one - and clamped
        # away from the ends, because a fraction at zero is a phase that is not there.
        scale = step / (step + 3.0)
        updated: list[_Phase] = []
        total = 0.0
        for k, phase in enumerate(phases):
            candidate = phase.fraction - scale * correction[k]
            fraction = min(max(candidate, FRACTION_FLOOR), 1.0 - FRACTION_FLOOR)
            updated.append(_Phase(fraction, phase.composition, phase.liquid))
            total += fraction
        for k in range(count):
            updated[k] = _Phase(
                updated[k].fraction / total, updated[k].composition, updated[k].liquid
            )

        # The compositions the material balance gives at those fractions, in the same closed
        # form and from the same frozen coefficients.
        for k in range(count):
            row = [feed[i] / (e[i] * phi[k][i]) for i in range(n)]
            row_total = sum(row)
            if row_total > 0.0:
                row = [value / row_total for value in row]
            updated[k] = _Phase(updated[k].fraction, row, updated[k].liquid)
        phases = updated

        if step >= MINIMUM_STEPS and residual <= tolerance and gradient_norm <= GRADIENT_TOLERANCE:
            converged = True
            break

    return phases, iterations, max(residual, gradient_norm), converged


def tp_multiflash(mixture: Mixture, T: Q, P: Q, z: list[float]) -> TpMultiflashResult:
    """The number of phases a feed splits into at a temperature and pressure, and their amounts.

    Raises:
        InvalidInputError: if ``z`` is the wrong length, negative, or does not sum to one.
        OutOfRangeError: if ``T`` or ``P`` is not positive.
        SolverNotConvergedError: if the fraction solve's Hessian is singular.

    Example:
        >>> import azoth
        >>> from azoth.eos import components
        >>> q = azoth.ureg.Quantity
        >>> mix, _ = components.mixture_of(["CO2", "methane", "nc10"])
        >>> r = tp_multiflash(mix, q(200.0, "K"), q(1.0e6, "Pa"), [0.4, 0.3, 0.3])
        >>> r.phase_count
        3
    """
    spec = _models_gen.model(MODEL_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    t_si = input_to_si(spec, "T", T)
    p_si = input_to_si(spec, "P", P)
    n = len(mixture.components)
    if len(z) != n:
        raise InvalidInputError("z", f"a feed for {n} components has {len(z)} entries")
    if any(value < 0.0 for value in z):
        raise InvalidInputError(
            "z",
            "a mole fraction cannot be negative",
        )
    total = sum(z)
    if abs(total - 1.0) > 1.0e-09:
        raise InvalidInputError(
            "z",
            f"the mole fractions sum to {total}, not to one. Renormalising them here "
            f"would make a composition error invisible in every number downstream, "
            f"so it is refused instead",
        )
    apply_checks(checks.on_input, {"T": t_si, "P": p_si}.get, warnings)

    algorithm = spec["algorithm"]
    tolerance = float(algorithm["tolerance"])
    cap = int(algorithm["max_iterations"])

    reduced = reduced_parameters(mixture, t_si, p_si)
    warnings.extend(reduced.warnings)
    kij = mixture.kij
    feed = list(z)

    flash = pt_flash(mixture, T=from_si(t_si, "K"), P=from_si(p_si, "Pa"), z=feed)
    warnings.extend(flash.warnings)
    phases = _phases_of(flash, reduced, kij, feed)
    base_count = len(phases)

    stability = stability_test(mixture, T=T, P=P, z=feed)
    warnings.extend(stability.warnings)

    for trial, distance in zip(stability.w, stability.tm, strict=True):
        if distance >= TM_LIMIT:
            continue
        if any(_distance(list(trial), phase.composition) < TRIVIAL_TOLERANCE for phase in phases):
            continue
        trial_row = list(trial)
        liquid = _is_liquid(reduced, kij, trial_row)
        # Upstream's starting fraction: the feed's mole fraction of the trial's largest
        # component. A share of the feed rather than a guess at the answer.
        #
        # The **last** largest, matching Rust's `max_by`, which keeps the last element among
        # equal maxima. A first-maximum rule picks a different component when two entries tie,
        # which is a different starting fraction and so a different iteration.
        dominant = 0
        for index in range(1, n):
            if trial_row[index] >= trial_row[dominant]:
                dominant = index
        phases.append(_Phase(feed[dominant], trial_row, liquid=liquid))
        # Upstream returns after the first phase it adds; a second seeding would start from a
        # set the solve has not seen.
        break

    # **Nothing seeded, nothing to solve.** When the trial added no phase the answer is the
    # two-phase flash's, and re-solving the fractions it already converged would replace a
    # converged split with an iterate of the same equations. Upstream's numbers say the same:
    # with the flag on and off, every state the sweep found no seeding at reports identical
    # betas to all seventeen digits.
    if len(phases) > base_count:
        total = sum(phase.fraction for phase in phases)
        phases = [
            _Phase(phase.fraction / total, phase.composition, phase.liquid) for phase in phases
        ]
        phases, iterations, residual, converged = _solve_phase_fractions(
            reduced, kij, feed, phases, tolerance, cap
        )
    else:
        iterations, residual, converged = 0, 0.0, True
    # The merge. A phase the solve pinned at the floor is a phase that is not there, and two
    # phases at the same composition are one phase however their roots were labelled. The side
    # is deliberately not part of that second test - see the Rust kernel for the measurement.
    kept: list[_Phase] = []
    for phase in phases:
        if phase.fraction < FRACTION_DROP:
            continue
        match = next(
            (
                index
                for index, other in enumerate(kept)
                if _distance(other.composition, phase.composition) < DUPLICATE_TOLERANCE
            ),
            None,
        )
        if match is not None:
            existing = kept[match]
            kept[match] = _Phase(
                existing.fraction + phase.fraction, existing.composition, existing.liquid
            )
            continue
        kept.append(phase)

    # A dropped phase's share returns to the others, which is upstream's
    # `removePhaseTotalComposition`: removing a phase says it is not there, not that part of
    # the feed has gone.
    total = sum(phase.fraction for phase in kept)
    if total > 0.0:
        kept = [_Phase(phase.fraction / total, phase.composition, phase.liquid) for phase in kept]

    # Whether the solve contributed to the answer, which is the same predicate `seeded` reports.
    # When the merge folded the seeded phase away the answer is the two-phase flash's and the
    # solve produced nothing, so its steps and residual are not properties of the answer.
    # Measured: at N2/CO2/n-octane 190 K / 1 bar the seeding fires and is then merged away, and
    # the discarded solve took 32 steps in one kernel and 34 in the other.
    grew = len(kept) > base_count
    if not grew:
        iterations, residual = 0, 0.0

    if grew and not converged:
        warnings.append(
            Warning(
                code=WarningCode.SOLVER_NOT_CONVERGED,
                message=(
                    f"the fraction solve left by its cap after {iterations} steps with a "
                    f"residual of {residual:.3e}, against the {tolerance} the spec declares. "
                    f"The fractions are the iterate it reached, not a state it settled at, "
                    f"and upstream accepts this rather than refusing the flash"
                ),
            )
        )

    min_t_over_tc = min(flash.min_t_over_tc, stability.min_t_over_tc)
    geometry = [phase_state(reduced, kij, phase.composition, liquid=phase.liquid) for phase in kept]

    return TpMultiflashResult(
        phase_count=len(kept),
        beta=tuple(phase.fraction for phase in kept),
        x=tuple(tuple(phase.composition) for phase in kept),
        z_factor=tuple(state.z for state in geometry),
        ln_phi=tuple(tuple(state.ln_phi) for state in geometry),
        seeded=(TpMultiflashSeed.STABILITY_SEEDED if grew else TpMultiflashSeed.TWO_PHASE_FLASH),
        tm=tuple(stability.tm),
        iterations=iterations,
        residual=residual,
        min_t_over_tc=min_t_over_tc,
        warnings=tuple(warnings),
    )
