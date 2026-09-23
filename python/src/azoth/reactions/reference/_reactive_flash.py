"""The reactive flash driver, the Python twin of
``crates/azoth-reactions/src/reactive_flash.rs``.

``run`` is a sequence of decisions around the two solvers: it builds the element matrix,
**short-circuits at ``NR = 0``** (no independent reaction means the element balance fully
determines the composition), initialises a phase split with a conventional VLE flash when
more than one phase is allowed, asks the stability analysis whether a second phase forms,
and then drives the modified-RAND solve in an outer loop that removes the phases that turn
out negligible.

**``SystemSrkEos`` constructs with two phase objects**, each holding the whole feed at
``beta = 1.0``, and ``init(0)`` restores that count - so the captured systems hand the
driver ``np = 2`` and ``maxPhases = 2`` and its ``skipStability`` rule takes it into the
outer loop **without any trial phase ever being added**.

**The weight ``computeGibbsEnergy`` gives each phase is stale.** ``updateSystem`` writes the
solver's converged fractions into the phase objects and then ``system.init(1)``, inside the
same method, re-initialises them from the *system's* array - refreshed only on the ionic
branch. So the converged split never reaches the system, and the measure uses whatever the
last system-level write left: ``(1 - V, V)`` from the VLE initialisation, or ``(1.0, 1.0)``
from the constructor. The same stale array is what ``removeNegligiblePhases`` *tests*, so a
phase the solve drove to nothing is not removed.

``removePhaseKeepTotalComposition`` is the class's removal and it **touches no moles**: it
shifts the phase index array and decrements the count. The name describes an intent the
method does not carry out.
"""

from __future__ import annotations

import math
from collections.abc import Callable
from typing import NamedTuple

from azoth.core.errors import InvalidInputError
from azoth.reactions.reference import _reactive_stability as stability
from azoth.reactions.reference._rand_solver import (
    PhaseFeed,
    RandSolution,
    solve,
    solve_single_phase,
)
from azoth.reactions.reference._reactive_stability import CriticalConstants

#: The floor ``initializeWithVLEFlash`` keeps a trial mole fraction above, which is its own
#: ``1e-30``.
VLE_MIN_MOLES = 1.0e-30

#: The bound ``initializeWithVLEFlash`` rejects a Rachford-Rice answer at.
VLE_SINGLE_PHASE_TOLERANCE = 1.0e-10

#: The ``|ln K|`` above which a component counts as volatile.
VLE_VOLATILE_LOG_K = 0.1

#: The outer loop's pass cap, from ``MAX_OUTER_ITER``.
MAX_OUTER_ITERATIONS = 20

#: The phase fraction below which ``removeNegligiblePhases`` drops a phase.
MIN_PHASE_FRACTION = 1.0e-12

#: The fraction ``addTrialPhases`` gives a new phase.
TRIAL_PHASE_BETA = 0.01

#: The residual ``run`` accepts a *non-converged* multiphase solve at.
NEAR_CONVERGED_RESIDUAL = 5.0e-3


class VleInitialisation(NamedTuple):
    """What ``initializeWithVLEFlash`` leaves behind."""

    vapour_fraction: float
    liquid: list[float]
    vapour: list[float]


class SinglePhaseOutcome(NamedTuple):
    """What the driver reports for a fluid it finds stable in one phase."""

    moles: list[float]
    total_moles: float
    gibbs_energy: float
    iterations: int
    converged: bool
    solution: RandSolution


class NonReactiveOutcome(NamedTuple):
    """What ``runNonReactiveFlash`` leaves: a two-phase split of a fluid with no reaction."""

    phases: list[PhaseFeed]
    converged: bool
    iterations: int
    all_supercritical: bool


class FlashOutcome(NamedTuple):
    """What the driver's ``run`` answers with."""

    phases: list[PhaseFeed]
    converged: bool
    total_iterations: int
    equilibrium_total_moles: float
    gibbs_energy: float
    solution: RandSolution | None


def _normalise(values: list[float]) -> list[float]:
    """``Phase.normalize()``: divide by the sum, leaving an all-zero vector alone."""
    total = sum(values)
    if total <= 0.0:
        return list(values)
    return [value / total for value in values]


def vle_initialization(
    fractions: list[float],
    constants: list[CriticalConstants],
    temperature: float,
    pressure: float,
) -> VleInitialisation | None:
    """``initializeWithVLEFlash``: Wilson K-values and a Rachford-Rice solve.

    ``None`` is the class's early return and it means one of two things: no component is
    volatile, or the root came back at a bound. The driver then adds no phase.

    The composition is ``getz()``, the *overall* one. The Rachford-Rice loop is the class's
    own: 100 passes, a Newton step clamped to ``[0, 1]`` after each, a ``1e-30`` floor on the
    denominator that *skips* the component rather than guarding it, and a ``1e-10`` bound.
    """
    nc = len(fractions)
    log_k = [0.0] * nc
    has_volatile = False
    for i in range(nc):
        entry = constants[i]
        if entry.pc > 0.0 and entry.tc > 0.0:
            log_k[i] = math.log(entry.pc / pressure) + 5.373 * (1.0 + entry.omega) * (
                1.0 - entry.tc / temperature
            )
            if abs(log_k[i]) > VLE_VOLATILE_LOG_K:
                has_volatile = True
    if not has_volatile:
        return None

    vapour = 0.5
    for _ in range(100):
        f = 0.0
        derivative = 0.0
        for i in range(nc):
            k = math.exp(log_k[i])
            denominator = 1.0 + vapour * (k - 1.0)
            if abs(denominator) < 1e-30:
                continue
            f += fractions[i] * (k - 1.0) / denominator
            derivative -= fractions[i] * (k - 1.0) * (k - 1.0) / (denominator * denominator)
        if abs(f) < 1e-10:
            break
        if abs(derivative) > 1e-30:
            vapour -= f / derivative
        vapour = max(0.0, min(1.0, vapour))

    if vapour < VLE_SINGLE_PHASE_TOLERANCE or vapour > 1.0 - VLE_SINGLE_PHASE_TOLERANCE:
        return None

    unfloored = [fractions[i] / (1.0 + vapour * (math.exp(log_k[i]) - 1.0)) for i in range(nc)]
    liquid = [max(x, VLE_MIN_MOLES) for x in unfloored]
    vapour_composition = [max(math.exp(log_k[i]) * unfloored[i], VLE_MIN_MOLES) for i in range(nc)]
    return VleInitialisation(
        vapour_fraction=vapour,
        liquid=_normalise(liquid),
        vapour=_normalise(vapour_composition),
    )


def solve_rachford_rice(feed: list[float], log_k: list[float], beta_guess: float) -> float:
    """``solveRachfordRice``: the class's own Newton solve, 50 passes, root clamped inside.

    **A different routine from ``initializeWithVLEFlash``'s inline loop**: 50 passes against
    100, a step test against a residual test, and a clamp that keeps the root inside the
    interval rather than rejecting it at the bound. Both are in the class.
    """
    beta = beta_guess
    for _ in range(50):
        f = 0.0
        derivative = 0.0
        for i in range(len(feed)):
            k = math.exp(log_k[i])
            denominator = 1.0 + beta * (k - 1.0)
            if abs(denominator) < 1e-30:
                continue
            f += feed[i] * (k - 1.0) / denominator
            derivative -= feed[i] * (k - 1.0) * (k - 1.0) / (denominator * denominator)
        if abs(derivative) < 1e-30:
            break
        step = f / derivative
        beta -= step
        beta = max(1e-15, min(1.0 - 1e-15, beta))
        if abs(step) < 1e-12:
            break
    return beta


def _compositions(
    feed: list[float], log_k: list[float], beta: float
) -> tuple[list[float], list[float]]:
    """The class's two composition updates, which it writes twice."""
    liquid: list[float] = []
    vapour: list[float] = []
    for i in range(len(feed)):
        k = math.exp(log_k[i])
        x = feed[i] / (1.0 + beta * (k - 1.0))
        liquid.append(max(x, VLE_MIN_MOLES))
        vapour.append(max(k * x, VLE_MIN_MOLES))
    return liquid, vapour


def non_reactive_flash(
    feed: list[float],
    constants: list[CriticalConstants],
    temperature: float,
    pressure: float,
    liquid_index: int,
    ln_phi: Callable[[int, list[float]], list[float]],
) -> NonReactiveOutcome:
    """``runNonReactiveFlash``: the conventional VLE flash a fluid with no reaction gets.

    Wilson K-values seed a split, then successive substitution replaces them with
    ``K_i = phi_liq,i / phi_vap,i`` until the log-K change comes in under ``1e-10``.

    Raises:
        InvalidInputError: where a shape disagrees.
    """
    nc = len(feed)
    if len(constants) != nc or liquid_index > 1:
        raise InvalidInputError(
            "constants",
            f"{len(constants)} constant set(s), {nc} component(s) and a liquid index of "
            f"{liquid_index}",
        )
    gas_index = 1 - liquid_index

    all_supercritical = all(temperature >= constants[i].tc for i in range(nc))
    if all_supercritical:
        return NonReactiveOutcome(
            phases=[
                PhaseFeed(fractions=list(feed), beta=1.0),
                PhaseFeed(fractions=list(feed), beta=1.0),
            ],
            converged=True,
            iterations=0,
            all_supercritical=True,
        )

    log_k = [
        math.log(constants[i].pc / pressure)
        + 5.373 * (1.0 + constants[i].omega) * (1.0 - constants[i].tc / temperature)
        if constants[i].pc > 0.0 and constants[i].tc > 0.0
        else 0.0
        for i in range(nc)
    ]

    beta = max(1e-15, min(1.0 - 1e-15, solve_rachford_rice(feed, log_k, 0.5)))
    liquid, vapour = _compositions(feed, log_k, beta)

    converged = False
    iterations = 0
    for iteration in range(200):
        iterations = iteration + 1
        ln_phi_liquid = ln_phi(liquid_index, _normalise(liquid))
        ln_phi_vapour = ln_phi(gas_index, _normalise(vapour))

        error = 0.0
        for i in range(nc):
            new_log_k = ln_phi_liquid[i] - ln_phi_vapour[i]
            error += abs(new_log_k - log_k[i])
            log_k[i] = new_log_k

        beta = max(1e-15, min(1.0 - 1e-15, solve_rachford_rice(feed, log_k, beta)))
        liquid, vapour = _compositions(feed, log_k, beta)
        if error < 1e-10:
            converged = True
            break

    phases = [PhaseFeed(fractions=[], beta=0.0), PhaseFeed(fractions=[], beta=0.0)]
    phases[liquid_index] = PhaseFeed(fractions=_normalise(liquid), beta=1.0 - beta)
    phases[gas_index] = PhaseFeed(fractions=_normalise(vapour), beta=beta)
    return NonReactiveOutcome(
        phases=phases, converged=converged, iterations=iterations, all_supercritical=False
    )


def normalise_betas(phases: list[PhaseFeed]) -> None:
    """``normalizeBeta``: divide every phase's fraction by their sum.

    Raises:
        InvalidInputError: where the fractions sum to zero.
    """
    total = sum(phase.beta for phase in phases)
    if total <= 0.0:
        raise InvalidInputError("phases", "the phase fractions sum to zero")
    for j, phase in enumerate(phases):
        phases[j] = PhaseFeed(fractions=phase.fractions, beta=phase.beta / total)


def remove_negligible_phases(phases: list[PhaseFeed]) -> bool:
    """``removeNegligiblePhases``: drop every phase under :data:`MIN_PHASE_FRACTION`.

    **The fraction it tests is the one the *system* holds, and on a neutral fluid that is the
    stale array** - so a phase the *solve* drove to nothing is not removed by this test.
    Reading the solve's own ``phase_amounts`` here would be a different algorithm.
    """
    removed = False
    j = len(phases)
    while j > 0:
        j -= 1
        if len(phases) <= 1:
            break
        if phases[j].beta < MIN_PHASE_FRACTION:
            phases.pop(j)
            removed = True
    if removed:
        normalise_betas(phases)
    return removed


def add_trial_phase(phases: list[PhaseFeed], trials: list[list[float]], max_phases: int) -> bool:
    """``addTrialPhases``: append **one** phase, the first of the trials, at
    :data:`TRIAL_PHASE_BETA`, then renormalise every fraction.

    Raises:
        InvalidInputError: where a trial's length disagrees with the phases'.
    """
    if not trials:
        return False
    if len(phases) >= max_phases:
        return False
    trial = trials[0]
    if phases and len(trial) != len(phases[0].fractions):
        raise InvalidInputError(
            "trials",
            f"a trial has {len(trial)} entries against {len(phases[0].fractions)} component(s)",
        )
    phases.append(
        PhaseFeed(
            fractions=_normalise([max(x, VLE_MIN_MOLES) for x in trial]),
            beta=TRIAL_PHASE_BETA,
        )
    )
    normalise_betas(phases)
    return True


def gibbs_energy(fractions: list[float], ln_phi: list[float]) -> float:
    """``sum_i x_i (ln x_i + ln phi_i)`` over one phase."""
    return sum(x * (math.log(x) + phi) for x, phi in zip(fractions, ln_phi, strict=True))


def total_gibbs_energy(phases: list[tuple[float, list[float], list[float]]]) -> float:
    """``computeGibbsEnergy``: the beta-weighted sum over the phase list.

    **The weights are the phase objects' own betas, and on a neutral fluid those are
    stale** - a port that substituted the solver's fractions would reproduce a number NeqSim
    does not report.
    """
    return sum(beta * gibbs_energy(fractions, ln_phi) for beta, fractions, ln_phi in phases)


def single_phase_equilibrium(
    a_matrix: list[list[float]],
    g0: list[float],
    b: list[float],
    feed_moles: list[float],
    ln_phi: Callable[[list[float]], list[float]],
) -> SinglePhaseOutcome:
    """``solveSinglePhaseChemicalEquilibrium``, and **this branch computes no Gibbs energy**:
    it returns before the driver's ``computeGibbsEnergy`` call, so the caller sees ``0.0``."""
    solution = solve_single_phase(a_matrix, g0, b, feed_moles, ln_phi)
    total = sum(solution.moles)
    fractions = [moles / total for moles in solution.moles]
    energy = gibbs_energy(fractions, ln_phi(fractions))
    return SinglePhaseOutcome(
        moles=list(solution.moles),
        total_moles=total,
        gibbs_energy=energy,
        iterations=solution.iterations,
        converged=solution.converged,
        solution=solution,
    )


def outer_loop(
    a_matrix: list[list[float]],
    g0: list[float],
    b: list[float],
    total_moles: float,
    initial: list[PhaseFeed],
    max_phases: int,
    ln_phi: Callable[[int, list[float]], list[float]],
    stability_check: Callable[[list[PhaseFeed]], list[list[float]]],
) -> FlashOutcome:
    """``run``'s outer loop, from the point the driver has a phase list to iterate on."""
    phases = list(initial)
    converged = False
    total_iterations = 0
    equilibrium_total_moles = total_moles
    outcome: RandSolution | None = None

    for _ in range(MAX_OUTER_ITERATIONS):
        solution = solve(a_matrix, g0, b, total_moles, phases, ln_phi)
        total_iterations += solution.iterations
        equilibrium_total_moles = solution.total_moles
        rand_converged = solution.converged
        outcome = solution

        if (
            not rand_converged
            and len(phases) > 1
            and solution.final_residual < NEAR_CONVERGED_RESIDUAL
        ):
            converged = True
            break

        phase_removed = remove_negligible_phases(phases)

        if rand_converged and not phase_removed:
            if len(phases) >= max_phases:
                converged = True
                break
            trials = stability_check(phases)
            phases_before = len(phases)
            add_trial_phase(phases, trials, max_phases)
            if len(phases) == phases_before:
                converged = True
                break

        if rand_converged and phase_removed:
            converged = True
            break

    assert outcome is not None
    return FlashOutcome(
        phases=phases,
        converged=converged,
        total_iterations=total_iterations,
        equilibrium_total_moles=equilibrium_total_moles,
        gibbs_energy=0.0,
        solution=outcome,
    )


def _render_gibbs(
    phases: list[PhaseFeed],
    solution: RandSolution,
    ln_phi: Callable[[int, list[float]], list[float]],
) -> float:
    """``computeGibbsEnergy`` over a solved phase list, at the fractions in the phase rows."""
    entries: list[tuple[float, list[float], list[float]]] = []
    for index, phase in enumerate(phases):
        moles = (
            solution.phase_moles[index]
            if index < len(solution.phase_moles)
            else list(phase.fractions)
        )
        total = sum(moles)
        fractions = [value / total for value in moles] if total > 0.0 else list(phase.fractions)
        entries.append((phase.beta, fractions, ln_phi(index, fractions)))
    return total_gibbs_energy(entries)


class DriverState(NamedTuple):
    """Everything the driver needs that is not a closure."""

    feed_moles: list[float]
    a_matrix: list[list[float]]
    g0: list[float]
    b: list[float]
    total_moles: float
    constants: list[CriticalConstants]
    charges: list[float]
    temperature: float
    pressure: float
    max_phases: int
    phases: list[PhaseFeed]


def run(
    state: DriverState,
    phase_ln_phi: Callable[[int, list[float]], list[float]],
    single_ln_phi: Callable[[list[float]], list[float]],
    ce: Callable[[list[float]], list[float]],
) -> FlashOutcome:
    """``ReactiveMultiphaseTPflash.run``: the driver, composed.

    A stable answer goes to :func:`single_phase_equilibrium` and returns - **with
    ``converged = True`` regardless of what that solve reported**, which is the class's own
    overwrite - while an unstable one adds a trial phase and runs :func:`outer_loop`.

    Raises:
        InvalidInputError: on a shape disagreement.
    """
    nc = len(state.feed_moles)
    feed_total = sum(state.feed_moles)
    feed_fractions = [moles / feed_total for moles in state.feed_moles]
    phases = list(state.phases)

    independent_reactions = nc - _rank(state.a_matrix)
    if independent_reactions == 0:
        if len(phases) <= 1 or any(charge != 0.0 for charge in state.charges):
            return FlashOutcome(phases, True, 0, state.total_moles, 0.0, None)
        flashed = non_reactive_flash(
            feed_fractions,
            state.constants,
            state.temperature,
            state.pressure,
            1,
            phase_ln_phi,
        )
        return FlashOutcome(flashed.phases, True, 0, state.total_moles, 0.0, None)

    if state.max_phases > 1 and len(phases) == 1:
        initialisation = vle_initialization(
            feed_fractions, state.constants, state.temperature, state.pressure
        )
        if initialisation is not None:
            phases = [
                PhaseFeed(
                    fractions=initialisation.liquid, beta=1.0 - initialisation.vapour_fraction
                ),
                PhaseFeed(fractions=initialisation.vapour, beta=initialisation.vapour_fraction),
            ]

    if state.max_phases == 1 and len(phases) > 1:
        phases = phases[:1]

    skip_stability = len(phases) >= 2 and state.max_phases >= 2
    unstable_trials: list[list[float]] = []
    if not skip_stability:
        outcome = stability.analyse(
            phases[0].fractions,
            state.constants,
            state.temperature,
            state.pressure,
            state.charges,
            ce,
            single_ln_phi,
        )
        if not outcome.unstable:
            single = single_phase_equilibrium(
                state.a_matrix, state.g0, state.b, state.feed_moles, single_ln_phi
            )
            total = sum(single.moles)
            return FlashOutcome(
                phases=[PhaseFeed(fractions=[moles / total for moles in single.moles], beta=1.0)],
                converged=True,
                total_iterations=single.iterations,
                equilibrium_total_moles=single.total_moles,
                gibbs_energy=0.0,
                solution=single.solution,
            )
        unstable_trials = outcome.unstable_trials

    if unstable_trials:
        add_trial_phase(phases, unstable_trials, state.max_phases)

    def stability_check(current: list[PhaseFeed]) -> list[list[float]]:
        reference = list(current[0].fractions) if current else []
        return stability.analyse(
            reference,
            state.constants,
            state.temperature,
            state.pressure,
            state.charges,
            ce,
            single_ln_phi,
        ).unstable_trials

    looped = outer_loop(
        state.a_matrix,
        state.g0,
        state.b,
        state.total_moles,
        phases,
        state.max_phases,
        phase_ln_phi,
        stability_check,
    )
    assert looped.solution is not None
    return FlashOutcome(
        phases=looped.phases,
        converged=looped.converged,
        total_iterations=looped.total_iterations,
        equilibrium_total_moles=looped.equilibrium_total_moles,
        gibbs_energy=_render_gibbs(looped.phases, looped.solution, phase_ln_phi),
        solution=looped.solution,
    )


def _rank(matrix: list[list[float]]) -> int:
    """The element matrix's rank, by the exact elimination the reactions namespace uses.

    The Rust twin calls the *ported* NeqSim elimination with its ``1e-12`` floor; the two
    agree on every matrix here, whose entries are the element table's integers.
    """
    from azoth.reactions.reference._linalg import rank_of_integer_matrix

    return rank_of_integer_matrix(matrix)


__all__ = [
    "DriverState",
    "FlashOutcome",
    "NonReactiveOutcome",
    "SinglePhaseOutcome",
    "VleInitialisation",
    "add_trial_phase",
    "gibbs_energy",
    "non_reactive_flash",
    "normalise_betas",
    "outer_loop",
    "remove_negligible_phases",
    "run",
    "single_phase_equilibrium",
    "solve_rachford_rice",
    "total_gibbs_energy",
    "vle_initialization",
]
