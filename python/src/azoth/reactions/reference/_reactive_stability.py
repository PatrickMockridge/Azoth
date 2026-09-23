"""The reactive tangent-plane stability analysis, the Python twin of
``crates/azoth-reactions/src/reactive_stability.rs``.

It decides whether a second phase forms, in four steps: bring the feed to **homogeneous**
chemical equilibrium, take the reference potentials ``d_i = ln x_i + ln phi_i`` from *that*
composition, seed trial phases with Wilson K-values and pure components, and run a
tangent-plane trial from each. A trial whose distance comes in under
:data:`TPD_THRESHOLD` is a phase the fluid is unstable with respect to.

The trial is a successive substitution on the log-composition:

.. code-block:: text

    logW_i = d_i - ln phi_i(W)          unnormalised, iterated to SS_TOL
    TPD    = 1 - sum_i exp(logW_i)      the tangent-plane distance at the stationary point

**Two answers are not distances and the class reports both as the same number.** A trial
that walks back to the reference's own composition is the *trivial* solution, and a trial
whose phase cannot be initialised is not a distance either; both come back as
:data:`STABLE_TPD` = ``10.0``, which is how the caller reads "this seed found nothing".

**The reference potentials have three branches and the middle one is easy to get wrong**: a
component the reference holds takes the logarithm, one at or below :data:`MIN_MOLES` takes
:data:`REFERENCE_ABSENT` (``-100.0``, *not* ``ln`` of the floor), and an ion takes
:data:`REFERENCE_ION`.
"""

from __future__ import annotations

import math
from collections.abc import Callable
from typing import NamedTuple

from azoth.core.errors import InvalidInputError

#: The successive-substitution tolerance, from ``SS_TOL``.
SS_TOL = 1.0e-9

#: The pass cap, from ``MAX_SS_ITER``.
MAX_SS_ITERATIONS = 100

#: The distance below which a trial counts as a new phase, from ``TPD_THRESHOLD``.
TPD_THRESHOLD = -1.0e-8

#: The floor a trial mole fraction is kept above, from ``MIN_MOLES``.
MIN_MOLES = 1.0e-30

#: What a trial reports when it is stable, trivial, or cannot be initialised.
STABLE_TPD = 10.0

#: The composition difference below which a trial is the reference's own solution.
TRIVIAL_TOLERANCE = 1.0e-4

#: The mole fraction above which the reference holds a component - **a looser test than
#: :data:`MIN_MOLES`** - and the one that decides whether a trial's log-composition is
#: updated for that component at all.
TRIAL_ABSENT_FLOOR = 1.0e-100

#: The value given a component the reference does not hold, and the one given an ion.
REFERENCE_ABSENT = -100.0
REFERENCE_ION = -1000.0


class CriticalConstants(NamedTuple):
    """One component's critical constants, which the Wilson seed reads."""

    #: The critical temperature, in K.
    tc: float
    #: The critical pressure, in **bara** - the unit the Wilson correlation is written in.
    pc: float
    #: The acentric factor.
    omega: float


def wilson_k(constants: CriticalConstants, temperature: float, pressure: float) -> float:
    """``(pc / P) exp(5.373 (1 + omega) (1 - tc / T))``, floored at ``1e-20``.

    **A component with no critical constants gets a K of 1**, which is the class's own
    fallback and the reason a fluid of pseudo-components still produces trials.
    """
    if constants.pc > 0.0 and constants.tc > 0.0:
        return max(
            (constants.pc / pressure)
            * math.exp(5.373 * (1.0 + constants.omega) * (1.0 - constants.tc / temperature)),
            1.0e-20,
        )
    return 1.0


def trial_seeds(
    fractions: list[float],
    constants: list[CriticalConstants],
    temperature: float,
    pressure: float,
) -> list[list[float]]:
    """``generateTrialPhases``: a liquid-like and a vapour-like trial, then one per component.

    The first two are built **only when some ``|ln K|`` exceeds ``0.01``**, and the pure
    component trials - ``1.0`` on one component and ``1e-12`` on the rest - are built for
    every component the feed holds. The vectors are what a trial *starts* from, not
    fractions.
    """
    nc = len(fractions)
    k = [wilson_k(entry, temperature, pressure) for entry in constants]
    all_near_one = all(abs(math.log(value)) <= 0.01 for value in k)

    seeds: list[list[float]] = []
    if not all_near_one:
        seeds.append(
            [fractions[i] / k[i] if fractions[i] > MIN_MOLES else MIN_MOLES for i in range(nc)]
        )
        seeds.append(
            [k[i] * fractions[i] if fractions[i] > MIN_MOLES else MIN_MOLES for i in range(nc)]
        )
    for j in range(nc):
        if fractions[j] > TRIAL_ABSENT_FLOOR:
            seeds.append([1.0 if i == j else 1.0e-12 for i in range(nc)])
    return seeds


def reference_potentials(
    fractions: list[float], ln_phi: list[float], charges: list[float]
) -> list[float]:
    """``computeReferencePotentials``, with all three of its branches.

    Raises:
        InvalidInputError: on a shape disagreement.
    """
    if len(fractions) != len(ln_phi) or len(fractions) != len(charges):
        raise InvalidInputError(
            "fractions",
            f"{len(fractions)} composition entr(ies) against {len(ln_phi)} fugacity "
            f"coefficient(s) and {len(charges)} charge(s)",
        )
    out: list[float] = []
    for i, x in enumerate(fractions):
        if charges[i] != 0.0:
            out.append(REFERENCE_ION)
        elif x > MIN_MOLES:
            out.append(math.log(x) + ln_phi[i])
        else:
            out.append(REFERENCE_ABSENT)
    return out


def run_trial(
    d: list[float],
    seed: list[float],
    fractions: list[float],
    ln_phi: Callable[[list[float]], list[float]],
) -> float:
    """One tangent-plane trial, returning its distance or :data:`STABLE_TPD`."""
    nc = len(d)
    log_w = [0.0] * nc
    sum_w = 0.0
    for i in range(nc):
        w = max(seed[i], MIN_MOLES)
        log_w[i] = math.log(w)
        sum_w += w
    fractions_w = [math.exp(log_w[i]) / sum_w for i in range(nc)]

    for _ in range(MAX_SS_ITERATIONS):
        previous = list(log_w)
        ln_phi_here = ln_phi(fractions_w)

        error = 0.0
        sum_w = 0.0
        for i in range(nc):
            # A component the *reference* does not hold keeps the log-composition it had:
            # the class updates only where `x > 1e-100`, and still counts the value it kept
            # in the error and the sum.
            if fractions[i] > TRIAL_ABSENT_FLOOR:
                log_w[i] = d[i] - ln_phi_here[i]
            error += abs(log_w[i] - previous[i])
            sum_w += math.exp(log_w[i])
        if not math.isfinite(sum_w) or sum_w <= 0.0:
            return STABLE_TPD
        fractions_w = [math.exp(log_w[i]) / sum_w for i in range(nc)]
        if error < SS_TOL:
            break

    tpd = 1.0
    for i in range(nc):
        if fractions[i] > TRIAL_ABSENT_FLOOR:
            tpd -= math.exp(log_w[i])

    trivial = sum(abs(fractions_w[i] - fractions[i]) for i in range(nc))
    if trivial < TRIVIAL_TOLERANCE:
        return STABLE_TPD
    return tpd


def is_unstable(tpd: float) -> bool:
    """Whether a distance counts as a new phase."""
    return tpd < TPD_THRESHOLD


class StabilityOutcome(NamedTuple):
    """What ``ReactiveStabilityAnalysis.run`` answers with."""

    #: The composition the reference potentials came from.
    reference: list[float]
    #: ``d_i``, the reference potentials.
    potentials: list[float]
    #: Each unstable trial **brought to chemical equilibrium**, in seed order.
    unstable_trials: list[list[float]]
    #: The distances of those trials, one each, in the same order.
    tpd_values: list[float]
    #: Whether any trial came in under :data:`TPD_THRESHOLD`.
    unstable: bool


def solve_trial_equilibrium(
    trial: list[float], ce: Callable[[list[float]], list[float]]
) -> list[float]:
    """``solveTrialChemicalEquilibrium``: normalise the trial, then bring it to equilibrium."""
    floored = [max(w, MIN_MOLES) for w in trial]
    total = sum(floored)
    return ce([w / total for w in floored])


def analyse(
    fractions: list[float],
    constants: list[CriticalConstants],
    temperature: float,
    pressure: float,
    charges: list[float],
    ce: Callable[[list[float]], list[float]],
    ln_phi: Callable[[list[float]], list[float]],
) -> StabilityOutcome:
    """``ReactiveStabilityAnalysis.run``: the four steps, composed.

    ``ce`` is the single-phase reactive solve, which **returns its input where it did not
    converge** - the class's own behaviour in both of the places it calls one - and
    ``ln_phi`` is the single-phase model the trials are evaluated against.

    Raises:
        InvalidInputError: where ``ce`` returns a composition of the wrong length.
    """
    reference = ce(fractions)
    if len(reference) != len(fractions) or len(charges) != len(fractions):
        raise InvalidInputError(
            "ce", f"{len(reference)} entries back against {len(fractions)} composition(s)"
        )

    ln_phi_reference = ln_phi(reference)
    potentials = reference_potentials(reference, ln_phi_reference, charges)

    unstable_trials: list[list[float]] = []
    tpd_values: list[float] = []
    for seed in trial_seeds(reference, constants, temperature, pressure):
        tpd = run_trial(potentials, seed, reference, ln_phi)
        if is_unstable(tpd):
            unstable_trials.append(solve_trial_equilibrium(seed, ce))
            tpd_values.append(tpd)

    return StabilityOutcome(
        reference=reference,
        potentials=potentials,
        unstable_trials=unstable_trials,
        tpd_values=tpd_values,
        unstable=bool(unstable_trials),
    )
