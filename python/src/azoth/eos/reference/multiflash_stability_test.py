"""``eos.multiflash_stability_test``: the Python reference.

The mirror of ``crates/azoth-eos/src/multiflash_stability_test.rs``. NeqSim has **two** methods
named ``stabilityAnalysis``: ``Flash``'s, which is a two-trial Michelsen test that reports
instability by setting K-values and which :mod:`azoth.eos.reference.stability_test` ports, and
``TPmultiflash``'s **override** of it, which ``TPmultiflash.run()`` actually calls. This is the
override.

Its second half - one trial per component, each seeded at that component alone - is what
``eos.tp_multiflash`` was missing, and it is measured rather than inferred: on NeqSim's own
water-dew-point gas the pure-water trial reaches ``x_water = 0.999993878160488``, which is the
composition NeqSim's own seeding reports.
"""

from __future__ import annotations

import math
from typing import NamedTuple

from azoth import _models_gen
from azoth.core.errors import InvalidInputError, PropertyUnavailableError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import MultiflashSeedPath, MultiflashStabilityTestResult
from azoth.core.units import Q, input_to_si
from azoth.core.warnings import Warning
from azoth.eos.mixture import Mixture
from azoth.eos.reference._mixture_state import (
    PhaseState,
    ReducedParameters,
    phase_state,
    reduced_parameters,
    wilson_k,
)
from azoth.eos.reference.stability_test import _feed_state, _ln

MODEL_ID = "eos.multiflash_stability_test"

#: A component at or below this feed fraction is not in the fluid at all. Upstream's.
PRESENT = 1.0e-100

#: The logarithm a component that is not present carries; ``exp`` of it underflows to zero.
SENTINEL = -10000.0

#: The Wilson K-value floor, upstream's.
WILSON_FLOOR = 1.0e-20

#: Below the largest ``|ln K|`` this, the Wilson pair is skipped.
WILSON_NEUTRAL = 0.01

#: The mole number every component but the trial's own is seeded at.
TRACE = 1.0e-12

#: A trial this close to the feed, in sum-absolute composition, *is* the feed.
TRIVIAL = 1.0e-4

#: Below this a trial is below the tangent plane, and the phase it found is real.
TM_LIMIT = -1.0e-8

#: The distance a trial that failed or went trivial is marked with. Far above any real one.
FAILED_TM = 10.0

#: Steps the Wilson pair may take, upstream's ``maxiter``.
WILSON_CAP = 50

#: Steps a pure-component trial may take, upstream's ``maxiter``.
PURE_CAP = 200

#: The change in ``ln W`` a trial must fall below to have settled - **a sum of absolute
#: changes, not their root mean square**, which is upstream's ``err``.
TOLERANCE = 1.0e-9


class _Settled:
    """One trial's outcome."""

    __slots__ = ("composition", "failed", "iterations", "totals")

    def __init__(
        self, composition: list[float], totals: float, iterations: int, failed: bool
    ) -> None:
        self.composition = composition
        self.totals = totals
        self.iterations = iterations
        self.failed = failed


def _trial(
    reduced: ReducedParameters,
    kij: tuple[tuple[float, ...], ...],
    d: list[float],
    valid: list[bool],
    seed: list[float],
    *,
    liquid: bool,
    cap: int,
) -> _Settled:
    """Iterate one trial to a stationary point of the tangent-plane distance.

    Unlike ``_trial`` in :mod:`azoth.eos.reference.stability_test` this **does not raise** on
    its cap: upstream computes the distance from whatever iterate it has and leaves the
    decision to the threshold.
    """
    w = list(seed)
    total = sum(w)
    w = [value / total for value in w]
    ln_w = [_ln(value) for value in w]
    previous = math.inf
    iterations = 0
    totals = 0.0

    for step in range(1, cap + 1):
        iterations = step
        try:
            state: PhaseState = phase_state(reduced, kij, w, liquid=liquid)
        except Exception:
            return _Settled(w, 0.0, iterations, True)
        fresh = [
            d[index] - state.ln_phi[index] if valid[index] else SENTINEL for index in range(len(w))
        ]
        error = sum(abs(a - b) for a, b in zip(fresh, ln_w, strict=True))
        mole_numbers = [math.exp(value) for value in fresh]
        totals = sum(mole_numbers)
        total = totals
        w = [value / total for value in mole_numbers]
        ln_w = fresh
        if error <= TOLERANCE and error <= previous:
            return _Settled(w, totals, iterations, False)
        previous = error

    return _Settled(w, totals, iterations, False)


def _is_hydrocarbon(mixture: Mixture, index: int) -> bool:
    """Whether a component counts as a hydrocarbon for the fallback's skip rule.

    NeqSim's ``isHydrocarbon`` is ``HC`` **or** a TBP or plus fraction, and this library has
    no pseudo-components, so the class is the whole of it. It is *not* the predicate
    :func:`azoth.eos.reference.hydrate_inhibitor_wt.phase_label` uses, which counts inerts on
    the same side - the two differ for carbon dioxide.
    """
    return mixture.components[index].component_class == "hc"


def _pure_trials(mixture: Mixture, z: list[float]) -> list[int]:
    """Which components the pure-component fallback runs, in the order it visits them."""
    count = len(mixture.components)
    heaviest = 0
    lightest = 0
    heaviest_mass = 0.0
    lightest_mass = math.inf
    for index in range(count):
        if not _is_hydrocarbon(mixture, index):
            continue
        molar_mass = mixture.components[index].molar_mass
        if molar_mass is None:
            raise PropertyUnavailableError(
                "component", "molar mass", "a card-added component needs its own molar mass"
            )
        mass = float(molar_mass.to("kg/mol").magnitude)
        if mass > heaviest_mass:
            heaviest_mass = mass
            heaviest = index
        if mass < lightest_mass:
            lightest_mass = mass
            lightest = index

    out: list[int] = []
    for index in reversed(range(count)):
        if z[index] < PRESENT:
            continue
        if _is_hydrocarbon(mixture, index) and index not in (heaviest, lightest):
            continue
        out.append(index)
    return out


def _seed_of(settled: _Settled, existing: list[list[float]]) -> float:
    """The distance a trial reports, or ``FAILED_TM`` where it failed or went trivial.

    Upstream's ``xTrivialCheck`` is the sum-absolute distance from **each existing phase**,
    and a trial that reproduces any of them is the plane through that phase rather than a
    new one. The single-phase model is the case with one entry.
    """
    if settled.failed:
        return FAILED_TM
    trivial = any(
        sum(abs(a - b) for a, b in zip(settled.composition, phase, strict=True)) < TRIVIAL
        for phase in existing
    )
    return FAILED_TM if trivial else 1.0 - settled.totals


class _EnhancedSeed(NamedTuple):
    """The phase an enhanced trial admits."""

    composition: list[float]
    liquid: bool
    fraction: float


#: Steps an enhanced trial may take, upstream's ``maxiter``.
ENHANCED_CAP = 150


def enhanced_against_phases(
    mixture: Mixture,
    t_si: float,
    p_si: float,
    z: list[float],
    phases: list[tuple[list[float], bool]],
) -> _EnhancedSeed | None:
    """``stabilityAnalysisEnhanced``, upstream's second analysis.

    For **each existing phase as the reference**, three trials are started - Wilson ``K`` on
    the vapour root, ``z * (hydrocarbon ? 0.5 : 2.0)`` on the liquid root, and ``1 / K`` on the
    liquid root, in that order - and the first whose tangent-plane distance is below ``-1e-8``
    and whose composition is not one of the existing phases is admitted. The starting fraction
    is the feed's mole fraction of the admitted composition's dominant component.

    **The Wegstein step every seventh iteration is not carried**, for the same reason the basic
    analysis' acceleration is not: it changes the path, not the fixed point, and the spec names
    the omission.
    """
    count = len(mixture.components)
    reduced = reduced_parameters(mixture, t_si, p_si)
    kij = mixture.effective_kij(reduced.t_kelvin)
    valid = [value > PRESENT for value in z]
    k = wilson_k(mixture, t_si, p_si)
    existing = [list(composition) for composition, _ in phases]

    for reference, liquid in phases:
        try:
            state = phase_state(reduced, kij, list(reference), liquid=liquid)
        except Exception:
            continue
        d = [
            _ln(fraction) + lp if fraction > PRESENT else 0.0
            for fraction, lp in zip(reference, state.ln_phi, strict=True)
        ]
        seeds = (
            (
                [max(k[i], WILSON_FLOOR) if valid[i] else 0.0 for i in range(count)],
                False,
            ),
            (
                [
                    z[i] * (0.5 if _is_hydrocarbon(mixture, i) else 2.0) if valid[i] else 0.0
                    for i in range(count)
                ],
                True,
            ),
            (
                [1.0 / max(k[i], WILSON_FLOOR) if valid[i] else 0.0 for i in range(count)],
                True,
            ),
        )
        for composition, trial_liquid in seeds:
            settled = _trial(
                reduced, kij, d, valid, composition, liquid=trial_liquid, cap=ENHANCED_CAP
            )
            if settled.failed:
                continue
            if any(
                sum(abs(a - b) for a, b in zip(settled.composition, phase, strict=True)) < TRIVIAL
                for phase in existing
            ):
                continue
            if 1.0 - settled.totals < TM_LIMIT:
                return _EnhancedSeed(
                    list(settled.composition),
                    trial_liquid,
                    z[_dominant(settled.composition)],
                )
    return None


def _dominant(values: list[float]) -> int:
    """The index of the largest entry, upstream's ``dominantComp``."""
    best = 0
    for index, value in enumerate(values):
        if value > values[best]:
            best = index
    return best


def multiflash_stability_test(
    mixture: Mixture, T: Q, P: Q, z: list[float]
) -> MultiflashStabilityTestResult:
    """The phase to add to make a fluid a Gibbs minimum, by the multiflash's own seeding.

    Raises:
        InvalidInputError: if ``z`` is the wrong length, has a negative entry, or does not sum
            to one.
        OutOfRangeError: if ``T`` or ``P`` is not positive.

    See :func:`azoth.eos.reference.multiflash_stability_test`.
    """
    spec = _models_gen.model(MODEL_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []
    t_si = float(input_to_si(spec, "T", T))
    p_si = float(input_to_si(spec, "P", P))
    apply_checks(checks.on_input, {"T": t_si, "P": p_si}.get, warnings)

    count = len(mixture.components)
    if len(z) != count:
        raise InvalidInputError("z", f"a feed for {count} components has {len(z)} entries")
    if any(value < 0.0 for value in z):
        raise InvalidInputError("z", "a mole fraction cannot be negative")
    if abs(sum(z) - 1.0) > 1.0e-9:
        raise InvalidInputError("z", f"the mole fractions sum to {sum(z)}, not to one")

    return seeding_against_reference(
        mixture,
        t_si,
        p_si,
        z,
        reference=list(z),
        reference_liquid=None,
        existing=[list(z)],
        warnings=warnings,
    )


def seeding_against_reference(
    mixture: Mixture,
    t_si: float,
    p_si: float,
    z: list[float],
    *,
    reference: list[float],
    reference_liquid: bool | None,
    existing: list[list[float]],
    warnings: list[Warning],
) -> MultiflashStabilityTestResult:
    """The analysis measured against a reference phase rather than the feed.

    NeqSim's ``TPmultiflash.stabilityAnalysis`` takes its reference potentials from **phase 0
    of the flashed system** - the least-dense existing phase, because ``TPflash`` orders by
    density - and marks a trial trivial when it reproduces **any** existing phase. The public
    model above is the single-phase case, where phase 0 is the feed; ``eos.tp_multiflash``
    calls this with the two-phase set and the reference the density order names.
    """
    spec = _models_gen.model(MODEL_ID)
    algorithm = spec["algorithm"]
    count = len(mixture.components)
    reduced = reduced_parameters(mixture, t_si, p_si)
    kij = mixture.effective_kij(reduced.t_kelvin)

    valid = [value > PRESENT for value in z]
    state = (
        _feed_state(reduced, kij, list(reference))
        if reference_liquid is None
        else phase_state(reduced, kij, list(reference), liquid=reference_liquid)
    )
    d = [
        _ln(fraction) + lp if fraction > PRESENT else 0.0
        for fraction, lp in zip(reference, state.ln_phi, strict=True)
    ]

    tm: list[float] = []
    iterations = 0
    k = wilson_k(mixture, t_si, p_si)
    neutral = max(
        (max(value, WILSON_FLOOR) for value, held in zip(k, valid, strict=True) if held),
        default=0.0,
    )
    neutral = abs(math.log(max(neutral, WILSON_FLOOR)))

    if neutral >= WILSON_NEUTRAL:
        liquid_seed = [
            z[index] / max(k[index], WILSON_FLOOR) if valid[index] else 0.0
            for index in range(count)
        ]
        vapour_seed = [
            z[index] * max(k[index], WILSON_FLOOR) if valid[index] else 0.0
            for index in range(count)
        ]
        pair = (
            (liquid_seed, True, MultiflashSeedPath.WILSON_LIQUID),
            (vapour_seed, False, MultiflashSeedPath.WILSON_VAPOUR),
        )
        cap = min(WILSON_CAP, int(algorithm["max_iterations"]))
        for seed, liquid, path in pair:
            settled = _trial(reduced, kij, d, valid, seed, liquid=liquid, cap=cap)
            iterations += settled.iterations
            distance = _seed_of(settled, existing)
            tm.append(distance)
            if distance < TM_LIMIT:
                return MultiflashStabilityTestResult(
                    path=path,
                    tm=tuple(tm),
                    seed_composition=tuple(settled.composition),
                    seed_fraction=z[_dominant(seed)],
                    seed_component_index=-1,
                    iterations=iterations,
                    min_t_over_tc=_min_t_over_tc(mixture, t_si),
                    warnings=tuple(warnings),
                )

    for index in _pure_trials(mixture, z):
        seed = [
            0.0 if not valid[other] else (1.0 if other == index else TRACE)
            for other in range(count)
        ]
        settled = _trial(reduced, kij, d, valid, seed, liquid=True, cap=PURE_CAP)
        iterations += settled.iterations
        distance = _seed_of(settled, existing)
        tm.append(distance)
        if distance < TM_LIMIT:
            return MultiflashStabilityTestResult(
                path=MultiflashSeedPath.PURE_COMPONENT,
                tm=tuple(tm),
                seed_composition=tuple(settled.composition),
                seed_fraction=z[index],
                seed_component_index=index,
                iterations=iterations,
                min_t_over_tc=_min_t_over_tc(mixture, t_si),
                warnings=tuple(warnings),
            )

    return MultiflashStabilityTestResult(
        path=MultiflashSeedPath.STABLE,
        tm=tuple(tm),
        seed_composition=(),
        seed_fraction=0.0,
        seed_component_index=-1,
        iterations=iterations,
        min_t_over_tc=_min_t_over_tc(mixture, t_si),
        warnings=tuple(warnings),
    )


def _min_t_over_tc(mixture: Mixture, temperature: float) -> float:
    """The smallest ``T / Tc`` over the components, in kelvin over kelvin."""
    return min(temperature / component.Tc.to("K").magnitude for component in mixture.components)
