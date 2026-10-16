"""``eos.co2_brine_phase_equilibrium`` - NeqSim's constrained CO2/brine phase selection.

Spec: ``specs/models/eos/co2_brine_phase_equilibrium.toml``. A *procedure*: a reference
phase, two normalised tangent-plane trials, and a split where a trial says the feed is
unstable.

This is the pure-Python reference: a second, independent expression of the same physics as
the Rust kernel, held to it by ``test_cross_impl.py``. The fluid is
``electrolyte_cpa_statoil_mixture_of``'s and the phase states are ``_mixture_state``'s, so
what is here is the class's own solve.

The Gibbs comparison is written as a *mixing* term - see the Rust module's docs, which state
why ordering candidates by ``sum_k beta_k sum_i x_ik ln(x_ik phi_ik)`` orders them exactly as
the absolute Gibbs energy upstream compares does, and which guard that leaves unported.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from dataclasses import dataclass

from azoth import _models_gen
from azoth.core.errors import InvalidInputError, SolverNotConvergedError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import Co2BrinePhaseEquilibriumResult
from azoth.core.units import Q, input_to_si
from azoth.core.warnings import Warning
from azoth.eos.components import electrolyte_cpa_statoil_mixture_of, entry
from azoth.eos.reference._mixture_state import phase_state, reduced_parameters

MODEL_ID = "eos.co2_brine_phase_equilibrium"

#: The class's `FLOOR`: what an ion's mole fraction is written as outside the aqueous phase.
_FLOOR = 1.0e-50

#: The class's `TOLERANCE`, which decides whether a trial is unstable and how far apart the
#: two phases' molecular fugacities may be.
_TOLERANCE = 1.0e-8

#: The composition change that ends the trial and the split.
_CHANGE_TOLERANCE = 1.0e-12

_MAX_TRIAL_STEPS = 100
_MAX_SPLIT_STEPS = 150
_BISECTION_STEPS = 80

#: The class's seed for the rich phase's fraction.
_SEED_BETA = 0.1

#: The upper end of the `beta` bracket: one is the single-aqueous answer and is excluded.
_BETA_CEILING = 1.0 - 1.0e-12

#: The worst `|z_i - sum_k beta_k x_ik|` an accepted answer may have.
_CONSERVATION_TOLERANCE = 1.0e-10


@dataclass(frozen=True, slots=True)
class _Layer:
    """One phase's state as the solve reads it: a composition and its log coefficients."""

    x: tuple[float, ...]
    ln_phi: tuple[float, ...]
    ln_x_phi: tuple[float, ...]

    @classmethod
    def of(
        cls,
        reduced: object,
        kij: tuple[tuple[float, ...], ...],
        x: Sequence[float],
        *,
        liquid: bool,
    ) -> _Layer:
        state = phase_state(reduced, kij, list(x), liquid=liquid)
        return cls(
            tuple(x),
            tuple(state.ln_phi),
            tuple(
                math.log(fraction) + ln_phi
                for fraction, ln_phi in zip(x, state.ln_phi, strict=True)
            ),
        )

    def mixing_term(self) -> float:
        """``sum_i x_i ln(x_i phi_i)``, one phase's part of the mixing Gibbs.

        A component the phase does not carry contributes nothing rather than ``0 * -inf``.
        """
        return sum(
            fraction * term
            for fraction, term in zip(self.x, self.ln_x_phi, strict=True)
            if fraction > 0.0
        )


def _mixing_gibbs(layers: Sequence[_Layer], fractions: Sequence[float]) -> float:
    """The feed's mixing Gibbs over ``R T``, per mole of feed."""
    return sum(beta * layer.mixing_term() for layer, beta in zip(layers, fractions, strict=True))


def _rachford_rice(z: Sequence[float], ratios: Sequence[float], beta: float) -> float:
    """The constrained Rachford-Rice residual."""
    return sum(
        overall * (ratio - 1.0) / (1.0 + beta * (ratio - 1.0))
        for overall, ratio in zip(z, ratios, strict=True)
    )


def _unqualified(iterations: int, residual: float, tolerance: float) -> SolverNotConvergedError:
    """The error a state that cannot be qualified gets."""
    return SolverNotConvergedError(iterations, residual, tolerance)


def co2_brine_phase_equilibrium(
    components: Sequence[str],
    T: Q,
    P: Q,
    moles: Sequence[float],
) -> Co2BrinePhaseEquilibriumResult:
    """The constrained equilibrium of a CO2/water/ion feed with the electrolyte-CPA EoS.

    Args:
        components: the feed's substances, by name: CO2, water and the ions.
        T: absolute temperature.
        P: absolute pressure.
        moles: the feed's mole numbers.

    Returns:
        The phase count, each phase's root, fraction and composition, and the trial distance
        the answer was decided on.

    Raises:
        InvalidInputError: if ``moles`` is not one entry per component or is not a feed, or
            the feed is outside the class's scope.
        PropertyUnavailableError: if a name is not in the databank.
        SolverNotConvergedError: if a trial or a split does not converge, or the answer fails
            its acceptance checks.

    See :func:`azoth.eos.co2_brine_phase_equilibrium`.
    """
    spec = _models_gen.model(MODEL_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    t_si = input_to_si(spec, "T", T)
    p_si = input_to_si(spec, "P", P)
    apply_checks(checks.on_input, {"T": t_si, "P": p_si}.get, warnings)

    moles = [input_to_si(spec, "moles", value) for value in moles]
    if len(components) != len(moles):
        raise InvalidInputError(
            "moles",
            f"{len(components)} components but {len(moles)} mole numbers; the two must match",
        )
    if any(not math.isfinite(m) or m < 0.0 for m in moles):
        raise InvalidInputError("moles", "every mole number must be finite and not negative")
    total = sum(moles)
    if not total > 0.0 or not math.isfinite(total):
        raise InvalidInputError(
            "moles", f"the feed's mole numbers total {total}, which is not a feed"
        )
    z = [m / total for m in moles]

    resolved = [name.strip().lower() for name in components]
    entries = [entry(name) for name in resolved]
    is_ion = [e.ionic_charge != 0.0 for e in entries]
    try:
        co2 = next(i for i, name in enumerate(resolved) if name == "co2")
    except StopIteration:
        raise _unsupported("it has no CO2") from None
    try:
        water = next(i for i, name in enumerate(resolved) if name == "water")
    except StopIteration:
        raise _unsupported("it has no water") from None
    if any(
        not is_ion[i] and name not in ("co2", "water") and moles[i] > 0.0
        for i, name in enumerate(resolved)
    ):
        raise _unsupported("it carries a substance that is neither CO2, water nor an ion")
    if not moles[co2] > 0.0 or not moles[water] > moles[co2]:
        raise _unsupported("its water is not in excess of its CO2")

    fluid, _ = electrolyte_cpa_statoil_mixture_of(list(components))
    reduced = reduced_parameters(fluid, t_si, p_si)
    warnings.extend(reduced.warnings)
    kij = fluid.effective_kij(t_si)

    reference = _Layer.of(reduced, kij, z, liquid=True)

    iterations = 0
    distance = math.inf
    best: tuple[str, float, _Layer, _Layer, float] | None = None
    for label, liquid in (("gas", False), ("oil", True)):
        trial_distance, seed, trial_steps = _normalised_trial(
            reduced, kij, z, is_ion, reference, liquid=liquid
        )
        iterations += trial_steps
        distance = min(distance, trial_distance)
        if trial_distance >= -_TOLERANCE:
            continue
        split = _solve_split(reduced, kij, z, is_ion, seed, liquid=liquid)
        if split is None:
            raise _unqualified(_MAX_SPLIT_STEPS, math.nan, _CHANGE_TOLERANCE)
        beta, rich, aqueous, split_steps = split
        iterations += split_steps
        gibbs = _mixing_gibbs([rich, aqueous], [beta, 1.0 - beta])
        if best is None or gibbs < best[1]:
            best = (label, gibbs, rich, aqueous, beta)

    if best is None:
        accepted_root = "aqueous"
        phase_fractions = [1.0]
        layers: list[_Layer] = [reference]
    else:
        label, _, rich, aqueous, beta = best
        accepted_root = label
        phase_fractions = [beta, 1.0 - beta]
        layers = [rich, aqueous]

    max_material_balance_residual = max(
        abs(
            z[i]
            - sum(beta * layer.x[i] for layer, beta in zip(layers, phase_fractions, strict=True))
        )
        for i in range(len(z))
    )
    if max_material_balance_residual > _CONSERVATION_TOLERANCE:
        raise _unqualified(iterations, max_material_balance_residual, _CONSERVATION_TOLERANCE)

    max_log_fugacity_residual = 0.0
    if len(layers) == 2:
        if layers[0].x[co2] <= 0.5:
            raise _unqualified(iterations, 0.5 - layers[0].x[co2], 0.5)
        for i in range(len(z)):
            if is_ion[i] or z[i] <= _FLOOR:
                continue
            ratio = layers[0].ln_x_phi[i] - layers[1].ln_x_phi[i]
            max_log_fugacity_residual = max(max_log_fugacity_residual, abs(ratio))
            if abs(ratio) > _TOLERANCE:
                raise _unqualified(iterations, max_log_fugacity_residual, _TOLERANCE)

    return Co2BrinePhaseEquilibriumResult(
        phase_count=len(phase_fractions),
        accepted_root=accepted_root,
        phase_fractions=tuple(phase_fractions),
        x=tuple(layer.x for layer in layers),
        ln_phi=tuple(layer.ln_phi for layer in layers),
        minimum_trial_distance=distance,
        iterations=iterations,
        max_material_balance_residual=max_material_balance_residual,
        max_log_fugacity_residual=max_log_fugacity_residual,
        warnings=tuple(warnings),
    )


def _normalised_trial(
    reduced: object,
    kij: tuple[tuple[float, ...], ...],
    z: Sequence[float],
    is_ion: Sequence[bool],
    reference: _Layer,
    *,
    liquid: bool,
) -> tuple[float, list[float], int]:
    """The normalised tangent-plane trial on one root.

    ``w_i = z_i phi_i^reference / phi_i^trial`` for a molecular component and zero for an
    ion, which is the class's ``f_i^aqueous / (P phi_i^trial)`` with the ``P`` cancelled.
    ``x_i`` is ``max(FLOOR, w_i / sum w)`` and the distance is ``-ln(sum w)``.
    """
    x = [_FLOOR if is_ion[i] else (0.9999 if i == 0 else 0.0001) for i in range(len(z))]
    distance = 0.0
    change = 0.0
    for step in range(1, _MAX_TRIAL_STEPS + 1):
        trial = _Layer.of(reduced, kij, x, liquid=liquid)
        weights = [0.0] * len(z)
        for i in range(len(z)):
            if is_ion[i] or z[i] <= _FLOOR:
                continue
            weight = z[i] * math.exp(reference.ln_phi[i]) / math.exp(trial.ln_phi[i])
            if not weight > 0.0 or not math.isfinite(weight):
                raise _unqualified(step, weight, _CHANGE_TOLERANCE)
            weights[i] = weight
        total = sum(weights)
        if not total > 0.0 or not math.isfinite(total):
            raise _unqualified(step, total, _CHANGE_TOLERANCE)
        change = 0.0
        for i in range(len(z)):
            updated = max(_FLOOR, weights[i] / total)
            change = max(change, abs(updated - x[i]))
            x[i] = updated
        distance = -math.log(total)
        if change < _CHANGE_TOLERANCE:
            return distance, x, step
    raise _unqualified(_MAX_TRIAL_STEPS, change, _CHANGE_TOLERANCE)


def _solve_split(
    reduced: object,
    kij: tuple[tuple[float, ...], ...],
    z: Sequence[float],
    is_ion: Sequence[bool],
    seed: Sequence[float],
    *,
    liquid: bool,
) -> tuple[float, _Layer, _Layer, int] | None:
    """The ion-confined split, or ``None`` where no bracket exists.

    ``K_i`` is ``phi_i^aqueous / phi_i^rich`` and an ion's is zero. ``beta`` is a bisection on
    the constrained Rachford-Rice - not a Newton, because the zero ``K`` puts the function's
    curvature on the wrong side for one.
    """
    n = len(z)
    rich = _Layer.of(reduced, kij, seed, liquid=liquid)
    aqueous = _Layer.of(reduced, kij, z, liquid=True)
    beta = _SEED_BETA
    for step in range(1, _MAX_SPLIT_STEPS + 1):
        ratios = [0.0] * n
        for i in range(n):
            if is_ion[i]:
                continue
            ratio = math.exp(aqueous.ln_phi[i]) / math.exp(rich.ln_phi[i])
            if not ratio > 0.0 or not math.isfinite(ratio):
                return None
            ratios[i] = ratio
        if (
            not _rachford_rice(z, ratios, 0.0) > 0.0
            or not _rachford_rice(z, ratios, _BETA_CEILING) < 0.0
        ):
            return None
        lower, upper = 0.0, _BETA_CEILING
        for _ in range(_BISECTION_STEPS):
            middle = 0.5 * (lower + upper)
            if _rachford_rice(z, ratios, middle) > 0.0:
                lower = middle
            else:
                upper = middle
        updated_beta = 0.5 * (lower + upper)
        change = abs(updated_beta - beta)
        beta = updated_beta

        aqueous_x: list[float] = []
        rich_x: list[float] = []
        for i in range(n):
            liquid_x = z[i] / (1.0 + beta * (ratios[i] - 1.0))
            rich_value = max(_FLOOR, ratios[i] * liquid_x)
            change = max(change, abs(liquid_x - aqueous.x[i]))
            change = max(change, abs(rich_value - rich.x[i]))
            aqueous_x.append(max(_FLOOR, liquid_x))
            rich_x.append(rich_value)
        aqueous = _Layer.of(reduced, kij, aqueous_x, liquid=True)
        rich = _Layer.of(reduced, kij, rich_x, liquid=liquid)
        if change < _CHANGE_TOLERANCE:
            return beta, rich, aqueous, step
    return None


def _unsupported(why: str) -> InvalidInputError:
    """The refusal a feed outside the class's scope gets."""
    return InvalidInputError(
        "components",
        "the feed is outside `CO2BrinePhaseEquilibrium`'s scope, which is CO2, water and "
        f"ions alone with water in excess: {why}",
    )


__all__ = ["co2_brine_phase_equilibrium"]
