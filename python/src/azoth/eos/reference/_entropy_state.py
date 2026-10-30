"""One equilibrium state at a fixed specific entropy.

The Python twin of ``crates/azoth-eos/src/entropy_state.rs``, written to mirror it.
``eos.equilibrium_sound_speed`` and the HEM release model both re-equilibrate a fluid at a
fixed specific entropy: a pure fluid through its saturation quality, a mixture through a
bracketed TP root. This module is that one arithmetic, with the two roots' stopping
policies as parameters rather than two copies.

The policies are the measured differences between the two classes.
``EquilibriumSoundSpeed`` returns early when its first trial already closes, stops a
bisection on a ``1e-12 max(1, T)`` bracket, and lets a non-finite residual fall through to
the "did not close" report; ``HomogeneousEquilibriumReleaseModel`` does none of the three
and throws on a non-finite residual. The acceptance thresholds stay with the callers.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any

from azoth.core.units import quantity
from azoth.eos.reference._mixture_state import phase_state, reduced_parameters
from azoth.eos.reference.molar_enthalpy_entropy import molar_enthalpy_entropy
from azoth.eos.reference.pr_molar_volume import MOLAR_GAS_CONSTANT
from azoth.eos.reference.pt_flash import pt_flash
from azoth.eos.reference.pure_saturation import pure_saturation

#: The root's target residual, J/(kg K).
ENTROPY_TOLERANCE = 1.0e-7
#: The bisection's relative temperature tolerance.
ROOT_TOLERANCE = 1.0e-12
#: Entropy-bracket widenings before the bracket is refused.
BRACKET_STEPS = 30
#: The bisection cap.
BISECTION_STEPS = 100


@dataclass(frozen=True, slots=True)
class EntropyState:
    """One equilibrium state at the target specific entropy."""

    pressure: float
    temperature: float
    density: float
    molar_enthalpy: float
    molar_entropy: float
    entropy_residual: float
    component_residual: float
    fugacity_residual: float
    phase_composition_residual: float
    phase: str
    vapour_fraction: float | None


@dataclass(frozen=True, slots=True)
class RootPolicy:
    """How one class stops its entropy root."""

    accept_initial: bool
    bracket_width_tolerance: float | None
    nonfinite_is_error: bool


#: ``EquilibriumSoundSpeed``'s root.
SOUND_SPEED = RootPolicy(True, ROOT_TOLERANCE, False)
#: ``HomogeneousEquilibriumReleaseModel``'s root.
HEM = RootPolicy(False, None, True)


class EntropyStateFailure(Exception):
    """The entropy state could not be obtained."""


class Evaluations:
    """The ``tp``/quality evaluations a calculation has spent."""

    __slots__ = ("count",)

    def __init__(self) -> None:
        self.count = 0


def isentropic_state(
    mixture: Any,
    ideal_gas: Any,
    z: list[float],
    molar_mass: float,
    target: float,
    pressure: float,
    start_temperature: float,
    evaluations: Evaluations,
    policy: RootPolicy,
) -> EntropyState:
    """The equilibrium state at ``pressure`` and the target specific entropy.

    ``evaluations`` counts every TP flash and every accepted saturation-quality state,
    including rejected trials, as both classes do.
    """
    if not pressure > 0.0 or not math.isfinite(pressure):
        raise EntropyStateFailure(f"the pressure sample {pressure} Pa is not positive")
    if len(mixture.components) == 1:
        state = _saturation_state(mixture, ideal_gas, z, molar_mass, target, pressure, evaluations)
        if state is not None:
            return state
    return _tp_root(
        mixture,
        ideal_gas,
        z,
        molar_mass,
        target,
        pressure,
        start_temperature,
        evaluations,
        policy,
    )


def _saturation_state(
    mixture: Any,
    ideal_gas: Any,
    z: list[float],
    molar_mass: float,
    target: float,
    pressure: float,
    evaluations: Evaluations,
) -> EntropyState | None:
    """The pure fluid's saturation-quality state, or ``None`` for the TP root."""
    component = mixture.components[0]
    if component.Tc is None or component.Pc is None:
        return None
    tc = component.Tc.to("K").magnitude
    low_t = 0.0
    p_low = math.nan
    for fraction in (0.3, 0.4, 0.5, 0.6):
        candidate = max(1.0, fraction * tc)
        try:
            result = pure_saturation(
                component.Tc, component.Pc, component.omega, quantity(candidate, "K")
            )
        except Exception:
            continue
        low_t = candidate
        p_low = result.p_sat.to("Pa").magnitude
        break
    high_t = 0.0
    p_high = math.nan
    for fraction in (0.94, 0.92, 0.90, 0.85, 0.80):
        try:
            result = pure_saturation(
                component.Tc, component.Pc, component.omega, quantity(fraction * tc, "K")
            )
        except Exception:
            continue
        if result.p_sat.to("Pa").magnitude > pressure:
            high_t = fraction * tc
            p_high = result.p_sat.to("Pa").magnitude
            break
    if not high_t > 0.0 or not low_t > 0.0 or not (p_low < pressure < p_high):
        return None
    lower, upper = low_t, high_t
    temperature = 0.5 * (lower + upper)
    for _ in range(BISECTION_STEPS):
        temperature = 0.5 * (lower + upper)
        try:
            trial = (
                pure_saturation(
                    component.Tc, component.Pc, component.omega, quantity(temperature, "K")
                )
                .p_sat.to("Pa")
                .magnitude
            )
        except Exception:
            return None
        if trial < pressure:
            lower = temperature
        else:
            upper = temperature
        if upper - lower <= ROOT_TOLERANCE * max(1.0, temperature):
            break
    temperature = 0.5 * (lower + upper)
    try:
        reduced = reduced_parameters(mixture, temperature, pressure)
    except Exception:
        return None
    kij = mixture.effective_kij(reduced.t_kelvin)
    try:
        liquid_root = phase_state(reduced, kij, z, liquid=True).z
        vapour_root = phase_state(reduced, kij, z, liquid=False).z
    except Exception:
        return None
    if abs(vapour_root - liquid_root) < 1.0e-12:
        return None
    try:
        liquid = molar_enthalpy_entropy(
            mixture,
            ideal_gas,
            quantity(temperature, "K"),
            quantity(pressure, "Pa"),
            z,
            liquid_root,
        )
        vapour = molar_enthalpy_entropy(
            mixture,
            ideal_gas,
            quantity(temperature, "K"),
            quantity(pressure, "Pa"),
            z,
            vapour_root,
        )
    except Exception:
        return None
    s_liquid = liquid.s.to("J/(mol*K)").magnitude / molar_mass
    s_vapour = vapour.s.to("J/(mol*K)").magnitude / molar_mass
    if not s_vapour > s_liquid or not (s_liquid <= target <= s_vapour):
        return None
    evaluations.count += 1
    quality = (target - s_liquid) / (s_vapour - s_liquid)
    v_liquid = liquid_root * MOLAR_GAS_CONSTANT * temperature / pressure
    v_vapour = vapour_root * MOLAR_GAS_CONSTANT * temperature / pressure
    volume = (1.0 - quality) * v_liquid + quality * v_vapour
    if not volume > 0.0 or not math.isfinite(volume):
        raise EntropyStateFailure(
            f"the saturation-quality volume at {pressure} Pa is {volume} m3/mol"
        )
    molar_entropy = (1.0 - quality) * liquid.s.to("J/(mol*K)").magnitude + quality * vapour.s.to(
        "J/(mol*K)"
    ).magnitude
    molar_enthalpy = (1.0 - quality) * liquid.h.to("J/mol").magnitude + quality * vapour.h.to(
        "J/mol"
    ).magnitude
    return EntropyState(
        pressure=pressure,
        temperature=temperature,
        density=molar_mass / volume,
        molar_enthalpy=molar_enthalpy,
        molar_entropy=molar_entropy,
        entropy_residual=(1.0 - quality) * s_liquid + quality * s_vapour - target,
        component_residual=0.0,
        fugacity_residual=0.0,
        phase_composition_residual=0.0,
        phase="two_phase",
        vapour_fraction=quality,
    )


def _tp_root(
    mixture: Any,
    ideal_gas: Any,
    z: list[float],
    molar_mass: float,
    target: float,
    pressure: float,
    start_temperature: float,
    evaluations: Evaluations,
    policy: RootPolicy,
) -> EntropyState:
    """The equilibrium state at a pressure by the bracketed TP entropy root."""
    guess = start_temperature
    initial = _tp_at(mixture, ideal_gas, z, molar_mass, target, pressure, guess, evaluations)
    if policy.accept_initial and abs(initial.entropy_residual) <= ENTROPY_TOLERANCE:
        return initial
    lower = upper = guess
    low_error = high_error = initial.entropy_residual
    span = max(1.0, guess * 0.005)
    for _ in range(BRACKET_STEPS):
        if not (low_error > 0.0 or high_error < 0.0):
            break
        if low_error > 0.0:
            lower = max(1.0, guess - span)
            low_error = _tp_at(
                mixture, ideal_gas, z, molar_mass, target, pressure, lower, evaluations
            ).entropy_residual
        if high_error < 0.0:
            upper = guess + span
            high_error = _tp_at(
                mixture, ideal_gas, z, molar_mass, target, pressure, upper, evaluations
            ).entropy_residual
        span *= 1.6
    if not (low_error <= 0.0 and high_error >= 0.0):
        raise EntropyStateFailure(f"could not bracket the specific entropy at {pressure:.6f} Pa")
    last_error = math.nan
    for _ in range(BISECTION_STEPS):
        temperature = 0.5 * (lower + upper)
        trial = _tp_at(
            mixture,
            ideal_gas,
            z,
            molar_mass,
            target,
            pressure,
            temperature,
            evaluations,
        )
        last_error = trial.entropy_residual
        if policy.nonfinite_is_error and not math.isfinite(last_error):
            raise EntropyStateFailure("nonfinite entropy root")
        if abs(last_error) <= ENTROPY_TOLERANCE:
            return trial
        if last_error > 0.0:
            upper = temperature
        else:
            lower = temperature
        if policy.bracket_width_tolerance is not None and upper - lower <= (
            policy.bracket_width_tolerance * max(1.0, temperature)
        ):
            break
    raise EntropyStateFailure(
        f"the entropy root did not close at {pressure:.6f} Pa; residual = {last_error} J/(kg K)"
    )


def _tp_at(
    mixture: Any,
    ideal_gas: Any,
    z: list[float],
    molar_mass: float,
    target: float,
    pressure: float,
    temperature: float,
    evaluations: Evaluations,
) -> EntropyState:
    """One TP flash and its assembled state."""
    evaluations.count += 1
    flash = pt_flash(mixture, quantity(temperature, "K"), quantity(pressure, "Pa"), z)
    return _assemble(mixture, ideal_gas, z, molar_mass, target, pressure, temperature, flash)


def _assemble(
    mixture: Any,
    ideal_gas: Any,
    z: list[float],
    molar_mass: float,
    target: float,
    pressure: float,
    temperature: float,
    flash: Any,
) -> EntropyState:
    """Turn a flash into the assembled state, or refuse it."""
    if flash.phase == "two_phase":
        beta = flash.vapour_fraction
        if beta is None:
            raise EntropyStateFailure("the flash split the feed and reported no vapour fraction")
        if not 0.0 <= beta <= 1.0:
            raise EntropyStateFailure(f"the vapour fraction {beta} is outside [0, 1]")
        liquid = molar_enthalpy_entropy(
            mixture,
            ideal_gas,
            quantity(temperature, "K"),
            quantity(pressure, "Pa"),
            list(flash.x),
            flash.z_liquid,
        )
        vapour = molar_enthalpy_entropy(
            mixture,
            ideal_gas,
            quantity(temperature, "K"),
            quantity(pressure, "Pa"),
            list(flash.y),
            flash.z_vapour,
        )
        molar_entropy = (1.0 - beta) * liquid.s.to("J/(mol*K)").magnitude + beta * vapour.s.to(
            "J/(mol*K)"
        ).magnitude
        molar_enthalpy = (1.0 - beta) * liquid.h.to("J/mol").magnitude + beta * vapour.h.to(
            "J/mol"
        ).magnitude
        v_liquid = flash.z_liquid * MOLAR_GAS_CONSTANT * temperature / pressure
        v_vapour = flash.z_vapour * MOLAR_GAS_CONSTANT * temperature / pressure
        volume = (1.0 - beta) * v_liquid + beta * v_vapour
        if not volume > 0.0 or not math.isfinite(volume):
            raise EntropyStateFailure(f"the two-phase volume at {pressure} Pa is {volume} m3/mol")
        balance = 0.0
        fugacity = 0.0
        for index, fraction in enumerate(z):
            recovered = (1.0 - beta) * flash.x[index] + beta * flash.y[index]
            balance = max(balance, abs(recovered - fraction))
            if fraction > 1.0e-12:
                if not flash.x[index] > 0.0 or not flash.y[index] > 0.0:
                    raise EntropyStateFailure("a phase composition is not positive")
                ln_ratio = (
                    math.log(flash.x[index])
                    + flash.ln_phi_liquid[index]
                    - math.log(flash.y[index])
                    - flash.ln_phi_vapour[index]
                )
                fugacity = max(fugacity, abs(ln_ratio))
        phase_residual = max(abs(sum(flash.x) - 1.0), abs(sum(flash.y) - 1.0))
        return EntropyState(
            pressure=pressure,
            temperature=temperature,
            density=molar_mass / volume,
            molar_enthalpy=molar_enthalpy,
            molar_entropy=molar_entropy,
            entropy_residual=molar_entropy / molar_mass - target,
            component_residual=balance,
            fugacity_residual=fugacity,
            phase_composition_residual=phase_residual,
            phase="two_phase",
            vapour_fraction=beta,
        )
    reduced = reduced_parameters(mixture, temperature, pressure)
    kij = mixture.effective_kij(reduced.t_kelvin)
    liquid_side = flash.phase != "all_vapour"
    root = phase_state(reduced, kij, z, liquid=liquid_side).z
    state = molar_enthalpy_entropy(
        mixture,
        ideal_gas,
        quantity(temperature, "K"),
        quantity(pressure, "Pa"),
        z,
        root,
    )
    molar_entropy = state.s.to("J/(mol*K)").magnitude
    molar_enthalpy = state.h.to("J/mol").magnitude
    volume = root * MOLAR_GAS_CONSTANT * temperature / pressure
    if not volume > 0.0 or not math.isfinite(volume):
        raise EntropyStateFailure(f"the single-phase volume at {pressure} Pa is {volume} m3/mol")
    return EntropyState(
        pressure=pressure,
        temperature=temperature,
        density=molar_mass / volume,
        molar_enthalpy=molar_enthalpy,
        molar_entropy=molar_entropy,
        entropy_residual=molar_entropy / molar_mass - target,
        component_residual=0.0,
        fugacity_residual=0.0,
        phase_composition_residual=0.0,
        phase=flash.phase,
        vapour_fraction=None,
    )
