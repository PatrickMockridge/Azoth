"""``eos.equilibrium_sound_speed`` - the homogeneous-equilibrium acoustic derivative.

Spec: ``specs/models/eos/equilibrium_sound_speed.toml``.

The Python twin of ``crates/azoth-eos/src/equilibrium_sound_speed.rs``, written to mirror it.
``c^2 = (dp/drho)`` at fixed specific entropy and total composition: each pressure sample is
re-equilibrated at fixed ``s`` - mixtures by a bracketed TP root, pure fluids by the
saturation-quality solve - and differentiated on a central or one-sided stencil that stays on
the centre state's phase branch.

This is not :func:`azoth.eos.sound_speed`: that is the frozen phase speed off one phase,
while this one lets the phase split relax with pressure.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any

from azoth.core.errors import InvalidInputError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import EquilibriumSoundSpeedResult
from azoth.core.units import Q, from_si, input_to_si, quantity
from azoth.core.warnings import Warning
from azoth.eos import components as _components
from azoth.eos.reference._mixture_state import phase_state, reduced_parameters
from azoth.eos.reference.molar_enthalpy_entropy import molar_enthalpy_entropy
from azoth.eos.reference.pr_molar_volume import MOLAR_GAS_CONSTANT
from azoth.eos.reference.pt_flash import check_composition, pt_flash
from azoth.eos.reference.pure_saturation import pure_saturation

MODEL_ID = "eos.equilibrium_sound_speed"

#: ``EquilibriumSoundSpeed.ENTROPY_TOLERANCE``, J/(kg K).
ENTROPY_TOLERANCE = 1.0e-7
#: ``BALANCE_TOLERANCE``.
BALANCE_TOLERANCE = 1.0e-8
#: ``FUGACITY_TOLERANCE``.
FUGACITY_TOLERANCE = 1.0e-6
#: The successive-derivative band that accepts a step.
DERIVATIVE_TOLERANCE = 2.0e-3
#: The step is halved this many times after the first estimate.
MAX_REFINEMENTS = 10
#: The temperature root's relative tolerance.
ROOT_TOLERANCE = 1.0e-12
#: Entropy-bracket widenings before the bracket is refused.
BRACKET_STEPS = 30
#: The bisection cap.
BISECTION_STEPS = 100
#: The class's relative-step band.
MIN_RELATIVE_STEP = 1.0e-6
MAX_RELATIVE_STEP = 0.05


@dataclass(frozen=True, slots=True)
class _Sample:
    """One accepted or trial pressure sample."""

    pressure: float
    temperature: float
    density: float
    entropy_residual: float
    component_residual: float
    fugacity_residual: float
    phase: str
    vapour_fraction: float | None


@dataclass(frozen=True, slots=True)
class _Estimate:
    stencil: str
    states: list[_Sample]
    derivative: float


class _FlashFailure(Exception):
    """The class's ``FLASH_FAILED`` inside the calculation."""


def equilibrium_sound_speed(
    components: list[str],
    T: Q,
    P: Q,
    z: list[float],
    relative_pressure_step: float | None = None,
    specific_entropy: Q | None = None,
) -> EquilibriumSoundSpeedResult:
    """The equilibrium acoustic speed of a mixture at a state.

    Args:
        components: the fluid's substances, by name.
        T: the state's absolute temperature.
        P: the state's absolute pressure.
        z: the mixture's mole fractions, checked rather than renormalised.
        relative_pressure_step: the first step over the centre pressure, in ``[1e-6, 0.05]``.
        specific_entropy: the fixed specific entropy; absent takes the flash's own value.

    Returns:
        The speed (absent on failure), the status and message, the stencil and its samples.

    Raises:
        InvalidInputError: for a step outside the band, a non-positive state, or a
            composition of the wrong length.

    See :func:`azoth.eos.equilibrium_sound_speed`.
    """
    spec = _spec()
    checks = checks_for(spec)
    warnings: list[Warning] = []

    t = input_to_si(spec, "T", T)
    p = input_to_si(spec, "P", P)
    apply_checks(
        checks.on_input,
        {
            "T": t,
            "P": p,
            "relative_pressure_step": relative_pressure_step,
        }.get,
        warnings,
    )
    step = 1.0e-3 if relative_pressure_step is None else float(relative_pressure_step)
    if not (MIN_RELATIVE_STEP <= step <= MAX_RELATIVE_STEP):
        raise InvalidInputError(
            "relative_pressure_step",
            f"the step must be in [1e-6, 0.05], and {step} is not",
        )
    target = (
        None
        if specific_entropy is None
        else input_to_si(spec, "specific_entropy", specific_entropy)
    )

    mixture, ideal_gas = _components.mixture_of(components, eos="pr")
    if len(z) != len(mixture.components):
        raise InvalidInputError(
            "z",
            f"a mixture of {len(mixture.components)} components needs that many mole "
            f"fractions, and {z!r} has {len(z)}",
        )
    check_composition(list(z), len(mixture.components), "z")
    mass = _molar_mass(mixture, list(z))
    if target is None:
        target = _specific_entropy_at(mixture, ideal_gas, t, p, list(z), mass)

    calculation = _Calculation(mixture, ideal_gas, list(z), mass, float(target), t)
    outcome = calculation.run(p, step)

    from azoth.core.result import EquilibriumSoundSpeedSample

    speed = outcome["sound_speed"]
    return EquilibriumSoundSpeedResult(
        sound_speed=None if speed is None else from_si(speed, "m/s"),
        status=outcome["status"],
        message=outcome["message"],
        pressure_step=from_si(outcome["pressure_step"], "Pa"),
        relative_step_error=outcome["relative_step_error"],
        specific_entropy=from_si(float(target), "J/(kg*K)"),
        stencil=outcome["stencil"],
        phase_boundary_encountered=outcome["phase_boundary"],
        flash_evaluations=outcome["evaluations"],
        samples=tuple(
            EquilibriumSoundSpeedSample(
                pressure=from_si(sample.pressure, "Pa"),
                temperature=from_si(sample.temperature, "K"),
                density=from_si(sample.density, "kg/m**3"),
                entropy_residual=from_si(sample.entropy_residual, "J/(kg*K)"),
                phase=sample.phase,
                vapour_fraction=sample.vapour_fraction,
            )
            for sample in outcome["samples"]
        ),
        warnings=tuple(warnings),
    )


class _Calculation:
    """One calculation's state, mirroring the class's inner ``Calculation``."""

    def __init__(
        self,
        mixture: Any,
        ideal_gas: Any,
        z: list[float],
        molar_mass: float,
        target: float,
        start_temperature: float,
    ) -> None:
        self.mixture = mixture
        self.ideal_gas = ideal_gas
        self.z = z
        self.molar_mass = molar_mass
        self.target = target
        self.start_temperature = start_temperature
        self.evaluations = 0
        self.boundary = False

    def run(self, pressure: float, relative_step: float) -> dict[str, Any]:
        """The class's ``run``: centre, stencils, refinement, status."""
        samples: list[_Sample] = []
        stencil = "unavailable"
        step = math.nan
        error: float | None = None
        try:
            centre = self.solve(pressure)
        except _FlashFailure as failure:
            return self.outcome("flash_failed", str(failure), None, step, error, stencil, samples)
        samples = [centre]
        previous: tuple[str, float] | None = None
        step = pressure * relative_step
        for refinement in range(MAX_REFINEMENTS + 1):
            if refinement > 0:
                step *= 0.5
            try:
                current = self.estimate(centre, step)
            except _FlashFailure as failure:
                return self.outcome(
                    "flash_failed", str(failure), None, step, error, stencil, samples
                )
            if current is None:
                previous = None
                continue
            samples = current.states
            stencil = current.stencil
            if previous is not None and previous[0] == current.stencil:
                last_derivative = previous[1]
                if (
                    last_derivative > 0.0
                    and current.derivative > 0.0
                    and math.isfinite(last_derivative)
                    and math.isfinite(current.derivative)
                ):
                    error = abs(current.derivative - last_derivative) / abs(current.derivative)
                    speed = math.sqrt(1.0 / current.derivative)
                    if error <= DERIVATIVE_TOLERANCE and math.isfinite(speed) and speed > 0.0:
                        return self.outcome(
                            "converged",
                            "entropy, inventory, fugacity and step convergence checks passed",
                            speed,
                            step,
                            error,
                            stencil,
                            samples,
                        )
            previous = (current.stencil, current.derivative)
        return self.outcome(
            "derivative_not_converged",
            "no positive step-converged equilibrium density derivative",
            None,
            step * 2.0,
            error,
            stencil,
            samples,
        )

    def outcome(
        self,
        status: str,
        message: str,
        sound_speed: float | None,
        pressure_step: float,
        relative_step_error: float | None,
        stencil: str,
        samples: list[_Sample],
    ) -> dict[str, Any]:
        return {
            "status": status,
            "message": message,
            "sound_speed": sound_speed,
            "pressure_step": pressure_step,
            "relative_step_error": relative_step_error,
            "stencil": stencil,
            "phase_boundary": self.boundary,
            "evaluations": self.evaluations,
            "samples": samples,
        }

    def estimate(self, centre: _Sample, step: float) -> _Estimate | None:
        pressure = centre.pressure
        low = self.solve(pressure - step)
        high = self.solve(pressure + step)
        same_low = centre.phase == low.phase
        same_high = centre.phase == high.phase
        self.boundary |= not same_low or not same_high
        if same_low and same_high:
            return _Estimate(
                "central",
                [centre, low, high],
                (high.density - low.density) / (2.0 * step),
            )
        if same_high:
            far = self.solve(pressure + 2.0 * step)
            self.boundary |= centre.phase != far.phase
            if centre.phase != far.phase:
                return None
            return _Estimate(
                "forward",
                [centre, high, far],
                (-3.0 * centre.density + 4.0 * high.density - far.density) / (2.0 * step),
            )
        if same_low:
            far = self.solve(pressure - 2.0 * step)
            self.boundary |= centre.phase != far.phase
            if centre.phase != far.phase:
                return None
            return _Estimate(
                "backward",
                [centre, low, far],
                (3.0 * centre.density - 4.0 * low.density + far.density) / (2.0 * step),
            )
        return None

    def solve(self, pressure: float) -> _Sample:
        if not pressure > 0.0 or not math.isfinite(pressure):
            raise _FlashFailure(f"the pressure sample {pressure} Pa is not positive")
        if len(self.mixture.components) == 1:
            state = self.saturation_state(pressure)
            if state is not None:
                return state
        return self.tp(pressure)

    def saturation_state(self, pressure: float) -> _Sample | None:
        component = self.mixture.components[0]
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
            reduced = reduced_parameters(self.mixture, temperature, pressure)
        except Exception:
            return None
        kij = self.mixture.effective_kij(reduced.t_kelvin)
        try:
            liquid_root = phase_state(reduced, kij, self.z, liquid=True).z
            vapour_root = phase_state(reduced, kij, self.z, liquid=False).z
        except Exception:
            return None
        if abs(vapour_root - liquid_root) < 1.0e-12:
            return None
        try:
            liquid = molar_enthalpy_entropy(
                self.mixture,
                self.ideal_gas,
                quantity(temperature, "K"),
                quantity(pressure, "Pa"),
                self.z,
                liquid_root,
            )
            vapour = molar_enthalpy_entropy(
                self.mixture,
                self.ideal_gas,
                quantity(temperature, "K"),
                quantity(pressure, "Pa"),
                self.z,
                vapour_root,
            )
        except Exception:
            return None
        s_liquid = liquid.s.to("J/(mol*K)").magnitude / self.molar_mass
        s_vapour = vapour.s.to("J/(mol*K)").magnitude / self.molar_mass
        if not s_vapour > s_liquid or not (s_liquid <= self.target <= s_vapour):
            return None
        self.evaluations += 1
        quality = (self.target - s_liquid) / (s_vapour - s_liquid)
        v_liquid = liquid_root * MOLAR_GAS_CONSTANT * temperature / pressure
        v_vapour = vapour_root * MOLAR_GAS_CONSTANT * temperature / pressure
        volume = (1.0 - quality) * v_liquid + quality * v_vapour
        if not volume > 0.0 or not math.isfinite(volume):
            raise _FlashFailure(
                f"the saturation-quality volume at {pressure} Pa is {volume} m3/mol"
            )
        return self.checked(
            _Sample(
                pressure=pressure,
                temperature=temperature,
                density=self.molar_mass / volume,
                entropy_residual=(1.0 - quality) * s_liquid + quality * s_vapour - self.target,
                component_residual=0.0,
                fugacity_residual=0.0,
                phase="two_phase",
                vapour_fraction=quality,
            )
        )

    def tp(self, pressure: float) -> _Sample:
        guess = self.start_temperature
        initial = self.tp_at(pressure, guess)
        if abs(initial.entropy_residual) <= ENTROPY_TOLERANCE:
            return self.checked(initial)
        lower = upper = guess
        low_error = high_error = initial.entropy_residual
        span = max(1.0, guess * 0.005)
        for _ in range(BRACKET_STEPS):
            if not (low_error > 0.0 or high_error < 0.0):
                break
            if low_error > 0.0:
                lower = max(1.0, guess - span)
                low_error = self.tp_at(pressure, lower).entropy_residual
            if high_error < 0.0:
                upper = guess + span
                high_error = self.tp_at(pressure, upper).entropy_residual
            span *= 1.6
        if not (low_error <= 0.0 and high_error >= 0.0):
            raise _FlashFailure(f"could not bracket the specific entropy at {pressure:.6f} Pa")
        last_error = math.nan
        for _ in range(BISECTION_STEPS):
            temperature = 0.5 * (lower + upper)
            trial = self.tp_at(pressure, temperature)
            last_error = trial.entropy_residual
            if abs(last_error) <= ENTROPY_TOLERANCE:
                return self.checked(trial)
            if last_error > 0.0:
                upper = temperature
            else:
                lower = temperature
            if upper - lower <= ROOT_TOLERANCE * max(1.0, temperature):
                break
        raise _FlashFailure(
            f"the entropy root did not close at {pressure:.6f} Pa; residual = {last_error} J/(kg K)"
        )

    def tp_at(self, pressure: float, temperature: float) -> _Sample:
        self.evaluations += 1
        flash = pt_flash(self.mixture, quantity(temperature, "K"), quantity(pressure, "Pa"), self.z)
        return self.assemble(pressure, temperature, flash)

    def assemble(self, pressure: float, temperature: float, flash: Any) -> _Sample:
        if flash.phase == "two_phase":
            beta = flash.vapour_fraction
            if beta is None:
                raise _FlashFailure("the flash split the feed and reported no vapour fraction")
            if not 0.0 <= beta <= 1.0:
                raise _FlashFailure(f"the vapour fraction {beta} is outside [0, 1]")
            liquid = molar_enthalpy_entropy(
                self.mixture,
                self.ideal_gas,
                quantity(temperature, "K"),
                quantity(pressure, "Pa"),
                list(flash.x),
                flash.z_liquid,
            )
            vapour = molar_enthalpy_entropy(
                self.mixture,
                self.ideal_gas,
                quantity(temperature, "K"),
                quantity(pressure, "Pa"),
                list(flash.y),
                flash.z_vapour,
            )
            entropy = (
                (1.0 - beta) * liquid.s.to("J/(mol*K)").magnitude
                + beta * vapour.s.to("J/(mol*K)").magnitude
            ) / self.molar_mass
            v_liquid = flash.z_liquid * MOLAR_GAS_CONSTANT * temperature / pressure
            v_vapour = flash.z_vapour * MOLAR_GAS_CONSTANT * temperature / pressure
            volume = (1.0 - beta) * v_liquid + beta * v_vapour
            if not volume > 0.0 or not math.isfinite(volume):
                raise _FlashFailure(f"the two-phase volume at {pressure} Pa is {volume} m3/mol")
            balance = 0.0
            fugacity = 0.0
            for index, fraction in enumerate(self.z):
                recovered = (1.0 - beta) * flash.x[index] + beta * flash.y[index]
                balance = max(balance, abs(recovered - fraction))
                if fraction > 1.0e-12:
                    if not flash.x[index] > 0.0 or not flash.y[index] > 0.0:
                        raise _FlashFailure("a phase composition is not positive")
                    ln_ratio = (
                        math.log(flash.x[index])
                        + flash.ln_phi_liquid[index]
                        - math.log(flash.y[index])
                        - flash.ln_phi_vapour[index]
                    )
                    fugacity = max(fugacity, abs(ln_ratio))
            density = self.molar_mass / volume
            phase = "two_phase"
            vapour_fraction = beta
        else:
            reduced = reduced_parameters(self.mixture, temperature, pressure)
            kij = self.mixture.effective_kij(reduced.t_kelvin)
            liquid_side = flash.phase != "all_vapour"
            root = phase_state(reduced, kij, self.z, liquid=liquid_side).z
            state = molar_enthalpy_entropy(
                self.mixture,
                self.ideal_gas,
                quantity(temperature, "K"),
                quantity(pressure, "Pa"),
                self.z,
                root,
            )
            entropy = state.s.to("J/(mol*K)").magnitude / self.molar_mass
            volume = root * MOLAR_GAS_CONSTANT * temperature / pressure
            if not volume > 0.0 or not math.isfinite(volume):
                raise _FlashFailure(f"the single-phase volume at {pressure} Pa is {volume} m3/mol")
            density = self.molar_mass / volume
            balance = 0.0
            fugacity = 0.0
            phase = flash.phase
            vapour_fraction = None
        return _Sample(
            pressure=pressure,
            temperature=temperature,
            density=density,
            entropy_residual=entropy - self.target,
            component_residual=balance,
            fugacity_residual=fugacity,
            phase=phase,
            vapour_fraction=vapour_fraction,
        )

    def checked(self, sample: _Sample) -> _Sample:
        closed = (
            abs(sample.entropy_residual) <= ENTROPY_TOLERANCE
            and sample.component_residual <= BALANCE_TOLERANCE
            and math.isfinite(sample.component_residual)
            and sample.fugacity_residual <= FUGACITY_TOLERANCE
            and math.isfinite(sample.fugacity_residual)
            and sample.density > 0.0
            and math.isfinite(sample.density)
            and sample.temperature > 0.0
            and math.isfinite(sample.temperature)
            and sample.pressure > 0.0
            and math.isfinite(sample.pressure)
        )
        if closed:
            return sample
        raise _FlashFailure(
            f"the equilibrium closure failed at {sample.pressure:.6f} Pa: entropy "
            f"residual = {sample.entropy_residual:.3e} J/(kg K), component residual = "
            f"{sample.component_residual:.3e}, log fugacity residual = "
            f"{sample.fugacity_residual:.3e}"
        )


def _molar_mass(mixture: Any, z: list[float]) -> float:
    """The mixture's molar mass, kg/mol."""
    total = 0.0
    for component, fraction in zip(mixture.components, z, strict=True):
        if component.molar_mass is None:
            raise InvalidInputError(
                "molar_mass",
                "a component carries no molar mass, and a specific entropy divides by it",
            )
        total += fraction * component.molar_mass.to("kg/mol").magnitude
    if not total > 0.0:
        raise InvalidInputError("molar_mass", f"the mixture's molar mass is {total} kg/mol")
    return total


def _specific_entropy_at(
    mixture: Any,
    ideal_gas: Any,
    t: float,
    p: float,
    z: list[float],
    molar_mass: float,
) -> float:
    """The flashed state's specific entropy, J/(kg K)."""
    flash = pt_flash(mixture, quantity(t, "K"), quantity(p, "Pa"), z)
    if flash.phase == "two_phase":
        beta = flash.vapour_fraction
        if beta is None:
            raise InvalidInputError(
                "state", "the flash split the feed and reported no vapour fraction"
            )
        liquid = molar_enthalpy_entropy(
            mixture,
            ideal_gas,
            quantity(t, "K"),
            quantity(p, "Pa"),
            list(flash.x),
            flash.z_liquid,
        )
        vapour = molar_enthalpy_entropy(
            mixture,
            ideal_gas,
            quantity(t, "K"),
            quantity(p, "Pa"),
            list(flash.y),
            flash.z_vapour,
        )
        molar = (1.0 - beta) * liquid.s.to("J/(mol*K)").magnitude + beta * vapour.s.to(
            "J/(mol*K)"
        ).magnitude
    else:
        reduced = reduced_parameters(mixture, t, p)
        kij = mixture.effective_kij(reduced.t_kelvin)
        liquid_side = flash.phase != "all_vapour"
        root = phase_state(reduced, kij, z, liquid=liquid_side).z
        molar = (
            molar_enthalpy_entropy(mixture, ideal_gas, quantity(t, "K"), quantity(p, "Pa"), z, root)
            .s.to("J/(mol*K)")
            .magnitude
        )
    return molar / molar_mass


def _spec() -> dict[str, object]:
    """The generated spec, imported here so the module has no import-time cycle."""
    from azoth._models_gen import model

    return model(MODEL_ID)
