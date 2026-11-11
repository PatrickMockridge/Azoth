"""``thermal.fouling_model`` - fixed, Kern-Seaton and Ebert-Panchal fouling.

Spec: ``specs/calcs/thermal/fouling_model.toml``

The Python twin of ``crates/azoth-thermal/src/fouling_model.rs``, and both mirror
``FoulingModel``: the fixed TEMA resistance, the Kern-Seaton asymptote, and the Ebert-Panchal
threshold model with its analytical resistance, rate, threshold temperature and threshold
velocity.

**The conditions are always applied** through the class's ``updateConditions``, so the
Reynolds number and wall shear stress are derived from the stated velocity, density,
viscosity, wall temperature and bore rather than from the class's cached fields. **The
threshold temperature answers the class's ``Double.MAX_VALUE`` sentinel** when the target
rate exceeds the deposition rate at zero fouling, and **the threshold velocity is the
class's 50-step bisection over ``[0.1, 10]`` m/s**.
"""

from __future__ import annotations

import math

from azoth._registry_gen import spec as _spec_for
from azoth.core.errors import InvalidInputError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import FoulingModelResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning
from azoth.thermal.reference.shell_and_tube_rating import darcy_friction

CALC_ID = "thermal.fouling_model"
R_GAS = 8.314
SECONDS_PER_HOUR = 3600.0


def conditions(
    velocity: float, density: float, viscosity: float, diameter: float
) -> tuple[float, float]:
    """The Reynolds number and wall shear stress, ``Re = rho v D/mu`` and ``f/8 rho v**2``."""
    if viscosity <= 0.0:
        return 0.0, 0.0
    reynolds = density * velocity * diameter / viscosity
    friction = darcy_friction(reynolds)
    shear = friction / 8.0 * density * velocity * velocity
    return reynolds, shear


def deposition(
    alpha: float, beta: float, activation_energy: float, reynolds: float, wall_k: float
) -> float:
    """``alpha Re**beta exp(-E/(R Tw))``."""
    if reynolds <= 0.0 or wall_k <= 0.0:
        return 0.0
    return alpha * math.pow(reynolds, beta) * math.exp(-activation_energy / (R_GAS * wall_k))


def ebert_panchal_fouling_rate(
    alpha: float,
    beta: float,
    activation_energy: float,
    gamma: float,
    wall_shear_stress: float,
    reynolds: float,
    wall_k: float,
    current_rf: float,
) -> float:
    """The instantaneous rate at a resistance, per hour."""
    if wall_k <= 0.0 or reynolds <= 0.0:
        return 0.0
    dep = deposition(alpha, beta, activation_energy, reynolds, wall_k)
    removal = gamma * wall_shear_stress * current_rf
    return (dep - removal) * SECONDS_PER_HOUR


def ebert_panchal_resistance(
    alpha: float,
    beta: float,
    activation_energy: float,
    gamma: float,
    wall_shear_stress: float,
    reynolds: float,
    wall_k: float,
    time_hours: float,
) -> float:
    """The analytical resistance at a time."""
    if time_hours <= 0.0:
        return 0.0
    dep = deposition(alpha, beta, activation_energy, reynolds, wall_k)
    removal = gamma * wall_shear_stress
    if removal <= 0.0:
        return dep * time_hours * SECONDS_PER_HOUR
    asymptotic = dep / removal
    tau_hours = 1.0 / removal / SECONDS_PER_HOUR
    return asymptotic * (1.0 - math.exp(-time_hours / tau_hours))


def kern_seaton_resistance(rf_max: float, time_constant: float, time_hours: float) -> float:
    """``Rf_max (1 - exp(-t/tau))``."""
    if time_hours <= 0.0 or rf_max <= 0.0 or time_constant <= 0.0:
        return 0.0
    return rf_max * (1.0 - math.exp(-time_hours / time_constant))


def threshold_temperature(
    alpha: float, beta: float, activation_energy: float, reynolds: float, threshold_rate: float
) -> float:
    """The wall temperature whose zero-resistance rate equals the target."""
    if alpha <= 0.0 or threshold_rate <= 0.0 or reynolds <= 0.0:
        return 0.0
    ln_term = math.log(threshold_rate / (alpha * math.pow(reynolds, beta)))
    if ln_term >= 0.0:
        return float.fromhex("0x1.fffffffffffffp+1023")
    return -activation_energy / (R_GAS * ln_term)


def threshold_velocity(
    alpha: float,
    beta: float,
    activation_energy: float,
    gamma: float,
    density: float,
    viscosity: float,
    wall_k: float,
    diameter: float,
    fixed_resistance: float,
) -> float:
    """The class's 50-step bisection for the velocity at which the asymptote crosses."""
    if wall_k <= 0.0 or density <= 0.0 or viscosity <= 0.0 or diameter <= 0.0:
        return 0.0
    low = 0.1
    high = 10.0
    for _ in range(50):
        mid = (low + high) / 2.0
        reynolds = density * mid * diameter / viscosity
        friction = darcy_friction(reynolds)
        shear = friction / 8.0 * density * mid * mid
        dep = deposition(alpha, beta, activation_energy, reynolds, wall_k)
        removal = gamma * shear
        asymptotic = dep / removal if removal > 0.0 else float.fromhex("0x1.fffffffffffffp+1023")
        if asymptotic > fixed_resistance:
            low = mid
        else:
            high = mid
    return (low + high) / 2.0


def fouling_model(
    model: str,
    velocity: Q,
    fluid_density: Q,
    fluid_viscosity: Q,
    wall_temperature: Q,
    tube_inner_diameter: Q,
    alpha: float,
    beta: float,
    activation_energy: Q,
    gamma: float,
    fixed_fouling_resistance: float,
    rf_max: float,
    time_constant: Q,
    time_hours: Q,
    target_fouling_resistance: float,
    threshold_rate: float,
    euler_steps: float | None = None,
    euler_step_hours: Q | None = None,
) -> FoulingModelResult:
    """The composed fouling calculation.

    Raises:
        InvalidInputError: for a non-integral or negative Euler step count.
        OutOfRangeError: for an input outside the spec's declared bounds.
    """
    spec = _spec_for(CALC_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []
    values = {
        "velocity": input_to_si(spec, "velocity", velocity),
        "fluid_density": input_to_si(spec, "fluid_density", fluid_density),
        "fluid_viscosity": input_to_si(spec, "fluid_viscosity", fluid_viscosity),
        "wall_temperature": input_to_si(spec, "wall_temperature", wall_temperature),
        "tube_inner_diameter": input_to_si(spec, "tube_inner_diameter", tube_inner_diameter),
        "alpha": alpha,
        "activation_energy": input_to_si(spec, "activation_energy", activation_energy),
        "gamma": gamma,
        "fixed_fouling_resistance": fixed_fouling_resistance,
        "rf_max": rf_max,
        "time_constant": input_to_si(spec, "time_constant", time_constant),
        "time_hours": input_to_si(spec, "time_hours", time_hours),
        "target_fouling_resistance": target_fouling_resistance,
        "threshold_rate": threshold_rate,
    }
    apply_checks(checks.on_input, values.get, warnings)

    if euler_steps is None:
        steps = 0
    else:
        if (
            not math.isfinite(euler_steps)
            or euler_steps != math.floor(euler_steps)
            or euler_steps < 0
        ):
            raise InvalidInputError("euler_steps", f"{euler_steps} is not a whole step count")
        steps = int(euler_steps)
    # The boundary carries SI seconds; the class's arithmetic is in hours.
    time_constant_h = values["time_constant"] / SECONDS_PER_HOUR
    time_hours_h = values["time_hours"] / SECONDS_PER_HOUR
    step_hours = (
        0.0
        if euler_step_hours is None
        else input_to_si(spec, "euler_step_hours", euler_step_hours) / SECONDS_PER_HOUR
    )

    reynolds, shear = conditions(
        values["velocity"],
        values["fluid_density"],
        values["fluid_viscosity"],
        values["tube_inner_diameter"],
    )
    alpha_value = values["alpha"]
    energy = values["activation_energy"]
    gamma_value = values["gamma"]
    wall_k = values["wall_temperature"]
    if model == "fixed":
        asymptotic = values["fixed_fouling_resistance"]
        resistance_at_time = values["fixed_fouling_resistance"]
    elif model == "ebert_panchal":
        removal = gamma_value * shear
        asymptotic = (
            deposition(alpha_value, beta, energy, reynolds, wall_k) / removal
            if removal > 0.0
            else float.fromhex("0x1.fffffffffffffp+1023")
        )
        resistance_at_time = ebert_panchal_resistance(
            alpha_value, beta, energy, gamma_value, shear, reynolds, wall_k, time_hours_h
        )
    elif model == "kern_seaton":
        asymptotic = values["rf_max"]
        resistance_at_time = kern_seaton_resistance(values["rf_max"], time_constant_h, time_hours_h)
    else:
        raise InvalidInputError(
            "model", f"{model} is not one of fixed, ebert_panchal or kern_seaton"
        )

    current = 0.0
    operating = 0.0
    for _ in range(steps):
        if model == "ebert_panchal":
            rate = ebert_panchal_fouling_rate(
                alpha_value, beta, energy, gamma_value, shear, reynolds, wall_k, current
            )
            current = max(0.0, current + rate * step_hours)
        operating += step_hours
    if model == "fixed":
        current_resistance = values["fixed_fouling_resistance"]
    elif model == "ebert_panchal":
        current_resistance = current
    else:
        current_resistance = kern_seaton_resistance(values["rf_max"], time_constant_h, operating)

    target = values["target_fouling_resistance"]
    if target >= asymptotic:
        time_to_target = -1.0
    elif model == "fixed":
        time_to_target = 0.0
    elif model == "ebert_panchal":
        removal = gamma_value * shear
        tau_hours = (
            1.0 / (removal * SECONDS_PER_HOUR)
            if removal > 0.0
            else float.fromhex("0x1.fffffffffffffp+1023")
        )
        time_to_target = -tau_hours * math.log(1.0 - target / asymptotic)
    else:
        time_to_target = -time_constant_h * math.log(1.0 - target / values["rf_max"])

    return FoulingModelResult(
        fouling_resistance_at_time=resistance_at_time,
        current_fouling_resistance=current_resistance,
        asymptotic_fouling_resistance=asymptotic,
        threshold_wall_temperature=from_si(
            threshold_temperature(alpha_value, beta, energy, reynolds, values["threshold_rate"]),
            "K",
        ),
        threshold_velocity=from_si(
            threshold_velocity(
                alpha_value,
                beta,
                energy,
                gamma_value,
                values["fluid_density"],
                values["fluid_viscosity"],
                wall_k,
                values["tube_inner_diameter"],
                values["fixed_fouling_resistance"],
            ),
            "m/s",
        ),
        time_to_target=from_si(time_to_target * SECONDS_PER_HOUR, "h"),
        operating_time=from_si(operating * SECONDS_PER_HOUR, "h"),
        fouling_rate=ebert_panchal_fouling_rate(
            alpha_value, beta, energy, gamma_value, shear, reynolds, wall_k, current_resistance
        ),
        reynolds_number=reynolds,
        wall_shear_stress=from_si(shear, "Pa"),
        warnings=tuple(warnings),
    )
