"""``process.leak_blowdown`` - NeqSim's legacy vessel-blowdown time loop.

Spec: ``specs/models/process/leak_blowdown.toml``

The Python twin of ``crates/azoth-process/src/models/leak_blowdown.rs``, written to mirror it.
The port source is ``LeakModel.calculateSourceTerm(duration, timeStep)``: the inventory from
the corrected density and the vessel volume, the record-then-update loop with its simplified
orifice, capped jet velocity, unclamped momentum and droplet SMD, and the isentropic pressure
and temperature update. The modern ``calculateReleaseFlow`` path is
:mod:`azoth.process.reference.hem_release`.

The class's silent fallbacks - an ideal-gas density, ``gamma = 1.3``, methane's molar mass,
``rho = 1``, a NaN rate as zero - are refused here, and the class's cached density is not
reproduced: every step flashes and uses the fresh corrected density.
"""

from __future__ import annotations

import math
from typing import Any, NamedTuple

from azoth.core.errors import InvalidInputError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import LeakBlowdownDiagnostic, LeakBlowdownPoint, LeakBlowdownResult
from azoth.core.units import Q, from_si, input_to_si, quantity
from azoth.core.warnings import Warning
from azoth.eos import components as _components
from azoth.eos.reference._entropy_state import molar_mass_of
from azoth.eos.reference._heat_capacity import phase_capacity
from azoth.eos.reference._mixture_state import phase_state, reduced_parameters
from azoth.eos.reference.pr_molar_volume import MOLAR_GAS_CONSTANT
from azoth.eos.reference.pt_flash import check_composition, pt_flash

MODEL_ID = "process.leak_blowdown"

#: The class's fixed gas constant, J/(kmol K).
R = 8314.0
#: The class's surrogate surface tension for the droplet correlation, N/m.
SURFACE_TENSION = 0.02
#: The class's Nukiyama-Tanasawa coefficient.
SMD_COEFFICIENT = 585.0
#: The class's liquid-fraction exponent.
SMD_LIQUID_EXPONENT = 0.45
#: The class's liquid-fraction gate.
SMD_LIQUID_FLOOR = 0.001
#: The inventory floor, relative to the initial mass.
MASS_FLOOR = 0.01
#: The declared-point cap.
MAX_POINTS = 100_000
#: The velocity cap, as a multiple of the ideal-gas sound speed.
VELOCITY_CAP = 1.2


class _Failure(Exception):
    """A state failure the class would throw out of the loop; the port reports it."""

    def __init__(self, code: str, message: str) -> None:
        super().__init__(message)
        self.code = code


class _Point(NamedTuple):
    time: float
    mass_flow_rate: float
    temperature: float
    pressure: float
    vapour_mole_fraction: float
    jet_velocity: float
    jet_momentum: float
    droplet_smd: float


class _Phase(NamedTuple):
    molar_mass: float
    molar_volume: float
    density: float
    cp: float
    cv: float


class _State(NamedTuple):
    flash: Any
    pressure: float
    temperature: float
    density: float
    gamma: float
    vapour_fraction: float | None
    vapour_density: float | None
    liquid_density: float | None

    def vapour_mole_fraction(self) -> float:
        if self.flash.phase == "two_phase":
            return 0.0 if self.vapour_fraction is None else self.vapour_fraction
        return 1.0 if self.flash.phase == "all_vapour" else 0.0


class _Outcome(NamedTuple):
    status: str
    scenario_name: str
    hole_diameter: float
    number_of_points: int
    points: tuple[_Point, ...]
    total_mass_released: float | None
    peak_mass_flow_rate: float | None
    time_to_empty: float | None
    diagnostics: tuple[tuple[str, str], ...]
    initial_mass: float
    initial_gamma: float
    gammas: tuple[float, ...]


def leak_blowdown(
    components: list[str],
    initial_z: list[float],
    initial_t: Q,
    initial_p: Q,
    duration: Q,
    vessel_volume: Q | None = None,
    hole_diameter: Q | None = None,
    discharge_coefficient: float | None = None,
    back_pressure: Q | None = None,
    time_step: Q | None = None,
    scenario_name: str | None = None,
) -> LeakBlowdownResult:
    """Blow a vessel down through a circular opening."""
    spec = _spec()
    checks = checks_for(spec)
    warnings: list[Warning] = []

    t = input_to_si(spec, "initial_t", initial_t)
    p = input_to_si(spec, "initial_p", initial_p)
    span = input_to_si(spec, "duration", duration)
    volume = 1.0 if vessel_volume is None else input_to_si(spec, "vessel_volume", vessel_volume)
    diameter = 0.01 if hole_diameter is None else input_to_si(spec, "hole_diameter", hole_diameter)
    coefficient = 0.62 if discharge_coefficient is None else float(discharge_coefficient)
    back = 101325.0 if back_pressure is None else input_to_si(spec, "back_pressure", back_pressure)
    step = 1.0 if time_step is None else input_to_si(spec, "time_step", time_step)
    scenario = "Release" if scenario_name is None else scenario_name
    apply_checks(
        checks.on_input,
        {
            "initial_t": t,
            "initial_p": p,
            "duration": span,
            "vessel_volume": volume,
            "hole_diameter": diameter,
            "discharge_coefficient": coefficient,
            "back_pressure": back,
            "time_step": step,
        }.get,
        warnings,
    )
    if not 0.0 < coefficient <= 1.0:
        raise InvalidInputError(
            "discharge_coefficient",
            f"the effective-area factor must be in (0, 1], and {coefficient} is not",
        )
    if not span > 0.0 or not step > 0.0:
        raise InvalidInputError(
            "duration", f"a positive duration and step are required, and they are {span} and {step}"
        )
    declared = math.ceil(span / step)
    if declared < 0 or declared > MAX_POINTS - 1:
        raise InvalidInputError(
            "duration",
            f"the declared point count is {declared}+1, and the port's cap is {MAX_POINTS}",
        )

    try:
        outcome = _states(
            components,
            list(initial_z),
            t,
            p,
            span,
            volume,
            diameter,
            coefficient,
            back,
            step,
            scenario,
        )
    except _Failure as failure:
        return LeakBlowdownResult(
            status="invalid",
            scenario_name=scenario,
            hole_diameter=from_si(diameter, "m"),
            number_of_points=declared + 1,
            recorded_points=0,
            points=(),
            total_mass_released=None,
            peak_mass_flow_rate=None,
            time_to_empty=None,
            diagnostics=(LeakBlowdownDiagnostic(code=failure.code, message=str(failure)),),
            diagnostic_codes=(failure.code,),
            warnings=tuple(warnings),
        )

    return LeakBlowdownResult(
        status=outcome.status,
        scenario_name=outcome.scenario_name,
        hole_diameter=from_si(outcome.hole_diameter, "m"),
        number_of_points=outcome.number_of_points,
        recorded_points=len(outcome.points),
        points=tuple(
            LeakBlowdownPoint(
                time=from_si(row.time, "s"),
                mass_flow_rate=from_si(row.mass_flow_rate, "kg/s"),
                temperature=from_si(row.temperature, "K"),
                pressure=from_si(row.pressure, "Pa"),
                vapour_mole_fraction=row.vapour_mole_fraction,
                jet_velocity=from_si(row.jet_velocity, "m/s"),
                jet_momentum=from_si(row.jet_momentum, "N"),
                droplet_smd=from_si(row.droplet_smd, "m"),
            )
            for row in outcome.points
        ),
        total_mass_released=(
            None
            if outcome.total_mass_released is None
            else from_si(outcome.total_mass_released, "kg")
        ),
        peak_mass_flow_rate=(
            None
            if outcome.peak_mass_flow_rate is None
            else from_si(outcome.peak_mass_flow_rate, "kg/s")
        ),
        time_to_empty=(
            None if outcome.time_to_empty is None else from_si(outcome.time_to_empty, "s")
        ),
        diagnostics=tuple(
            LeakBlowdownDiagnostic(code=code, message=message)
            for code, message in outcome.diagnostics
        ),
        diagnostic_codes=tuple(code for code, _ in outcome.diagnostics),
        warnings=tuple(warnings),
    )


def _states(
    components: list[str],
    z: list[float],
    t: float,
    p: float,
    duration: float,
    volume: float,
    diameter: float,
    coefficient: float,
    back: float,
    step: float,
    scenario: str,
) -> _Outcome:
    """The whole calculation in SI magnitudes, for the public function and the layer diff."""
    mixture, ideal_gas = _components.mixture_of(components, eos="pr")
    if len(z) != len(mixture.components):
        raise InvalidInputError(
            "initial_z",
            f"a mixture of {len(mixture.components)} components needs that many mole "
            f"fractions, and {z!r} has {len(z)}",
        )
    check_composition(list(z), len(mixture.components), "initial_z")
    if any(not math.isfinite(value) for value in z):
        raise InvalidInputError("initial_z", "a mole fraction must be finite")

    diagnostics: list[tuple[str, str]] = [
        (
            "MODEL_ASSUMPTIONS",
            "Rigid adiabatic vessel; equilibrium TP flashes; simplified orifice; no slip "
            "or friction in the jet",
        )
    ]
    if p <= back:
        diagnostics.append(
            ("NO_FORWARD_FLOW", "The initial pressure is at or below the receiving pressure")
        )

    state = _flash_state(mixture, ideal_gas, t, p, z)
    initial_mass = state.density * volume
    if not initial_mass > 0.0 or not math.isfinite(initial_mass):
        raise _Failure(
            "DENSITY_UNAVAILABLE",
            f"the initial density is {state.density} kg/m3 and the volume {volume} m3",
        )
    initial_gamma = state.gamma

    area = math.pi * diameter * diameter / 4.0
    number_of_points = math.ceil(duration / step) + 1
    mass = initial_mass
    total = 0.0
    peak = 0.0
    points: list[_Point] = []
    gammas: list[float] = []
    floor_reached = False
    mass_floor = False

    for index in range(number_of_points):
        if not mass > MASS_FLOOR * initial_mass:
            mass_floor = True
            break
        time = index * step
        row, mdot = _record_point(mixture, state, z, back, diameter, coefficient, area, time)
        peak = max(peak, mdot)
        points.append(row)
        gammas.append(state.gamma)

        mass -= mdot * step
        total += mdot * step
        if not mass > MASS_FLOOR * initial_mass:
            mass_floor = True
            break

        ratio = (mass / initial_mass) ** state.gamma
        new_pressure = state.pressure * ratio
        floored = max(new_pressure, back)
        if floored > new_pressure:
            floor_reached = True
        new_temperature = state.temperature * (floored / state.pressure) ** (
            (state.gamma - 1.0) / state.gamma
        )
        state = _flash_state(mixture, ideal_gas, new_temperature, floored, z)

    if floor_reached:
        diagnostics.append(
            (
                "PRESSURE_FLOOR_REACHED",
                "The isentropic update was floored at the receiving pressure",
            )
        )
    if mass_floor:
        diagnostics.append(
            ("MASS_FLOOR_REACHED", "The loop ended below 1% of the initial inventory")
        )
    time_to_empty = total / (0.5 * peak) if peak > 0.0 else 0.0
    return _Outcome(
        status="valid",
        scenario_name=scenario,
        hole_diameter=diameter,
        number_of_points=number_of_points,
        points=tuple(points),
        total_mass_released=total,
        peak_mass_flow_rate=peak,
        time_to_empty=time_to_empty,
        diagnostics=tuple(diagnostics),
        initial_mass=initial_mass,
        initial_gamma=initial_gamma,
        gammas=tuple(gammas),
    )


def _record_point(
    mixture: Any,
    state: _State,
    z: list[float],
    back: float,
    diameter: float,
    coefficient: float,
    area: float,
    time: float,
) -> tuple[_Point, float]:
    molar_mass = molar_mass_of(mixture, z)
    mw = 1000.0 * molar_mass
    if state.pressure <= back:
        mdot = 0.0
    else:
        critical = (2.0 / (state.gamma + 1.0)) ** (state.gamma / (state.gamma - 1.0))
        ratio = back / state.pressure
        if ratio <= critical:
            factor = (2.0 / (state.gamma + 1.0)) ** (
                (state.gamma + 1.0) / (2.0 * (state.gamma - 1.0))
            )
            mdot = (
                coefficient
                * area
                * state.pressure
                * math.sqrt(state.gamma * mw / (R * state.temperature))
                * factor
            )
        else:
            drop = state.pressure - back
            mdot = coefficient * area * math.sqrt(2.0 * state.density * drop)
    if not math.isfinite(mdot) or mdot < 0.0:
        raise _Failure("FLASH_FAILED", f"the rate at {state.pressure} Pa is {mdot} kg/s")
    if mdot <= 0.0:
        velocity = 0.0
    else:
        raw = mdot / (state.density * area * coefficient)
        sound = math.sqrt(state.gamma * R * state.temperature / mw)
        velocity = min(raw, VELOCITY_CAP * sound)
    momentum = mdot * velocity + (state.pressure - back) * area
    smd = _droplet_smd(state, diameter, velocity)
    return (
        _Point(
            time=time,
            mass_flow_rate=mdot,
            temperature=state.temperature,
            pressure=state.pressure,
            vapour_mole_fraction=state.vapour_mole_fraction(),
            jet_velocity=velocity,
            jet_momentum=momentum,
            droplet_smd=smd,
        ),
        mdot,
    )


def _droplet_smd(state: _State, diameter: float, velocity: float) -> float:
    if state.vapour_fraction is None:
        return 0.0
    liquid_fraction = 1.0 - state.vapour_fraction
    if liquid_fraction < SMD_LIQUID_FLOOR:
        return 0.0
    vapour = state.vapour_density
    liquid = state.liquid_density
    if vapour is None or liquid is None or not vapour > 0.0 or not liquid > 0.0:
        return 0.0
    weber = vapour * velocity * velocity * diameter / SURFACE_TENSION
    if math.isnan(weber):
        raise _Failure("DROPLET_SMD_NOT_FINITE", "the droplet Weber number is NaN")
    smd = (
        SMD_COEFFICIENT
        * diameter
        / math.sqrt(weber)
        * math.sqrt(liquid / vapour)
        * liquid_fraction**SMD_LIQUID_EXPONENT
    )
    if math.isnan(smd):
        raise _Failure("DROPLET_SMD_NOT_FINITE", "the droplet Sauter mean diameter is NaN")
    return float(min(smd, diameter / 2.0))


def _flash_state(
    mixture: Any, ideal_gas: Any, temperature: float, pressure: float, z: list[float]
) -> _State:
    reduced = reduced_parameters(mixture, temperature, pressure)
    kij = mixture.effective_kij(reduced.t_kelvin)
    flash = pt_flash(mixture, quantity(temperature, "K"), quantity(pressure, "Pa"), z)
    if flash.phase == "two_phase":
        beta = flash.vapour_fraction
        if beta is None:
            raise _Failure(
                "FLASH_FAILED", "the flash split the feed and reported no vapour fraction"
            )
        vapour = _phase_of(
            mixture, ideal_gas, reduced, kij, temperature, pressure, list(flash.y), flash.z_vapour
        )
        liquid = _phase_of(
            mixture, ideal_gas, reduced, kij, temperature, pressure, list(flash.x), flash.z_liquid
        )
        mass = beta * vapour.molar_mass + (1.0 - beta) * liquid.molar_mass
        volume = beta * vapour.molar_volume + (1.0 - beta) * liquid.molar_volume
        density = mass / volume
        cp = beta * vapour.cp + (1.0 - beta) * liquid.cp
        cv = beta * vapour.cv + (1.0 - beta) * liquid.cv
        gamma = cp / cv
        vapour_fraction: float | None = beta
        vapour_density: float | None = vapour.density
        liquid_density: float | None = liquid.density
    else:
        if flash.phase == "trivial":
            raise _Failure("FLASH_FAILED", f"the flash found a trivial solution at {pressure} Pa")
        liquid_side = flash.phase != "all_vapour"
        root = phase_state(reduced, kij, z, liquid=liquid_side).z
        phase = _phase_of(mixture, ideal_gas, reduced, kij, temperature, pressure, z, root)
        density = phase.density
        gamma = phase.cp / phase.cv
        vapour_fraction = None
        vapour_density = None
        liquid_density = None
    if not density > 0.0 or not math.isfinite(density):
        raise _Failure("DENSITY_UNAVAILABLE", f"the density at {pressure} Pa is {density} kg/m3")
    if not math.isfinite(gamma) or not 1.0 < gamma <= 2.0:
        raise _Failure(
            "GAMMA_OUT_OF_RANGE",
            f"the capacity ratio at {pressure} Pa is {gamma}, and the port refuses one outside "
            "(1, 2]",
        )
    return _State(
        flash=flash,
        pressure=pressure,
        temperature=temperature,
        density=density,
        gamma=gamma,
        vapour_fraction=vapour_fraction,
        vapour_density=vapour_density,
        liquid_density=liquid_density,
    )


def _phase_of(
    mixture: Any,
    ideal_gas: Any,
    reduced: Any,
    kij: Any,
    t: float,
    p: float,
    composition: list[float],
    root: float,
) -> _Phase:
    molar_mass = molar_mass_of(mixture, composition)
    shift = mixture.volume_shift(composition)
    molar_volume = root * MOLAR_GAS_CONSTANT * t / p - shift
    if not molar_volume > 0.0 or not math.isfinite(molar_volume):
        raise _Failure(
            "DENSITY_UNAVAILABLE",
            f"the corrected molar volume at {p} Pa is {molar_volume} m3/mol",
        )
    capacity = phase_capacity(mixture, ideal_gas, reduced, kij, t, p, composition, root)
    return _Phase(
        molar_mass=molar_mass,
        molar_volume=molar_volume,
        density=molar_mass / molar_volume,
        cp=capacity.cp,
        cv=capacity.cv,
    )


def _spec() -> dict[str, object]:
    """The generated spec, imported here so the module has no import-time cycle."""
    from azoth._models_gen import model

    return model(MODEL_ID)
