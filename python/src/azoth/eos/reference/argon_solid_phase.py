"""``eos.argon_solid_phase`` - the solid argon phase state at a temperature and pressure.

Spec: ``specs/models/eos/argon_solid_phase.toml``. The Maltby-Hammer-Wilhelmsen Helmholtz
equation for solid argon: a three-body-corrected Buckingham static lattice, a zero-point
term, Debye and Einstein vibrations and the anharmonic corrections. The volume is solved
by bracketing and bisection/Newton in logarithmic volume.

This is the pure-Python reference: a second, independent expression of the same physics
as the Rust kernel, held to it by ``test_cross_impl.py``.
"""

from __future__ import annotations

import math

from azoth.core.range import apply_checks, checks_for
from azoth.core.result import ArgonSolidPhaseResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning
from azoth.eos.reference._solid_helmholtz import Dual

MODEL_ID = "eos.argon_solid_phase"

R = 8.314462618

_AVOGADRO = 6.02214076e23
_BOLTZMANN = R / _AVOGADRO
_EPSILON = _BOLTZMANN * 134.7
_REPULSIVE_EXPONENT = 14.19
_MINIMUM_POTENTIAL_DISTANCE = 3.802e-10
_THREE_BODY_COEFFICIENT = _BOLTZMANN * 3.202e5 * 1.0e-90

_A_D = 10.4
_B_D = 0.02503
_C_D = 0.0001568
_V0 = 22.56e-6
_THETA_D0 = 92.0
_GAMMA_D0 = 2.563
_Q_D = 0.2874

_B1 = -0.0004475
_B2 = 2.041e-6
_B3 = 5.75e-7
_C1 = 0.7204
_C2 = -1.614
_C3 = -0.01943
_C4 = -27.64

_EINSTEIN_WEIGHTS = (0.0261, 0.03784, 0.04512)
_EINSTEIN_TEMPERATURES = (77.81, 550.0, 45.36)
_EINSTEIN_GRUNEISEN = (6.221, 1.617e-6, 3.1278)

_Z1 = 140.0
_Z2 = 2.34
_Z3 = 0.683
_Z4 = 19.7e-6

_MINIMUM_VOLUME = 1.0e-7
_MAXIMUM_VOLUME = 1.0e-3


def _buckingham_zero_distance() -> float:
    lower = 0.5 * _MINIMUM_POTENTIAL_DISTANCE
    upper = _MINIMUM_POTENTIAL_DISTANCE
    for _ in range(100):
        middle = 0.5 * (lower + upper)
        if _buckingham_pair_potential_scalar(middle) > 0.0:
            lower = middle
        else:
            upper = middle
    return 0.5 * (lower + upper)


def _buckingham_pair_potential_scalar(distance: float) -> float:
    repulsive = (
        _EPSILON
        * 6.0
        / (_REPULSIVE_EXPONENT - 6.0)
        * math.exp(_REPULSIVE_EXPONENT * (1.0 - distance / _MINIMUM_POTENTIAL_DISTANCE))
    )
    attractive = (
        _EPSILON
        * _REPULSIVE_EXPONENT
        / (_REPULSIVE_EXPONENT - 6.0)
        * math.pow(_MINIMUM_POTENTIAL_DISTANCE / distance, 6.0)
    )
    return repulsive - attractive


def _buckingham_pair_potential(distance: Dual) -> Dual:
    repulsive = (
        distance * (-_REPULSIVE_EXPONENT / _MINIMUM_POTENTIAL_DISTANCE) + _REPULSIVE_EXPONENT
    ).exp() * (_EPSILON * 6.0 / (_REPULSIVE_EXPONENT - 6.0))
    attractive = (distance.reciprocal() * _MINIMUM_POTENTIAL_DISTANCE).powf(6.0) * (
        _EPSILON * _REPULSIVE_EXPONENT / (_REPULSIVE_EXPONENT - 6.0)
    )
    return repulsive - attractive


def _fcc_shells() -> list[tuple[int, int]]:
    counts: dict[int, int] = {}
    for h in range(-24, 25):
        for k in range(-24, 25):
            for l in range(-24, 25):
                squared = h * h + k * k + l * l
                if squared == 0 or (h + k + l) % 2 != 0:
                    continue
                counts[squared] = counts.get(squared, 0) + 1
    return sorted(counts.items())[:20]


def _static_lattice_energy(volume: Dual) -> Dual:
    lattice_constant = (volume * (4.0 / _AVOGADRO)).powf(1.0 / 3.0)
    three_body_factor = 1.0 - volume.reciprocal() * (
        _THREE_BODY_COEFFICIENT * _AVOGADRO / (_EPSILON * _buckingham_zero_distance() ** 6)
    )
    shell_energy = Dual.constant(0.0)
    for squared_distance, coordination in _fcc_shells():
        distance = lattice_constant * (0.5 * math.sqrt(squared_distance))
        shell_energy = (
            shell_energy + _buckingham_pair_potential(distance) * three_body_factor * coordination
        )
    return shell_energy * (0.5 * _AVOGADRO)


def _helmholtz(temperature: float, molar_volume: float) -> Dual:
    temp = Dual.temperature(temperature)
    volume = Dual.volume(molar_volume)
    reduced_volume = volume / _V0

    static_lattice = _static_lattice_energy(volume)
    zero_point = ((_Z2 / _Z3) * (1.0 - (reduced_volume * (_V0 / _Z4)).powf(_Z3))).exp() * (_Z1 * R)

    theta_d_temperature = _THETA_D0 + _A_D * (
        (-_B_D * temp.powf(2.0) - _C_D * temp.powf(3.0)).exp() - 1.0
    )
    theta_d = ((_GAMMA_D0 / _Q_D) * (1.0 - reduced_volume.powf(_Q_D))).exp() * theta_d_temperature
    weight_sum = sum(_EINSTEIN_WEIGHTS)
    debye = (theta_d / temp).debye_free_energy() * temp * (3.0 * R * (1.0 - weight_sum))

    einstein = Dual.constant(0.0)
    for weight, temperature, gruneisen in zip(
        _EINSTEIN_WEIGHTS, _EINSTEIN_TEMPERATURES, _EINSTEIN_GRUNEISEN, strict=True
    ):
        theta = (gruneisen * (1.0 - reduced_volume)).exp() * temperature
        einstein = einstein + (theta / temp).log_one_minus_exp_negative() * weight
    einstein = einstein * temp * (3.0 * R)

    temperature_ratio = temp / _THETA_D0
    anharmonic_temperature = (
        temperature_ratio.powf(4.0)
        / (temperature_ratio.powf(2.0) * _B2 + 1.0)
        * theta_d_temperature
        * _B1
        * (_B3 * (reduced_volume - 1.0)).exp()
        * R
    )
    anharmonic_cold = (
        _C1 * (_C2 * (1.0 - reduced_volume)).exp()
        + reduced_volume.reciprocal() * _C3 * (_C4 * (1.0 - reduced_volume)).exp()
    ) * R

    return static_lattice + zero_point + debye + einstein + anharmonic_temperature + anharmonic_cold


def _pressure(temperature: float, molar_volume: float) -> float:
    return -_helmholtz(temperature, molar_volume).dv


def _pressure_residual(temperature: float, log_volume: float, pressure_pa: float) -> float:
    return _pressure(temperature, math.exp(log_volume)) - pressure_pa


def _solve_volume(temperature: float, pressure_pa: float) -> float:
    lower = math.log(_V0)
    upper = lower
    lower_residual = _pressure_residual(temperature, lower, pressure_pa)
    upper_residual = lower_residual
    expansion = math.log(1.02)
    if lower_residual > 0.0:
        iteration = 0
        while upper_residual > 0.0 and iteration < 600:
            lower = upper
            lower_residual = upper_residual
            upper = min(math.log(_MAXIMUM_VOLUME), upper + expansion)
            upper_residual = _pressure_residual(temperature, upper, pressure_pa)
            iteration += 1
    else:
        iteration = 0
        while lower_residual < 0.0 and iteration < 600:
            upper = lower
            upper_residual = lower_residual
            lower = max(math.log(_MINIMUM_VOLUME), lower - expansion)
            lower_residual = _pressure_residual(temperature, lower, pressure_pa)
            iteration += 1

    current = 0.5 * (lower + upper)
    for _ in range(100):
        volume = math.exp(current)
        helmholtz = _helmholtz(temperature, volume)
        residual = -helmholtz.dv - pressure_pa
        if abs(residual) / max(pressure_pa, 1000.0) < 1.0e-9:
            return volume
        derivative = -helmholtz.dvv * volume
        candidate = current - residual / derivative
        if derivative >= 0.0 or candidate <= lower or candidate >= upper:
            candidate = 0.5 * (lower + upper)
        candidate_residual = _pressure_residual(temperature, candidate, pressure_pa)
        if lower_residual * candidate_residual <= 0.0:
            upper = candidate
        else:
            lower = candidate
            lower_residual = candidate_residual
        current = candidate
    raise RuntimeError("solid argon volume solver did not converge")


def argon_solid_phase(T: Q, P: Q) -> ArgonSolidPhaseResult:
    """The solid argon phase state at a temperature and pressure.

    Args:
        T: absolute temperature, at most 300 K.
        P: absolute pressure, at most 16 GPa.

    Raises:
        OutOfRangeError: if ``T`` or ``P`` is outside the published range.
    """
    from azoth import _models_gen

    spec = _models_gen.model(MODEL_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    t_si = input_to_si(spec, "T", T)
    p_si = input_to_si(spec, "P", P)
    apply_checks(checks.on_input, {"T": t_si, "P": p_si}.get, warnings)

    molar_volume = _solve_volume(t_si, p_si)
    helmholtz = _helmholtz(t_si, molar_volume)
    entropy = -helmholtz.dt
    internal_energy = helmholtz.value + t_si * entropy

    return ArgonSolidPhaseResult(
        z_factor=p_si * molar_volume / (R * t_si),
        u=from_si(internal_energy, "J/mol"),
        h=from_si(internal_energy + p_si * molar_volume, "J/mol"),
        s=from_si(entropy, "J/(mol*K)"),
        cv=from_si(-t_si * helmholtz.dtt, "J/(mol*K)"),
        cp=from_si(
            -t_si * helmholtz.dtt + t_si * helmholtz.dtv * helmholtz.dtv / helmholtz.dvv,
            "J/(mol*K)",
        ),
        g=from_si(helmholtz.value + p_si * molar_volume, "J/mol"),
        warnings=tuple(warnings),
    )
