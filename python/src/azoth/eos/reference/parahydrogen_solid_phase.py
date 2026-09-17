"""``eos.parahydrogen_solid_phase`` - the solid para-hydrogen phase state at a temperature
and pressure.

Spec: ``specs/models/eos/parahydrogen_solid_phase.toml``. The Sannerhaugen Helmholtz
equation for hcp phase-I solid para-hydrogen: a Vinet cold curve, Debye and three Einstein
external modes, an internal-vibration mode and the anharmonic terms. The volume is solved
by bracketing and bisection/Newton in logarithmic volume.

This is the pure-Python reference: a second, independent expression of the same physics
as the Rust kernel, held to it by ``test_cross_impl.py``.
"""

from __future__ import annotations

import math

from azoth.core.range import apply_checks, checks_for
from azoth.core.result import ParahydrogenSolidPhaseResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning
from azoth.eos.reference._solid_helmholtz import Dual

MODEL_ID = "eos.parahydrogen_solid_phase"

R = 8.3144621

_TRIPLE_POINT_TEMPERATURE = 13.8033

_MINIMUM_VOLUME = 1.0e-7
_MAXIMUM_VOLUME = 1.0e-3

_V0 = 23.14e-6
_B0 = 180.5e6
_B0_PRIME = 7.131

_THETA_D0 = 127.0
_GAMMA_D0 = 2.843
_Q_D = 1.647
_A_D = 0.0008734
_B_D = 0.9295
_C_D = 1.020

_EINSTEIN_WEIGHTS = (0.09616, 0.1597, 0.004324)
_EINSTEIN_TEMPERATURES = (51.62, 247.0, 187.9)
_EINSTEIN_GRUNEISEN = (2.290, 0.9952, 1.367)

_THETA_INTERNAL0 = 1086.0
_GAMMA_INTERNAL = 2.775
_EXTERNAL_MODES = 5.0

_B1 = 0.08337
_B2 = 22.97
_B3 = -5.5412
_C1 = -1.704
_C2 = -0.08174
_C3 = -0.03986
_C4 = 3.828


def _vinet_energy(reduced_volume: Dual) -> Dual:
    x = reduced_volume.powf(1.0 / 3.0)
    exponent_coefficient = 1.5 * (_B0_PRIME - 1.0)
    vinet_exponent = (Dual.constant(1.0) - x) * exponent_coefficient
    return ((vinet_exponent - 1.0) * vinet_exponent.exp() + 1.0) * (
        4.0 * _B0 * _V0 / (_B0_PRIME - 1.0) ** 2
    )


def _helmholtz(temperature: float, molar_volume: float) -> Dual:
    temp = Dual.temperature(temperature)
    volume = Dual.volume(molar_volume)
    reduced_volume = volume / _V0

    vinet = _vinet_energy(reduced_volume)

    theta_d_temperature = _THETA_D0 + _A_D * (
        (-_B_D * temp.powf(2.0) - _C_D * temp.powf(3.0)).exp() - 1.0
    )
    theta_d = ((_GAMMA_D0 / _Q_D) * (1.0 - reduced_volume.powf(_Q_D))).exp() * theta_d_temperature
    weight_sum = sum(_EINSTEIN_WEIGHTS)
    debye = (theta_d / temp).debye_free_energy() * temp * (R * _EXTERNAL_MODES * (1.0 - weight_sum))

    einstein = Dual.constant(0.0)
    for weight, temperature, gruneisen in zip(
        _EINSTEIN_WEIGHTS, _EINSTEIN_TEMPERATURES, _EINSTEIN_GRUNEISEN, strict=True
    ):
        theta = (gruneisen * (1.0 - reduced_volume)).exp() * temperature
        einstein = einstein + (theta / temp).log_one_minus_exp_negative() * weight
    einstein = einstein * temp * (R * _EXTERNAL_MODES)

    theta_internal = (_GAMMA_INTERNAL * (1.0 - reduced_volume)).exp() * _THETA_INTERNAL0
    internal = (theta_internal / temp).log_one_minus_exp_negative() * temp * R

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

    return vinet + debye + einstein + internal + anharmonic_temperature + anharmonic_cold


def _pressure(temperature: float, molar_volume: float) -> float:
    return -_helmholtz(temperature, molar_volume).dv


def _pressure_residual(temperature: float, log_volume: float, pressure_pa: float) -> float:
    return _pressure(temperature, math.exp(log_volume)) - pressure_pa


def _solve_volume(temperature: float, pressure_pa: float) -> float:
    guessed_volume = _V0
    if temperature >= _TRIPLE_POINT_TEMPERATURE:
        melting_volume = (27.1788 - 0.283044 * temperature) * 1.0e-6
        if melting_volume > _MINIMUM_VOLUME:
            guessed_volume = melting_volume

    lower = math.log(guessed_volume)
    upper = lower
    lower_residual = _pressure_residual(temperature, lower, pressure_pa)
    upper_residual = lower_residual
    expansion = math.log(1.1)
    for _ in range(160):
        if lower_residual * upper_residual < 0.0:
            break
        lower = max(math.log(_MINIMUM_VOLUME), lower - expansion)
        upper = min(math.log(_MAXIMUM_VOLUME), upper + expansion)
        lower_residual = _pressure_residual(temperature, lower, pressure_pa)
        upper_residual = _pressure_residual(temperature, upper, pressure_pa)

    current = min(upper, max(lower, math.log(guessed_volume)))
    for _ in range(100):
        volume = math.exp(current)
        helmholtz = _helmholtz(temperature, volume)
        residual = -helmholtz.dv - pressure_pa
        if abs(residual) / max(pressure_pa, 1000.0) < 1.0e-10:
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
    raise RuntimeError("solid para-hydrogen volume solver did not converge")


def parahydrogen_solid_phase(T: Q, P: Q) -> ParahydrogenSolidPhaseResult:
    """The solid para-hydrogen phase state at a temperature and pressure.

    Args:
        T: absolute temperature, at most 200 K.
        P: absolute pressure, at most 10 GPa.

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

    return ParahydrogenSolidPhaseResult(
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
