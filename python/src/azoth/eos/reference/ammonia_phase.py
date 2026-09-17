"""``eos.ammonia_phase`` - the ammonia reference phase state at a temperature and pressure.

Spec: ``specs/models/eos/ammonia_phase.toml``. A multiparameter Helmholtz equation with a
Planck-Einstein ideal-gas part and polynomial, exponential, Gaussian and Gao-B residual
terms, in reduced variables ``tau = Tc/T`` and ``delta = rho/rho_c``.

This is the pure-Python reference: a second, independent expression of the same physics
as the Rust kernel, held to it by ``test_cross_impl.py``.
"""

from __future__ import annotations

import math

from azoth.core.range import apply_checks, checks_for
from azoth.core.result import AmmoniaPhaseResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning

MODEL_ID = "eos.ammonia_phase"

#: The gas constant, J/(mol*K).
R = 8.31446261815324
#: The molar mass, kg/mol.
MOLAR_MASS = 0.01703052
#: The critical temperature, K.
T_CRIT = 405.56
#: The critical molar density, mol/m^3.
RHO_CRIT = 13696.0

_A1 = -6.59406093943886
_A2 = 5.601011519879
_C0 = 3.0
_IDEAL_N = (2.224, 3.148, 0.9579)
_IDEAL_T = (4.0585856593352405, 9.776605187888352, 17.829667620080876)

_N = (
    0.006132232,
    1.7395866,
    -2.2261792,
    -0.30127553,
    0.08967023,
    -0.076387037,
    -0.84063963,
    -0.27026327,
    6.212578,
    -5.7844357,
    2.4817542,
    -2.3739168,
    0.01493697,
    -3.7749264,
    0.0006254348,
    -1.7359e-05,
    -0.13462033,
    0.07749072839,
)
_T = (
    1.0,
    0.382,
    1.0,
    1.0,
    0.677,
    2.915,
    3.51,
    1.063,
    0.655,
    1.3,
    3.1,
    1.4395,
    1.623,
    0.643,
    1.13,
    4.5,
    1.0,
    4.0,
)
_D = (4, 1, 1, 2, 3, 3, 2, 3, 1, 1, 1, 2, 2, 1, 3, 3, 1, 1)
_L = (2, 2, 1)
_G = (1, 1, 1)
_ETA = (0.42776, 0.6424, 0.8175, 0.7995, 0.91, 0.3574, 1.21, 4.14, 22.56, 22.68)
_BETA = (1.708, 1.4865, 2.0915, 2.43, 0.488, 1.1, 0.85, 1.14, 945.64, 993.85)
_GAMMA = (1.036, 1.2777, 1.083, 1.2906, 0.928, 0.934, 0.919, 1.852, 1.05897, 1.05277)
_EPS = (-0.0726, -0.1274, 0.7527, 0.57, 2.2, -0.243, 2.96, 3.02, 0.9574, 0.9576)
_GAOB_N = (-1.6909858, 0.93739074)
_GAOB_T = (4.3315, 4.015)
_GAOB_D = (1.0, 1.0)
_GAOB_ETA = (-2.8452, -2.8342)
_GAOB_BETA = (0.3696, 0.2962)
_GAOB_GAMMA = (1.108, 1.313)
_GAOB_EPS = (0.4478, 0.44689)
_GAOB_B = (1.244, 0.6826)


def _ideal(delta: float, tau: float) -> tuple[float, float, float]:
    alpha = math.log(delta) + _A1 + _A2 * tau + _C0 * math.log(tau)
    dalpha_dtau = _A2 + _C0 / tau
    d2alpha_dtau2 = -_C0 / (tau * tau)
    for n, t in zip(_IDEAL_N, _IDEAL_T, strict=True):
        exp_t = math.exp(-t * tau)
        denom = 1.0 - exp_t
        alpha += n * math.log(denom)
        dalpha_dtau += n * t * exp_t / denom
        d2alpha_dtau2 -= n * t * t * exp_t / (denom * denom)
    return alpha, dalpha_dtau, d2alpha_dtau2


def _residual(delta: float, tau: float) -> dict[str, float]:
    alpha = 0.0
    da_dd = 0.0
    d2a_dd2 = 0.0
    da_dt = 0.0
    d2a_dt2 = 0.0
    d2a_dd_dt = 0.0

    for i in range(5):
        term = _N[i] * delta ** _D[i] * tau ** _T[i]
        alpha += term
        da_dd += term * _D[i] / delta
        d2a_dd2 += term * _D[i] * (_D[i] - 1.0) / (delta * delta)
        da_dt += term * _T[i] / tau
        d2a_dt2 += term * _T[i] * (_T[i] - 1.0) / (tau * tau)
        d2a_dd_dt += term * _D[i] * _T[i] / (delta * tau)

    for i in range(5, 8):
        li = i - 5
        term = _N[i] * delta ** _D[i] * tau ** _T[i] * math.exp(-_G[li] * delta ** _L[li])
        b = _D[i] / delta - _G[li] * _L[li] * delta ** (_L[li] - 1.0)
        alpha += term
        da_dd += term * b
        d2a_dd2 += term * (
            b * b
            - _D[i] / (delta * delta)
            - _G[li] * _L[li] * (_L[li] - 1.0) * delta ** (_L[li] - 2.0)
        )
        da_dt += term * _T[i] / tau
        d2a_dt2 += term * _T[i] * (_T[i] - 1.0) / (tau * tau)
        d2a_dd_dt += term * b * _T[i] / tau

    for i in range(8, 18):
        gi = i - 8
        a = delta - _EPS[gi]
        b = tau - _GAMMA[gi]
        term = (
            _N[i] * delta ** _D[i] * tau ** _T[i] * math.exp(-_ETA[gi] * a * a - _BETA[gi] * b * b)
        )
        b_delta = _D[i] / delta - 2.0 * _ETA[gi] * a
        b_tau = _T[i] / tau - 2.0 * _BETA[gi] * b
        alpha += term
        da_dd += term * b_delta
        d2a_dd2 += term * (b_delta * b_delta - _D[i] / (delta * delta) - 2.0 * _ETA[gi])
        da_dt += term * b_tau
        d2a_dt2 += term * (b_tau * b_tau - _T[i] / (tau * tau) - 2.0 * _BETA[gi])
        d2a_dd_dt += term * b_delta * b_tau

    for i in range(2):
        a = delta - _GAOB_EPS[i]
        y = _GAOB_BETA[i] * (tau - _GAOB_GAMMA[i]) ** 2 + _GAOB_B[i]
        term = (
            _GAOB_N[i]
            * delta ** _GAOB_D[i]
            * tau ** _GAOB_T[i]
            * math.exp(_GAOB_ETA[i] * a * a + 1.0 / y)
        )
        b_delta = _GAOB_D[i] / delta + 2.0 * _GAOB_ETA[i] * a
        b_t = _GAOB_T[i] / tau - (2.0 * _GAOB_BETA[i] * (tau - _GAOB_GAMMA[i])) / (y * y)
        alpha += term
        da_dd += term * b_delta
        d2a_dd2 += term * (b_delta * b_delta - _GAOB_D[i] / (delta * delta) + 2.0 * _GAOB_ETA[i])
        da_dt += term * b_t
        tau_offset = tau - _GAOB_GAMMA[i]
        d2a_dt2 += term * (
            b_t * b_t
            - _GAOB_T[i] / (tau * tau)
            - 2.0 * _GAOB_BETA[i] / (y * y)
            + 8.0 * _GAOB_BETA[i] * _GAOB_BETA[i] * tau_offset * tau_offset / (y * y * y)
        )
        d2a_dd_dt += term * b_delta * b_t

    return {
        "alpha": alpha,
        "alpha_delta": da_dd,
        "alpha_delta_delta": d2a_dd2,
        "alpha_tau": da_dt,
        "alpha_tau_tau": d2a_dt2,
        "alpha_delta_tau": d2a_dd_dt,
    }


def _properties(t: float, rho: float) -> dict[str, float]:
    delta = rho / RHO_CRIT
    tau = T_CRIT / t
    id0, id_t, id_tt = _ideal(delta, tau)
    res = _residual(delta, tau)

    cv = -R * tau * tau * (id_tt + res["alpha_tau_tau"])
    numer = 1.0 + delta * res["alpha_delta"] - delta * tau * res["alpha_delta_tau"]
    denom = 1.0 + 2.0 * delta * res["alpha_delta"] + delta * delta * res["alpha_delta_delta"]
    cp = cv + R * numer * numer / denom

    u = R * t * tau * (id_t + res["alpha_tau"])
    h = R * t * (1.0 + tau * (id_t + res["alpha_tau"]) + delta * res["alpha_delta"])
    s = R * (tau * (id_t + res["alpha_tau"]) - id0 - res["alpha"])
    g = R * t * (1.0 + id0 + res["alpha"] + delta * res["alpha_delta"])

    return {
        "z": 1.0 + delta * res["alpha_delta"],
        "u": u,
        "h": h,
        "s": s,
        "cv": cv,
        "cp": cp,
        "g": g,
    }


def _solve_density(t: float, p: float) -> float:
    rho = p / (R * t)
    for _ in range(100):
        delta = rho / RHO_CRIT
        tau = T_CRIT / t
        res = _residual(delta, tau)
        p_calc = rho * R * t * (1.0 + delta * res["alpha_delta"])
        dpdrho = (
            R
            * t
            * (1.0 + 2.0 * delta * res["alpha_delta"] + delta * delta * res["alpha_delta_delta"])
        )
        step = (p_calc - p) / dpdrho
        rho -= step
        if abs(step) < 1e-8 * abs(rho):
            break
    return rho


def ammonia_phase(T: Q, P: Q) -> AmmoniaPhaseResult:
    """The ammonia phase state at a temperature and pressure.

    Args:
        T: absolute temperature.
        P: absolute pressure.

    Raises:
        OutOfRangeError: if ``T`` or ``P`` is not positive.
    """
    from azoth import _models_gen

    spec = _models_gen.model(MODEL_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    t_si = input_to_si(spec, "T", T)
    p_si = input_to_si(spec, "P", P)
    apply_checks(checks.on_input, {"T": t_si, "P": p_si}.get, warnings)

    rho = _solve_density(t_si, p_si)
    props = _properties(t_si, rho)

    return AmmoniaPhaseResult(
        z_factor=props["z"],
        u=from_si(props["u"], "J/mol"),
        h=from_si(props["h"], "J/mol"),
        s=from_si(props["s"], "J/(mol*K)"),
        cv=from_si(props["cv"], "J/(mol*K)"),
        cp=from_si(props["cp"], "J/(mol*K)"),
        g=from_si(props["g"], "J/mol"),
        warnings=tuple(warnings),
    )
