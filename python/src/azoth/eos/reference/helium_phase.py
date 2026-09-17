"""``eos.helium_phase`` - the Vega helium phase state at a temperature and pressure.

Spec: ``specs/models/eos/helium_phase.toml``. A multiparameter Helmholtz equation with a
monatomic ideal-gas part and polynomial, exponential and Gaussian residual terms, in
reduced variables ``tau = Tc/T`` and ``delta = rho/rho_c``.

This is the pure-Python reference: a second, independent expression of the same physics
as the Rust kernel, held to it by ``test_cross_impl.py``.
"""

from __future__ import annotations

import math

from azoth.core.range import apply_checks, checks_for
from azoth.core.result import HeliumPhaseResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning

MODEL_ID = "eos.helium_phase"

#: The gas constant, J/(mol*K) - the NIST helium EOS value.
R = 8.314472
#: The critical temperature, K.
T_CRIT = 5.1953
#: The critical molar density, mol/m^3.
RHO_CRIT = 17383.7

_A1 = 0.1733487932835764
_A2 = 0.4674522201550815

_N = (
    0.015559018,
    3.0638932,
    -4.2420844,
    0.054418088,
    -0.18971904,
    0.087856262,
    2.2833566,
    -0.53331595,
    -0.53296502,
    0.99444915,
    -0.30078896,
    -1.6432563,
    0.8029102,
    0.026838669,
    0.04687678,
    -0.14832766,
    0.03016211,
    -0.019986041,
    0.14283514,
    0.007418269,
    -0.22989793,
    0.79224829,
    -0.049386338,
)
_T = (
    1.0,
    0.425,
    0.63,
    0.69,
    1.83,
    0.575,
    0.925,
    1.585,
    1.69,
    1.51,
    2.9,
    0.8,
    1.26,
    3.51,
    2.785,
    1.0,
    4.22,
    0.83,
    1.575,
    3.447,
    0.73,
    1.634,
    6.13,
)
_D = (4, 1, 1, 2, 2, 3, 1, 1, 3, 2, 2, 1, 2, 1, 2, 1, 1, 3, 2, 2, 3, 2, 2)
_L = (0, 0, 0, 0, 0, 0, 1, 2, 2, 1, 2, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
_ETA = (
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    1.5497,
    9.245,
    4.76323,
    6.3826,
    8.7023,
    0.255,
    0.3523,
    0.1492,
    0.05,
    0.1668,
    42.2358,
)
_BETA = (
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0.2471,
    0.0983,
    0.1556,
    2.6782,
    2.7077,
    0.6621,
    0.1775,
    0.4821,
    0.3069,
    0.1758,
    1357.6577,
)
_GAMMA = (
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    3.15,
    2.54505,
    1.2513,
    1.9416,
    0.5984,
    2.2282,
    1.606,
    3.815,
    1.61958,
    0.6407,
    1.076,
)
_EPS = (
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0.596,
    0.3423,
    0.761,
    0.9747,
    0.5868,
    0.5627,
    2.5346,
    3.6763,
    4.5245,
    5.039,
    0.959,
)


def _ideal(delta: float, tau: float) -> tuple[float, float, float]:
    return (
        _A1 + _A2 * tau + math.log(delta) + 1.5 * math.log(tau),
        _A2 + 1.5 / tau,
        -1.5 / (tau * tau),
    )


def _residual(delta: float, tau: float) -> dict[str, float]:
    alpha = 0.0
    da_dd = 0.0
    d2a_dd2 = 0.0
    da_dt = 0.0
    d2a_dt2 = 0.0
    d2a_dd_dt = 0.0

    for i in range(6):
        term = _N[i] * delta ** _D[i] * tau ** _T[i]
        alpha += term
        da_dd += term * _D[i] / delta
        d2a_dd2 += term * _D[i] * (_D[i] - 1) / (delta * delta)
        da_dt += term * _T[i] / tau
        d2a_dt2 += term * _T[i] * (_T[i] - 1) / (tau * tau)
        d2a_dd_dt += term * _D[i] * _T[i] / (delta * tau)

    for i in range(6, 12):
        term = _N[i] * delta ** _D[i] * tau ** _T[i] * math.exp(-(delta ** _L[i]))
        u = _D[i] / delta - _L[i] * delta ** (_L[i] - 1)
        u_prime = -_D[i] / (delta * delta) - _L[i] * (_L[i] - 1) * delta ** (_L[i] - 2)
        alpha += term
        da_dd += term * u
        d2a_dd2 += term * (u * u + u_prime)
        da_dt += term * _T[i] / tau
        d2a_dt2 += term * _T[i] * (_T[i] - 1) / (tau * tau)
        d2a_dd_dt += term * u * _T[i] / tau

    for i in range(12, 23):
        dr = delta - _EPS[i]
        tr = tau - _GAMMA[i]
        term = (
            _N[i]
            * delta ** _D[i]
            * tau ** _T[i]
            * math.exp(-_ETA[i] * dr * dr - _BETA[i] * tr * tr)
        )
        u = _D[i] / delta - 2.0 * _ETA[i] * dr
        u_prime = -_D[i] / (delta * delta) - 2.0 * _ETA[i]
        v = _T[i] / tau - 2.0 * _BETA[i] * tr
        v_prime = -_T[i] / (tau * tau) - 2.0 * _BETA[i]
        alpha += term
        da_dd += term * u
        d2a_dd2 += term * (u * u + u_prime)
        da_dt += term * v
        d2a_dt2 += term * (v * v + v_prime)
        d2a_dd_dt += term * u * v

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

    at = id_t + res["alpha_tau"]
    numer = 1.0 + delta * res["alpha_delta"] - delta * tau * res["alpha_delta_tau"]
    denom = 1.0 + 2.0 * delta * res["alpha_delta"] + delta * delta * res["alpha_delta_delta"]

    z = 1.0 + delta * res["alpha_delta"]
    u = R * t * tau * at
    h = R * t * (1.0 + delta * res["alpha_delta"] + tau * at)
    s = R * (tau * at - id0 - res["alpha"])
    cv = -R * tau * tau * (id_tt + res["alpha_tau_tau"])
    cp = cv + R * numer * numer / denom
    g = R * t * (1.0 + delta * res["alpha_delta"] + id0 + res["alpha"])

    return {"z": z, "u": u, "h": h, "s": s, "cv": cv, "cp": cp, "g": g}


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
        if abs(step) < 1e-12 * abs(rho):
            break
    return rho


def helium_phase(T: Q, P: Q) -> HeliumPhaseResult:
    """The Vega helium phase state at a temperature and pressure.

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

    return HeliumPhaseResult(
        z_factor=props["z"],
        u=from_si(props["u"], "J/mol"),
        h=from_si(props["h"], "J/mol"),
        s=from_si(props["s"], "J/(mol*K)"),
        cv=from_si(props["cv"], "J/(mol*K)"),
        cp=from_si(props["cp"], "J/(mol*K)"),
        g=from_si(props["g"], "J/mol"),
        warnings=tuple(warnings),
    )
