"""``eos.co2_phase`` - the Span-Wagner CO2 phase state at a temperature and pressure.

Spec: ``specs/models/eos/co2_phase.toml``. A multiparameter Helmholtz equation with a
Planck-Einstein ideal-gas part and power-exponential and Gaussian residual terms, in
reduced variables ``tau = Tc/T`` and ``delta = rho/rho_c``.

This is the pure-Python reference: a second, independent expression of the same physics
as the Rust kernel, held to it by ``test_cross_impl.py``.
"""

from __future__ import annotations

import math

from azoth.core.range import apply_checks, checks_for
from azoth.core.result import Co2PhaseResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning

MODEL_ID = "eos.co2_phase"

#: The gas constant, J/(mol*K) - Span-Wagner's published value.
R = 8.31451
#: The critical temperature, K.
T_CRIT = 304.1282
#: The critical molar density, mol/m^3.
RHO_CRIT = 10624.9063

_A1 = -6.1248710624319
_A2 = 5.11559631801453
_C0 = 2.5
_IDEAL_N = (1.99427042, 0.62105248, 0.41195293, 1.04028922, 0.08327678)
_IDEAL_T = (3.15163, 6.1119, 6.77708, 11.32384, 27.08792)

_N = (
    0.388568232032,
    2.93854759427,
    -5.5867188535,
    -0.767531995925,
    0.317290055804,
    0.548033158978,
    0.122794112203,
    2.16589615432,
    1.58417351097,
    -0.231327054055,
    0.0581169164314,
    -0.553691372054,
    0.489466159094,
    -0.0242757398435,
    0.0624947905017,
    -0.121758602252,
    -0.370556852701,
    -0.0167758797004,
    -0.11960736638,
    -0.0456193625088,
    0.0356127892703,
    -0.00744277271321,
    -0.00173957049024,
    -0.0218101212895,
    0.0243321665592,
    -0.0374401334235,
    0.143387157569,
    -0.134919690833,
    -0.0231512250535,
    0.0123631254929,
    0.00210583219729,
    -0.000339585190264,
    0.00559936517716,
    -0.000303351180556,
)
_D = (
    1,
    1,
    1,
    1,
    2,
    2,
    3,
    1,
    2,
    4,
    5,
    5,
    5,
    6,
    6,
    6,
    1,
    1,
    4,
    4,
    4,
    7,
    8,
    2,
    3,
    3,
    5,
    5,
    6,
    7,
    8,
    10,
    4,
    8,
)
_T = (
    0,
    0.75,
    1,
    2,
    0.75,
    2,
    0.75,
    1.5,
    1.5,
    2.5,
    0,
    1.5,
    2,
    0,
    1,
    2,
    3,
    6,
    3,
    6,
    8,
    6,
    0,
    7,
    12,
    16,
    22,
    24,
    16,
    24,
    8,
    2,
    28,
    14,
)
_L = (
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    2,
    2,
    2,
    2,
    2,
    2,
    2,
    3,
    3,
    3,
    4,
    4,
    4,
    4,
    4,
    4,
    5,
    6,
)
_GN = (-213.654886883, 26641.5691493, -24027.2122046, -283.41603424, 212.472844002)
_GD = (2, 2, 2, 3, 3)
_GT = (1, 0, 1, 3, 3)
_GBETA = (325, 300, 300, 275, 275)
_GGAMMA = (1.16, 1.19, 1.19, 1.25, 1.22)
_GEPS = (1, 1, 1, 1, 1)
_GETA = (25, 25, 25, 15, 20)


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

    for i in range(34):
        expc = 1.0 if _L[i] == 0 else math.exp(-(delta ** _L[i]))
        common = _N[i] * delta ** _D[i] * tau ** _T[i] * expc
        dterm = _D[i] / delta - (0.0 if _L[i] == 0 else _L[i] * delta ** (_L[i] - 1))
        alpha += common
        da_dd += common * dterm
        d2a_dd2 += common * (
            dterm * dterm
            - _D[i] / (delta * delta)
            - (0.0 if _L[i] == 0 else _L[i] * (_L[i] - 1) * delta ** (_L[i] - 2))
        )
        tterm = _T[i] / tau
        da_dt += common * tterm
        d2a_dt2 += common * tterm * (tterm - 1.0 / tau)
        d2a_dd_dt += common * dterm * tterm

    for i in range(5):
        dr = delta - _GEPS[i]
        tr = tau - _GGAMMA[i]
        common = (
            _GN[i]
            * delta ** _GD[i]
            * tau ** _GT[i]
            * math.exp(-_GETA[i] * dr * dr - _GBETA[i] * tr * tr)
        )
        dterm = _GD[i] / delta - 2 * _GETA[i] * dr
        alpha += common
        da_dd += common * dterm
        d2a_dd2 += common * (dterm * dterm - _GD[i] / (delta * delta) - 2 * _GETA[i])
        tterm = _GT[i] / tau - 2 * _GBETA[i] * tr
        da_dt += common * tterm
        d2a_dt2 += common * (tterm * tterm - _GT[i] / (tau * tau) - 2 * _GBETA[i])
        d2a_dd_dt += common * dterm * tterm

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

    z = 1.0 + delta * res["alpha_delta"]
    at = id_t + res["alpha_tau"]
    h = R * t * (1.0 + tau * at + delta * res["alpha_delta"])
    s = R * (tau * at - id0 - res["alpha"])
    cv = -R * tau * tau * (id_tt + res["alpha_tau_tau"])
    numer = 1.0 + delta * res["alpha_delta"] - delta * tau * res["alpha_delta_tau"]
    denom = 1.0 + 2.0 * delta * res["alpha_delta"] + delta * delta * res["alpha_delta_delta"]
    cp = cv + R * numer * numer / denom
    u = R * t * tau * at
    g = R * t * (id0 + res["alpha"] + 1.0 + delta * res["alpha_delta"])

    return {"z": z, "u": u, "h": h, "s": s, "cv": cv, "cp": cp, "g": g}


def _solve_density(t: float, p: float) -> float:
    delta = p / (R * t) / RHO_CRIT
    tau = T_CRIT / t
    for _ in range(100):
        res = _residual(delta, tau)
        f = R * t * RHO_CRIT * delta * (1.0 + delta * res["alpha_delta"]) - p
        df = (
            R
            * t
            * RHO_CRIT
            * (1.0 + 2.0 * delta * res["alpha_delta"] + delta * delta * res["alpha_delta_delta"])
        )
        new_delta = delta - f / df
        if abs(new_delta - delta) < 1e-12:
            return new_delta * RHO_CRIT
        delta = new_delta
    return delta * RHO_CRIT


def co2_phase(T: Q, P: Q) -> Co2PhaseResult:
    """The Span-Wagner CO2 phase state at a temperature and pressure.

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

    return Co2PhaseResult(
        z_factor=props["z"],
        u=from_si(props["u"], "J/mol"),
        h=from_si(props["h"], "J/mol"),
        s=from_si(props["s"], "J/(mol*K)"),
        cv=from_si(props["cv"], "J/(mol*K)"),
        cp=from_si(props["cp"], "J/(mol*K)"),
        g=from_si(props["g"], "J/mol"),
        warnings=tuple(warnings),
    )
