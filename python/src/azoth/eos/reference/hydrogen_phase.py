"""``eos.hydrogen_phase`` - the Leachman hydrogen phase state at a temperature and pressure.

Spec: ``specs/models/eos/hydrogen_phase.toml``. A multiparameter Helmholtz equation with a
Planck-Einstein ideal-gas part and polynomial, exponential and Gaussian residual terms, for
the three hydrogen spin-isomers (normal, para, ortho), in reduced variables ``tau = Tc/T``
and ``delta = rho/rho_c``.

This is the pure-Python reference: a second, independent expression of the same physics
as the Rust kernel, held to it by ``test_cross_impl.py``.
"""

from __future__ import annotations

import math

from azoth.core.errors import InvalidInputError, OutOfRangeError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import HydrogenPhaseResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning

MODEL_ID = "eos.hydrogen_phase"

#: The gas constant, J/(mol*K).
R = 8.31451

#: The densest reduced density the dense-root solve considers, and the dilute end of the same
#: search. Guards rather than physical bounds: the residual's fitted terms run out before the
#: densest state reached here, which is `delta = 2.46`.
_MAXIMUM_REDUCED_DENSITY = 4.0
_MINIMUM_REDUCED_DENSITY = 1.0e-3

#: The bracket expansion in log density - ten per cent a step, as the solid volume solve uses
#: - and the two caps.
_DENSITY_EXPANSION = 0.0953102
_MAXIMUM_BRACKET_STEPS = 200
_MAXIMUM_BISECTION_STEPS = 200

_TYPES = {
    "normal": (
        33.145,
        15508.0,
        (-1.4579856475, 1.888076782, 1.616, -0.4117, -0.792, 0.758, 1.217),
        (0, 0, -16.0205159149, -22.6580178006, -60.0090511389, -74.9434303817, -206.9392065168),
    ),
    "para": (
        32.938,
        15538.0,
        (
            -1.4485891134,
            1.884521239,
            4.30256,
            13.0289,
            -47.7365,
            50.0013,
            -18.6261,
            0.993973,
            0.536078,
        ),
        (
            0,
            0,
            -15.1496751472,
            -25.0925982148,
            -29.4735563787,
            -35.4059141417,
            -40.724998482,
            -163.7925799988,
            -309.2173173842,
        ),
    ),
    "ortho": (
        33.22,
        15445.0,
        (-1.4675442336, 1.8845068862, 2.54151, -2.3661, 1.00365, 1.22447),
        (0, 0, -25.7676098736, -43.4677904877, -66.0445514750, -209.7531607465),
    ),
}

_N = {
    "normal": (
        -6.93643,
        0.01,
        2.1101,
        4.52059,
        0.732564,
        -1.34086,
        0.130985,
        -0.777414,
        0.351944,
        -0.0211716,
        0.0226312,
        0.032187,
        -0.0231752,
        0.0557346,
    ),
    "para": (
        -7.33375,
        0.01,
        2.60375,
        4.66279,
        0.68239,
        -1.47078,
        0.135801,
        -1.05327,
        0.328239,
        -0.057783,
        0.044974,
        0.070346,
        -0.040176,
        0.11951,
    ),
    "ortho": (
        -6.83148,
        0.01,
        2.11505,
        4.38353,
        0.211292,
        -1.00939,
        0.142086,
        -0.87696,
        0.804927,
        -0.710775,
        0.0639688,
        0.0710858,
        -0.087654,
        0.647088,
    ),
}
_T = {
    "normal": (
        0.6844,
        1,
        0.989,
        0.489,
        0.803,
        1.1444,
        1.409,
        1.754,
        1.311,
        4.187,
        5.646,
        0.791,
        7.249,
        2.986,
    ),
    "para": (
        0.6855,
        1,
        1,
        0.489,
        0.774,
        1.133,
        1.386,
        1.619,
        1.162,
        3.96,
        5.276,
        0.99,
        6.791,
        3.19,
    ),
    "ortho": (
        0.7333,
        1,
        1.1372,
        0.5136,
        0.5638,
        1.6248,
        1.829,
        2.404,
        2.105,
        4.1,
        7.658,
        1.259,
        7.589,
        3.946,
    ),
}
_D = (1, 4, 1, 1, 2, 2, 3, 1, 3, 2, 1, 3, 1, 1)
_PHI = {
    "normal": (0, 0, 0, 0, 0, 0, 0, 0, 0, -1.685, -0.489, -0.103, -2.506, -1.607),
    "para": (0, 0, 0, 0, 0, 0, 0, 0, 0, -1.7437, -0.5516, -0.0634, -2.1341, -1.777),
    "ortho": (0, 0, 0, 0, 0, 0, 0, 0, 0, -1.169, -0.894, -0.04, -2.072, -1.306),
}
_BETA = {
    "normal": (0, 0, 0, 0, 0, 0, 0, 0, 0, -0.171, -0.2245, -0.1304, -0.2785, -0.3967),
    "para": (0, 0, 0, 0, 0, 0, 0, 0, 0, -0.194, -0.2019, -0.0301, -0.2383, -0.3253),
    "ortho": (0, 0, 0, 0, 0, 0, 0, 0, 0, -0.4555, -0.4046, -0.0869, -0.4415, -0.5743),
}
_GAMMA = {
    "normal": (0, 0, 0, 0, 0, 0, 0, 0, 0, 0.7164, 1.3444, 1.4517, 0.7204, 1.5445),
    "para": (0, 0, 0, 0, 0, 0, 0, 0, 0, 0.8048, 1.5248, 0.6648, 0.6832, 1.493),
    "ortho": (0, 0, 0, 0, 0, 0, 0, 0, 0, 1.5444, 0.6627, 0.763, 0.6587, 1.4327),
}
_EPS = {
    "normal": (0, 0, 0, 0, 0, 0, 0, 0, 0, 1.506, 0.156, 1.736, 0.67, 1.662),
    "para": (0, 0, 0, 0, 0, 0, 0, 0, 0, 1.5487, 0.1785, 1.28, 0.6319, 1.7104),
    "ortho": (0, 0, 0, 0, 0, 0, 0, 0, 0, 0.6366, 0.3876, 0.9437, 0.3976, 0.9626),
}


def _ideal(delta: float, tau: float, ht: str) -> tuple[float, float, float]:
    _, _, a0, b0 = _TYPES[ht]
    alpha = math.log(delta) + 1.5 * math.log(tau) + a0[0] + a0[1] * tau
    da_dt = 1.5 / tau + a0[1]
    d2a_dt2 = -1.5 / (tau * tau)
    for k in range(2, len(a0)):
        e = math.exp(b0[k] * tau)
        denom = 1.0 - e
        alpha += a0[k] * math.log(denom)
        da_dt += -a0[k] * b0[k] * e / denom
        d2a_dt2 += -a0[k] * b0[k] * b0[k] * e / (denom * denom)
    return alpha, da_dt, d2a_dt2


def _residual(delta: float, tau: float, ht: str) -> dict[str, float]:
    n = _N[ht]
    t = _T[ht]
    phi = _PHI[ht]
    beta = _BETA[ht]
    gamma = _GAMMA[ht]
    eps = _EPS[ht]

    alpha = 0.0
    da_dd = 0.0
    d2a_dd2 = 0.0
    da_dt = 0.0
    d2a_dt2 = 0.0
    d2a_dd_dt = 0.0

    for i in range(7):
        term = n[i] * delta ** _D[i] * tau ** t[i]
        alpha += term
        da_dd += term * _D[i] / delta
        d2a_dd2 += term * _D[i] * (_D[i] - 1) / (delta * delta)
        da_dt += term * t[i] / tau
        d2a_dt2 += term * t[i] * (t[i] - 1) / (tau * tau)
        d2a_dd_dt += term * _D[i] * t[i] / (delta * tau)

    for i in range(7, 9):
        term = n[i] * delta ** _D[i] * tau ** t[i] * math.exp(-delta)
        u = _D[i] / delta - 1.0
        u_prime = -_D[i] / (delta * delta)
        alpha += term
        da_dd += term * u
        d2a_dd2 += term * (u * u + u_prime)
        da_dt += term * t[i] / tau
        d2a_dt2 += term * t[i] * (t[i] - 1) / (tau * tau)
        d2a_dd_dt += term * u * t[i] / tau

    for i in range(9, 14):
        dr = delta - eps[i]
        tr = tau - gamma[i]
        term = n[i] * delta ** _D[i] * tau ** t[i] * math.exp(phi[i] * dr * dr + beta[i] * tr * tr)
        u = _D[i] / delta + 2.0 * phi[i] * dr
        u_prime = -_D[i] / (delta * delta) + 2.0 * phi[i]
        v = t[i] / tau + 2.0 * beta[i] * tr
        v_prime = -t[i] / (tau * tau) + 2.0 * beta[i]
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


def _properties(t: float, rho: float, ht: str) -> dict[str, float]:
    tc, rhoc, _, _ = _TYPES[ht]
    delta = rho / rhoc
    tau = tc / t
    id0, id_t, id_tt = _ideal(delta, tau, ht)
    res = _residual(delta, tau, ht)

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


def _pressure_residual(t: float, rho: float, p: float, ht: str) -> float:
    """The pressure a state carries at a molar density, minus the one asked for."""
    tc, rhoc, _, _ = _TYPES[ht]
    delta = rho / rhoc
    res = _residual(delta, tc / t, ht)
    return rho * R * t * (1.0 + delta * res["alpha_delta"]) - p


def _solve_density(t: float, p: float, ht: str) -> float:
    """The **dilute** root: Newton from the ideal-gas guess, which follows whichever branch
    that guess is on. For a compressed liquid that is not the liquid's - see
    :func:`_solve_density_dense`."""
    tc, rhoc, _, _ = _TYPES[ht]
    rho = p / (R * t)
    for _ in range(100):
        delta = rho / rhoc
        tau = tc / t
        res = _residual(delta, tau, ht)
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


def _solve_density_dense(t: float, p: float, ht: str) -> float | None:
    """The **dense** root, or ``None`` where the state has none in the fitted range.

    Bracketed from the dense side **down**, which is the whole difference from
    :func:`_solve_density`: below the critical temperature the isotherm crosses a pressure
    three times, so walking up from a dilute start meets the vapour's crossing first and
    stops there. The bracket is then bisected, which needs no derivative and cannot leave
    the interval.
    """
    _tc, rhoc, _, _ = _TYPES[ht]
    upper = math.log(_MAXIMUM_REDUCED_DENSITY * rhoc)
    upper_residual = _pressure_residual(t, math.exp(upper), p, ht)
    if not math.isfinite(upper_residual) or upper_residual < 0.0:
        return None

    lower = upper - _DENSITY_EXPANSION
    lower_residual = _pressure_residual(t, math.exp(lower), p, ht)
    bracketed = False
    for _ in range(_MAXIMUM_BRACKET_STEPS):
        if math.isfinite(lower_residual) and lower_residual * upper_residual < 0.0:
            bracketed = True
            break
        upper, upper_residual = lower, lower_residual
        lower -= _DENSITY_EXPANSION
        if math.exp(lower) <= _MINIMUM_REDUCED_DENSITY * rhoc:
            break
        lower_residual = _pressure_residual(t, math.exp(lower), p, ht)
    if not bracketed:
        return None

    # The residual rises with density on the dense branch, so a positive one means the trial
    # is above the root. Bisection in log density, which is monotone in the same direction.
    lo, hi = lower, upper
    for _ in range(_MAXIMUM_BISECTION_STEPS):
        if hi - lo < 1.0e-14:
            break
        mid = 0.5 * (lo + hi)
        trial = _pressure_residual(t, math.exp(mid), p, ht)
        if not math.isfinite(trial):
            return None
        if trial > 0.0:
            hi = mid
        else:
            lo = mid
    return math.exp(0.5 * (lo + hi))


def hydrogen_phase(
    T: Q, P: Q, hydrogen_type: str = "normal", compressed_phase: str = "vapour"
) -> HydrogenPhaseResult:
    """The Leachman hydrogen phase state at a temperature and pressure.

    Args:
        T: absolute temperature.
        P: absolute pressure.
        hydrogen_type: the spin-isomer, one of ``normal``, ``para`` or ``ortho``.
        compressed_phase: which root of the isotherm is wanted. Below the critical
            temperature the two are different states at the same temperature and pressure -
            at 13.8 K and 7042 Pa they are ``Z = 0.985`` and ``Z = 0.0016``.

    Raises:
        OutOfRangeError: if ``T`` or ``P`` is not positive, or the dense root was asked for
            at a state whose dense root is outside the range the equation is fitted to.
        InvalidInputError: if ``hydrogen_type`` is not a known isomer, or
            ``compressed_phase`` is neither ``liquid`` nor ``vapour``.
    """
    from azoth import _models_gen

    spec = _models_gen.model(MODEL_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    t_si = input_to_si(spec, "T", T)
    p_si = input_to_si(spec, "P", P)
    apply_checks(checks.on_input, {"T": t_si, "P": p_si}.get, warnings)

    ht = hydrogen_type.strip().lower()
    if ht not in _TYPES:
        raise InvalidInputError(
            "hydrogen_type",
            f"unknown hydrogen type {hydrogen_type!r}; expected normal, para or ortho",
        )

    if compressed_phase == "vapour":
        rho = _solve_density(t_si, p_si, ht)
    elif compressed_phase == "liquid":
        dense = _solve_density_dense(t_si, p_si, ht)
        if dense is None:
            raise OutOfRangeError(
                "compressed_phase",
                p_si,
                "the dense root was asked for, and this state has none within the density "
                "range the equation is fitted to: the isotherm does not cross the pressure "
                "there below the ceiling the solve brackets from. `vapour` is the root such a "
                "state is on.",
            )
        rho = dense
    else:
        raise InvalidInputError(
            "compressed_phase",
            f"{compressed_phase!r}; expected `liquid` or `vapour`",
        )
    props = _properties(t_si, rho, ht)

    return HydrogenPhaseResult(
        z_factor=props["z"],
        u=from_si(props["u"], "J/mol"),
        h=from_si(props["h"], "J/mol"),
        s=from_si(props["s"], "J/(mol*K)"),
        cv=from_si(props["cv"], "J/(mol*K)"),
        cp=from_si(props["cp"], "J/(mol*K)"),
        g=from_si(props["g"], "J/mol"),
        warnings=tuple(warnings),
    )
