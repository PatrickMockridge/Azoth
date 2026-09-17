"""``eos.bwrs_phase`` - the BWRS (MBWR-32) phase state at a temperature and pressure.

Spec: ``specs/models/eos/bwrs_phase.toml``. The 32-term modified Benedict-Webb-Rubin
equation of Younglove and Ely, integrated in closed form to the residual Helmholtz
free energy and differentiated for the pressure, the fugacity and the departure
functions, in the coefficients' native mol/L, MPa convention.

This is the pure-Python reference: a second, independent expression of the same
physics as the Rust kernel, held to it by ``test_cross_impl.py``.
"""

from __future__ import annotations

import math
from collections.abc import Sequence

from azoth.core.errors import InvalidInputError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import BwrsPhaseResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning
from azoth.eos.components import BwrsCoefficients

MODEL_ID = "eos.bwrs_phase"

#: The gas constant in the coefficients' native units, L.MPa/(mol.K).
R_MPA = 8.3144621e-3
#: The gas constant in SI, J/(mol.K), for the residual properties.
R_SI = 8.3144621

#: Factorials for the six exponential terms.
_FACT = (1.0, 1.0, 2.0, 6.0, 24.0, 120.0)


def _bp(t: float, a: Sequence[float]) -> list[float]:
    sr = t**0.5
    t2 = t * t
    return [
        R_MPA * t,
        a[0] * t + a[1] * sr + a[2] + a[3] / t + a[4] / t2,
        a[5] * t + a[6] + a[7] / t + a[8] / t2,
        a[9] * t + a[10] + a[11] / t,
        a[12],
        a[13] / t + a[14] / t2,
        a[15] / t,
        a[16] / t + a[17] / t2,
        a[18] / t2,
    ]


def _be(t: float, a: Sequence[float]) -> list[float]:
    t2 = t * t
    t3 = t2 * t
    t4 = t2 * t2
    return [
        a[19] / t2 + a[20] / t3,
        a[21] / t2 + a[22] / t4,
        a[23] / t2 + a[24] / t3,
        a[25] / t2 + a[26] / t4,
        a[27] / t2 + a[28] / t3,
        a[29] / t2 + a[30] / t3 + a[31] / t4,
    ]


def _bp_dt(t: float, a: Sequence[float]) -> list[float]:
    sr = t**0.5
    t2 = t * t
    t3 = t2 * t
    return [
        R_MPA,
        a[0] + a[1] / (2.0 * sr) - a[3] / t2 - 2.0 * a[4] / t3,
        a[5] - a[7] / t2 - 2.0 * a[8] / t3,
        a[9] - a[11] / t2,
        0.0,
        -a[13] / t2 - 2.0 * a[14] / t3,
        -a[15] / t2,
        -a[16] / t2 - 2.0 * a[17] / t3,
        -2.0 * a[18] / t3,
    ]


def _be_dt(t: float, a: Sequence[float]) -> list[float]:
    t2 = t * t
    t3 = t2 * t
    t4 = t2 * t2
    t5 = t4 * t
    return [
        -2.0 * a[19] / t3 - 3.0 * a[20] / t4,
        -2.0 * a[21] / t3 - 4.0 * a[22] / t5,
        -2.0 * a[23] / t3 - 3.0 * a[24] / t4,
        -2.0 * a[25] / t3 - 4.0 * a[26] / t5,
        -2.0 * a[27] / t3 - 3.0 * a[28] / t4,
        -2.0 * a[29] / t3 - 3.0 * a[30] / t4 - 4.0 * a[31] / t5,
    ]


def _bp_dt_dt(t: float, a: Sequence[float]) -> list[float]:
    t15 = t * t**0.5
    t3 = t * t * t
    t4 = t3 * t
    return [
        0.0,
        -a[1] / (4.0 * t15) + 2.0 * a[3] / t3 + 6.0 * a[4] / t4,
        2.0 * a[7] / t3 + 6.0 * a[8] / t4,
        2.0 * a[11] / t3,
        0.0,
        2.0 * a[13] / t3 + 6.0 * a[14] / t4,
        2.0 * a[15] / t3,
        2.0 * a[16] / t3 + 6.0 * a[17] / t4,
        6.0 * a[18] / t4,
    ]


def _be_dt_dt(t: float, a: Sequence[float]) -> list[float]:
    t2 = t * t
    t4 = t2 * t2
    t5 = t4 * t
    t6 = t5 * t
    return [
        6.0 * a[19] / t4 + 12.0 * a[20] / t5,
        6.0 * a[21] / t4 + 20.0 * a[22] / t6,
        6.0 * a[23] / t4 + 12.0 * a[24] / t5,
        6.0 * a[25] / t4 + 20.0 * a[26] / t6,
        6.0 * a[27] / t4 + 12.0 * a[28] / t5,
        6.0 * a[29] / t4 + 12.0 * a[30] / t5 + 20.0 * a[31] / t6,
    ]


def _pressure(rho: float, b: Sequence[float], e: Sequence[float], gamma: float) -> float:
    p = 0.0
    rp = rho
    for bi in b:
        p += bi * rp
        rp *= rho
    el = math.exp(-gamma * rho * rho)
    tail = 0.0
    rp = rho * rho * rho
    for ei in e:
        tail += ei * rp
        rp *= rho * rho
    return p + el * tail


def _d_pressure_drho(rho: float, b: Sequence[float], e: Sequence[float], gamma: float) -> float:
    dp = 0.0
    rp = 1.0
    for i, bi in enumerate(b):
        dp += (i + 1) * bi * rp
        rp *= rho
    el = math.exp(-gamma * rho * rho)
    tail = 0.0
    rp = rho * rho
    for i, ei in enumerate(e):
        n = 3 + 2 * i
        tail += ei * rp * (n - 2.0 * gamma * rho * rho)
        rp *= rho * rho
    return dp + el * tail


def _solve_density(
    t: float, p_target: float, b: Sequence[float], e: Sequence[float], gamma: float
) -> float:
    rho = p_target / (R_MPA * t)
    for _ in range(100):
        step = (_pressure(rho, b, e, gamma) - p_target) / _d_pressure_drho(rho, b, e, gamma)
        rho -= step
        if abs(step) < 1e-12 * abs(rho):
            break
    return rho


def _exp_integrals(rho: float, gamma: float) -> list[float]:
    g2 = gamma * rho * rho
    el = math.exp(-g2)
    out = []
    gpow = 1.0
    for i in range(6):
        gpow *= gamma
        s = 0.0
        pk = 1.0
        for k in range(i + 1):
            s += pk / _FACT[k]
            pk *= g2
        out.append(_FACT[i] / (2.0 * gpow) * (1.0 - el * s))
    return out


def _g_sum(rho: float, e: Sequence[float]) -> tuple[float, float, float]:
    g = 0.0
    gp = 0.0
    gpp = 0.0
    for i, ei in enumerate(e):
        n = 2 * i + 1
        g += ei * rho ** (2 * i + 1)
        gp += ei * n * rho ** (2 * i)
        if i >= 1:
            gpp += ei * n * (n - 1) * rho ** (2 * i - 1)
    return g, gp, gpp


def _helmholtz(t: float, rho: float, b: Sequence[float], e: Sequence[float], gamma: float) -> float:
    pol = 0.0
    rp = rho
    for i in range(1, 9):
        pol += b[i] / i * rp
        rp *= rho
    exp = sum(ti * ei for ti, ei in zip(_exp_integrals(rho, gamma), e, strict=True))
    return (pol + exp) / (R_MPA * t)


def _d_helmholtz_drho(
    t: float, rho: float, b: Sequence[float], e: Sequence[float], gamma: float
) -> float:
    pol = 0.0
    rp = 1.0
    for i in range(1, 9):
        pol += b[i] * rp
        rp *= rho
    g, _, _ = _g_sum(rho, e)
    el = math.exp(-gamma * rho * rho)
    return (pol + el * g) / (R_MPA * t)


def _d_helmholtz_dt(
    t: float,
    rho: float,
    b: Sequence[float],
    e: Sequence[float],
    bt: Sequence[float],
    et: Sequence[float],
    gamma: float,
) -> float:
    pol = 0.0
    rp = rho
    for i in range(1, 9):
        pol += (bt[i] - b[i] / t) / i * rp
        rp *= rho
    integrals = _exp_integrals(rho, gamma)
    exp = sum((et[i] - e[i] / t) * integrals[i] for i in range(6))
    return (pol + exp) / (R_MPA * t)


def _d2_helmholtz_drho2(
    t: float, rho: float, b: Sequence[float], e: Sequence[float], gamma: float
) -> float:
    pol = 0.0
    rp = 1.0
    for i in range(2, 9):
        pol += b[i] * (i - 1) * rp
        rp *= rho
    g, gp, _ = _g_sum(rho, e)
    el = math.exp(-gamma * rho * rho)
    return (pol + el * (gp - 2.0 * gamma * rho * g)) / (R_MPA * t)


def _d2_helmholtz_dt2(
    t: float,
    rho: float,
    b: Sequence[float],
    e: Sequence[float],
    bt: Sequence[float],
    et: Sequence[float],
    btt: Sequence[float],
    ett: Sequence[float],
    gamma: float,
) -> float:
    t2 = t * t
    pol = 0.0
    rp = rho
    for i in range(1, 9):
        pol += (btt[i] - 2.0 * bt[i] / t + 2.0 * b[i] / t2) / i * rp
        rp *= rho
    integrals = _exp_integrals(rho, gamma)
    exp = sum((ett[i] - 2.0 * et[i] / t + 2.0 * e[i] / t2) * integrals[i] for i in range(6))
    return (pol + exp) / (R_MPA * t)


def _d2_helmholtz_dtdrho(
    t: float,
    rho: float,
    b: Sequence[float],
    e: Sequence[float],
    bt: Sequence[float],
    et: Sequence[float],
    gamma: float,
) -> float:
    pol = 0.0
    rp = 1.0
    for i in range(1, 9):
        pol += (bt[i] - b[i] / t) * rp
        rp *= rho
    g_t = 0.0
    rp = rho
    for i in range(6):
        g_t += (et[i] - e[i] / t) * rp
        rp *= rho * rho
    el = math.exp(-gamma * rho * rho)
    return (pol + el * g_t) / (R_MPA * t)


def _departure(
    t: float,
    rho: float,
    b: Sequence[float],
    bt: Sequence[float],
    btt: Sequence[float],
    e: Sequence[float],
    et: Sequence[float],
    ett: Sequence[float],
    gamma: float,
) -> tuple[float, float, float, float, float, float]:
    f = _helmholtz(t, rho, b, e, gamma)
    phi_rho = rho * _d_helmholtz_drho(t, rho, b, e, gamma)
    phi_t = t * _d_helmholtz_dt(t, rho, b, e, bt, et, gamma)
    phi_rho_rho = rho * rho * _d2_helmholtz_drho2(t, rho, b, e, gamma)
    phi_t_t = t * t * _d2_helmholtz_dt2(t, rho, b, e, bt, et, btt, ett, gamma)
    phi_rho_t = rho * t * _d2_helmholtz_dtdrho(t, rho, b, e, bt, et, gamma)

    z = 1.0 + phi_rho
    a_res = R_SI * t * f
    s_res = R_SI * (math.log(z) - f - phi_t)
    h_res = R_SI * t * (phi_rho - phi_t)
    g_res = R_SI * t * (f + phi_rho - math.log(z))
    cv_res = -R_SI * (2.0 * phi_t + phi_t_t)
    cp_res = R_SI * (
        -(2.0 * phi_t + phi_t_t) + (z + phi_rho_t) ** 2 / (1.0 + 2.0 * phi_rho + phi_rho_rho) - 1.0
    )
    return a_res, s_res, h_res, g_res, cv_res, cp_res


def _mix(
    x: Sequence[float],
    bps: Sequence[Sequence[float]],
    bes: Sequence[Sequence[float]],
    rhocs: Sequence[float],
) -> tuple[list[float], list[float], float]:
    mb = [sum(xj * bps[j][k] for j, xj in enumerate(x)) for k in range(9)]
    me = [sum(xj * bes[j][k] for j, xj in enumerate(x)) for k in range(6)]
    rhoc_mix = sum(xj * rhocs[j] for j, xj in enumerate(x))
    return mb, me, 1.0 / (rhoc_mix * rhoc_mix)


def _extensive_helmholtz(
    t: float,
    rho: float,
    x: Sequence[float],
    i: int,
    delta: float,
    bps: Sequence[Sequence[float]],
    bes: Sequence[Sequence[float]],
    rhocs: Sequence[float],
) -> float:
    n = sum(x)
    n_new = n + delta
    x_new = [(xj * n + (delta if j == i else 0.0)) / n_new for j, xj in enumerate(x)]
    rho_new = rho * n_new / n
    mb, me, gamma = _mix(x_new, bps, bes, rhocs)
    return n_new * _helmholtz(t, rho_new, mb, me, gamma)


def _d_f_dn(
    t: float,
    rho: float,
    x: Sequence[float],
    bps: Sequence[Sequence[float]],
    bes: Sequence[Sequence[float]],
    rhocs: Sequence[float],
) -> list[float]:
    n = sum(x)
    dn = n / 10_000.0
    out = []
    for i in range(len(x)):
        f_plus = _extensive_helmholtz(t, rho, x, i, dn, bps, bes, rhocs)
        f_minus = _extensive_helmholtz(t, rho, x, i, -dn, bps, bes, rhocs)
        out.append((f_plus - f_minus) / (2.0 * dn))
    return out


def _ln_fugacity(
    t: float,
    rho: float,
    x: Sequence[float],
    bps: Sequence[Sequence[float]],
    bes: Sequence[Sequence[float]],
    rhocs: Sequence[float],
) -> list[float]:
    mb, me, gamma = _mix(x, bps, bes, rhocs)
    ln_z = math.log(1.0 + rho * _d_helmholtz_drho(t, rho, mb, me, gamma))
    return [d - ln_z for d in _d_f_dn(t, rho, x, bps, bes, rhocs)]


def bwrs_phase(
    coeffs: Sequence[BwrsCoefficients], T: Q, P: Q, z: Sequence[float]
) -> BwrsPhaseResult:
    """The BWRS phase state of a mixture of MBWR-32 substances.

    Args:
        coeffs: the per-component coefficient sets, resolved by name against
            :func:`azoth.eos.components.bwrs_coefficients`.
        T: absolute temperature.
        P: absolute pressure.
        z: the mixture's mole fractions.

    Raises:
        InvalidInputError: if ``z`` is not a composition of ``coeffs``'s length.
        OutOfRangeError: if ``T`` or ``P`` is not positive.
    """
    from azoth import _models_gen

    spec = _models_gen.model(MODEL_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    t_si = input_to_si(spec, "T", T)
    p_si = input_to_si(spec, "P", P)

    apply_checks(checks.on_input, {"T": t_si, "P": p_si}.get, warnings)

    n = len(coeffs)
    if len(z) != n:
        raise InvalidInputError(
            "z",
            f"a mixture of {n} components needs {n} mole fractions, but z has {len(z)}",
        )
    if any(value < 0.0 for value in z):
        raise InvalidInputError("z", "a mole fraction cannot be negative")
    total = sum(z)
    if abs(total - 1.0) > 1.0e-9:
        raise InvalidInputError(
            "z",
            f"the composition sums to {total}, not to one; renormalising it here "
            "would hide a caller's error",
        )

    t = t_si
    p_mpa = p_si / 1.0e6

    bps = [_bp(t, c.a) for c in coeffs]
    bes = [_be(t, c.a) for c in coeffs]
    bts = [_bp_dt(t, c.a) for c in coeffs]
    ets = [_be_dt(t, c.a) for c in coeffs]
    btts = [_bp_dt_dt(t, c.a) for c in coeffs]
    etts = [_be_dt_dt(t, c.a) for c in coeffs]
    rhocs = [c.rhoc for c in coeffs]

    mb, me, gamma = _mix(z, bps, bes, rhocs)
    mbt, met, _ = _mix(z, bts, ets, rhocs)
    mbtt, mett, _ = _mix(z, btts, etts, rhocs)

    rho = _solve_density(t, p_mpa, mb, me, gamma)
    z_factor = 1.0 + rho * _d_helmholtz_drho(t, rho, mb, me, gamma)
    ln_phi = _ln_fugacity(t, rho, z, bps, bes, rhocs)
    _, s_res, h_res, _, _, cp_res = _departure(t, rho, mb, mbt, mbtt, me, met, mett, gamma)

    return BwrsPhaseResult(
        z_factor=z_factor,
        ln_phi=tuple(ln_phi),
        h_res=from_si(h_res, "J/mol"),
        s_res=from_si(s_res, "J/(mol*K)"),
        cp_res=from_si(cp_res, "J/(mol*K)"),
        warnings=tuple(warnings),
    )
