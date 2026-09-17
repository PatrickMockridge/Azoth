"""``eos.eos_cg_phase`` - the EOS-CG phase state at a temperature, pressure and composition.

Spec: ``specs/models/eos/eos_cg_phase.toml``. The 28-component combustion-gas Helmholtz
model: composition-averaged reducing parameters, an ideal-gas and a residual Helmholtz
part, and a pairwise departure contribution. The density is solved in logarithmic volume.

This is the pure-Python reference: a second, independent expression of the same physics
as the Rust kernel, held to it by ``test_cross_impl.py``.
"""

from __future__ import annotations

import math

from azoth.core.errors import InvalidInputError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import EosCgPhaseResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning
from azoth.eos.reference import _eos_cg_data as data

MODEL_ID = "eos.eos_cg_phase"

R = data.R
NCOMP = data.NCOMP

_COMPONENT_NAMES = [
    "",
    "methane",
    "nitrogen",
    "CO2",
    "ethane",
    "propane",
    "i-butane",
    "n-butane",
    "i-pentane",
    "n-pentane",
    "n-hexane",
    "n-heptane",
    "n-octane",
    "n-nonane",
    "nC10",
    "hydrogen",
    "oxygen",
    "CO",
    "water",
    "H2S",
    "helium",
    "argon",
    "SO2",
    "ammonia",
    "chlorine",
    "HCl",
    "MEA",
    "DEA",
    "MDEA",
]
_NAME_TO_INDEX = {name: i for i, name in enumerate(_COMPONENT_NAMES) if name}

_EPSILON = 1e-15


def _composition(components: list[str], z: list[float]) -> list[float]:
    x = [0.0] * (NCOMP + 1)
    for name, fraction in zip(components, z, strict=True):
        index = _NAME_TO_INDEX.get(name)
        if index is None:
            raise InvalidInputError("components", f"unknown EOS-CG component {name!r}")
        x[index] += fraction
    return x


def _reducing_parameters(x: list[float]) -> tuple[float, float]:
    tr = 0.0
    vr = 0.0
    for i in range(1, NCOMP + 1):
        if x[i] > _EPSILON:
            f = 1.0
            for j in range(i, NCOMP + 1):
                if x[j] > _EPSILON:
                    xij = f * (x[i] * x[j]) * (x[i] + x[j])
                    vr += xij * data.GVIJ[i][j] / (data.BVIJ[i][j] * x[i] + x[j])
                    tr += xij * data.GTIJ[i][j] / (data.BTIJ[i][j] * x[i] + x[j])
                    f = 2.0
    dr = 1.0 / vr if vr > _EPSILON else 0.0
    return tr, dr


def _alpha0(t: float, d: float, x: list[float]) -> list[float]:
    log_d = math.log(d) if d > _EPSILON else math.log(_EPSILON)
    log_t = math.log(t)
    a0 = [0.0, 0.0, 0.0]
    for i in range(1, NCOMP + 1):
        if x[i] > _EPSILON:
            log_xd = log_d + math.log(x[i])
            sum_hyp0 = 0.0
            sum_hyp1 = 0.0
            sum_hyp2 = 0.0
            for j in range(4, 8):
                if data.TH0I[i][j] < -0.5:
                    term = data.N0I[i][j] * t / data.TC[i]
                    a0[0] += x[i] * term
                    a0[1] += x[i] * (-term)
                    a0[2] += x[i] * (2.0 * term)
                elif data.TH0I[i][j] > _EPSILON:
                    th0t = data.TH0I[i][j] / t
                    ep = math.exp(th0t)
                    em = 1.0 / ep
                    hsn = (ep - em) / 2.0
                    hcn = (ep + em) / 2.0
                    if j == 4 or j == 6:
                        log_hyp = math.log(abs(hsn))
                        sum_hyp0 += data.N0I[i][j] * log_hyp
                        sum_hyp1 += data.N0I[i][j] * th0t * hcn / hsn
                        sum_hyp2 += data.N0I[i][j] * (th0t / hsn) * (th0t / hsn)
                    else:
                        log_hyp = math.log(abs(hcn))
                        sum_hyp0 -= data.N0I[i][j] * log_hyp
                        sum_hyp1 -= data.N0I[i][j] * th0t * hsn / hcn
                        sum_hyp2 += data.N0I[i][j] * (th0t / hcn) * (th0t / hcn)
            a0[0] += x[i] * (
                log_xd + data.N0I[i][1] + data.N0I[i][2] / t - data.N0I[i][3] * log_t + sum_hyp0
            )
            a0[1] += x[i] * (data.N0I[i][3] + data.N0I[i][2] / t + sum_hyp1)
            a0[2] += -x[i] * (data.N0I[i][3] + sum_hyp2)
    return a0


def _alphar(itau: bool, t: float, d: float, x: list[float]) -> list[list[float]]:
    ar = [[0.0] * 4 for _ in range(4)]
    tr, dr = _reducing_parameters(x)
    del_ = d / dr
    tau = tr / t
    lntau = math.log(tau)

    delp = [0.0] * 8
    expd = [0.0] * 8
    delp[1] = del_
    expd[1] = math.exp(-delp[1])
    for i in range(2, 8):
        delp[i] = delp[i - 1] * del_
        expd[i] = math.exp(-delp[i])

    # The tau-dependent parts, with the short-form sharing of the propane exponents.
    taup0 = [0.0] * 13
    for k in range(1, data.KPOL[5] + data.KEXP[5] + 1):
        taup0[k] = math.exp(data.TOIK[5][k] * lntau)
    taup = [[0.0] * 25 for _ in range(NCOMP + 1)]
    for i in range(1, NCOMP + 1):
        if x[i] > _EPSILON:
            short = i > 4 and i != 15 and i != 18 and i != 20 and i != 22 and i < 23
            for k in range(1, data.KPOL[i] + data.KEXP[i] + data.KGAUSS[i] + 1):
                if short:
                    taup[i][k] = data.NOIK[i][k] * taup0[k]
                else:
                    taup[i][k] = data.NOIK[i][k] * math.exp(data.TOIK[i][k] * lntau)

    # Pure-fluid contributions.
    for i in range(1, NCOMP + 1):
        if x[i] <= _EPSILON:
            continue
        for k in range(1, data.KPOL[i] + 1):
            ndt = x[i] * delp[data.DOIK[i][k]] * taup[i][k]
            ndtd = ndt * data.DOIK[i][k]
            ar[0][1] += ndtd
            ar[0][2] += ndtd * (data.DOIK[i][k] - 1)
            if itau:
                ndtt = ndt * data.TOIK[i][k]
                ar[0][0] += ndt
                ar[1][0] += ndtt
                ar[2][0] += ndtt * (data.TOIK[i][k] - 1.0)
                ar[1][1] += ndtt * data.DOIK[i][k]
                ar[1][2] += ndtt * data.DOIK[i][k] * (data.DOIK[i][k] - 1)
                ar[0][3] += ndtd * (data.DOIK[i][k] - 1) * (data.DOIK[i][k] - 2)
        for k in range(1 + data.KPOL[i], data.KPOL[i] + data.KEXP[i] + 1):
            ndt = x[i] * delp[data.DOIK[i][k]] * taup[i][k] * expd[data.COIK[i][k]]
            ex = data.COIK[i][k] * delp[data.COIK[i][k]]
            ex2 = data.DOIK[i][k] - ex
            ex3 = ex2 * (ex2 - 1.0)
            ar[0][1] += ndt * ex2
            ar[0][2] += ndt * (ex3 - data.COIK[i][k] * ex)
            if itau:
                ndtt = ndt * data.TOIK[i][k]
                ar[0][0] += ndt
                ar[1][0] += ndtt
                ar[2][0] += ndtt * (data.TOIK[i][k] - 1.0)
                ar[1][1] += ndtt * ex2
                ar[1][2] += ndtt * (ex3 - data.COIK[i][k] * ex)
                ar[0][3] += ndt * (
                    ex3 * (ex2 - 2.0) - ex * (3.0 * ex2 - 3.0 + data.COIK[i][k]) * data.COIK[i][k]
                )
        for k in range(
            1 + data.KPOL[i] + data.KEXP[i], data.KPOL[i] + data.KEXP[i] + data.KGAUSS[i] + 1
        ):
            d_val = data.DOIK[i][k]
            t_val = data.TOIK[i][k]
            eta = data.ETA_PURE[i][k]
            eps = data.EPSILON_PURE[i][k]
            beta = data.BETA_PURE[i][k]
            gam = data.GAMMA_PURE[i][k]
            delta = del_ - eps
            theta = tau - gam
            exp_val = math.exp(-eta * delta * delta - beta * theta * theta)
            ndt = x[i] * delp[d_val] * taup[i][k] * exp_val

            term_d = d_val - 2.0 * eta * del_ * delta
            ar[0][1] += ndt * term_d
            delta_term_d_prime = -2.0 * eta * del_ * (2.0 * del_ - eps)
            term_d2 = term_d * term_d - term_d + delta_term_d_prime
            ar[0][2] += ndt * term_d2

            if itau:
                ar[0][0] += ndt
                term_t = t_val - 2.0 * beta * tau * theta
                ar[1][0] += ndt * term_t
                tau_term_t_prime = -2.0 * beta * tau * (2.0 * tau - gam)
                term_t2 = term_t * term_t - term_t + tau_term_t_prime
                ar[2][0] += ndt * term_t2
                ar[1][1] += ndt * term_t * term_d
                ar[1][2] += ndt * term_t * term_d2
                delta2_term_d_double_prime = -4.0 * eta * del_ * del_
                delta_b_prime = 2.0 * term_d * delta_term_d_prime + delta2_term_d_double_prime
                term_d3 = term_d * term_d2 + delta_b_prime - 2.0 * term_d2
                ar[0][3] += ndt * term_d3

    # Mixture departure contributions.
    for i in range(1, NCOMP):
        if x[i] <= _EPSILON:
            continue
        for j in range(i + 1, NCOMP + 1):
            if x[j] <= _EPSILON:
                continue
            mn = data.MNUMB[i][j]
            if mn < 0:
                continue
            xijf = x[i] * x[j] * data.FIJ[i][j]
            for k in range(1, data.KPOLIJ[mn] + 1):
                taupijk = data.NIJK[mn][k] * math.exp(data.TIJK[mn][k] * lntau)
                ndt = xijf * delp[data.DIJK[mn][k]] * taupijk
                ndtd = ndt * data.DIJK[mn][k]
                ar[0][1] += ndtd
                ar[0][2] += ndtd * (data.DIJK[mn][k] - 1)
                if itau:
                    ndtt = ndt * data.TIJK[mn][k]
                    ar[0][0] += ndt
                    ar[1][0] += ndtt
                    ar[2][0] += ndtt * (data.TIJK[mn][k] - 1.0)
                    ar[1][1] += ndtt * data.DIJK[mn][k]
                    ar[1][2] += ndtt * data.DIJK[mn][k] * (data.DIJK[mn][k] - 1)
                    ar[0][3] += ndtd * (data.DIJK[mn][k] - 1) * (data.DIJK[mn][k] - 2)
            for k in range(1 + data.KPOLIJ[mn], data.KPOLIJ[mn] + data.KEXPIJ[mn] + 1):
                cij0 = data.CIJK[mn][k] * delp[2]
                eij0 = data.EIJK[mn][k] * del_
                ndt = (
                    xijf
                    * data.NIJK[mn][k]
                    * delp[data.DIJK[mn][k]]
                    * math.exp(cij0 + eij0 + data.GIJK[mn][k] + data.TIJK[mn][k] * lntau)
                )
                ex = data.DIJK[mn][k] + 2.0 * cij0 + eij0
                ex2 = ex * ex - data.DIJK[mn][k] + 2.0 * cij0
                ar[0][1] += ndt * ex
                ar[0][2] += ndt * ex2
                if itau:
                    ndtt = ndt * data.TIJK[mn][k]
                    ar[0][0] += ndt
                    ar[1][0] += ndtt
                    ar[2][0] += ndtt * (data.TIJK[mn][k] - 1.0)
                    ar[1][1] += ndtt * ex
                    ar[1][2] += ndtt * ex2
                    ar[0][3] += ndt * (
                        ex * (ex2 - 2.0 * (data.DIJK[mn][k] - 2.0 * cij0)) + 2.0 * data.DIJK[mn][k]
                    )

    return ar


def _molar_mass(x: list[float]) -> float:
    return sum(x[i] * data.MM[i] for i in range(1, NCOMP + 1))


def _properties(t: float, d: float, x: list[float]) -> dict[str, float]:
    mm = _molar_mass(x)
    a0 = _alpha0(t, d, x)
    ar = _alphar(True, t, d, x)

    rt = R * t
    z = 1.0 + ar[0][1]
    pressure_kpa = d * rt * z
    dpdd = rt * (1.0 + 2.0 * ar[0][1] + ar[0][2])
    dpdt = d * R * (1.0 + ar[0][1] - ar[1][1])

    u = rt * (a0[1] + ar[1][0])
    h = rt * (1.0 + ar[0][1] + a0[1] + ar[1][0])
    s = R * (a0[1] + ar[1][0] - a0[0] - ar[0][0])
    cv = -R * (a0[2] + ar[2][0])
    g = rt * (1.0 + ar[0][1] + a0[0] + ar[0][0])

    cp = cv + t * (dpdt / d) * (dpdt / d) / dpdd if d > _EPSILON else cv + R
    w = math.sqrt(max(0.0, 1000.0 * cp / cv * dpdd / mm)) if cv > _EPSILON else 0.0

    return {"z": z, "p": pressure_kpa, "u": u, "h": h, "s": s, "cv": cv, "cp": cp, "g": g, "w": w}


def _pressure(t: float, d: float, x: list[float]) -> tuple[float, float]:
    ar = _alphar(False, t, d, x)
    z = 1.0 + ar[0][1]
    return d * R * t * z, R * t * (1.0 + 2.0 * ar[0][1] + ar[0][2])


def _pseudo_critical(x: list[float]) -> tuple[float, float]:
    tcx = sum(x[i] * data.TC[i] for i in range(1, NCOMP + 1))
    vcx = sum(x[i] / data.DC[i] for i in range(1, NCOMP + 1))
    dcx = 1.0 / vcx if vcx > _EPSILON else 0.0
    return tcx, dcx


def _solve_density(t: float, pressure_kpa: float, x: list[float]) -> float:
    _, dcx = _pseudo_critical(x)
    d = pressure_kpa / R / t
    plog = math.log(pressure_kpa)
    vlog = -math.log(d)
    n_fail = 0
    i_fail = 0

    for it in range(1, 51):
        if not -7.0 <= vlog <= 100.0 or it in (20, 30, 40) or i_fail == 1:
            i_fail = 0
            if n_fail > 2:
                return pressure_kpa / R / t
            n_fail += 1
            d = {1: dcx * 3.0, 2: dcx * 2.5}.get(n_fail, dcx * 2.0)
            vlog = -math.log(d)
        d = math.exp(-vlog)
        p2, dpdd = _pressure(t, d, x)
        if dpdd < _EPSILON or p2 < _EPSILON:
            vinc = -0.1 if d > dcx else 0.1
            if it > 5:
                vinc /= 2.0
            if 10 < it < 20:
                vinc /= 5.0
            vlog += vinc
        else:
            dpdlv = -d * dpdd
            vdiff = (math.log(p2) - plog) * p2 / dpdlv
            vlog += -vdiff
            if abs(vdiff) < 1.0e-7:
                if dpdd < 0.0:
                    i_fail = 1
                else:
                    return math.exp(-vlog)
    return pressure_kpa / R / t


def eos_cg_phase(T: Q, P: Q, components: list[str], z: list[float]) -> EosCgPhaseResult:
    """The EOS-CG phase state at a temperature, pressure and composition.

    Args:
        T: absolute temperature.
        P: absolute pressure.
        components: the EOS-CG component names, one per entry.
        z: the mole fractions, one per component, summing to one.

    Raises:
        OutOfRangeError: if ``T`` or ``P`` is not positive.
        InvalidInputError: if a component name is not an EOS-CG component.
    """
    from azoth import _models_gen

    spec = _models_gen.model(MODEL_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    t_si = input_to_si(spec, "T", T)
    p_si = input_to_si(spec, "P", P)
    apply_checks(checks.on_input, {"T": t_si, "P": p_si}.get, warnings)

    if len(components) != len(z):
        raise InvalidInputError(
            "z", f"{len(components)} components but {len(z)} fractions; the two must match"
        )
    x = _composition(components, z)
    molar_density = _solve_density(t_si, p_si / 1000.0, x)
    props = _properties(t_si, molar_density, x)

    return EosCgPhaseResult(
        z_factor=props["z"],
        u=from_si(props["u"], "J/mol"),
        h=from_si(props["h"], "J/mol"),
        s=from_si(props["s"], "J/(mol*K)"),
        cv=from_si(props["cv"], "J/(mol*K)"),
        cp=from_si(props["cp"], "J/(mol*K)"),
        g=from_si(props["g"], "J/mol"),
        warnings=tuple(warnings),
    )
