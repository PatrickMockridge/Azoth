"""``eos.water_phase`` - the IAPWS-IF97 water phase state at a temperature and pressure.

Spec: ``specs/models/eos/water_phase.toml``. The Gibbs-energy steam tables: Region 1
(subcooled liquid) and Region 2 (superheated vapour) as dimensionless Gibbs functions in
``pi = p/p*`` and ``tau = T*/T``, selected by the Region 4 saturation curve. IF97 is
mass-based, so this reference converts the specific property set to molar SI.

This is the pure-Python reference: a second, independent expression of the same physics
as the Rust kernel, held to it by ``test_cross_impl.py``.
"""

from __future__ import annotations

import math

from azoth.core.range import apply_checks, checks_for
from azoth.core.result import WaterPhaseResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning

MODEL_ID = "eos.water_phase"

#: The specific gas constant for water, kJ/(kg*K).
R = 0.461526
#: The molar mass of water, kg/mol.
MOLAR_MASS = 0.01801528

_I1 = [
    0,
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
    5,
    8,
    8,
    21,
    23,
    29,
    30,
    31,
    32,
]
_J1 = [
    -2,
    -1,
    0,
    1,
    2,
    3,
    4,
    5,
    -9,
    -7,
    -1,
    0,
    1,
    3,
    -3,
    0,
    1,
    3,
    17,
    -4,
    0,
    6,
    -5,
    -2,
    10,
    -8,
    -11,
    -6,
    -29,
    -31,
    -38,
    -39,
    -40,
    -41,
]
_N1 = [
    0.14632971213167,
    -0.84548187169114,
    -3.756360367204,
    3.3855169168385,
    -0.95791963387872,
    0.15772038513228,
    -0.016616417199501,
    8.1214629983568e-4,
    2.8319080123804e-4,
    -6.0706301565874e-4,
    -0.018990068218419,
    -0.032529748770505,
    -0.021841717175414,
    -5.283835796993e-5,
    -4.7184321073267e-4,
    -3.0001780793026e-4,
    4.7661393906987e-5,
    -4.4141845330846e-6,
    -7.2694996297594e-16,
    -3.1679644845054e-5,
    -2.8270797985312e-6,
    -8.5205128120103e-10,
    -2.2425281908e-6,
    -6.5171222895601e-7,
    -1.4341729937924e-13,
    -4.0516996860117e-7,
    -1.2734301741641e-9,
    -1.7424871230634e-10,
    -6.8762131295531e-19,
    1.4478307828521e-20,
    2.6335781662795e-23,
    -1.1947622640071e-23,
    1.8228094581404e-24,
    -9.3537087292458e-26,
]
_J0 = [0, 1, -5, -4, -3, -2, -1, 2, 3]
_N0 = [
    -9.6927686500217,
    10.086655968018,
    -0.005608791128302,
    0.071452738081455,
    -0.40710498223928,
    1.4240819171444,
    -4.383951131945,
    -0.28408632460772,
    0.021268463753307,
]
_IR = [
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
    3,
    3,
    3,
    3,
    3,
    4,
    4,
    4,
    5,
    6,
    6,
    6,
    7,
    7,
    7,
    8,
    8,
    9,
    10,
    10,
    10,
    16,
    16,
    18,
    20,
    20,
    20,
    21,
    22,
    23,
    24,
    24,
    24,
]
_JR = [
    0,
    1,
    2,
    3,
    6,
    1,
    2,
    4,
    7,
    36,
    0,
    1,
    3,
    6,
    35,
    1,
    2,
    3,
    7,
    3,
    16,
    35,
    0,
    11,
    25,
    8,
    36,
    13,
    4,
    10,
    14,
    29,
    50,
    57,
    20,
    35,
    48,
    21,
    53,
    39,
    26,
    40,
    58,
]
_NR = [
    -1.7731742473213e-3,
    -0.017834862292358,
    -0.045996013696365,
    -0.057581259083432,
    -0.05032527872793,
    -3.3032641670203e-5,
    -1.8948987516315e-4,
    -3.9392777243355e-3,
    -0.043797295650573,
    -2.6674547914087e-5,
    2.0481737692309e-8,
    4.3870667284435e-7,
    -3.227767723857e-5,
    -1.5033924542148e-3,
    -0.040668253562649,
    -7.8847309559367e-10,
    1.2790717852285e-8,
    4.8225372718507e-7,
    2.2922076337661e-6,
    -1.6714766451061e-11,
    -2.1171472321355e-3,
    -23.895741934104,
    -5.905956432427e-18,
    -1.2621808899101e-6,
    -0.038946842435739,
    1.1256211360459e-11,
    -8.2311340897998,
    1.9809712802088e-8,
    1.0406965210174e-19,
    -1.0234747095929e-13,
    -1.0018179379511e-9,
    -8.0882908646985e-11,
    0.10693031879409,
    -0.33662250574171,
    8.9185845355421e-25,
    3.0629316876232e-13,
    -4.2002467698208e-6,
    -5.9056029685639e-26,
    3.7826947613457e-6,
    -1.2768608934681e-15,
    7.3087610595061e-29,
    5.5414715350778e-17,
    -9.436970724121e-7,
]


def _gamma1(pi: float, tau: float) -> dict[str, float]:
    g = gp = gt = gpp = gpt = gtt = 0.0
    for i in range(34):
        p1 = (7.1 - pi) ** _I1[i]
        t1 = (tau - 1.222) ** _J1[i]
        g += _N1[i] * p1 * t1
        gp += -_N1[i] * _I1[i] * (7.1 - pi) ** (_I1[i] - 1) * t1
        gt += _N1[i] * p1 * _J1[i] * (tau - 1.222) ** (_J1[i] - 1)
        gpp += _N1[i] * _I1[i] * (_I1[i] - 1) * (7.1 - pi) ** (_I1[i] - 2) * t1
        gpt += (
            -_N1[i] * _I1[i] * (7.1 - pi) ** (_I1[i] - 1) * _J1[i] * (tau - 1.222) ** (_J1[i] - 1)
        )
        gtt += _N1[i] * p1 * _J1[i] * (_J1[i] - 1) * (tau - 1.222) ** (_J1[i] - 2)
    return {"g": g, "gp": gp, "gt": gt, "gpp": gpp, "gpt": gpt, "gtt": gtt}


def _gamma0(pi: float, tau: float) -> dict[str, float]:
    g = math.log(pi)
    gt = gtt = 0.0
    for i in range(9):
        g += _N0[i] * tau ** _J0[i]
        gt += _N0[i] * _J0[i] * tau ** (_J0[i] - 1)
        gtt += _N0[i] * _J0[i] * (_J0[i] - 1) * tau ** (_J0[i] - 2)
    return {"g": g, "gt": gt, "gtt": gtt}


def _gammar(pi: float, tau: float) -> dict[str, float]:
    g = gp = gt = gpp = gpt = gtt = 0.0
    for i in range(43):
        g += _NR[i] * pi ** _IR[i] * (tau - 0.5) ** _JR[i]
        gp += _NR[i] * _IR[i] * pi ** (_IR[i] - 1) * (tau - 0.5) ** _JR[i]
        gt += _NR[i] * pi ** _IR[i] * _JR[i] * (tau - 0.5) ** (_JR[i] - 1)
        gpp += _NR[i] * _IR[i] * (_IR[i] - 1) * pi ** (_IR[i] - 2) * (tau - 0.5) ** _JR[i]
        gpt += _NR[i] * _IR[i] * pi ** (_IR[i] - 1) * _JR[i] * (tau - 0.5) ** (_JR[i] - 1)
        gtt += _NR[i] * pi ** _IR[i] * _JR[i] * (_JR[i] - 1) * (tau - 0.5) ** (_JR[i] - 2)
    return {"g": g, "gp": gp, "gt": gt, "gpp": gpp, "gpt": gpt, "gtt": gtt}


def _saturation_temperature(p: float) -> float:
    beta = math.pow(p, 0.25)
    e = beta * beta - 17.073846940092 * beta + 14.91510861353
    f = 1167.0521452767 * beta * beta + 12020.82470247 * beta - 4823.2657361591
    g = -724213.16703206 * beta * beta - 3232555.0322333 * beta + 405113.40542057
    d = 2.0 * g / (-f - math.sqrt(f * f - 4.0 * e * g))
    s = 650.17534844798 + d
    return (s - math.sqrt(s * s - 4.0 * (-0.23855557567849 + 650.17534844798 * d))) / 2.0


def _properties(p_mpa: float, t: float) -> dict[str, float]:
    if t <= _saturation_temperature(p_mpa):
        pi = p_mpa / 16.53
        tau = 1386.0 / t
        d = _gamma1(pi, tau)
    else:
        pi = p_mpa
        tau = 540.0 / t
        d0 = _gamma0(pi, tau)
        dr = _gammar(pi, tau)
        d = {
            "g": d0["g"] + dr["g"],
            "gp": 1.0 / pi + dr["gp"],
            "gt": d0["gt"] + dr["gt"],
            "gpp": -1.0 / (pi * pi) + dr["gpp"],
            "gpt": dr["gpt"],
            "gtt": d0["gtt"] + dr["gtt"],
        }

    z = pi * d["gp"]
    h = R * t * tau * d["gt"]
    s = R * (tau * d["gt"] - d["g"])
    cp = -R * tau * tau * d["gtt"]
    cv = cp + R * (d["gp"] - tau * d["gpt"]) ** 2 / d["gpp"]
    g = R * t * d["g"]
    u = h - R * t * z

    return {"z": z, "u": u, "h": h, "s": s, "cv": cv, "cp": cp, "g": g}


def water_phase(T: Q, P: Q) -> WaterPhaseResult:
    """The IAPWS-IF97 water phase state at a temperature and pressure.

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

    props = _properties(p_si / 1.0e6, t_si)

    molar = MOLAR_MASS * 1.0e3

    return WaterPhaseResult(
        z_factor=props["z"],
        u=from_si(props["u"] * molar, "J/mol"),
        h=from_si(props["h"] * molar, "J/mol"),
        s=from_si(props["s"] * molar, "J/(mol*K)"),
        cv=from_si(props["cv"] * molar, "J/(mol*K)"),
        cp=from_si(props["cp"] * molar, "J/(mol*K)"),
        g=from_si(props["g"] * molar, "J/mol"),
        warnings=tuple(warnings),
    )
