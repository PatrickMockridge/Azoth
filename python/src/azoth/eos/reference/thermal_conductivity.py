"""``eos.thermal_conductivity`` - the liquid thermal conductivity from the Pedersen (PFCT)
corresponding-states correlation.

Spec: ``specs/models/eos/thermal_conductivity.toml``. A methane SRK reference flash, the
Tc/Pc/molar-mass mixing rule, the dilute-gas methane viscosity, the ideal-gas heat
capacity, the reference-conductivity correlation and the corresponding-states scaling -
the same arithmetic as the Rust kernel, so the two implementations can be compared.
"""

from __future__ import annotations

import math

from azoth.core.errors import InvalidInputError, PropertyUnavailableError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import ThermalConductivityResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning
from azoth.eos.mixture import Mixture
from azoth.eos.reference.molar_enthalpy_entropy import IdealGasModel
from azoth.eos.reference.srk_alpha_ab import srk_alpha_ab
from azoth.eos.reference.srk_departure import srk_departure
from azoth.eos.reference.srk_kappa import srk_kappa
from azoth.eos.reference.srk_z_factor import srk_z_factor

MODEL_ID = "eos.thermal_conductivity"

_METHANE_TC = 190.56
_METHANE_PC = 4.599e6
_METHANE_M = 0.016043
_METHANE_OMEGA = 0.0115
_NEQSIM_GAS_CONSTANT = 8.3144621
_ATM = 101325.0
_CRIT_MOL_DENS = 10.1521197

_METHANE_CP = (37.978352, -0.07461815, 0.000301881, -2.83e-7, 9.070574e-11)

_COND_GV = (
    -2.147621e5,
    2.190461e5,
    -8.618097e4,
    1.496099e4,
    -4.730660e2,
    -2.331178e2,
    3.778439e1,
    -2.320481,
    5.311764e-2,
)
_COND_REF_A = -0.25276292
_COND_REF_B = 0.33432859
_COND_REF_C = 1.12
_COND_REF_F = 168.0
_COND_REF_J = (
    -7.04036339907,
    12.319512908,
    -8.8525979933e2,
    72.835897919,
    0.74421462902,
    -2.9706914540,
    2.2209758501e3,
)
_COND_REF_K = (-8.55109, 12.5539, -1020.85, 238.394, 1.31563, -72.5759, 1411.6)

_VISC_GV = (
    -2.090975e5,
    2.647269e5,
    -1.472818e5,
    4.716740e4,
    -9.491872e3,
    1.219979e3,
    -9.627993e1,
    4.274152,
    -8.141531e-2,
)
_VISC_REF_J = (
    -1.035060586e1,
    1.7571599671e1,
    -3.0193918656e3,
    1.8873011594e2,
    4.2903609488e-2,
    1.4529023444e2,
    6.1276818706e3,
)


def _methane_srk_density(t: float, p_pa: float) -> float:
    kappa = srk_kappa(_METHANE_OMEGA).kappa
    tr = t / _METHANE_TC
    pr = p_pa / _METHANE_PC
    ab = srk_alpha_ab(kappa, tr, pr)
    z = srk_z_factor(ab.a_reduced, ab.b_reduced)
    phi_min = srk_departure(ab.a_reduced, ab.b_reduced, z.z_min, kappa, tr).ln_phi
    phi_max = srk_departure(ab.a_reduced, ab.b_reduced, z.z_max, kappa, tr).ln_phi
    z_stable = z.z_min if phi_min <= phi_max else z.z_max
    return _METHANE_M * p_pa / (z_stable * _NEQSIM_GAS_CONSTANT * t)


def _methane_cp0(t: float) -> float:
    cp = _METHANE_CP
    return cp[0] + cp[1] * t + cp[2] * t**2 + cp[3] * t**3 + cp[4] * t**4


def _polynomial(coef: tuple[float, ...], t: float) -> float:
    exponents = (-1.0, -2.0 / 3.0, -1.0 / 3.0, 0.0, 1.0 / 3.0, 2.0 / 3.0, 1.0, 4.0 / 3.0, 5.0 / 3.0)
    total = 0.0
    for c, e in zip(coef, exponents, strict=True):
        total += c * t**e
    return total


def _dilute_gas_viscosity(t: float) -> float:
    mol_dens = _ATM / _NEQSIM_GAS_CONSTANT / t / 1e3
    red = (mol_dens - 10.15) / 10.15
    visc_ref_o = _polynomial(_VISC_GV, t)
    temp_1 = mol_dens**0.1 * (_VISC_REF_J[1] + _VISC_REF_J[2] / t**1.5)
    temp_2 = red * mol_dens**0.5 * (_VISC_REF_J[4] + _VISC_REF_J[5] / t + _VISC_REF_J[6] / t**2)
    temp_3 = math.exp(temp_1 + temp_2)
    htan = math.tanh(t - 90.69)
    visc_ref_2 = (htan + 1.0) / 2.0 * math.exp(_VISC_REF_J[0] + _VISC_REF_J[3] / t) * (temp_3 - 1.0)
    if math.isnan(visc_ref_2):
        visc_ref_2 = 0.0
    return (visc_ref_o + visc_ref_2) / 1.0e7


def _reference_conductivity(t: float, p_pa: float) -> float:
    rho = _methane_srk_density(t, p_pa)
    mol_dens_molar = rho / _METHANE_M * 1e-3
    red = (mol_dens_molar - 10.15) / 10.15
    mol_dens = rho * 1e-3

    visc_ref_o = _polynomial(_COND_GV, t)
    visc_ref_1 = (
        _COND_REF_A + _COND_REF_B * (_COND_REF_C - math.log(t / _COND_REF_F)) ** 2
    ) * mol_dens

    temp_1 = mol_dens**0.1 * (_COND_REF_J[1] + _COND_REF_J[2] / t**1.5)
    temp_2 = red * mol_dens**0.5 * (_COND_REF_J[4] + _COND_REF_J[5] / t + _COND_REF_J[6] / t**2)
    temp_3 = math.exp(temp_1 + temp_2)
    htan = math.tanh(t - 90.69)
    visc_ref_2 = (htan + 1.0) / 2.0 * math.exp(_COND_REF_J[0] + _COND_REF_J[3] / t) * (temp_3 - 1.0)
    if math.isnan(visc_ref_2):
        visc_ref_2 = 0.0

    temp_4 = mol_dens**0.1 * (_COND_REF_K[1] + _COND_REF_K[2] / t**1.5)
    temp_5 = red * mol_dens**0.5 * (_COND_REF_K[4] + _COND_REF_K[5] / t + _COND_REF_K[6] / t**2)
    temp_6 = math.exp(temp_4 + temp_5)
    visc_ref_3 = (1.0 - htan) / 2.0 * math.exp(_COND_REF_K[0] + _COND_REF_K[3] / t) * (temp_6 - 1.0)
    if math.isnan(visc_ref_3):
        visc_ref_3 = 0.0

    ref_cond = (visc_ref_o + visc_ref_1 + visc_ref_2 + visc_ref_3) * 1e-3
    if t > 400.0:
        ref_cond *= max(1.0 + (-9.0e-4) * (t - 400.0), 0.70)
    return ref_cond


def thermal_conductivity(
    mixture: Mixture, ideal_gas: IdealGasModel, T: Q, P: Q, z: list[float]
) -> ThermalConductivityResult:
    """The liquid thermal conductivity of a mixture, from the Pedersen (PFCT) correlation.

    Args:
        mixture: the components and their interaction parameters.
        ideal_gas: the heat-capacity coefficients, one vector per component.
        T: absolute temperature of the state.
        P: absolute pressure of the state.
        z: the mixture's mole fractions, checked rather than renormalised.

    Returns:
        The liquid thermal conductivity, in ``W/(m*K)``.

    Raises:
        InvalidInputError: if the composition is the wrong length, has a negative
            entry, or does not sum to one.
        PropertyUnavailableError: if a component carries no molar mass.
        OutOfRangeError: if ``T`` or ``P`` is not positive.
    """
    from azoth import _models_gen

    spec = _models_gen.model(MODEL_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    t_si = input_to_si(spec, "T", T)
    p_si = input_to_si(spec, "P", P)
    apply_checks(checks.on_input, {"T": t_si, "P": p_si}.get, warnings)

    n = len(mixture)
    if len(z) != n:
        raise InvalidInputError(
            "z", f"a mixture of {n} components needs {n} mole fractions, got {len(z)}"
        )
    if any(v < 0.0 for v in z) or abs(sum(z) - 1.0) > 1e-9:
        raise InvalidInputError("z", "mole fractions must be non-negative and sum to one")

    tcs = [c.Tc.to_base_units().magnitude for c in mixture.components]
    pcs = [c.Pc.to_base_units().magnitude for c in mixture.components]
    molar_mass = []
    for c in mixture.components:
        if c.molar_mass is None:
            raise PropertyUnavailableError(
                "component", "molar mass", "a card-added component needs its own molar mass"
            )
        molar_mass.append(c.molar_mass.to_base_units().magnitude)

    temp_tc1 = temp_tc2 = temp_pc1 = temp_pc2 = mw_temp = mm_temp = 0.0
    for i in range(n):
        for j in range(n):
            temp_var = (
                z[i]
                * z[j]
                * ((tcs[i] / pcs[i]) ** (1.0 / 3.0) + (tcs[j] / pcs[j]) ** (1.0 / 3.0)) ** 3
            )
            temp_tc1 += temp_var * math.sqrt(tcs[i] * tcs[j])
            temp_tc2 += temp_var
            temp_pc1 += temp_var * math.sqrt(tcs[i] * tcs[j])
            temp_pc2 += temp_var
        mw_temp += z[i] * molar_mass[i] ** 2
        mm_temp += z[i] * molar_mass[i]

    pc_mix = 8.0 * temp_pc1 / (temp_pc2 * temp_pc2)
    tc_mix = temp_tc1 / temp_tc2
    m_mix = (mm_temp + 1.304e-4 * ((mw_temp / mm_temp) ** 2.303 - mm_temp**2.303)) * 1e3

    t_o_ref = t_si * _METHANE_TC / tc_mix
    p_o_ref = p_si * _METHANE_PC / pc_mix
    rho = _methane_srk_density(t_o_ref, p_o_ref)
    red_dens = (rho / _METHANE_M * 1e-3) / _CRIT_MOL_DENS

    def alpha(m: float) -> float:
        return 1.0 + 6.004e-4 * math.pow(red_dens, 2.043) * math.pow(m * 1e3, 1.086)

    alfa_mix = 0.0
    for i in range(n):
        for j in range(n):
            alfa_mix += z[i] * z[j] * math.sqrt(alpha(molar_mass[i]) * alpha(molar_mass[j]))
    alfa0 = 1.0 + 6.004e-4 * red_dens**2.043 * (_METHANE_M * 1e3) ** 1.086
    if alfa_mix < 1e-10:
        return ThermalConductivityResult(k=from_si(0.0, "W/(m*K)"), warnings=tuple(warnings))

    t0 = t_o_ref * alfa0 / alfa_mix
    p0 = p_o_ref * alfa0 / alfa_mix

    nstar_ref = _dilute_gas_viscosity(t0)
    f_func = 1.0 + 0.053432 * red_dens - 0.030182 * red_dens**2 - 0.029725 * red_dens**3
    cond_int_ref = (
        1.18653 * nstar_ref * (_methane_cp0(t0) - 2.5 * _NEQSIM_GAS_CONSTANT) * f_func / _METHANE_M
    )

    red_dens_2 = _ATM / _NEQSIM_GAS_CONSTANT / t_si / 1e3 / 10.15

    def alpha_2(m: float) -> float:
        return 1.0 + 6.004e-4 * math.pow(red_dens_2, 2.043) * math.pow(m * 1e3, 1.086)

    alfa_mix_2 = 0.0
    for i in range(n):
        for j in range(n):
            alfa_mix_2 += z[i] * z[j] * math.sqrt(alpha_2(molar_mass[i]) * alpha_2(molar_mass[j]))
    alfa0_2 = 1.0 + 6.004e-4 * red_dens_2**2.043 * (_METHANE_M * 1e3) ** 1.086
    t0_b = t_si * _METHANE_TC / tc_mix * alfa_mix_2 / alfa0_2
    nstar_mix = (
        _dilute_gas_viscosity(t0_b)
        * math.pow(tc_mix / _METHANE_TC, -1.0 / 6.0)
        * math.pow(pc_mix / _METHANE_PC, 2.0 / 3.0)
        * math.pow(m_mix / (_METHANE_M * 1e3), 0.5)
        * alfa_mix_2
        / alfa0_2
    )

    cp_id_mix = 0.0
    for i in range(n):
        cp = (
            ideal_gas.cp_a[i],
            ideal_gas.cp_b[i],
            ideal_gas.cp_c[i],
            ideal_gas.cp_d[i],
            ideal_gas.cp_e[i],
        )
        cp_id_mix += z[i] * (
            cp[0] + cp[1] * t_si + cp[2] * t_si**2 + cp[3] * t_si**3 + cp[4] * t_si**4
        )
    cond_int_mix = (
        1.18653 * nstar_mix * (cp_id_mix - 2.5 * _NEQSIM_GAS_CONSTANT) * f_func / (m_mix / 1e3)
    )

    ref_conductivity = _reference_conductivity(t0, p0)
    conductivity = (
        math.pow(tc_mix / _METHANE_TC, -1.0 / 6.0)
        * math.pow(pc_mix / _METHANE_PC, 2.0 / 3.0)
        * math.pow(m_mix / (_METHANE_M * 1e3), -0.5)
        * alfa_mix
        / alfa0
        * (ref_conductivity - cond_int_ref)
        + cond_int_mix
    )

    return ThermalConductivityResult(k=from_si(conductivity, "W/(m*K)"), warnings=tuple(warnings))
