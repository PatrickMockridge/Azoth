"""``eos.viscosity`` - the liquid viscosity from the Pedersen (PFCT) heavy-oil
corresponding-states correlation.

Spec: ``specs/models/eos/viscosity.toml``. A methane SRK reference flash, a
Tc/Pc/molar-mass mixing rule, the reference-viscosity correlation and the
corresponding-states scaling - the same arithmetic as the Rust kernel, so the two
implementations can be compared case by case.
"""

from __future__ import annotations

import math

from azoth.core.errors import InvalidInputError, PropertyUnavailableError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import ViscosityResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning
from azoth.eos.mixture import Mixture
from azoth.eos.reference.srk_alpha_ab import srk_alpha_ab
from azoth.eos.reference.srk_departure import srk_departure
from azoth.eos.reference.srk_kappa import srk_kappa
from azoth.eos.reference.srk_z_factor import srk_z_factor

MODEL_ID = "eos.viscosity"

#: Methane's critical constants, the reference component's.
_METHANE_TC = 190.56
_METHANE_PC = 4.599e6
_METHANE_M = 0.016043
_METHANE_OMEGA = 0.0115

#: NeqSim's gas constant, the one the reference density uses.
_NEQSIM_GAS_CONSTANT = 8.3144621

_GVCOEF = (
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
_VIS_REF_A = 1.696985927
_VIS_REF_B = -0.133372346
_VIS_REF_C = 1.4
_VIS_REF_F = 168.0
_VISC_REF_J = (
    -1.035060586e1,
    1.7571599671e1,
    -3.0193918656e3,
    1.8873011594e2,
    4.2903609488e-2,
    1.4529023444e2,
    6.1276818706e3,
)
_VISC_REF_K = (-9.74602, 18.0834, -4126.66, 44.6055, 0.976544, 81.8134, 15649.9)
_CRIT_MOL_DENS = 10.15


def _methane_srk_density(t: float, p_pa: float) -> float:
    """Methane's stable-phase SRK density in kg/m^3, the lower-fugacity root."""
    kappa = srk_kappa(_METHANE_OMEGA).kappa
    tr = t / _METHANE_TC
    pr = p_pa / _METHANE_PC
    ab = srk_alpha_ab(kappa, tr, pr)
    z = srk_z_factor(ab.a_reduced, ab.b_reduced)
    phi_min = srk_departure(ab.a_reduced, ab.b_reduced, z.z_min, kappa, tr).ln_phi
    phi_max = srk_departure(ab.a_reduced, ab.b_reduced, z.z_max, kappa, tr).ln_phi
    z_stable = z.z_min if phi_min <= phi_max else z.z_max
    return _METHANE_M * p_pa / (z_stable * _NEQSIM_GAS_CONSTANT * t)


def _get_ref_viscosity(t0: float, p0: float) -> float:
    """The methane reference viscosity, `getRefComponentViscosity`."""
    rho = _methane_srk_density(t0, p0)
    mol_dens_molar = rho / _METHANE_M * 1e-3
    red_mol_dens = (mol_dens_molar - _CRIT_MOL_DENS) / _CRIT_MOL_DENS
    mol_dens = rho * 1e-3

    exponents = (-1.0, -2.0 / 3.0, -1.0 / 3.0, 0.0, 1.0 / 3.0, 2.0 / 3.0, 1.0, 4.0 / 3.0, 5.0 / 3.0)
    visc_ref_o = 0.0
    for coef, exp in zip(_GVCOEF, exponents, strict=True):
        visc_ref_o += coef * t0**exp
    visc_ref_1 = (
        _VIS_REF_A + _VIS_REF_B * (_VIS_REF_C - math.log(t0 / _VIS_REF_F)) ** 2
    ) * mol_dens

    temp_1 = mol_dens**0.1 * (_VISC_REF_J[1] + _VISC_REF_J[2] / t0**1.5)
    temp_2 = (
        red_mol_dens
        * mol_dens**0.5
        * (_VISC_REF_J[4] + _VISC_REF_J[5] / t0 + _VISC_REF_J[6] / t0**2)
    )
    temp_3 = math.exp(temp_1 + temp_2)

    htan = math.tanh(t0 - 90.69)
    visc_ref_2 = (
        (htan + 1.0) / 2.0 * math.exp(_VISC_REF_J[0] + _VISC_REF_J[3] / t0) * (temp_3 - 1.0)
    )
    if math.isnan(visc_ref_2):
        visc_ref_2 = 0.0

    temp_4 = mol_dens**0.1 * (_VISC_REF_K[1] + _VISC_REF_K[2] / t0**1.5)
    temp_5 = (
        red_mol_dens
        * mol_dens**0.5
        * (_VISC_REF_K[4] + _VISC_REF_K[5] / t0 + _VISC_REF_K[6] / t0**2)
    )
    temp_6 = math.exp(temp_4 + temp_5)
    visc_ref_3 = (
        (1.0 - htan) / 2.0 * math.exp(_VISC_REF_K[0] + _VISC_REF_K[3] / t0) * (temp_6 - 1.0)
    )
    if math.isnan(visc_ref_3):
        visc_ref_3 = 0.0

    return (visc_ref_o + visc_ref_1 + visc_ref_2 + visc_ref_3) / 1.0e7


def _csp(
    ref_visc: float, tc_mix: float, pc_mix: float, m_mix: float, alfa_mix: float, alfa0: float
) -> float:
    """The corresponding-states scaling, four correction factors at their default 1.0."""
    return (
        ref_visc
        * math.pow(tc_mix / _METHANE_TC, -1.0 / 6.0)
        * math.pow(pc_mix / _METHANE_PC, 2.0 / 3.0)
        * math.pow(m_mix / (_METHANE_M * 1e3), 0.5)
        * (alfa_mix / alfa0)
    )


def viscosity(mixture: Mixture, T: Q, P: Q, z: list[float]) -> ViscosityResult:
    """The liquid viscosity of a mixture, from the Pedersen (PFCT) correlation.

    Args:
        mixture: the components and their interaction parameters.
        T: absolute temperature of the state.
        P: absolute pressure of the state.
        z: the mixture's mole fractions, checked rather than renormalised.

    Returns:
        The liquid dynamic viscosity, in ``Pa*s``.

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

    if temp_tc2 < 1e-10:
        return ViscosityResult(mu=from_si(0.0, "Pa*s"), warnings=tuple(warnings))

    pc_mix = 8.0 * temp_pc1 / (temp_pc2 * temp_pc2)
    tc_mix = temp_tc1 / temp_tc2
    m_mix = (mm_temp + 1.304e-4 * ((mw_temp / mm_temp) ** 2.303 - mm_temp**2.303)) * 1e3

    t_scaled = t_si * _METHANE_TC / tc_mix
    p_scaled = p_si * _METHANE_PC / pc_mix
    rho = _methane_srk_density(t_scaled, p_scaled)
    red_dens = (rho / _METHANE_M * 1e-3) / _CRIT_MOL_DENS
    alfa_mix = 1.0 + 7.378e-3 * red_dens**1.847 * m_mix**0.5173
    alfa0 = 1.0 + 7.378e-3 * red_dens**1.847 * (_METHANE_M * 1e3) ** 0.5173
    t0 = t_scaled * alfa0 / alfa_mix
    p0 = p_scaled * alfa0 / alfa_mix

    if t0 < 75.0:
        mol_m = (
            mm_temp * 1e3
            if mw_temp / mm_temp / mm_temp <= 1.5
            else mm_temp * (mw_temp / mm_temp / (1.5 * mm_temp)) ** 0.5 * 1e3
        )
        sign = -1.0 if t_si > 564.49 else 1.0
        termm = -0.07955 - sign * 0.01101 * mol_m - 371.8 / t_si + 6.215 * mol_m / t_si
        ho_viscosity = 10.0**termm * 1e-3
        ho_viscosity += ho_viscosity * 0.008 * (p_si / 1e5 - 1.0)
        if t0 < 65.0:
            mu = ho_viscosity
        else:
            lo_viscosity = _csp(_get_ref_viscosity(t0, p0), tc_mix, pc_mix, m_mix, alfa_mix, alfa0)
            mu = lo_viscosity * (1.0 - (75.0 - t0) / 10.0) + ho_viscosity * (75.0 - t0) / 10.0
    else:
        mu = _csp(_get_ref_viscosity(t0, p0), tc_mix, pc_mix, m_mix, alfa_mix, alfa0)

    return ViscosityResult(mu=from_si(mu, "Pa*s"), warnings=tuple(warnings))
