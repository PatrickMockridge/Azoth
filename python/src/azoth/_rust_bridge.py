"""The Rust backend, presented with the same interface as the reference.

Each function here takes and returns exactly what
:mod:`azoth.hydraulics.reference` does, so the dispatcher can swap one for the
other and a caller cannot tell which answered. Two jobs sit between the two
representations:

* **Units in.** Quantities are converted to SI magnitudes before crossing into
  Rust. That conversion lives here, once, rather than in the extension - one
  conversion site used by the reference and the Rust path alike cannot disagree
  with itself about what a number is in.
* **Results out.** The extension's transport objects become the same frozen
  dataclasses the reference returns, so ``result.dp`` is a pint quantity either
  way and ``isinstance`` checks behave identically.

If this file is missing or out of step with the extension, the dispatcher raises
rather than falling back - a silent fallback would mean the cross-language
guarantee was being verified by nothing.
"""

from __future__ import annotations

from collections.abc import Callable, Sequence
from typing import Any

from azoth import _core, _models_gen
from azoth._registry_gen import spec as _spec_for
from azoth.core.errors import PropertyUnavailableError
from azoth.core.result import (
    AmmoniaPhaseResult,
    AntoineVaporPressureResult,
    ArgonSolidPhaseResult,
    BubblePressureResult,
    BubbleTemperatureResult,
    BwrsPhaseResult,
    CapillaryDewPointResult,
    ChokedFlowAreaResult,
    ChungConductivityResult,
    ChungViscosityResult,
    Co2PhaseResult,
    Co2WaterDiffusivityResult,
    ColebrookResult,
    ConductionPlaneWallResult,
    ControlValveCvResult,
    CostaldMolarVolumeResult,
    CriticalPointResult,
    DarcyWeisbachResult,
    DewPressureResult,
    DewTemperatureResult,
    EosCgPhaseResult,
    FlowRegime,
    GeNrtlFlashResult,
    GeNrtlPhaseResult,
    Gerg2008PhaseResult,
    GeUnifacPhaseResult,
    GeUniquacPhaseResult,
    GeVanLaarAcidPhaseResult,
    GeWilsonPhaseResult,
    HaalandResult,
    HaydukMinhasDiffusivityResult,
    HeatOfVaporizationResult,
    HeliumPhaseResult,
    HydrogenPhaseResult,
    IdealGasCpResult,
    KComponent,
    KFactorsResult,
    LiquidHeatCapacityResult,
    MasonSaxenaConductivityResult,
    Matcop5PrumrAlphaResult,
    MatcopAlphaResult,
    MatcopPrAlphaResult,
    MatcopPrumrAlphaResult,
    MatcopPrumrNewAlphaResult,
    MolarEnthalpyEntropyResult,
    MollerupAlphaResult,
    NitricSulfuricAcidVaporPressureResult,
    NrtlActivityCoefficientsResult,
    OrificeFlowResult,
    ParachorSurfaceTensionResult,
    ParahydrogenSolidPhaseResult,
    PhFlashResult,
    Pr78KappaResult,
    PrAlphaAbResult,
    PrDaneshAlphaResult,
    PrDelft1998AlphaResult,
    PrDepartureResult,
    PrGassem2001AlphaResult,
    PrKappaResult,
    PrLeeKeslerAlphaResult,
    PrMassDensityResult,
    PrMolarVolumeResult,
    PrPenelouxShiftResult,
    PrsvKappaResult,
    PrZFactorResult,
    PsFlashResult,
    PtFlashResult,
    PtPhaseEnvelopeResult,
    PuFlashResult,
    PumpPowerResult,
    PureSaturationResult,
    PvfFlashResult,
    PvFlashResult,
    PvRefluxFlashResult,
    RachfordRiceBinaryResult,
    RachfordRiceResult,
    RackettMolarVolumeResult,
    ReynoldsNumberResult,
    RkAlphaAbResult,
    RkDepartureResult,
    RootStructure,
    SchwartzentruberAlphaResult,
    SiddiqiLucasDiffusivityResult,
    SoreideWhitsonAlphaResult,
    SrkAlphaAbResult,
    SrkDepartureResult,
    SrkKappaResult,
    SrkPenelouxShiftResult,
    SrkZFactorResult,
    StabilityTestResult,
    SwameeJainResult,
    ThermalConductivityResult,
    ThFlashResult,
    TpMultiflashResult,
    TsFlashResult,
    TuFlashResult,
    TvFlashResult,
    TvFractionFlashResult,
    TwucoonAlphaResult,
    TwucoonParamAlphaResult,
    TwucoonStatoilAlphaResult,
    TwuKappaResult,
    TynCalusDiffusivityResult,
    UmrprAlphaResult,
    UnifacActivityCoefficientsResult,
    UnifacPsrkActivityCoefficientsResult,
    UnifacUmrpruActivityCoefficientsResult,
    UniquacActivityCoefficientsResult,
    VanLaarAcidActivityCoefficientsResult,
    Vdw1fMixBinaryResult,
    VhFlashResult,
    ViscosityResult,
    VsFlashResult,
    VuFlashResult,
    VuFlashSingleCompResult,
    WaterPhaseResult,
    WilkeChangDiffusivityResult,
    WilkeViscosityResult,
    WilsonActivityCoefficientsResult,
)
from azoth.core.result import Phase as _Phase
from azoth.core.result import StabilityVerdict as _StabilityVerdict
from azoth.core.result import TpMultiflashSeed as _TpMultiflashSeed
from azoth.core.units import Q, from_si, input_to_si, to_si
from azoth.core.warnings import Warning, WarningCode


def _warnings(raw: Sequence[_core.Warning]) -> tuple[Warning, ...]:
    """Convert transported warnings to the real ones.

    Codes become the enum rather than staying strings, so a caller can compare
    ``warning.code == WarningCode.TRANSITIONAL_FLOW`` regardless of which backend
    produced it. A code the Python side does not know is impossible - the two
    sets are asserted equal by a test - but the ValueError would surface here if
    it ever happened, which is better than a string quietly failing to match.
    """
    return tuple(Warning(WarningCode(w.code), w.message, w.field) for w in raw)


def reynolds_number(rho: Q, v: Q, D: Q, mu: Q) -> ReynoldsNumberResult:
    """Reynolds number, computed in Rust."""
    result = _core.reynolds_number(
        to_si(rho, "kg/m**3", "rho"),
        to_si(v, "m/s", "v"),
        to_si(D, "m", "D"),
        to_si(mu, "Pa*s", "mu"),
    )
    return ReynoldsNumberResult(
        re=result.re,
        regime=FlowRegime(result.regime),
        warnings=_warnings(result.warnings),
    )


def friction_factor_colebrook(re: float, relative_roughness: float) -> ColebrookResult:
    """Colebrook friction factor, solved in Rust."""
    result = _core.friction_factor_colebrook(re, relative_roughness)
    return ColebrookResult(
        f=result.f,
        iterations=result.iterations,
        converged=result.converged,
        residual=result.residual,
        warnings=_warnings(result.warnings),
    )


def friction_factor_swamee_jain(re: float, relative_roughness: float) -> SwameeJainResult:
    """Swamee-Jain friction factor, evaluated in Rust."""
    result = _core.friction_factor_swamee_jain(re, relative_roughness)
    return SwameeJainResult(f=result.f, warnings=_warnings(result.warnings))


def friction_factor_haaland(re: float, relative_roughness: float) -> HaalandResult:
    """Haaland explicit friction factor, evaluated in Rust."""
    result = _core.friction_factor_haaland(re, relative_roughness)
    return HaalandResult(f=result.f, warnings=_warnings(result.warnings))


def crane_k_factors(fittings: Sequence[str], f_t: float) -> KFactorsResult:
    """Fitting resistance coefficients, resolved in Rust."""
    result = _core.crane_k_factors(list(fittings), f_t)
    return KFactorsResult(
        k_total=result.k_total,
        f_t=result.f_t,
        components=tuple(
            KComponent(fitting_id=c.fitting_id, n_ld=c.n_ld, k=c.k) for c in result.components
        ),
        warnings=_warnings(result.warnings),
    )


def darcy_weisbach(
    f: float,
    L: Q,
    D: Q,
    rho: Q,
    v: Q,
    mu: Q | None = None,
) -> DarcyWeisbachResult:
    """Darcy-Weisbach pressure drop, computed in Rust."""
    result = _core.darcy_weisbach(
        f,
        to_si(L, "m", "L"),
        to_si(D, "m", "D"),
        to_si(rho, "kg/m**3", "rho"),
        to_si(v, "m/s", "v"),
        None if mu is None else to_si(mu, "Pa*s", "mu"),
    )
    return DarcyWeisbachResult(
        # Rebuilt as a real pint quantity from the SI magnitude and the unit the
        # Rust side reported, so the field has the same type as the reference's.
        # `magnitude_si` is SI base, which is the `from_si` direction, not a number
        # already stated in pascals - the two coincide for pascals alone.
        dp=from_si(result.dp.magnitude_si, result.dp.unit),
        f=result.f,
        re=result.re,
        regime=None if result.regime is None else FlowRegime(result.regime),
        warnings=_warnings(result.warnings),
    )


def conduction_plane_wall(k: Q, A: Q, dT: Q, L: Q) -> ConductionPlaneWallResult:
    """Plane-wall conduction, computed in Rust."""
    spec = _spec_for("thermal.conduction_plane_wall")
    result = _core.conduction_plane_wall(
        input_to_si(spec, "k", k),
        input_to_si(spec, "A", A),
        # `interval: true` in the spec, so an absolute `Q(30, "degC")` is refused
        # here rather than converted to 303.15 K. Both backends must make the same
        # choice, and this is the only place the Rust path's units are decided.
        input_to_si(spec, "dT", dT),
        input_to_si(spec, "L", L),
    )
    return ConductionPlaneWallResult(
        # `q` is an SI base magnitude from the extension, so it is rebuilt as a real
        # pint quantity with `from_si` - the same direction the reference uses.
        q=from_si(result.q.magnitude_si, result.q.unit),
        warnings=_warnings(result.warnings),
    )


def pr_kappa(omega: float) -> PrKappaResult:
    """The Peng-Robinson attraction-parameter coefficient, computed in Rust.

    No conversion in either direction: both this input and this output are
    genuinely dimensionless, so the extension carries a bare float and there is no
    unit string for the two implementations to disagree about. The acentric factor
    crosses this boundary as the same number the callers on both sides used.
    """
    result = _core.pr_kappa(omega)
    return PrKappaResult(kappa=result.kappa, warnings=_warnings(result.warnings))


def matcop_alpha(mc1: float, mc2: float, mc3: float, Tr: float) -> MatcopAlphaResult:
    """The Mathias-Copeman alpha function, computed in Rust."""
    result = _core.matcop_alpha(mc1, mc2, mc3, Tr)
    return MatcopAlphaResult(alpha=result.alpha, warnings=_warnings(result.warnings))


def matcop_pr_alpha(
    omega: float, mc1: float, mc2: float, mc3: float, Tr: float
) -> MatcopPrAlphaResult:
    """The Mathias-Copeman alpha with a Peng-Robinson fallback, computed in Rust."""
    result = _core.matcop_pr_alpha(omega, mc1, mc2, mc3, Tr)
    return MatcopPrAlphaResult(alpha=result.alpha, warnings=_warnings(result.warnings))


def matcop_prumr_alpha(
    omega: float, mc1: float, mc2: float, mc3: float, Tr: float
) -> MatcopPrumrAlphaResult:
    """The Mathias-Copeman alpha with the UMR-PR fallback, computed in Rust."""
    result = _core.matcop_prumr_alpha(omega, mc1, mc2, mc3, Tr)
    return MatcopPrumrAlphaResult(alpha=result.alpha, warnings=_warnings(result.warnings))


def matcop_prumr_new_alpha(
    omega: float, mc1: float, mc2: float, mc3: float, mc4: float, mc5: float, Tr: float
) -> MatcopPrumrNewAlphaResult:
    """The five-parameter Mathias-Copeman alpha, UMR-PR new variant, in Rust."""
    result = _core.matcop_prumr_new_alpha(omega, mc1, mc2, mc3, mc4, mc5, Tr)
    return MatcopPrumrNewAlphaResult(alpha=result.alpha, warnings=_warnings(result.warnings))


def matcop5_prumr_alpha(
    omega: float, mc1: float, mc2: float, mc3: float, mc4: float, mc5: float, Tr: float
) -> Matcop5PrumrAlphaResult:
    """The five-parameter Mathias-Copeman alpha function, computed in Rust."""
    result = _core.matcop5_prumr_alpha(omega, mc1, mc2, mc3, mc4, mc5, Tr)
    return Matcop5PrumrAlphaResult(alpha=result.alpha, warnings=_warnings(result.warnings))


def mollerup_alpha(p1: float, p2: float, p3: float, Tr: float) -> MollerupAlphaResult:
    """The Mollerup alpha function, computed in Rust."""
    result = _core.mollerup_alpha(p1, p2, p3, Tr)
    return MollerupAlphaResult(alpha=result.alpha, warnings=_warnings(result.warnings))


def pr_danesh_alpha(omega: float, Tr: float) -> PrDaneshAlphaResult:
    """The Danesh alpha function, computed in Rust."""
    result = _core.pr_danesh_alpha(omega, Tr)
    return PrDaneshAlphaResult(alpha=result.alpha, warnings=_warnings(result.warnings))


def pr_delft1998_alpha(omega: float, Tr: float) -> PrDelft1998AlphaResult:
    """The Peng-Robinson alpha, Delft (1998), computed in Rust."""
    result = _core.pr_delft1998_alpha(omega, Tr)
    return PrDelft1998AlphaResult(alpha=result.alpha, warnings=_warnings(result.warnings))


def pr_gassem2001_alpha(omega: float, Tr: float) -> PrGassem2001AlphaResult:
    """The Gassem (2001) alpha function, computed in Rust."""
    result = _core.pr_gassem2001_alpha(omega, Tr)
    return PrGassem2001AlphaResult(alpha=result.alpha, warnings=_warnings(result.warnings))


def pr_lee_kesler_alpha(omega: float, Tr: float) -> PrLeeKeslerAlphaResult:
    """The Peng-Robinson alpha with a Soave-form m-factor, computed in Rust."""
    result = _core.pr_lee_kesler_alpha(omega, Tr)
    return PrLeeKeslerAlphaResult(alpha=result.alpha, warnings=_warnings(result.warnings))


def pr_alpha_ab(kappa: float, Tr: float, Pr: float) -> PrAlphaAbResult:
    """The Peng-Robinson alpha function and reduced parameters, computed in Rust.

    Dimensionless end to end, so like `pr_kappa` there is no conversion in either
    direction - the numbers that cross are the numbers the callers used.
    """
    result = _core.pr_alpha_ab(kappa, Tr, Pr)
    return PrAlphaAbResult(
        alpha=result.alpha,
        a_reduced=result.a_reduced,
        b_reduced=result.b_reduced,
        warnings=_warnings(result.warnings),
    )


def pr_z_factor(a_reduced: float, b_reduced: float) -> PrZFactorResult:
    """The Peng-Robinson compressibility factor, computed in Rust.

    `root_structure` crosses as the spec's string and is rebuilt into the enum, the
    same way `regime` is for `reynolds_number` - so a caller cannot tell which
    backend answered, which is the point of the adapter.
    """
    result = _core.pr_z_factor(a_reduced, b_reduced)
    return PrZFactorResult(
        z_min=result.z_min,
        z_max=result.z_max,
        root_structure=RootStructure(result.root_structure),
        iterations=result.iterations,
        converged=result.converged,
        residual=result.residual,
        warnings=_warnings(result.warnings),
    )


def prsv_kappa(omega: float, Tr: float, kappa1: float) -> PrsvKappaResult:
    """The PRSV alpha-function coefficient, computed in Rust.

    Dimensionless end to end, so nothing is converted in either direction.
    """
    result = _core.prsv_kappa(omega, Tr, kappa1)
    return PrsvKappaResult(kappa=result.kappa, warnings=_warnings(result.warnings))


def pr78_kappa(omega: float) -> Pr78KappaResult:
    """The 1978 Peng-Robinson alpha-function coefficient, computed in Rust."""
    result = _core.pr78_kappa(omega)
    return Pr78KappaResult(kappa=result.kappa, warnings=_warnings(result.warnings))


def twu_kappa(omega: float) -> TwuKappaResult:
    """Twu's alpha-function coefficient, computed in Rust."""
    result = _core.twu_kappa(omega)
    return TwuKappaResult(kappa=result.kappa, warnings=_warnings(result.warnings))


def twucoon_alpha(omega: float, Tr: float) -> TwucoonAlphaResult:
    """The Twu-Coon alpha function, computed in Rust."""
    result = _core.twucoon_alpha(omega, Tr)
    return TwucoonAlphaResult(alpha=result.alpha, warnings=_warnings(result.warnings))


def twucoon_param_alpha(a: float, b: float, c: float, Tr: float) -> TwucoonParamAlphaResult:
    """The Twu-Coon parameter alpha function, computed in Rust."""
    result = _core.twucoon_param_alpha(a, b, c, Tr)
    return TwucoonParamAlphaResult(alpha=result.alpha, warnings=_warnings(result.warnings))


def twucoon_statoil_alpha(a: float, b: float, c: float, Tr: float) -> TwucoonStatoilAlphaResult:
    """The Twu-Coon Statoil alpha function, computed in Rust."""
    result = _core.twucoon_statoil_alpha(a, b, c, Tr)
    return TwucoonStatoilAlphaResult(alpha=result.alpha, warnings=_warnings(result.warnings))


def pr_departure(
    a_reduced: float, b_reduced: float, z: float, kappa: float, Tr: float
) -> PrDepartureResult:
    """The Peng-Robinson fugacity coefficient and departures, computed in Rust.

    Five dimensionless arguments and four dimensionless outputs, so there is
    nothing to convert and no unit string to keep in step.
    """
    result = _core.pr_departure(a_reduced, b_reduced, z, kappa, Tr)
    return PrDepartureResult(
        ln_phi=result.ln_phi,
        h_dep_rt=result.h_dep_rt,
        s_dep_r=result.s_dep_r,
        cp_dep_r=result.cp_dep_r,
        warnings=_warnings(result.warnings),
    )


def srk_kappa(omega: float) -> SrkKappaResult:
    """The Soave-Redlich-Kwong attraction-parameter coefficient, computed in Rust."""
    result = _core.srk_kappa(omega)
    return SrkKappaResult(kappa=result.kappa, warnings=_warnings(result.warnings))


def srk_alpha_ab(kappa: float, Tr: float, Pr: float) -> SrkAlphaAbResult:
    """The Soave-Redlich-Kwong alpha function and reduced parameters, in Rust."""
    result = _core.srk_alpha_ab(kappa, Tr, Pr)
    return SrkAlphaAbResult(
        alpha=result.alpha,
        a_reduced=result.a_reduced,
        b_reduced=result.b_reduced,
        warnings=_warnings(result.warnings),
    )


def srk_z_factor(a_reduced: float, b_reduced: float) -> SrkZFactorResult:
    """The Soave-Redlich-Kwong compressibility factor, computed in Rust."""
    result = _core.srk_z_factor(a_reduced, b_reduced)
    return SrkZFactorResult(
        z_min=result.z_min,
        z_max=result.z_max,
        root_structure=RootStructure(result.root_structure),
        iterations=result.iterations,
        converged=result.converged,
        residual=result.residual,
        warnings=_warnings(result.warnings),
    )


def srk_departure(
    a_reduced: float, b_reduced: float, z: float, kappa: float, Tr: float
) -> SrkDepartureResult:
    """The Soave-Redlich-Kwong fugacity coefficient and departures, in Rust."""
    result = _core.srk_departure(a_reduced, b_reduced, z, kappa, Tr)
    return SrkDepartureResult(
        ln_phi=result.ln_phi,
        h_dep_rt=result.h_dep_rt,
        s_dep_r=result.s_dep_r,
        cp_dep_r=result.cp_dep_r,
        warnings=_warnings(result.warnings),
    )


def rk_alpha_ab(Tr: float, Pr: float) -> RkAlphaAbResult:
    """The Redlich-Kwong alpha function and reduced parameters, in Rust."""
    result = _core.rk_alpha_ab(Tr, Pr)
    return RkAlphaAbResult(
        alpha=result.alpha,
        a_reduced=result.a_reduced,
        b_reduced=result.b_reduced,
        warnings=_warnings(result.warnings),
    )


def rk_departure(a_reduced: float, b_reduced: float, z: float) -> RkDepartureResult:
    """The Redlich-Kwong fugacity coefficient and departures, in Rust."""
    result = _core.rk_departure(a_reduced, b_reduced, z)
    return RkDepartureResult(
        ln_phi=result.ln_phi,
        h_dep_rt=result.h_dep_rt,
        s_dep_r=result.s_dep_r,
        cp_dep_r=result.cp_dep_r,
        warnings=_warnings(result.warnings),
    )


def vdw1f_mix_binary(
    z1: float, a1: float, a2: float, b1: float, b2: float, k12: float
) -> Vdw1fMixBinaryResult:
    """Van der Waals one-fluid mixing, computed in Rust. Dimensionless throughout."""
    result = _core.vdw1f_mix_binary(z1, a1, a2, b1, b2, k12)
    return Vdw1fMixBinaryResult(
        a_mix=result.a_mix, b_mix=result.b_mix, warnings=_warnings(result.warnings)
    )


def rachford_rice_binary(z1: float, K1: float, K2: float) -> RachfordRiceBinaryResult:
    """The binary Rachford-Rice vapour fraction, computed in Rust."""
    result = _core.rachford_rice_binary(z1, K1, K2)
    return RachfordRiceBinaryResult(beta=result.beta, warnings=_warnings(result.warnings))


def rachford_rice(z: Sequence[float], K: Sequence[float]) -> RachfordRiceResult:
    """The Rachford-Rice vapour fraction, computed in Rust."""
    result = _core.rachford_rice(list(z), list(K))
    return RachfordRiceResult(beta=result.beta, warnings=_warnings(result.warnings))


def pr_molar_volume(z: float, T: Q, P: Q) -> PrMolarVolumeResult:
    """Molar volume, computed in Rust.

    One of the two dimensional calcs in this namespace, so unlike its neighbours it
    converts: `T` and `P` to SI magnitudes on the way in, and the volume back to a
    real pint quantity on the way out.
    """
    spec = _spec_for("eos.pr_molar_volume")
    result = _core.pr_molar_volume(
        z,
        input_to_si(spec, "T", T),
        input_to_si(spec, "P", P),
    )
    return PrMolarVolumeResult(
        v=from_si(result.v.magnitude_si, result.v.unit),
        warnings=_warnings(result.warnings),
    )


def pure_saturation(Tc: Q, Pc: Q, omega: float, T: Q) -> PureSaturationResult:
    """The saturation pressure of a pure component, computed in Rust.

    A model rather than a calculation - its spec fixes a procedure and the Rust side
    composes the same kernels this bridge's scalar functions call, so the two
    implementations run the same search over the same arithmetic.
    """
    spec = _models_gen.model("eos.pure_saturation")
    result = _core.pure_saturation(
        # Not `input_to_si`: that reads the unit out of the spec's declaration, and the
        # spec no longer declares `Tc` or `Pc` - a component's constants come from the
        # databank by name, and a caller holding them holds quantities already. `T` is
        # still declared and still goes through the spec.
        Tc.to("K").magnitude,
        Pc.to("Pa").magnitude,
        # A plain float: `omega` is genuinely dimensionless, so it crosses as
        # the number the caller used - the same rule the other eos calcs follow.
        omega,
        input_to_si(spec, "T", T),
    )
    return PureSaturationResult(
        p_sat=from_si(result.p_sat.magnitude_si, result.p_sat.unit),
        ln_phi=result.ln_phi,
        iterations=result.iterations,
        residual=result.residual,
        warnings=_warnings(result.warnings),
    )


def pr_mass_density(M: Q, v: Q) -> PrMassDensityResult:
    """Mass density, computed in Rust.

    `M` crosses as kg/mol, which is what the spec declares - a caller supplying the
    g/mol a table quotes has already been refused by `to_si`'s dimensionality check
    or has converted themselves, and there is no third case.
    """
    spec = _spec_for("eos.pr_mass_density")
    result = _core.pr_mass_density(
        input_to_si(spec, "M", M),
        input_to_si(spec, "v", v),
    )
    return PrMassDensityResult(
        rho=from_si(result.rho.magnitude_si, result.rho.unit),
        warnings=_warnings(result.warnings),
    )


def pr_peneloux_shift(omega: float, Tc: Q, Pc: Q) -> PrPenelouxShiftResult:
    """The Peng-Robinson Peneloux volume-translation parameter, computed in Rust."""
    spec = _spec_for("eos.pr_peneloux_shift")
    result = _core.pr_peneloux_shift(
        omega,
        input_to_si(spec, "Tc", Tc),
        input_to_si(spec, "Pc", Pc),
    )
    return PrPenelouxShiftResult(
        c=from_si(result.c.magnitude_si, result.c.unit),
        warnings=_warnings(result.warnings),
    )


def srk_peneloux_shift(omega: float, Tc: Q, Pc: Q) -> SrkPenelouxShiftResult:
    """The Soave-Redlich-Kwong Peneloux volume-translation parameter, in Rust."""
    spec = _spec_for("eos.srk_peneloux_shift")
    result = _core.srk_peneloux_shift(
        omega,
        input_to_si(spec, "Tc", Tc),
        input_to_si(spec, "Pc", Pc),
    )
    return SrkPenelouxShiftResult(
        c=from_si(result.c.magnitude_si, result.c.unit),
        warnings=_warnings(result.warnings),
    )


def heat_of_vaporization(
    c0: float, c1: float, c2: float, c3: float, Tc: Q, T: Q
) -> HeatOfVaporizationResult:
    """The pure-component heat of vaporisation, computed in Rust."""
    spec = _spec_for("eos.heat_of_vaporization")
    result = _core.heat_of_vaporization(
        c0,
        c1,
        c2,
        c3,
        input_to_si(spec, "Tc", Tc),
        input_to_si(spec, "T", T),
    )
    return HeatOfVaporizationResult(
        hov=from_si(result.hov.magnitude_si, result.hov.unit),
        warnings=_warnings(result.warnings),
    )


def liquid_heat_capacity(
    c0: float, c1: float, c2: float, c3: float, c4: float, T: Q
) -> LiquidHeatCapacityResult:
    """The pure-component liquid heat capacity, computed in Rust."""
    spec = _spec_for("eos.liquid_heat_capacity")
    result = _core.liquid_heat_capacity(
        c0,
        c1,
        c2,
        c3,
        c4,
        input_to_si(spec, "T", T),
    )
    return LiquidHeatCapacityResult(
        cp=from_si(result.cp.magnitude_si, result.cp.unit),
        warnings=_warnings(result.warnings),
    )


def antoine_vapor_pressure(
    A: float, B: float, C: float, D: float, E: float, form: str, Tc: Q, Pc: Q, T: Q
) -> AntoineVaporPressureResult:
    """The pure-component vapour pressure, computed in Rust."""
    spec = _spec_for("eos.antoine_vapor_pressure")
    result = _core.antoine_vapor_pressure(
        A,
        B,
        C,
        D,
        E,
        form,
        input_to_si(spec, "Tc", Tc),
        input_to_si(spec, "Pc", Pc),
        input_to_si(spec, "T", T),
    )
    return AntoineVaporPressureResult(
        p_sat=from_si(result.p_sat.magnitude_si, result.p_sat.unit),
        warnings=_warnings(result.warnings),
    )


def rackett_molar_volume(omega: float, Tc: Q, Pc: Q, T: Q) -> RackettMolarVolumeResult:
    """The saturated liquid molar volume, computed in Rust."""
    spec = _spec_for("eos.rackett_molar_volume")
    result = _core.rackett_molar_volume(
        omega,
        input_to_si(spec, "Tc", Tc),
        input_to_si(spec, "Pc", Pc),
        input_to_si(spec, "T", T),
    )
    return RackettMolarVolumeResult(
        v=from_si(result.v.magnitude_si, result.v.unit),
        warnings=_warnings(result.warnings),
    )


def costald_molar_volume(
    omega: float, Tc: Q, Vc: Q, M: Q, rho_normal: Q, T: Q
) -> CostaldMolarVolumeResult:
    """The saturated liquid molar volume, computed in Rust."""
    spec = _spec_for("eos.costald_molar_volume")
    result = _core.costald_molar_volume(
        omega,
        input_to_si(spec, "Tc", Tc),
        input_to_si(spec, "Vc", Vc),
        input_to_si(spec, "M", M),
        input_to_si(spec, "rho_normal", rho_normal),
        input_to_si(spec, "T", T),
    )
    return CostaldMolarVolumeResult(
        v=from_si(result.v.magnitude_si, result.v.unit),
        warnings=_warnings(result.warnings),
    )


def chung_viscosity(
    omega: float, Tc: Q, Vc: Q, M: Q, dipole: float, kappa: float, T: Q, V: Q
) -> ChungViscosityResult:
    """The gas dynamic viscosity, computed in Rust."""
    spec = _spec_for("eos.chung_viscosity")
    result = _core.chung_viscosity(
        omega,
        input_to_si(spec, "Tc", Tc),
        input_to_si(spec, "Vc", Vc),
        input_to_si(spec, "M", M),
        dipole,
        kappa,
        input_to_si(spec, "T", T),
        input_to_si(spec, "V", V),
    )
    return ChungViscosityResult(
        mu=from_si(result.mu.magnitude_si, result.mu.unit),
        warnings=_warnings(result.warnings),
    )


def chung_conductivity(
    Cv0: Q, M: Q, omega: float, Tc: Q, Vc: Q, dipole: float, kappa: float, T: Q
) -> ChungConductivityResult:
    """The gas thermal conductivity, computed in Rust."""
    spec = _spec_for("eos.chung_conductivity")
    result = _core.chung_conductivity(
        input_to_si(spec, "Cv0", Cv0),
        input_to_si(spec, "M", M),
        omega,
        input_to_si(spec, "Tc", Tc),
        input_to_si(spec, "Vc", Vc),
        dipole,
        kappa,
        input_to_si(spec, "T", T),
    )
    return ChungConductivityResult(
        k=from_si(result.k.magnitude_si, result.k.unit),
        warnings=_warnings(result.warnings),
    )


def wilke_viscosity(
    Tc: Sequence[Q],
    Vc: Sequence[Q],
    M: Sequence[Q],
    omega: Sequence[float],
    dipole: Sequence[float],
    kappa: Sequence[float],
    T: Q,
    V: Q,
    z: Sequence[float],
) -> WilkeViscosityResult:
    """The gas mixture dynamic viscosity, computed in Rust."""
    spec = _models_gen.model("eos.wilke_viscosity")
    result = _core.wilke_viscosity(
        [input_to_si(spec, "Tc", value) for value in Tc],
        [input_to_si(spec, "Vc", value) for value in Vc],
        [input_to_si(spec, "M", value) for value in M],
        list(omega),
        list(dipole),
        list(kappa),
        input_to_si(spec, "T", T),
        input_to_si(spec, "V", V),
        list(z),
    )
    return WilkeViscosityResult(
        mu=from_si(result.mu.magnitude_si, result.mu.unit),
        warnings=_warnings(result.warnings),
    )


def mason_saxena_conductivity(
    Cv0: Sequence[Q],
    M: Sequence[Q],
    omega: Sequence[float],
    Tc: Sequence[Q],
    Vc: Sequence[Q],
    dipole: Sequence[float],
    kappa: Sequence[float],
    T: Q,
    z: Sequence[float],
) -> MasonSaxenaConductivityResult:
    """The gas mixture thermal conductivity, computed in Rust."""
    spec = _models_gen.model("eos.mason_saxena_conductivity")
    result = _core.mason_saxena_conductivity(
        [input_to_si(spec, "Cv0", value) for value in Cv0],
        [input_to_si(spec, "M", value) for value in M],
        list(omega),
        [input_to_si(spec, "Tc", value) for value in Tc],
        [input_to_si(spec, "Vc", value) for value in Vc],
        list(dipole),
        list(kappa),
        input_to_si(spec, "T", T),
        list(z),
    )
    return MasonSaxenaConductivityResult(
        k=from_si(result.k.magnitude_si, result.k.unit),
        warnings=_warnings(result.warnings),
    )


def nitric_sulfuric_acid_vapor_pressure(T: Q) -> NitricSulfuricAcidVaporPressureResult:
    """The three acid-system vapour pressures, computed in Rust."""
    spec = _spec_for("eos.nitric_sulfuric_acid_vapor_pressure")
    result = _core.nitric_sulfuric_acid_vapor_pressure(input_to_si(spec, "T", T))
    return NitricSulfuricAcidVaporPressureResult(
        p_water=from_si(result.p_water.magnitude_si, result.p_water.unit),
        p_nitric_acid=from_si(result.p_nitric_acid.magnitude_si, result.p_nitric_acid.unit),
        p_sulfuric_acid=from_si(result.p_sulfuric_acid.magnitude_si, result.p_sulfuric_acid.unit),
        warnings=_warnings(result.warnings),
    )


def nrtl_activity_coefficients(
    params: Any,
    T: Q,
    x: Sequence[float],
) -> NrtlActivityCoefficientsResult:
    """The activity coefficients of a mixture, computed in Rust.

    The resolved matrices cross the boundary flattened row-major, the same two vectors
    the Rust kernel takes.
    """
    spec = _models_gen.model("eos.nrtl_activity_coefficients")
    result = _core.nrtl_activity_coefficients(
        list(params.alpha),
        list(params.dij),
        input_to_si(spec, "T", T),
        list(x),
    )
    return NrtlActivityCoefficientsResult(
        ln_gamma=tuple(result.ln_gamma),
        gamma=tuple(result.gamma),
        warnings=_warnings(result.warnings),
    )


def unifac_activity_coefficients(
    params: Any,
    T: Q,
    x: Sequence[float],
) -> UnifacActivityCoefficientsResult:
    """The activity coefficients of a mixture, computed in Rust.

    The resolved group tables cross the boundary flattened, in the dataclass's own
    field order, the same vectors the Rust kernel takes.
    """
    spec = _models_gen.model("eos.unifac_activity_coefficients")
    result = _core.unifac_activity_coefficients(
        list(params.groups),
        list(params.group_r),
        list(params.group_q),
        list(params.aij),
        input_to_si(spec, "T", T),
        list(x),
    )
    return UnifacActivityCoefficientsResult(
        ln_gamma=tuple(result.ln_gamma),
        gamma=tuple(result.gamma),
        warnings=_warnings(result.warnings),
    )


def van_laar_acid_activity_coefficients(
    params: Any,
    T: Q,
    x: Sequence[float],
) -> VanLaarAcidActivityCoefficientsResult:
    """The activity coefficients of a mixture, computed in Rust.

    The resolved acid identities cross the boundary as the flattened integers they are.
    """
    spec = _models_gen.model("eos.van_laar_acid_activity_coefficients")
    result = _core.van_laar_acid_activity_coefficients(
        list(params.acid_index),
        input_to_si(spec, "T", T),
        list(x),
    )
    return VanLaarAcidActivityCoefficientsResult(
        ln_gamma=tuple(result.ln_gamma),
        gamma=tuple(result.gamma),
        warnings=_warnings(result.warnings),
    )


def unifac_psrk_activity_coefficients(
    params: Any,
    T: Q,
    x: Sequence[float],
) -> UnifacPsrkActivityCoefficientsResult:
    """The activity coefficients of a mixture, computed in Rust.

    The resolved basis and the three interaction matrices cross flattened, in the
    dataclass's own field order.
    """
    spec = _models_gen.model("eos.unifac_psrk_activity_coefficients")
    result = _core.unifac_psrk_activity_coefficients(
        list(params.groups),
        list(params.group_r),
        list(params.group_q),
        list(params.aij),
        list(params.bij),
        list(params.cij),
        input_to_si(spec, "T", T),
        list(x),
    )
    return UnifacPsrkActivityCoefficientsResult(
        ln_gamma=tuple(result.ln_gamma),
        gamma=tuple(result.gamma),
        warnings=_warnings(result.warnings),
    )


def unifac_umrpru_activity_coefficients(
    params: Any,
    T: Q,
    x: Sequence[float],
) -> UnifacUmrpruActivityCoefficientsResult:
    """The activity coefficients of a mixture, computed in Rust.

    The resolved basis and the chosen set's three matrices cross flattened, in the
    dataclass's own field order.
    """
    spec = _models_gen.model("eos.unifac_umrpru_activity_coefficients")
    result = _core.unifac_umrpru_activity_coefficients(
        list(params.groups),
        list(params.group_r),
        list(params.group_q),
        list(params.aij),
        list(params.bij),
        list(params.cij),
        input_to_si(spec, "T", T),
        list(x),
    )
    return UnifacUmrpruActivityCoefficientsResult(
        ln_gamma=tuple(result.ln_gamma),
        gamma=tuple(result.gamma),
        warnings=_warnings(result.warnings),
    )


def uniquac_activity_coefficients(
    params: Any,
    T: Q,
    x: Sequence[float],
    aij: Sequence[Sequence[Q]],
) -> UniquacActivityCoefficientsResult:
    """The activity coefficients of a mixture, computed in Rust.

    The resolved `r` and `q` cross the boundary flattened, in the dataclass's own field
    order; `aij` crosses as the caller supplied it.
    """
    spec = _models_gen.model("eos.uniquac_activity_coefficients")
    result = _core.uniquac_activity_coefficients(
        list(params.r),
        list(params.q),
        input_to_si(spec, "T", T),
        list(x),
        [[input_to_si(spec, "aij", value) for value in row] for row in aij],
    )
    return UniquacActivityCoefficientsResult(
        ln_gamma=tuple(result.ln_gamma),
        gamma=tuple(result.gamma),
        warnings=_warnings(result.warnings),
    )


def wilson_activity_coefficients(
    mixture: Any, T: Q, x: Sequence[float]
) -> WilsonActivityCoefficientsResult:
    """The activity coefficients of a mixture, computed in Rust.

    The mixture is flattened into the critical constants and the molar mass the boundary
    carries, the same prefix `viscosity` sends.
    """
    spec = _models_gen.model("eos.wilson_activity_coefficients")
    molar_mass = []
    for c in mixture.components:
        if c.molar_mass is None:
            raise PropertyUnavailableError(
                "component",
                "molar mass",
                "the paraffin-wax Wilson correlation needs a molar mass for the carbon "
                "number, and a card-added component carries none",
            )
        molar_mass.append(c.molar_mass.to_base_units().magnitude)
    result = _core.wilson_activity_coefficients(
        [c.Tc.to_base_units().magnitude for c in mixture.components],
        [c.Pc.to_base_units().magnitude for c in mixture.components],
        [c.omega for c in mixture.components],
        mixture.flattened_kij(),
        molar_mass,
        input_to_si(spec, "T", T),
        list(x),
    )
    return WilsonActivityCoefficientsResult(
        ln_gamma=tuple(result.ln_gamma),
        gamma=tuple(result.gamma),
        warnings=_warnings(result.warnings),
    )


def umrpr_alpha(omega: float, Tr: float) -> UmrprAlphaResult:
    """The UMR-PR alpha function, computed in Rust."""
    result = _core.umrpr_alpha(omega, Tr)
    return UmrprAlphaResult(alpha=result.alpha, warnings=_warnings(result.warnings))


def tyn_calus_diffusivity(VA: Q, VB: Q, T: Q, eta: Q) -> TynCalusDiffusivityResult:
    """The liquid binary diffusivity, computed in Rust."""
    spec = _spec_for("eos.tyn_calus_diffusivity")
    result = _core.tyn_calus_diffusivity(
        input_to_si(spec, "VA", VA),
        input_to_si(spec, "VB", VB),
        input_to_si(spec, "T", T),
        input_to_si(spec, "eta", eta),
    )
    return TynCalusDiffusivityResult(
        d=from_si(result.d.magnitude_si, result.d.unit),
        warnings=_warnings(result.warnings),
    )


def wilke_chang_diffusivity(phi: float, M: Q, T: Q, eta: Q, VA: Q) -> WilkeChangDiffusivityResult:
    """The liquid binary diffusivity, computed in Rust."""
    spec = _spec_for("eos.wilke_chang_diffusivity")
    result = _core.wilke_chang_diffusivity(
        phi,
        input_to_si(spec, "M", M),
        input_to_si(spec, "T", T),
        input_to_si(spec, "eta", eta),
        input_to_si(spec, "VA", VA),
    )
    return WilkeChangDiffusivityResult(
        d=from_si(result.d.magnitude_si, result.d.unit),
        warnings=_warnings(result.warnings),
    )


def hayduk_minhas_diffusivity(form: str, VA: Q, T: Q, eta: Q) -> HaydukMinhasDiffusivityResult:
    """The liquid binary diffusivity, computed in Rust."""
    spec = _spec_for("eos.hayduk_minhas_diffusivity")
    result = _core.hayduk_minhas_diffusivity(
        form,
        input_to_si(spec, "VA", VA),
        input_to_si(spec, "T", T),
        input_to_si(spec, "eta", eta),
    )
    return HaydukMinhasDiffusivityResult(
        d=from_si(result.d.magnitude_si, result.d.unit),
        warnings=_warnings(result.warnings),
    )


def schwartzentruber_alpha(
    omega: float, p1: float, p2: float, p3: float, Tr: float
) -> SchwartzentruberAlphaResult:
    """The Schwartzentruber-Renon alpha function, computed in Rust."""
    result = _core.schwartzentruber_alpha(omega, p1, p2, p3, Tr)
    return SchwartzentruberAlphaResult(alpha=result.alpha, warnings=_warnings(result.warnings))


def soreide_whitson_alpha(salinity: float, Tr: float) -> SoreideWhitsonAlphaResult:
    """The Soreide-Whitson alpha function for water, computed in Rust."""
    result = _core.soreide_whitson_alpha(salinity, Tr)
    return SoreideWhitsonAlphaResult(alpha=result.alpha, warnings=_warnings(result.warnings))


def siddiqi_lucas_diffusivity(
    form: str, VA: Q, VB: Q, T: Q, eta: Q
) -> SiddiqiLucasDiffusivityResult:
    """The liquid binary diffusivity, computed in Rust."""
    spec = _spec_for("eos.siddiqi_lucas_diffusivity")
    result = _core.siddiqi_lucas_diffusivity(
        form,
        input_to_si(spec, "VA", VA),
        input_to_si(spec, "VB", VB),
        input_to_si(spec, "T", T),
        input_to_si(spec, "eta", eta),
    )
    return SiddiqiLucasDiffusivityResult(
        d=from_si(result.d.magnitude_si, result.d.unit),
        warnings=_warnings(result.warnings),
    )


def co2_water_diffusivity(T: Q) -> Co2WaterDiffusivityResult:
    """The CO2-in-water binary diffusivity, computed in Rust."""
    spec = _spec_for("eos.co2_water_diffusivity")
    result = _core.co2_water_diffusivity(input_to_si(spec, "T", T))
    return Co2WaterDiffusivityResult(
        d=from_si(result.d.magnitude_si, result.d.unit),
        warnings=_warnings(result.warnings),
    )


def parachor_surface_tension(
    parachor: float, rho_l: Q, rho_v: Q, M: Q
) -> ParachorSurfaceTensionResult:
    """The surface tension from the parachor correlation, computed in Rust."""
    spec = _spec_for("eos.parachor_surface_tension")
    result = _core.parachor_surface_tension(
        parachor,
        input_to_si(spec, "rho_l", rho_l),
        input_to_si(spec, "rho_v", rho_v),
        input_to_si(spec, "M", M),
    )
    return ParachorSurfaceTensionResult(
        sigma=from_si(result.sigma.magnitude_si, result.sigma.unit),
        warnings=_warnings(result.warnings),
    )


def viscosity(mixture: Any, T: Q, P: Q, z: Sequence[float]) -> ViscosityResult:
    """The liquid viscosity from the Pedersen correlation, computed in Rust."""
    spec = _models_gen.model("eos.viscosity")
    molar_mass = []
    for c in mixture.components:
        if c.molar_mass is None:
            raise PropertyUnavailableError(
                "component", "molar mass", "a card-added component needs its own molar mass"
            )
        molar_mass.append(c.molar_mass.to_base_units().magnitude)
    result = _core.viscosity(
        [c.Tc.to_base_units().magnitude for c in mixture.components],
        [c.Pc.to_base_units().magnitude for c in mixture.components],
        [c.omega for c in mixture.components],
        mixture.flattened_kij(),
        molar_mass,
        input_to_si(spec, "T", T),
        input_to_si(spec, "P", P),
        list(z),
        mixture.cubic.name,
        mixture.alpha,
        [list(c.alpha_params) for c in mixture.components],
    )
    return ViscosityResult(
        mu=from_si(result.mu.magnitude_si, result.mu.unit),
        warnings=_warnings(result.warnings),
    )


def thermal_conductivity(
    mixture: Any, ideal_gas: Any, T: Q, P: Q, z: Sequence[float]
) -> ThermalConductivityResult:
    """The liquid thermal conductivity from the Pedersen correlation, computed in Rust."""
    spec = _models_gen.model("eos.thermal_conductivity")
    molar_mass = []
    for c in mixture.components:
        if c.molar_mass is None:
            raise PropertyUnavailableError(
                "component", "molar mass", "a card-added component needs its own molar mass"
            )
        molar_mass.append(c.molar_mass.to_base_units().magnitude)
    result = _core.thermal_conductivity(
        [c.Tc.to_base_units().magnitude for c in mixture.components],
        [c.Pc.to_base_units().magnitude for c in mixture.components],
        [c.omega for c in mixture.components],
        mixture.flattened_kij(),
        molar_mass,
        list(ideal_gas.cp_a),
        list(ideal_gas.cp_b),
        list(ideal_gas.cp_c),
        list(ideal_gas.cp_d),
        list(ideal_gas.cp_e),
        input_to_si(spec, "T", T),
        input_to_si(spec, "P", P),
        list(z),
        mixture.cubic.name,
        mixture.alpha,
        [list(c.alpha_params) for c in mixture.components],
    )
    return ThermalConductivityResult(
        k=from_si(result.k.magnitude_si, result.k.unit),
        warnings=_warnings(result.warnings),
    )


def pump_power(rho: Q, q: Q, H: Q, eta: float) -> PumpPowerResult:
    """Pump shaft power, computed in Rust."""
    result = _core.pump_power(
        to_si(rho, "kg/m**3", "rho"),
        to_si(q, "m**3/s", "q"),
        to_si(H, "m", "H"),
        # Dimensionless: no unit to convert, so it crosses as a plain float.
        eta,
    )
    return PumpPowerResult(
        power=from_si(result.power.magnitude_si, result.power.unit),
        warnings=_warnings(result.warnings),
    )


def orifice_flow(d: Q, dP: Q, rho: Q, Cd: float) -> OrificeFlowResult:
    """Orifice flow, computed in Rust."""
    result = _core.orifice_flow(
        to_si(d, "mm", "d"),
        to_si(dP, "Pa", "dP"),
        to_si(rho, "kg/m**3", "rho"),
        # Dimensionless: no unit to convert, so it crosses as a plain float.
        Cd,
    )
    return OrificeFlowResult(
        q=from_si(result.q.magnitude_si, result.q.unit),
        warnings=_warnings(result.warnings),
    )


def control_valve_cv(Cv: float, dP: Q, SG: float) -> ControlValveCvResult:
    """Control-valve flow, computed in Rust."""
    result = _core.control_valve_cv(
        Cv,
        to_si(dP, "Pa", "dP"),
        # Dimensionless: no unit to convert, so it crosses as a plain float.
        SG,
    )
    return ControlValveCvResult(
        q=from_si(result.q.magnitude_si, result.q.unit),
        warnings=_warnings(result.warnings),
    )


def choked_flow_area(m_dot: Q, P0: Q, rho0: Q, k: float) -> ChokedFlowAreaResult:
    """Choked-flow throat area, computed in Rust."""
    result = _core.choked_flow_area(
        to_si(m_dot, "kg/s", "m_dot"),
        to_si(P0, "Pa", "P0"),
        to_si(rho0, "kg/m**3", "rho0"),
        # Dimensionless: no unit to convert, so it crosses as a plain float.
        k,
    )
    return ChokedFlowAreaResult(
        a=from_si(result.a.magnitude_si, result.a.unit),
        warnings=_warnings(result.warnings),
    )


def ideal_gas_cp(
    cp_a: float, cp_b: float, cp_c: float, cp_d: float, cp_e: float, T: Q
) -> IdealGasCpResult:
    """The ideal-gas heat capacity, computed in Rust.

    The five coefficients cross as plain floats, in the units the spec declares them:
    a heat capacity and one per kelvin per degree. Only `T` needs a conversion, and
    only the result needs a unit put back on it.
    """
    spec = _spec_for("eos.ideal_gas_cp")
    result = _core.ideal_gas_cp(
        # Each coefficient carries its own power of temperature, so each has its own
        # conversion. Passing them through unstripped would hand Rust a quantity where
        # it wants an SI magnitude.
        input_to_si(spec, "cp_a", cp_a),
        input_to_si(spec, "cp_b", cp_b),
        input_to_si(spec, "cp_c", cp_c),
        input_to_si(spec, "cp_d", cp_d),
        input_to_si(spec, "cp_e", cp_e),
        input_to_si(spec, "T", T),
    )
    return IdealGasCpResult(
        cp=from_si(result.cp.magnitude_si, result.cp.unit),
        warnings=_warnings(result.warnings),
    )


#: The same table for *models*, kept separate because the calc table is asserted to
#: be exactly the calc registry - ``test_registration_completeness`` compares it by
#: equality in both directions, so a model id in it would break that contract rather
#: than extend it. Models have their own id list for the same reason, and
#: ``test_model_contract.py`` holds this table to it.
def pt_flash(mixture: Any, T: Q, P: Q, z: Sequence[float]) -> PtFlashResult:
    """The two-phase flash of a mixture, computed in Rust.

    The first function here whose arguments are not all scalars: the mixture's
    components cross as three parallel lists and its interaction matrix flattened
    row-major, and the feed as a list. The component order is the one thing the two
    sides have to agree about, and it is the caller's - both implementations are
    handed the same order and produce the same order back.

    `beta` crosses as `Option<f64>` and becomes `None`, not a sentinel. That is the
    design, and flattening it here would undo it one layer above where it was made.
    """
    spec = _models_gen.model("eos.pt_flash")
    result = _core.pt_flash(
        [c.Tc.to_base_units().magnitude for c in mixture.components],
        [c.Pc.to_base_units().magnitude for c in mixture.components],
        # Plain floats: an acentric factor is genuinely dimensionless, so it
        # crosses as the number the caller used.
        [c.omega for c in mixture.components],
        mixture.flattened_kij(),
        input_to_si(spec, "T", T),
        input_to_si(spec, "P", P),
        list(z),
        mixture.cubic.name,
        mixture.alpha,
        [list(c.alpha_params) for c in mixture.components],
    )
    return PtFlashResult(
        beta=result.beta,
        x=tuple(result.x),
        y=tuple(result.y),
        k=tuple(result.k),
        ln_phi_liquid=tuple(result.ln_phi_liquid),
        ln_phi_vapour=tuple(result.ln_phi_vapour),
        z_liquid=result.z_liquid,
        z_vapour=result.z_vapour,
        min_t_over_tc=result.min_t_over_tc,
        phase=_Phase(result.phase),
        iterations=result.iterations,
        residual=result.residual,
        warnings=_warnings(result.warnings),
    )


def pt_phase_envelope(mixture: Any, P: Q, z: Sequence[float]) -> PtPhaseEnvelopeResult:
    """The PT phase envelope of a mixture, computed in Rust.

    The two branches cross as parallel temperature and pressure lists; the scalar
    characteristic points cross as pint quantities.
    """
    spec = _models_gen.model("eos.pt_phase_envelope")
    raw = _core.pt_phase_envelope(
        [c.Tc.to_base_units().magnitude for c in mixture.components],
        [c.Pc.to_base_units().magnitude for c in mixture.components],
        [c.omega for c in mixture.components],
        mixture.flattened_kij(),
        input_to_si(spec, "P", P),
        list(z),
        mixture.cubic.name,
        mixture.alpha,
        [list(c.alpha_params) for c in mixture.components],
    )
    return PtPhaseEnvelopeResult(
        dew_temperature=tuple(raw.dew_temperature),
        dew_pressure=tuple(raw.dew_pressure),
        bubble_temperature=tuple(raw.bubble_temperature),
        bubble_pressure=tuple(raw.bubble_pressure),
        cricondenbar_temperature=from_si(raw.cricondenbar_temperature.magnitude_si, "K"),
        cricondenbar_pressure=from_si(raw.cricondenbar_pressure.magnitude_si, "Pa"),
        cricondentherm_temperature=from_si(raw.cricondentherm_temperature.magnitude_si, "K"),
        cricondentherm_pressure=from_si(raw.cricondentherm_pressure.magnitude_si, "Pa"),
        critical_temperature=from_si(raw.critical_temperature.magnitude_si, "K"),
        critical_pressure=from_si(raw.critical_pressure.magnitude_si, "Pa"),
        iterations=raw.iterations,
        residual=raw.residual,
        warnings=_warnings(raw.warnings),
    )


def ph_flash(mixture: Any, ideal_gas: Any, P: Q, H: Q, z: Sequence[float]) -> PhFlashResult:
    """The pressure-enthalpy flash of a mixture, solved in Rust.

    The ideal-gas vectors cross as six parallel lists plus the two reference states,
    exactly as `molar_enthalpy_entropy` sends them. They are the datum the requested
    enthalpy is a difference from, and a coefficient set without a reference state is
    not a thermodynamic model - which is why they are required here even though this
    model uses only four of the six.

    `beta` crosses as `Option<f64>` and becomes `None`, not a sentinel: a single-phase
    feed has no vapour fraction, and the flash's own extrapolated value would be a
    number a caller could use by mistake.
    """
    spec = _models_gen.model("eos.ph_flash")
    result = _core.ph_flash(
        [c.Tc.to_base_units().magnitude for c in mixture.components],
        [c.Pc.to_base_units().magnitude for c in mixture.components],
        [c.omega for c in mixture.components],
        mixture.flattened_kij(),
        list(ideal_gas.cp_a),
        list(ideal_gas.cp_b),
        list(ideal_gas.cp_c),
        list(ideal_gas.cp_d),
        list(ideal_gas.cp_e),
        input_to_si(spec, "P", P),
        input_to_si(spec, "H", H),
        list(z),
        mixture.cubic.name,
        mixture.alpha,
        [list(c.alpha_params) for c in mixture.components],
    )
    return PhFlashResult(
        T=from_si(result.T.magnitude_si, result.T.unit),
        beta=result.beta,
        x=tuple(result.x),
        y=tuple(result.y),
        k=tuple(result.k),
        phase=_Phase(result.phase),
        z_liquid=result.z_liquid,
        z_vapour=result.z_vapour,
        iterations=result.iterations,
        residual=result.residual,
        warnings=_warnings(result.warnings),
    )


def ps_flash(mixture: Any, ideal_gas: Any, P: Q, S: Q, z: Sequence[float]) -> PsFlashResult:
    """The pressure-entropy flash of a mixture, solved in Rust.

    The isentropic companion to `ph_flash`, and the same arguments: the ideal-gas
    vectors are the datum the requested entropy is a difference from.
    """
    spec = _models_gen.model("eos.ps_flash")
    result = _core.ps_flash(
        [c.Tc.to_base_units().magnitude for c in mixture.components],
        [c.Pc.to_base_units().magnitude for c in mixture.components],
        [c.omega for c in mixture.components],
        mixture.flattened_kij(),
        list(ideal_gas.cp_a),
        list(ideal_gas.cp_b),
        list(ideal_gas.cp_c),
        list(ideal_gas.cp_d),
        list(ideal_gas.cp_e),
        input_to_si(spec, "P", P),
        input_to_si(spec, "S", S),
        list(z),
        mixture.cubic.name,
        mixture.alpha,
        [list(c.alpha_params) for c in mixture.components],
    )
    return PsFlashResult(
        T=from_si(result.T.magnitude_si, result.T.unit),
        beta=result.beta,
        x=tuple(result.x),
        y=tuple(result.y),
        k=tuple(result.k),
        phase=_Phase(result.phase),
        z_liquid=result.z_liquid,
        z_vapour=result.z_vapour,
        iterations=result.iterations,
        residual=result.residual,
        warnings=_warnings(result.warnings),
    )


def tv_flash(mixture: Any, ideal_gas: Any, T: Q, V: Q, z: Sequence[float]) -> TvFlashResult:
    """The temperature-volume flash of a mixture, solved in Rust."""
    spec = _models_gen.model("eos.tv_flash")
    result = _core.tv_flash(
        [c.Tc.to_base_units().magnitude for c in mixture.components],
        [c.Pc.to_base_units().magnitude for c in mixture.components],
        [c.omega for c in mixture.components],
        mixture.flattened_kij(),
        list(ideal_gas.cp_a),
        list(ideal_gas.cp_b),
        list(ideal_gas.cp_c),
        list(ideal_gas.cp_d),
        list(ideal_gas.cp_e),
        input_to_si(spec, "T", T),
        input_to_si(spec, "V", V),
        list(z),
        mixture.cubic.name,
        mixture.alpha,
        [list(c.alpha_params) for c in mixture.components],
    )
    return TvFlashResult(
        P=from_si(result.P.magnitude_si, result.P.unit),
        beta=result.beta,
        x=tuple(result.x),
        y=tuple(result.y),
        k=tuple(result.k),
        phase=_Phase(result.phase),
        z_liquid=result.z_liquid,
        z_vapour=result.z_vapour,
        iterations=result.iterations,
        residual=result.residual,
        warnings=_warnings(result.warnings),
    )


def pvf_flash(
    mixture: Any, P: Q, beta: float, temperature: Q, z: Sequence[float]
) -> PvfFlashResult:
    """The pressure/vapour-fraction flash, solved in Rust.

    No ideal-gas model: nothing here needs an enthalpy or an entropy, and a signature
    that demanded one would oblige a caller to supply heat capacities for a calculation
    that never reads them.
    """
    spec = _models_gen.model("eos.pvf_flash")
    result = _core.pvf_flash(
        [c.Tc.to_base_units().magnitude for c in mixture.components],
        [c.Pc.to_base_units().magnitude for c in mixture.components],
        [c.omega for c in mixture.components],
        mixture.flattened_kij(),
        input_to_si(spec, "P", P),
        float(beta),
        input_to_si(spec, "temperature", temperature),
        list(z),
        mixture.cubic.name,
        mixture.alpha,
        [list(c.alpha_params) for c in mixture.components],
    )
    return PvfFlashResult(
        T=from_si(result.T.magnitude_si, result.T.unit),
        beta=result.beta,
        phase=_Phase(result.phase),
        x=tuple(result.x),
        y=tuple(result.y),
        k=tuple(result.k),
        z_liquid=result.z_liquid,
        z_vapour=result.z_vapour,
        iterations=result.iterations,
        residual=result.residual,
        warnings=_warnings(result.warnings),
    )


def pv_reflux_flash(
    mixture: Any, P: Q, reflux: float, phase: str, temperature: Q, z: Sequence[float]
) -> PvRefluxFlashResult:
    """The pressure/reflux-ratio flash, solved in Rust.

    No ideal-gas model: nothing here needs an enthalpy or an entropy.
    """
    spec = _models_gen.model("eos.pv_reflux_flash")
    result = _core.pv_reflux_flash(
        [c.Tc.to_base_units().magnitude for c in mixture.components],
        [c.Pc.to_base_units().magnitude for c in mixture.components],
        [c.omega for c in mixture.components],
        mixture.flattened_kij(),
        input_to_si(spec, "P", P),
        float(reflux),
        str(phase),
        input_to_si(spec, "temperature", temperature),
        list(z),
        mixture.cubic.name,
        mixture.alpha,
        [list(c.alpha_params) for c in mixture.components],
    )
    return PvRefluxFlashResult(
        T=from_si(result.T.magnitude_si, result.T.unit),
        beta=result.beta,
        phase=_Phase(result.phase),
        x=tuple(result.x),
        y=tuple(result.y),
        k=tuple(result.k),
        z_liquid=result.z_liquid,
        z_vapour=result.z_vapour,
        iterations=result.iterations,
        residual=result.residual,
        warnings=_warnings(result.warnings),
    )


def pv_flash(mixture: Any, ideal_gas: Any, P: Q, V: Q, z: Sequence[float]) -> PvFlashResult:
    """The pressure-volume flash of a mixture, solved in Rust."""
    spec = _models_gen.model("eos.pv_flash")
    result = _core.pv_flash(
        [c.Tc.to_base_units().magnitude for c in mixture.components],
        [c.Pc.to_base_units().magnitude for c in mixture.components],
        [c.omega for c in mixture.components],
        mixture.flattened_kij(),
        list(ideal_gas.cp_a),
        list(ideal_gas.cp_b),
        list(ideal_gas.cp_c),
        list(ideal_gas.cp_d),
        list(ideal_gas.cp_e),
        input_to_si(spec, "P", P),
        input_to_si(spec, "V", V),
        list(z),
        mixture.cubic.name,
        mixture.alpha,
        [list(c.alpha_params) for c in mixture.components],
    )
    return PvFlashResult(
        T=from_si(result.T.magnitude_si, result.T.unit),
        beta=result.beta,
        x=tuple(result.x),
        y=tuple(result.y),
        k=tuple(result.k),
        phase=_Phase(result.phase),
        z_liquid=result.z_liquid,
        z_vapour=result.z_vapour,
        iterations=result.iterations,
        residual=result.residual,
        warnings=_warnings(result.warnings),
    )


def tp_multiflash(mixture: Any, T: Q, P: Q, z: Sequence[float]) -> TpMultiflashResult:
    """How many phases a feed splits into, computed in Rust.

    The two-phase flash plus the stability seeding and a fraction solve on the whole phase
    set. `seeded` crosses as the spec's spelling and is rebuilt here as the enum, so
    `result.seeded is TpMultiflashSeed.STABILITY_SEEDED` holds whichever backend answered.
    """
    spec = _models_gen.model("eos.tp_multiflash")
    result = _core.tp_multiflash(
        [c.Tc.to_base_units().magnitude for c in mixture.components],
        [c.Pc.to_base_units().magnitude for c in mixture.components],
        [c.omega for c in mixture.components],
        mixture.flattened_kij(),
        input_to_si(spec, "T", T),
        input_to_si(spec, "P", P),
        list(z),
        mixture.cubic.name,
        mixture.alpha,
        [list(c.alpha_params) for c in mixture.components],
    )
    return TpMultiflashResult(
        phase_count=result.phase_count,
        beta=tuple(result.beta),
        x=tuple(tuple(row) for row in result.x),
        z_factor=tuple(result.z_factor),
        ln_phi=tuple(tuple(row) for row in result.ln_phi),
        seeded=_TpMultiflashSeed(result.seeded),
        tm=tuple(result.tm),
        iterations=result.iterations,
        residual=result.residual,
        min_t_over_tc=result.min_t_over_tc,
        warnings=_warnings(result.warnings),
    )


def stability_test(mixture: Any, T: Q, P: Q, z: Sequence[float]) -> StabilityTestResult:
    """Whether a feed is stable as a single phase, computed in Rust.

    The same seven arguments as `pt_flash` and the same unpacking, because it is the
    same state asked a different question - which is what makes the two results
    comparable, and what the cross-model test compares.

    `verdict` crosses as the spec's spelling and is rebuilt here as the enum, so
    `result.verdict is StabilityVerdict.UNSTABLE` holds whichever backend answered.
    The `tm` and `w` rows come back in trial order - vapour-like first - and are
    tuples because every other sequence in a result dataclass is one.
    """
    spec = _models_gen.model("eos.stability_test")
    result = _core.stability_test(
        [c.Tc.to_base_units().magnitude for c in mixture.components],
        [c.Pc.to_base_units().magnitude for c in mixture.components],
        [c.omega for c in mixture.components],
        mixture.flattened_kij(),
        input_to_si(spec, "T", T),
        input_to_si(spec, "P", P),
        list(z),
        mixture.cubic.name,
        mixture.alpha,
        [list(c.alpha_params) for c in mixture.components],
    )
    return StabilityTestResult(
        verdict=_StabilityVerdict(result.verdict),
        tm=tuple(result.tm),
        w=tuple(tuple(row) for row in result.w),
        iterations=tuple(result.iterations),
        min_t_over_tc=result.min_t_over_tc,
        warnings=_warnings(result.warnings),
    )


def th_flash(mixture: Any, ideal_gas: Any, T: Q, H: Q, z: Sequence[float]) -> ThFlashResult:
    """The temperature-enthalpy flash of a mixture, solved in Rust."""
    spec = _models_gen.model("eos.th_flash")
    result = _core.th_flash(
        [c.Tc.to_base_units().magnitude for c in mixture.components],
        [c.Pc.to_base_units().magnitude for c in mixture.components],
        [c.omega for c in mixture.components],
        mixture.flattened_kij(),
        list(ideal_gas.cp_a),
        list(ideal_gas.cp_b),
        list(ideal_gas.cp_c),
        list(ideal_gas.cp_d),
        list(ideal_gas.cp_e),
        input_to_si(spec, "T", T),
        input_to_si(spec, "H", H),
        list(z),
        mixture.cubic.name,
        mixture.alpha,
        [list(c.alpha_params) for c in mixture.components],
    )
    return ThFlashResult(
        P=from_si(result.P.magnitude_si, result.P.unit),
        beta=result.beta,
        x=tuple(result.x),
        y=tuple(result.y),
        k=tuple(result.k),
        phase=_Phase(result.phase),
        z_liquid=result.z_liquid,
        z_vapour=result.z_vapour,
        iterations=result.iterations,
        residual=result.residual,
        warnings=_warnings(result.warnings),
    )


def ts_flash(mixture: Any, ideal_gas: Any, T: Q, S: Q, z: Sequence[float]) -> TsFlashResult:
    """The temperature-entropy flash of a mixture, solved in Rust."""
    spec = _models_gen.model("eos.ts_flash")
    result = _core.ts_flash(
        [c.Tc.to_base_units().magnitude for c in mixture.components],
        [c.Pc.to_base_units().magnitude for c in mixture.components],
        [c.omega for c in mixture.components],
        mixture.flattened_kij(),
        list(ideal_gas.cp_a),
        list(ideal_gas.cp_b),
        list(ideal_gas.cp_c),
        list(ideal_gas.cp_d),
        list(ideal_gas.cp_e),
        input_to_si(spec, "T", T),
        input_to_si(spec, "S", S),
        list(z),
        mixture.cubic.name,
        mixture.alpha,
        [list(c.alpha_params) for c in mixture.components],
    )
    return TsFlashResult(
        P=from_si(result.P.magnitude_si, result.P.unit),
        beta=result.beta,
        x=tuple(result.x),
        y=tuple(result.y),
        k=tuple(result.k),
        phase=_Phase(result.phase),
        z_liquid=result.z_liquid,
        z_vapour=result.z_vapour,
        iterations=result.iterations,
        residual=result.residual,
        warnings=_warnings(result.warnings),
    )


def tu_flash(mixture: Any, ideal_gas: Any, T: Q, U: Q, z: Sequence[float]) -> TuFlashResult:
    """The temperature-internal-energy flash of a mixture, solved in Rust."""
    spec = _models_gen.model("eos.tu_flash")
    result = _core.tu_flash(
        [c.Tc.to_base_units().magnitude for c in mixture.components],
        [c.Pc.to_base_units().magnitude for c in mixture.components],
        [c.omega for c in mixture.components],
        mixture.flattened_kij(),
        list(ideal_gas.cp_a),
        list(ideal_gas.cp_b),
        list(ideal_gas.cp_c),
        list(ideal_gas.cp_d),
        list(ideal_gas.cp_e),
        input_to_si(spec, "T", T),
        input_to_si(spec, "U", U),
        list(z),
        mixture.cubic.name,
        mixture.alpha,
        [list(c.alpha_params) for c in mixture.components],
    )
    return TuFlashResult(
        P=from_si(result.P.magnitude_si, result.P.unit),
        beta=result.beta,
        x=tuple(result.x),
        y=tuple(result.y),
        k=tuple(result.k),
        phase=_Phase(result.phase),
        z_liquid=result.z_liquid,
        z_vapour=result.z_vapour,
        iterations=result.iterations,
        residual=result.residual,
        warnings=_warnings(result.warnings),
    )


def pu_flash(mixture: Any, ideal_gas: Any, P: Q, U: Q, z: Sequence[float]) -> PuFlashResult:
    """The pressure-internal-energy flash of a mixture, solved in Rust."""
    spec = _models_gen.model("eos.pu_flash")
    result = _core.pu_flash(
        [c.Tc.to_base_units().magnitude for c in mixture.components],
        [c.Pc.to_base_units().magnitude for c in mixture.components],
        [c.omega for c in mixture.components],
        mixture.flattened_kij(),
        list(ideal_gas.cp_a),
        list(ideal_gas.cp_b),
        list(ideal_gas.cp_c),
        list(ideal_gas.cp_d),
        list(ideal_gas.cp_e),
        input_to_si(spec, "P", P),
        input_to_si(spec, "U", U),
        list(z),
        mixture.cubic.name,
        mixture.alpha,
        [list(c.alpha_params) for c in mixture.components],
    )
    return PuFlashResult(
        T=from_si(result.T.magnitude_si, result.T.unit),
        beta=result.beta,
        x=tuple(result.x),
        y=tuple(result.y),
        k=tuple(result.k),
        phase=_Phase(result.phase),
        z_liquid=result.z_liquid,
        z_vapour=result.z_vapour,
        iterations=result.iterations,
        residual=result.residual,
        warnings=_warnings(result.warnings),
    )


def tv_fraction_flash(
    mixture: Any, T: Q, fraction: float, P: Q, z: Sequence[float]
) -> TvFractionFlashResult:
    """The temperature and vapour-volume-fraction flash, solved in Rust."""
    spec = _models_gen.model("eos.tv_fraction_flash")
    result = _core.tv_fraction_flash(
        [c.Tc.to_base_units().magnitude for c in mixture.components],
        [c.Pc.to_base_units().magnitude for c in mixture.components],
        [c.omega for c in mixture.components],
        mixture.flattened_kij(),
        input_to_si(spec, "T", T),
        float(fraction),
        input_to_si(spec, "P", P),
        list(z),
        mixture.cubic.name,
        mixture.alpha,
        [list(c.alpha_params) for c in mixture.components],
    )
    return TvFractionFlashResult(
        P=from_si(result.P.magnitude_si, result.P.unit),
        T=from_si(result.T.magnitude_si, result.T.unit),
        beta=result.beta,
        volume_fraction=result.volume_fraction,
        phase=_Phase(result.phase),
        x=tuple(result.x),
        y=tuple(result.y),
        k=tuple(result.k),
        z_liquid=result.z_liquid,
        z_vapour=result.z_vapour,
        iterations=result.iterations,
        residual=result.residual,
        warnings=_warnings(result.warnings),
    )


def vu_flash(mixture: Any, ideal_gas: Any, V: Q, U: Q, z: Sequence[float]) -> VuFlashResult:
    """The volume-internal-energy flash of a mixture, solved in Rust."""
    spec = _models_gen.model("eos.vu_flash")
    result = _core.vu_flash(
        [c.Tc.to_base_units().magnitude for c in mixture.components],
        [c.Pc.to_base_units().magnitude for c in mixture.components],
        [c.omega for c in mixture.components],
        mixture.flattened_kij(),
        list(ideal_gas.cp_a),
        list(ideal_gas.cp_b),
        list(ideal_gas.cp_c),
        list(ideal_gas.cp_d),
        list(ideal_gas.cp_e),
        input_to_si(spec, "V", V),
        input_to_si(spec, "U", U),
        list(z),
        mixture.cubic.name,
        mixture.alpha,
        [list(c.alpha_params) for c in mixture.components],
    )
    return VuFlashResult(
        P=from_si(result.P.magnitude_si, result.P.unit),
        T=from_si(result.T.magnitude_si, result.T.unit),
        beta=result.beta,
        x=tuple(result.x),
        y=tuple(result.y),
        k=tuple(result.k),
        phase=_Phase(result.phase),
        z_liquid=result.z_liquid,
        z_vapour=result.z_vapour,
        iterations=result.iterations,
        residual=result.residual,
        warnings=_warnings(result.warnings),
    )


def vs_flash(mixture: Any, ideal_gas: Any, V: Q, S: Q, z: Sequence[float]) -> VsFlashResult:
    """The volume-entropy flash of a mixture, solved in Rust."""
    spec = _models_gen.model("eos.vs_flash")
    result = _core.vs_flash(
        [c.Tc.to_base_units().magnitude for c in mixture.components],
        [c.Pc.to_base_units().magnitude for c in mixture.components],
        [c.omega for c in mixture.components],
        mixture.flattened_kij(),
        list(ideal_gas.cp_a),
        list(ideal_gas.cp_b),
        list(ideal_gas.cp_c),
        list(ideal_gas.cp_d),
        list(ideal_gas.cp_e),
        input_to_si(spec, "V", V),
        input_to_si(spec, "S", S),
        list(z),
        mixture.cubic.name,
        mixture.alpha,
        [list(c.alpha_params) for c in mixture.components],
    )
    return VsFlashResult(
        P=from_si(result.P.magnitude_si, result.P.unit),
        T=from_si(result.T.magnitude_si, result.T.unit),
        beta=result.beta,
        x=tuple(result.x),
        y=tuple(result.y),
        k=tuple(result.k),
        phase=_Phase(result.phase),
        z_liquid=result.z_liquid,
        z_vapour=result.z_vapour,
        iterations=result.iterations,
        residual=result.residual,
        warnings=_warnings(result.warnings),
    )


def vh_flash(mixture: Any, ideal_gas: Any, V: Q, H: Q, z: Sequence[float]) -> VhFlashResult:
    """The volume-enthalpy flash of a mixture, solved in Rust."""
    spec = _models_gen.model("eos.vh_flash")
    result = _core.vh_flash(
        [c.Tc.to_base_units().magnitude for c in mixture.components],
        [c.Pc.to_base_units().magnitude for c in mixture.components],
        [c.omega for c in mixture.components],
        mixture.flattened_kij(),
        list(ideal_gas.cp_a),
        list(ideal_gas.cp_b),
        list(ideal_gas.cp_c),
        list(ideal_gas.cp_d),
        list(ideal_gas.cp_e),
        input_to_si(spec, "V", V),
        input_to_si(spec, "H", H),
        list(z),
        mixture.cubic.name,
        mixture.alpha,
        [list(c.alpha_params) for c in mixture.components],
    )
    return VhFlashResult(
        P=from_si(result.P.magnitude_si, result.P.unit),
        T=from_si(result.T.magnitude_si, result.T.unit),
        beta=result.beta,
        x=tuple(result.x),
        y=tuple(result.y),
        k=tuple(result.k),
        phase=_Phase(result.phase),
        z_liquid=result.z_liquid,
        z_vapour=result.z_vapour,
        iterations=result.iterations,
        residual=result.residual,
        warnings=_warnings(result.warnings),
    )


def vu_flash_single_comp(mixture: Any, ideal_gas: Any, P: Q, V: Q, U: Q) -> VuFlashSingleCompResult:
    """The pure-component VU state, computed in Rust."""
    spec = _models_gen.model("eos.vu_flash_single_comp")
    result = _core.vu_flash_single_comp(
        [c.Tc.to_base_units().magnitude for c in mixture.components],
        [c.Pc.to_base_units().magnitude for c in mixture.components],
        [c.omega for c in mixture.components],
        mixture.flattened_kij(),
        list(ideal_gas.cp_a),
        list(ideal_gas.cp_b),
        list(ideal_gas.cp_c),
        list(ideal_gas.cp_d),
        list(ideal_gas.cp_e),
        input_to_si(spec, "P", P),
        input_to_si(spec, "V", V),
        input_to_si(spec, "U", U),
        mixture.cubic.name,
        mixture.alpha,
        [list(c.alpha_params) for c in mixture.components],
    )
    return VuFlashSingleCompResult(
        T=from_si(result.T.magnitude_si, result.T.unit),
        beta=result.beta,
        V=from_si(result.V.magnitude_si, result.V.unit),
        phase=_Phase(result.phase),
        warnings=_warnings(result.warnings),
    )


def _boundary_result(raw: Any, result_type: Any, *, liquid_first: bool) -> Any:
    """Shared unpacking for the two phase-boundary models.

    The two differ only in whether the held phase is the liquid, so the `z_liquid`
    and `z_vapour` fields land in the opposite order. Doing that here rather than
    twice is the same argument the Rust side makes for one shared iteration.
    """
    held, incipient = (raw.z_liquid, raw.z_vapour) if liquid_first else (raw.z_vapour, raw.z_liquid)
    return result_type(
        pressure=from_si(raw.pressure.magnitude_si, raw.pressure.unit),
        incipient=tuple(raw.incipient),
        k=tuple(raw.k),
        z_liquid=held if liquid_first else incipient,
        z_vapour=incipient if liquid_first else held,
        min_t_over_tc=raw.min_t_over_tc,
        iterations=raw.iterations,
        residual=raw.residual,
        warnings=_warnings(raw.warnings),
    )


def _boundary_temperature_result(raw: Any, result_type: Any) -> Any:
    """Unpack a temperature boundary, whose scalar is a temperature rather than a
    pressure. The `z_liquid` and `z_vapour` fields already name the phases they say
    they do, so there is no swap here."""
    return result_type(
        temperature=from_si(raw.temperature.magnitude_si, raw.temperature.unit),
        incipient=tuple(raw.incipient),
        k=tuple(raw.k),
        z_liquid=raw.z_liquid,
        z_vapour=raw.z_vapour,
        min_t_over_tc=raw.min_t_over_tc,
        iterations=raw.iterations,
        residual=raw.residual,
        warnings=_warnings(raw.warnings),
    )


def bubble_pressure(mixture: Any, T: Q, x: Sequence[float]) -> BubblePressureResult:
    """The bubble-point pressure of a mixture, computed in Rust.

    The held composition crosses as a list, as the flash's feed does. The Rust side
    rebuilds the mixture from the same three vectors and flattened matrix the flash
    uses, so the component order is the one thing the two sides agree about.
    """
    spec = _models_gen.model("eos.bubble_pressure")
    raw = _core.bubble_pressure(
        [c.Tc.to_base_units().magnitude for c in mixture.components],
        [c.Pc.to_base_units().magnitude for c in mixture.components],
        [c.omega for c in mixture.components],
        mixture.flattened_kij(),
        input_to_si(spec, "T", T),
        list(x),
        mixture.cubic.name,
        mixture.alpha,
        [list(c.alpha_params) for c in mixture.components],
    )
    return _boundary_result(raw, BubblePressureResult, liquid_first=True)  # type: ignore[no-any-return]


def bubble_temperature(mixture: Any, P: Q, x: Sequence[float]) -> BubbleTemperatureResult:
    """The bubble-point temperature of a mixture, computed in Rust.

    The mirror of :func:`bubble_pressure`: the held composition crosses as a list and
    the pressure is the fixed variable instead of the temperature.
    """
    spec = _models_gen.model("eos.bubble_temperature")
    raw = _core.bubble_temperature(
        [c.Tc.to_base_units().magnitude for c in mixture.components],
        [c.Pc.to_base_units().magnitude for c in mixture.components],
        [c.omega for c in mixture.components],
        mixture.flattened_kij(),
        input_to_si(spec, "P", P),
        list(x),
        mixture.cubic.name,
        mixture.alpha,
        [list(c.alpha_params) for c in mixture.components],
    )
    return _boundary_temperature_result(raw, BubbleTemperatureResult)  # type: ignore[no-any-return]


def critical_point(mixture: Any, z: Sequence[float]) -> CriticalPointResult:
    """The critical point of a mixture of composition ``z``, computed in Rust.

    The composition crosses as a list, as the flash's feed does. The Rust side rebuilds
    the mixture from the same three vectors and flattened matrix the other models use,
    so the component order is the one thing the two sides agree about.
    """
    raw = _core.critical_point(
        [c.Tc.to_base_units().magnitude for c in mixture.components],
        [c.Pc.to_base_units().magnitude for c in mixture.components],
        [c.omega for c in mixture.components],
        mixture.flattened_kij(),
        list(z),
        mixture.cubic.name,
        mixture.alpha,
        [list(c.alpha_params) for c in mixture.components],
    )
    return CriticalPointResult(
        tc=from_si(raw.tc.magnitude_si, "K"),
        pc=from_si(raw.pc.magnitude_si, "Pa"),
        vc=from_si(raw.vc.magnitude_si, "m**3/mol"),
        z_c=raw.z_c,
        iterations=raw.iterations,
        residual=raw.residual,
        warnings=_warnings(raw.warnings),
    )


def dew_pressure(mixture: Any, T: Q, y: Sequence[float]) -> DewPressureResult:
    """The dew-point pressure of a mixture, computed in Rust.

    The mirror of the bubble point: the held phase is the vapour, so the Rust
    result's `z_liquid` and `z_vapour` are swapped on the way back to keep the
    Python field names naming the phases they say they do.
    """
    spec = _models_gen.model("eos.dew_pressure")
    raw = _core.dew_pressure(
        [c.Tc.to_base_units().magnitude for c in mixture.components],
        [c.Pc.to_base_units().magnitude for c in mixture.components],
        [c.omega for c in mixture.components],
        mixture.flattened_kij(),
        input_to_si(spec, "T", T),
        list(y),
        mixture.cubic.name,
        mixture.alpha,
        [list(c.alpha_params) for c in mixture.components],
    )
    return _boundary_result(raw, DewPressureResult, liquid_first=False)  # type: ignore[no-any-return]


def capillary_dew_point(
    mixture: Any,
    P: Q,
    y: Sequence[float],
    pore_radius: Q,
    contact_angle: Q,
    surface_tension: Q,
) -> CapillaryDewPointResult:
    """The dew point of a vapour held in a pore, computed in Rust.

    The same unpacking as `dew_temperature` plus the three curvature arguments. The surface
    tension is the caller's, not the mixture's: it is a fitted quantity with its own provenance,
    and reading one inside a model whose other inputs are all stated would hide which was used.
    """
    spec = _models_gen.model("eos.capillary_dew_point")
    result = _core.capillary_dew_point(
        [c.Tc.to_base_units().magnitude for c in mixture.components],
        [c.Pc.to_base_units().magnitude for c in mixture.components],
        [c.omega for c in mixture.components],
        mixture.flattened_kij(),
        input_to_si(spec, "P", P),
        list(y),
        input_to_si(spec, "pore_radius", pore_radius),
        # The angle is dimensionless and a bare number by this library's convention for
        # dimensionless scalars; a caller who wrote it as a quantity is taken at its magnitude.
        (
            contact_angle.to_base_units().magnitude
            if hasattr(contact_angle, "to_base_units")
            else float(contact_angle)
        ),
        input_to_si(spec, "surface_tension", surface_tension),
        mixture.cubic.name,
        mixture.alpha,
        [list(c.alpha_params) for c in mixture.components],
    )
    return CapillaryDewPointResult(
        temperature=from_si(result.temperature.magnitude_si, result.temperature.unit),
        incipient=tuple(result.incipient),
        k=tuple(result.k),
        z_liquid=result.z_liquid,
        z_vapour=result.z_vapour,
        capillary_pressure=from_si(
            result.capillary_pressure.magnitude_si, result.capillary_pressure.unit
        ),
        min_t_over_tc=result.min_t_over_tc,
        iterations=result.iterations,
        residual=result.residual,
        warnings=_warnings(result.warnings),
    )


def dew_temperature(mixture: Any, P: Q, y: Sequence[float]) -> DewTemperatureResult:
    """The dew-point temperature of a mixture, computed in Rust.

    The mirror of the bubble point: the held phase is the vapour, so the answer is the
    temperature at which that vapour first gives off liquid at a fixed pressure.
    """
    spec = _models_gen.model("eos.dew_temperature")
    raw = _core.dew_temperature(
        [c.Tc.to_base_units().magnitude for c in mixture.components],
        [c.Pc.to_base_units().magnitude for c in mixture.components],
        [c.omega for c in mixture.components],
        mixture.flattened_kij(),
        input_to_si(spec, "P", P),
        list(y),
        mixture.cubic.name,
        mixture.alpha,
        [list(c.alpha_params) for c in mixture.components],
    )
    return _boundary_temperature_result(raw, DewTemperatureResult)  # type: ignore[no-any-return]


def molar_enthalpy_entropy(
    mixture: Any, ideal_gas: Any, T: Q, P: Q, z: Sequence[float], compressibility: float
) -> MolarEnthalpyEntropyResult:
    """The absolute enthalpy and entropy of a mixture, computed in Rust.

    The `IdealGasModel` and the `Mixture` are unpacked into the flat vectors the
    boundary carries, so the component order is the one thing the two sides agree
    about - the same arrangement `pt_flash` uses.
    """
    spec = _models_gen.model("eos.molar_enthalpy_entropy")
    result = _core.molar_enthalpy_entropy(
        [c.Tc.to_base_units().magnitude for c in mixture.components],
        [c.Pc.to_base_units().magnitude for c in mixture.components],
        [c.omega for c in mixture.components],
        mixture.flattened_kij(),
        list(ideal_gas.cp_a),
        list(ideal_gas.cp_b),
        list(ideal_gas.cp_c),
        list(ideal_gas.cp_d),
        list(ideal_gas.cp_e),
        input_to_si(spec, "T", T),
        input_to_si(spec, "P", P),
        list(z),
        compressibility,
        mixture.cubic.name,
        mixture.alpha,
        [list(c.alpha_params) for c in mixture.components],
    )
    return MolarEnthalpyEntropyResult(
        h=from_si(result.h.magnitude_si, result.h.unit),
        s=from_si(result.s.magnitude_si, result.s.unit),
        h_ideal=from_si(result.h_ideal.magnitude_si, result.h_ideal.unit),
        s_ideal=from_si(result.s_ideal.magnitude_si, result.s_ideal.unit),
        h_departure=from_si(result.h_departure.magnitude_si, result.h_departure.unit),
        s_departure=from_si(result.s_departure.magnitude_si, result.s_departure.unit),
        psi_bar=result.psi_bar,
        cp=from_si(result.cp.magnitude_si, result.cp.unit),
        cp_ideal=from_si(result.cp_ideal.magnitude_si, result.cp_ideal.unit),
        cp_departure=from_si(result.cp_departure.magnitude_si, result.cp_departure.unit),
        warnings=_warnings(result.warnings),
    )


def bwrs_phase(coeffs: Any, T: Q, P: Q, z: Sequence[float]) -> BwrsPhaseResult:
    """The BWRS (MBWR-32) phase state, computed in Rust.

    The per-component coefficient sets are unpacked into the flattened ``a`` vector
    (32 per component, component-major) and the ``rhoc`` vector the boundary carries.
    """
    spec = _models_gen.model("eos.bwrs_phase")
    result = _core.bwrs_phase(
        [v for c in coeffs for v in c.a],
        [c.rhoc for c in coeffs],
        input_to_si(spec, "T", T),
        input_to_si(spec, "P", P),
        list(z),
    )
    return BwrsPhaseResult(
        z_factor=result.z_factor,
        ln_phi=tuple(result.ln_phi),
        h_res=from_si(result.h_res.magnitude_si, result.h_res.unit),
        s_res=from_si(result.s_res.magnitude_si, result.s_res.unit),
        cp_res=from_si(result.cp_res.magnitude_si, result.cp_res.unit),
        warnings=_warnings(result.warnings),
    )


def ammonia_phase(T: Q, P: Q) -> AmmoniaPhaseResult:
    """The ammonia reference phase state, computed in Rust."""
    spec = _models_gen.model("eos.ammonia_phase")
    result = _core.ammonia_phase(input_to_si(spec, "T", T), input_to_si(spec, "P", P))
    return AmmoniaPhaseResult(
        z_factor=result.z_factor,
        u=from_si(result.u.magnitude_si, result.u.unit),
        h=from_si(result.h.magnitude_si, result.h.unit),
        s=from_si(result.s.magnitude_si, result.s.unit),
        cv=from_si(result.cv.magnitude_si, result.cv.unit),
        cp=from_si(result.cp.magnitude_si, result.cp.unit),
        g=from_si(result.g.magnitude_si, result.g.unit),
        warnings=_warnings(result.warnings),
    )


def co2_phase(T: Q, P: Q) -> Co2PhaseResult:
    """The Span-Wagner CO2 phase state, computed in Rust."""
    spec = _models_gen.model("eos.co2_phase")
    result = _core.co2_phase(input_to_si(spec, "T", T), input_to_si(spec, "P", P))
    return Co2PhaseResult(
        z_factor=result.z_factor,
        u=from_si(result.u.magnitude_si, result.u.unit),
        h=from_si(result.h.magnitude_si, result.h.unit),
        s=from_si(result.s.magnitude_si, result.s.unit),
        cv=from_si(result.cv.magnitude_si, result.cv.unit),
        cp=from_si(result.cp.magnitude_si, result.cp.unit),
        g=from_si(result.g.magnitude_si, result.g.unit),
        warnings=_warnings(result.warnings),
    )


def helium_phase(T: Q, P: Q) -> HeliumPhaseResult:
    """The Vega helium phase state, computed in Rust."""
    spec = _models_gen.model("eos.helium_phase")
    result = _core.helium_phase(input_to_si(spec, "T", T), input_to_si(spec, "P", P))
    return HeliumPhaseResult(
        z_factor=result.z_factor,
        u=from_si(result.u.magnitude_si, result.u.unit),
        h=from_si(result.h.magnitude_si, result.h.unit),
        s=from_si(result.s.magnitude_si, result.s.unit),
        cv=from_si(result.cv.magnitude_si, result.cv.unit),
        cp=from_si(result.cp.magnitude_si, result.cp.unit),
        g=from_si(result.g.magnitude_si, result.g.unit),
        warnings=_warnings(result.warnings),
    )


def water_phase(T: Q, P: Q) -> WaterPhaseResult:
    """The IAPWS-IF97 water phase state, computed in Rust."""
    spec = _models_gen.model("eos.water_phase")
    result = _core.water_phase(input_to_si(spec, "T", T), input_to_si(spec, "P", P))
    return WaterPhaseResult(
        z_factor=result.z_factor,
        u=from_si(result.u.magnitude_si, result.u.unit),
        h=from_si(result.h.magnitude_si, result.h.unit),
        s=from_si(result.s.magnitude_si, result.s.unit),
        cv=from_si(result.cv.magnitude_si, result.cv.unit),
        cp=from_si(result.cp.magnitude_si, result.cp.unit),
        g=from_si(result.g.magnitude_si, result.g.unit),
        warnings=_warnings(result.warnings),
    )


def argon_solid_phase(T: Q, P: Q) -> ArgonSolidPhaseResult:
    """The solid argon phase state, computed in Rust."""
    spec = _models_gen.model("eos.argon_solid_phase")
    result = _core.argon_solid_phase(input_to_si(spec, "T", T), input_to_si(spec, "P", P))
    return ArgonSolidPhaseResult(
        z_factor=result.z_factor,
        u=from_si(result.u.magnitude_si, result.u.unit),
        h=from_si(result.h.magnitude_si, result.h.unit),
        s=from_si(result.s.magnitude_si, result.s.unit),
        cv=from_si(result.cv.magnitude_si, result.cv.unit),
        cp=from_si(result.cp.magnitude_si, result.cp.unit),
        g=from_si(result.g.magnitude_si, result.g.unit),
        warnings=_warnings(result.warnings),
    )


def parahydrogen_solid_phase(T: Q, P: Q) -> ParahydrogenSolidPhaseResult:
    """The solid para-hydrogen phase state, computed in Rust."""
    spec = _models_gen.model("eos.parahydrogen_solid_phase")
    result = _core.parahydrogen_solid_phase(input_to_si(spec, "T", T), input_to_si(spec, "P", P))
    return ParahydrogenSolidPhaseResult(
        z_factor=result.z_factor,
        u=from_si(result.u.magnitude_si, result.u.unit),
        h=from_si(result.h.magnitude_si, result.h.unit),
        s=from_si(result.s.magnitude_si, result.s.unit),
        cv=from_si(result.cv.magnitude_si, result.cv.unit),
        cp=from_si(result.cp.magnitude_si, result.cp.unit),
        g=from_si(result.g.magnitude_si, result.g.unit),
        warnings=_warnings(result.warnings),
    )


def eos_cg_phase(components: Sequence[str], T: Q, P: Q, z: Sequence[float]) -> EosCgPhaseResult:
    """The EOS-CG phase state, computed in Rust."""
    spec = _models_gen.model("eos.eos_cg_phase")
    result = _core.eos_cg_phase(
        list(components), input_to_si(spec, "T", T), input_to_si(spec, "P", P), list(z)
    )
    return EosCgPhaseResult(
        z_factor=result.z_factor,
        u=from_si(result.u.magnitude_si, result.u.unit),
        h=from_si(result.h.magnitude_si, result.h.unit),
        s=from_si(result.s.magnitude_si, result.s.unit),
        cv=from_si(result.cv.magnitude_si, result.cv.unit),
        cp=from_si(result.cp.magnitude_si, result.cp.unit),
        g=from_si(result.g.magnitude_si, result.g.unit),
        warnings=_warnings(result.warnings),
    )


def ge_nrtl_phase(params: Any, T: Q, P: Q, x: Sequence[float]) -> GeNrtlPhaseResult:
    """The fugacity coefficients of an NRTL liquid, computed in Rust.

    The resolved record crosses flattened, one list per field, in the dataclass's own
    order.
    """
    spec = _models_gen.model("eos.ge_nrtl_phase")
    result = _core.ge_nrtl_phase(
        list(params.alpha),
        list(params.dij),
        list(params.antoine_type),
        list(params.antoine_coefficients),
        list(params.antoine_tc),
        list(params.antoine_pc),
        input_to_si(spec, "T", T),
        input_to_si(spec, "P", P),
        list(x),
    )
    return GeNrtlPhaseResult(
        gamma=tuple(result.gamma),
        ln_gamma=tuple(result.ln_gamma),
        ln_phi=tuple(result.ln_phi),
        p_sat=tuple(from_si(value.magnitude_si, value.unit) for value in result.p_sat),
        warnings=_warnings(result.warnings),
    )


def ge_wilson_phase(
    params: Any, mixture: Any, T: Q, P: Q, x: Sequence[float]
) -> GeWilsonPhaseResult:
    """The fugacity coefficients of a Wilson liquid, computed in Rust.

    Both objects cross: the mixture as three parallel lists - the Wilson correlation
    reads the molar mass and the critical temperature, and nothing else - then the
    record's own vapour-pressure fields.
    """
    spec = _models_gen.model("eos.ge_wilson_phase")
    result = _core.ge_wilson_phase(
        [c.Tc.to_base_units().magnitude for c in mixture.components],
        [c.Pc.to_base_units().magnitude for c in mixture.components],
        [c.omega for c in mixture.components],
        mixture.flattened_kij(),
        [c.molar_mass.to_base_units().magnitude for c in mixture.components],
        list(params.antoine_type),
        list(params.antoine_coefficients),
        list(params.antoine_tc),
        list(params.antoine_pc),
        input_to_si(spec, "T", T),
        input_to_si(spec, "P", P),
        list(x),
        mixture.cubic.name,
        mixture.alpha,
        [list(c.alpha_params) for c in mixture.components],
    )
    return GeWilsonPhaseResult(
        gamma=tuple(result.gamma),
        ln_gamma=tuple(result.ln_gamma),
        ln_phi=tuple(result.ln_phi),
        p_sat=tuple(from_si(value.magnitude_si, value.unit) for value in result.p_sat),
        warnings=_warnings(result.warnings),
    )


def ge_uniquac_phase(
    params: Any, T: Q, P: Q, x: Sequence[float], aij: Sequence[Sequence[Q]]
) -> GeUniquacPhaseResult:
    """The fugacity coefficients of a UNIQUAC liquid, computed in Rust.

    Both the resolved record and `aij` cross flattened: the record's fields in the
    dataclass's own order, then the interaction matrix row by row, row-major.
    """
    spec = _models_gen.model("eos.ge_uniquac_phase")
    result = _core.ge_uniquac_phase(
        list(params.r),
        list(params.q),
        list(params.antoine_type),
        list(params.antoine_coefficients),
        list(params.antoine_tc),
        list(params.antoine_pc),
        input_to_si(spec, "T", T),
        input_to_si(spec, "P", P),
        list(x),
        [[input_to_si(spec, "aij", value) for value in row] for row in aij],
    )
    return GeUniquacPhaseResult(
        gamma=tuple(result.gamma),
        ln_gamma=tuple(result.ln_gamma),
        ln_phi=tuple(result.ln_phi),
        p_sat=tuple(from_si(value.magnitude_si, value.unit) for value in result.p_sat),
        warnings=_warnings(result.warnings),
    )


def ge_van_laar_acid_phase(params: Any, T: Q, P: Q, x: Sequence[float]) -> GeVanLaarAcidPhaseResult:
    """The acid liquid's fugacity coefficients, computed in Rust.

    The resolved record crosses flattened, one list per field, in the dataclass's own
    order - the acid identities first, then the vapour-pressure columns beside them.
    """
    spec = _models_gen.model("eos.ge_van_laar_acid_phase")
    result = _core.ge_van_laar_acid_phase(
        list(params.acid_index),
        list(params.antoine_type),
        list(params.antoine_coefficients),
        list(params.antoine_tc),
        list(params.antoine_pc),
        input_to_si(spec, "T", T),
        input_to_si(spec, "P", P),
        list(x),
    )
    return GeVanLaarAcidPhaseResult(
        gamma=tuple(result.gamma),
        ln_gamma=tuple(result.ln_gamma),
        ln_phi=tuple(result.ln_phi),
        p_sat=tuple(from_si(value.magnitude_si, value.unit) for value in result.p_sat),
        warnings=_warnings(result.warnings),
    )


def ge_unifac_phase(params: Any, T: Q, P: Q, x: Sequence[float]) -> GeUnifacPhaseResult:
    """The fugacity coefficients of a UNIFAC liquid, computed in Rust.

    The resolved record crosses flattened, one list per field, in the dataclass's own
    order.
    """
    spec = _models_gen.model("eos.ge_unifac_phase")
    result = _core.ge_unifac_phase(
        list(params.groups),
        list(params.group_r),
        list(params.group_q),
        list(params.aij),
        list(params.antoine_type),
        list(params.antoine_coefficients),
        list(params.antoine_tc),
        list(params.antoine_pc),
        input_to_si(spec, "T", T),
        input_to_si(spec, "P", P),
        list(x),
    )
    return GeUnifacPhaseResult(
        gamma=tuple(result.gamma),
        ln_gamma=tuple(result.ln_gamma),
        ln_phi=tuple(result.ln_phi),
        p_sat=tuple(from_si(value.magnitude_si, value.unit) for value in result.p_sat),
        warnings=_warnings(result.warnings),
    )


def ge_nrtl_flash(params: Any, mixture: Any, T: Q, P: Q, z: Sequence[float]) -> GeNrtlFlashResult:
    """The gamma-phi flash of an SRK vapour over an NRTL liquid, computed in Rust.

    Both objects cross, in the order the stub declares them: the mixture as four
    parallel lists plus the cubic it is evaluated under, then the resolved parameter
    record's own fields. The record's first field is NRTL's ``alpha`` matrix, so the
    cubic's alpha *correlation* crosses beside it as ``cubic_alpha``.

    `beta` crosses as `Option<f64>` and becomes `None`, not a sentinel.
    """
    spec = _models_gen.model("eos.ge_nrtl_flash")
    result = _core.ge_nrtl_flash(
        [c.Tc.to_base_units().magnitude for c in mixture.components],
        [c.Pc.to_base_units().magnitude for c in mixture.components],
        [c.omega for c in mixture.components],
        mixture.flattened_kij(),
        list(params.alpha),
        list(params.dij),
        list(params.antoine_type),
        list(params.antoine_coefficients),
        list(params.antoine_tc),
        list(params.antoine_pc),
        input_to_si(spec, "T", T),
        input_to_si(spec, "P", P),
        list(z),
        mixture.cubic.name,
        mixture.alpha,
        [list(c.alpha_params) for c in mixture.components],
    )
    return GeNrtlFlashResult(
        beta=result.beta,
        x=tuple(result.x),
        y=tuple(result.y),
        k=tuple(result.k),
        ln_phi_liquid=tuple(result.ln_phi_liquid),
        ln_phi_vapour=tuple(result.ln_phi_vapour),
        z_vapour=result.z_vapour,
        min_t_over_tc=result.min_t_over_tc,
        phase=_Phase(result.phase),
        iterations=result.iterations,
        residual=result.residual,
        warnings=_warnings(result.warnings),
    )


def gerg2008_phase(
    components: Sequence[str], T: Q, P: Q, z: Sequence[float]
) -> Gerg2008PhaseResult:
    """The GERG-2008 phase state, computed in Rust."""
    spec = _models_gen.model("eos.gerg2008_phase")
    result = _core.gerg2008_phase(
        list(components), input_to_si(spec, "T", T), input_to_si(spec, "P", P), list(z)
    )
    return Gerg2008PhaseResult(
        z_factor=result.z_factor,
        u=from_si(result.u.magnitude_si, result.u.unit),
        h=from_si(result.h.magnitude_si, result.h.unit),
        s=from_si(result.s.magnitude_si, result.s.unit),
        cv=from_si(result.cv.magnitude_si, result.cv.unit),
        cp=from_si(result.cp.magnitude_si, result.cp.unit),
        g=from_si(result.g.magnitude_si, result.g.unit),
        warnings=_warnings(result.warnings),
    )


def hydrogen_phase(T: Q, P: Q, hydrogen_type: str = "normal") -> HydrogenPhaseResult:
    """The Leachman hydrogen phase state, computed in Rust."""
    spec = _models_gen.model("eos.hydrogen_phase")
    result = _core.hydrogen_phase(
        input_to_si(spec, "T", T), input_to_si(spec, "P", P), hydrogen_type
    )
    return HydrogenPhaseResult(
        z_factor=result.z_factor,
        u=from_si(result.u.magnitude_si, result.u.unit),
        h=from_si(result.h.magnitude_si, result.h.unit),
        s=from_si(result.s.magnitude_si, result.s.unit),
        cv=from_si(result.cv.magnitude_si, result.cv.unit),
        cp=from_si(result.cp.magnitude_si, result.cp.unit),
        g=from_si(result.g.magnitude_si, result.g.unit),
        warnings=_warnings(result.warnings),
    )


def overlay_from(card: Any) -> Any:
    """A keycard as the extension reads it.

    The card's ``pint`` quantities become SI magnitudes here, the way every other
    dimensioned value crosses - the extension takes numbers, and one conversion site is
    what keeps the two languages working in the same units.

    **Nothing is resolved here.** This hands the card over as data; what a name resolves
    to is `azoth_eos::databank`'s answer, and `python/tests/test_data_agreement.py`
    compares it against `azoth.eos.components`' answer for the same card. Two
    implementations of one merge rule, held to each other the way the two kernels are.
    """
    components = [
        (
            name,
            _si_of(parameters.get("Tc")),
            _si_of(parameters.get("Pc")),
            _plain(parameters.get("omega")),
        )
        for name, parameters in sorted(card.components.items())
    ]
    kij = [(a, b, value) for (a, b), value in sorted(card.kij.items())]
    return _core.overlay(components, kij)


def _si_of(value: Q | None) -> float | None:
    """A dimensioned card value as an SI magnitude, or ``None`` if the card omits it."""
    return None if value is None else float(value.to_base_units().magnitude)


def _plain(value: Q | None) -> float | None:
    """A dimensionless card value as a bare float, or ``None`` if the card omits it."""
    return None if value is None else float(value.to("dimensionless").magnitude)


def resolve(calc_id: str) -> Callable[..., Any]:
    """The bridge function for a calc id, derived from the id rather than listed.

    A calculation's id is its address here exactly as it is on the reference side:
    `hydraulics.orifice_flow` names `orifice_flow` in this module. So there is no
    table, and adding a calculation adds no entry.

    The table this replaces was hand-maintained and asserted equal to the registry by
    a test - which is a list kept in step with a list, so the test could only ever
    tell you that you had forgotten something. Deriving it means there is nothing to
    forget, and a calc whose bridge function is missing fails the same completeness
    test one step earlier.

    Raises:
        KeyError: if the id is not in the registry at all. That is a programming
            error rather than a runtime condition - the reference and the extension
            are supposed to cover the same calcs, and a test asserts they do.
    """
    from azoth._models_gen import MODELS
    from azoth._registry_gen import CALCS

    if calc_id not in {entry["id"] for entry in [*CALCS, *MODELS]}:
        raise KeyError(
            f"no Rust implementation for {calc_id!r}: it is not in the registry. "
            f"The extension covers {len(CALCS)} calc(s) and {len(MODELS)} model(s)."
        ) from None

    resolved: Callable[..., Any] = globals()[calc_id.rpartition(".")[2]]
    return resolved
