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

from collections.abc import Callable, Mapping, Sequence
from typing import Any

from azoth import _core, _models_gen
from azoth._registry_gen import spec as _spec_for
from azoth.core.errors import PropertyUnavailableError
from azoth.core.result import (
    AbsorptionColumnResult,
    AmmoniaPhaseResult,
    AntoineVaporPressureResult,
    AqueousViscosityResult,
    ArgonSolidPhaseResult,
    BubblePressureResult,
    BubbleTemperatureResult,
    BwrsPhaseResult,
    CapillaryDewPointResult,
    ChemicalEquilibriumResult,
    ChokedFlowAreaResult,
    ChungConductivityResult,
    ChungViscosityResult,
    Co2PhaseResult,
    Co2WaterDiffusivityResult,
    ColebrookResult,
    ComponentSplitterResult,
    CompressorResult,
    ConductionPlaneWallResult,
    ControlValveCvResult,
    CoolerResult,
    CostaldMolarVolumeResult,
    CriticalPointResult,
    DarcyWeisbachResult,
    DesmukhMatherPhaseResult,
    DewPressureResult,
    DewTemperatureResult,
    DistillationColumnResult,
    EjectorResult,
    FlareResult,
    StirredTankReactorResult,
    Iso6976Result,
    EffectiveDiffusionResult,
    EosCgPhaseResult,
    EquilibriumConstantResult,
    ExpanderResult,
    FilterResult,
    FlowRegime,
    FreezingPointResult,
    FurstElectrolyteMod2004PhaseResult,
    FurstElectrolytePhaseResult,
    GasScrubberResult,
    GeFlashResult,
    GeNrtlFlashResult,
    GeNrtlPhaseResult,
    Gerg2008PhaseResult,
    GeUnifacPhaseResult,
    GeUniquacPhaseResult,
    GeVanLaarAcidPhaseResult,
    GeWilsonPhaseResult,
    HaalandResult,
    HaydukMinhasDiffusivityResult,
    HeaterResult,
    HeatExchangerResult,
    HeatOfVaporizationResult,
    HeliumPhaseResult,
    HybridEosGeFlashResult,
    HydrateEquilibriumLineResult,
    HydrateFormationPressureResult,
    HydrateFormationTemperatureResult,
    HydrateFractionResult,
    HydrateInhibitorConcentrationResult,
    HydrateInhibitorWtResult,
    HydrateStructure,
    HydrogenPhaseResult,
    IapwsHenryLawResult,
    IdealGasCpResult,
    KComponent,
    KentEisenbergPhaseResult,
    KFactorsResult,
    KineticRateLawResult,
    KineticsResult,
    LiquidHeatCapacityResult,
    ManifoldResult,
    MasonSaxenaConductivityResult,
    Matcop5PrumrAlphaResult,
    MatcopAlphaResult,
    MatcopPrAlphaResult,
    MatcopPrumrAlphaResult,
    MatcopPrumrNewAlphaResult,
    MixerResult,
    MolarEnthalpyEntropyResult,
    MollerupAlphaResult,
    NitricSulfuricAcidVaporPressureResult,
    NrtlActivityCoefficientsResult,
    OrificeFlowResult,
    ParachorSurfaceTensionResult,
    ParahydrogenSolidPhaseResult,
    PcsaftRahmatPhaseResult,
    PhFlashResult,
    PipeResult,
    PitzerPhaseResult,
    Pr78KappaResult,
    PrAlphaAbResult,
    PrCpaPhaseResult,
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
    PumpResult,
    PureSaturationResult,
    PvfFlashResult,
    PvFlashResult,
    PvRefluxFlashResult,
    RachfordRiceBinaryResult,
    RachfordRiceResult,
    RackettMolarVolumeResult,
    ReactiveHybridEosGeFlashResult,
    ReactivePhaseEquilibriumResult,
    ReactivePhFlashResult,
    ReactiveTpFlashResult,
    ReferencePotentialsResult,
    ReynoldsNumberResult,
    RkAlphaAbResult,
    RkDepartureResult,
    RootStructure,
    SaftVrMiePhaseResult,
    SaltPrecipitationResult,
    ScaleSaturationRatioResult,
    SchwartzentruberAlphaResult,
    SeparatorResult,
    ShortcutDistillationColumnResult,
    SiddiqiLucasDiffusivityResult,
    SolidFugacityResult,
    SoreideWhitsonAlphaResult,
    SoreideWhitsonPhaseResult,
    SplitterResult,
    SrkAlphaAbResult,
    SrkCpaPhaseResult,
    SrkDepartureResult,
    SrkKappaResult,
    SrkPenelouxShiftResult,
    SrkZFactorResult,
    StabilityTestResult,
    SwameeJainResult,
    TankResult,
    TbpFractionPropertiesResult,
    ThermalConductivityResult,
    ThFlashResult,
    ThreePhaseSeparatorResult,
    ThrottlingValveResult,
    TpFlashSaftResult,
    TpMultiflashResult,
    TpMultiflashWaxResult,
    TpSolidFlashResult,
    TsFlashResult,
    TuFlashResult,
    TvFlashResult,
    TvFractionFlashResult,
    TwucoonAlphaResult,
    TwucoonParamAlphaResult,
    TwucoonStatoilAlphaResult,
    TwuKappaResult,
    TynCalusDiffusivityResult,
    UmrCpaPhaseResult,
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
    WaxSolidFugacityResult,
    WilkeChangDiffusivityResult,
    WilkeViscosityResult,
    WilsonActivityCoefficientsResult,
)
from azoth.core.result import HenryStatus as _HenryStatus
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


def _association_spec(mixture: Any) -> Any:
    """The mixture's association, in the form the Rust boundary takes.

    **Every model whose Python side takes a ``Mixture`` crosses this**, and none of them
    may default it away. A mixture whose association does not cross is a *different
    fluid* that converges: ``eos.pt_flash`` sent nine arguments and no association, so
    the Rust backend ran a classical SRK flash on a fluid carrying the CPA interaction
    column and returned ``all_liquid`` where the associating model splits at
    ``beta = 0.208383589``. The twelve numbers below are the fields of
    :class:`AssociationParameters` after the scheme, in the order Rust's
    ``AssociationRecord`` declares them, in the internal scale the table states them in.
    """
    schemes: list[str] = []
    values: list[list[float]] = []
    for component in mixture.components:
        record = component.association
        if record is None:
            schemes.append("")
            values.append([0.0] * 12)
            continue
        schemes.append(record.scheme)
        values.append(
            [
                float(record.sites),
                record.energy,
                record.volume_srk,
                record.a_srk,
                record.b_srk,
                record.m_srk,
                record.volume_pr,
                record.a_pr,
                record.b_pr,
                record.m_pr,
                record.racket_z,
                record.volume_correction,
            ]
        )
    return _core.AssociationSpec(bool(mixture.associating), schemes, values)


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


def equilibrium_constant(source: str, reaction: str, T: Q) -> EquilibriumConstantResult:
    """One reaction's equilibrium constant, computed in Rust.

    The reaction's name and its source cross unresolved, so the Rust side reads the
    table itself - which source, which row, which coefficients. No coefficient reaches
    Python, and the two languages cannot disagree about which row answered.
    """
    spec = _spec_for("reactions.equilibrium_constant")
    result = _core.equilibrium_constant(source, reaction, input_to_si(spec, "T", T))
    return EquilibriumConstantResult(
        ln_k=result.ln_k,
        k=result.k,
        ln_k_derivative=from_si(result.ln_k_derivative.magnitude_si, result.ln_k_derivative.unit),
        reaction_heat=from_si(result.reaction_heat.magnitude_si, result.reaction_heat.unit),
        reference=result.reference,
        warnings=_warnings(result.warnings),
    )


def reference_potentials(components: list[str], source: str, T: Q) -> ReferencePotentialsResult:
    """The standard-state reference potentials, computed in Rust.

    The component names cross unresolved and in the caller's order, which is the order the
    potentials come back in: the Rust side reads the reaction tables, chooses the basis
    and propagates, so no stoichiometric coefficient and no rank decision reaches Python.
    """
    # A model, not a calc: its spec is in the model registry, and `input_to_si` reads the
    # same declarations either way.
    spec = _models_gen.model("reactions.reference_potentials")
    result = _core.reference_potentials(list(components), source, input_to_si(spec, "T", T))
    return ReferencePotentialsResult(
        potentials=tuple(from_si(value.magnitude_si, value.unit) for value in result.potentials),
        independent=tuple(result.independent),
        survivors=tuple(result.survivors),
        rank=result.rank,
        warnings=_warnings(result.warnings),
    )


def _si(spec: dict[str, object], name: str, value: float | Q) -> float:
    """A declared input as its SI magnitude, quantity or not.

    **The boundary is mixed by design.** A vector whose spec declares a unit arrives as a
    pint quantity and has to be converted; a dimensionless one - ``chem_ref``,
    ``log_activity``, the element matrix - arrives as the bare number it is, and
    ``input_to_si`` refuses that rather than passing it through.
    """
    # A quantity or a bare number. `Q` is a type alias rather than a class, so the check
    # is made the other way round: the numeric branch is the one `isinstance` can name.
    if isinstance(value, int | float):
        return float(value)
    return input_to_si(spec, name, value)


def chemical_equilibrium(
    a_matrix: list[list[float]],
    b: list[float],
    whole_system: bool,
    moles: list[float],
    chem_ref: list[float],
    log_activity: list[float],
    T: Q,
    max_iterations: float,
    tolerance: float,
    concentration_basis: str,
    solvent_weight: Q,
    solvent_mask: Sequence[float],
    phase_moles: Q,
) -> ChemicalEquilibriumResult:
    """The reactive equilibrium solve, computed in Rust.

    The matrix crosses nested, so every argument here is one the spec declares.
    """
    spec = _models_gen.model("reactions.chemical_equilibrium")
    # The unit-carrying vectors cross as SI magnitudes, element by element, for the same
    # reason the reference converts them: a model's declared units are the boundary's.
    a_matrix = [[_si(spec, "a_matrix", value) for value in row] for row in a_matrix]
    b = [_si(spec, "b", value) for value in b]
    moles = [_si(spec, "moles", value) for value in moles]
    chem_ref = [_si(spec, "chem_ref", value) for value in chem_ref]
    log_activity = [_si(spec, "log_activity", value) for value in log_activity]
    result = _core.chemical_equilibrium(
        a_matrix,
        list(b),
        whole_system,
        list(moles),
        list(chem_ref),
        list(log_activity),
        input_to_si(spec, "T", T),
        int(max_iterations),
        float(tolerance),
        concentration_basis,
        input_to_si(spec, "solvent_weight", solvent_weight),
        [_si(spec, "solvent_mask", value) for value in solvent_mask],
        input_to_si(spec, "phase_moles", phase_moles),
    )
    return ChemicalEquilibriumResult(
        moles=tuple(from_si(value.magnitude_si, value.unit) for value in result.moles),
        iterations=result.iterations,
        error=result.error,
        converged=result.converged,
        warnings=_warnings(result.warnings),
    )


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


def pr_peneloux_shift(
    omega: float, Tc: Q, Pc: Q, z_ra: float | None = None
) -> PrPenelouxShiftResult:
    """The Peng-Robinson Peneloux volume-translation parameter, computed in Rust."""
    spec = _spec_for("eos.pr_peneloux_shift")
    result = _core.pr_peneloux_shift(
        omega,
        input_to_si(spec, "Tc", Tc),
        input_to_si(spec, "Pc", Pc),
        # Dimensionless, and `None` means the fallback correlation - which is what the
        # spec's optional input documents.
        z_ra,
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
        _association_spec(mixture),
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


def aqueous_viscosity(mixture: Any, T: Q, P: Q, z: Sequence[float]) -> AqueousViscosityResult:
    """The liquid viscosity NeqSim gives an aqueous phase, computed in Rust.

    The mixture's per-component data crosses as vectors, as `eos.viscosity`'s does, plus the
    two the correlation reads: `liqvisc` flattened and the model number beside it.
    """
    spec = _models_gen.model("eos.aqueous_viscosity")
    molar_mass = []
    for c in mixture.components:
        if c.molar_mass is None:
            raise PropertyUnavailableError(
                "component", "molar mass", "a card-added component needs its own molar mass"
            )
        molar_mass.append(c.molar_mass.to_base_units().magnitude)
    result = _core.aqueous_viscosity(
        [c.Tc.to_base_units().magnitude for c in mixture.components],
        [c.Pc.to_base_units().magnitude for c in mixture.components],
        [c.omega for c in mixture.components],
        mixture.flattened_kij(),
        _association_spec(mixture),
        molar_mass,
        [value for c in mixture.components for value in c.liqvisc],
        [c.liqvisc_model for c in mixture.components],
        input_to_si(spec, "T", T),
        input_to_si(spec, "P", P),
        list(z),
        mixture.cubic.name,
        mixture.alpha,
        [list(c.alpha_params) for c in mixture.components],
    )
    return AqueousViscosityResult(
        viscosity=from_si(result.viscosity.magnitude_si, result.viscosity.unit),
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
        _association_spec(mixture),
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
        _association_spec(mixture),
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
        _association_spec(mixture),
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
        _association_spec(mixture),
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
        _association_spec(mixture),
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
        _association_spec(mixture),
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
        _association_spec(mixture),
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
        _association_spec(mixture),
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
        _association_spec(mixture),
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
        _association_spec(mixture),
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
        _association_spec(mixture),
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
        _association_spec(mixture),
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
        _association_spec(mixture),
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
        _association_spec(mixture),
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
        _association_spec(mixture),
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
        _association_spec(mixture),
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
        _association_spec(mixture),
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
        _association_spec(mixture),
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
        _association_spec(mixture),
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
        _association_spec(mixture),
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
        _association_spec(mixture),
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
        _association_spec(mixture),
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
        _association_spec(mixture),
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
        _association_spec(mixture),
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
        _association_spec(mixture),
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
        _association_spec(mixture),
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
        _association_spec(mixture),
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
        _association_spec(mixture),
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
        _association_spec(mixture),
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


def desmukh_mather_phase(
    components: Sequence[str], T: Q, P: Q, x: Sequence[float]
) -> DesmukhMatherPhaseResult:
    """The activity coefficients of a Desmukh-Mather phase, computed in Rust.

    Only the names and the state cross: which component is the solvent, which pairs carry
    parameters and which branch each component's fugacity coefficient takes are all read
    from the Rust side's own databank.
    """
    spec = _models_gen.model("eos.desmukh_mather_phase")
    result = _core.desmukh_mather_phase(
        list(components),
        input_to_si(spec, "T", T),
        input_to_si(spec, "P", P),
        list(x),
    )
    return DesmukhMatherPhaseResult(
        gamma=tuple(result.gamma),
        ln_gamma=tuple(result.ln_gamma),
        molality=tuple(result.molality),
        ionic_strength=result.ionic_strength,
        solvent_molar_mass=result.solvent_molar_mass,
        ln_phi=tuple(result.ln_phi),
        warnings=_warnings(result.warnings),
    )


def kent_eisenberg_phase(
    components: Sequence[str], T: Q, P: Q, x: Sequence[float]
) -> KentEisenbergPhaseResult:
    """The fugacity coefficients of a Kent-Eisenberg phase, computed in Rust.

    Only the names and the state cross: the Rust side resolves each component's reference
    state, charge and correlations against its own copy of the databank, because which
    branch a component takes is what the model is.
    """
    spec = _models_gen.model("eos.kent_eisenberg_phase")
    result = _core.kent_eisenberg_phase(
        list(components),
        input_to_si(spec, "T", T),
        input_to_si(spec, "P", P),
        list(x),
    )
    return KentEisenbergPhaseResult(
        gamma=tuple(result.gamma),
        ln_gamma=tuple(result.ln_gamma),
        ln_phi=tuple(result.ln_phi),
        warnings=_warnings(result.warnings),
    )


def pitzer_phase(components: Sequence[str], T: Q, P: Q, x: Sequence[float]) -> PitzerPhaseResult:
    """The activity and fugacity coefficients of a Pitzer electrolyte phase, computed in Rust.

    Only the names, the state and the composition cross: the Rust side resolves each
    component's charge, molar mass and reference state against its own copy of the
    databank, because the parameter datasets are keyed by the name and the model reads
    them itself.
    """
    spec = _models_gen.model("eos.pitzer_phase")
    result = _core.pitzer_phase(
        list(components),
        input_to_si(spec, "T", T),
        input_to_si(spec, "P", P),
        list(x),
    )
    return PitzerPhaseResult(
        gamma=tuple(result.gamma),
        ln_gamma=tuple(result.ln_gamma),
        ln_phi=tuple(result.ln_phi),
        henry=tuple(from_si(value.magnitude_si, value.unit) for value in result.henry),
        gamma_inf=tuple(result.gamma_inf),
        molality=tuple(result.molality),
        ionic_strength=result.ionic_strength,
        osmotic_coefficient=result.osmotic_coefficient,
        water_activity=result.water_activity,
        dataset=result.dataset,
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
        _association_spec(mixture),
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


def hydrate_formation_temperature(
    components: Sequence[str],
    P: Q,
    z: Sequence[float],
    eos: str = "srk",
    hydrate_model: str = "pvtsim",
) -> HydrateFormationTemperatureResult:
    """The hydrate formation temperature of a fluid, computed in Rust.

    **The component names cross unresolved**, and the Rust side resolves them: the hydrate's
    guest tables are keyed by name, so a mixture built from constants alone could not carry
    them.
    """
    spec = _models_gen.model("eos.hydrate_formation_temperature")
    result = _core.hydrate_formation_temperature(
        list(components), input_to_si(spec, "P", P), list(z), eos, hydrate_model
    )
    return HydrateFormationTemperatureResult(
        temperature=from_si(result.temperature.magnitude_si, result.temperature.unit),
        structure=HydrateStructure(result.structure),
        iterations=result.iterations,
        residual=result.residual,
        warnings=_warnings(result.warnings),
    )


def hydrate_fraction(
    components: Sequence[str],
    T: Q,
    P: Q,
    z: Sequence[float],
    eos: str = "srk",
    hydrate_model: str = "pvtsim",
) -> HydrateFractionResult:
    """The fraction of a feed that is hydrate at a state, computed in Rust.

    **The component names cross unresolved**, and the Rust side resolves them: the hydrate's
    guest tables are keyed by name, so a mixture built from constants alone could not carry
    them.
    """
    spec = _models_gen.model("eos.hydrate_fraction")
    result = _core.hydrate_fraction(
        list(components),
        input_to_si(spec, "T", T),
        input_to_si(spec, "P", P),
        list(z),
        eos,
        hydrate_model,
    )
    return HydrateFractionResult(
        beta=result.beta,
        structure=HydrateStructure(result.structure),
        balance_error=result.balance_error,
        iterations=result.iterations,
        residual=result.residual,
        warnings=_warnings(result.warnings),
    )


def hydrate_equilibrium_line(
    components: Sequence[str],
    P_min: Q,
    P_max: Q,
    z: Sequence[float],
    eos: str = "srk",
    hydrate_model: str = "pvtsim",
) -> HydrateEquilibriumLineResult:
    """A hydrate curve over a pressure grid, computed in Rust.

    **The component names cross unresolved**, and the Rust side resolves them: the hydrate's
    guest tables are keyed by name, so a mixture built from constants alone could not carry
    them.
    """
    spec = _models_gen.model("eos.hydrate_equilibrium_line")
    result = _core.hydrate_equilibrium_line(
        list(components),
        input_to_si(spec, "P_min", P_min),
        input_to_si(spec, "P_max", P_max),
        list(z),
        eos,
        hydrate_model,
    )
    return HydrateEquilibriumLineResult(
        temperature=tuple(result.temperature),
        pressure=tuple(result.pressure),
        warnings=_warnings(result.warnings),
    )


def hydrate_inhibitor_wt(
    components: Sequence[str],
    moles: Sequence[float],
    inhibitor: str,
    wt_target: float,
    T: Q,
    P: Q,
    eos: str = "srk",
) -> HydrateInhibitorWtResult:
    """The inhibitor dose that reaches a target aqueous mass fraction, computed in Rust.

    **The feed crosses in moles**, as for its sibling: the secant adds an absolute amount to
    the inhibitor's entry.
    """
    spec = _models_gen.model("eos.hydrate_inhibitor_wt")
    result = _core.hydrate_inhibitor_wt(
        list(components),
        # A unit-bearing vector crosses as quantities; the kernel takes SI moles.
        [to_si(value, "mol", "moles") for value in moles],
        inhibitor,
        wt_target,
        input_to_si(spec, "T", T),
        input_to_si(spec, "P", P),
        eos,
    )
    return HydrateInhibitorWtResult(
        inhibitor_moles=result.inhibitor_moles,
        weight_fraction=result.weight_fraction,
        phases=result.phases,
        iterations=result.iterations,
        residual=result.residual,
        warnings=_warnings(result.warnings),
    )


def hydrate_inhibitor_concentration(
    components: Sequence[str],
    moles: Sequence[float],
    inhibitor: str,
    T_target: Q,
    P: Q,
    eos: str = "srk",
    hydrate_model: str = "pvtsim",
) -> HydrateInhibitorConcentrationResult:
    """The moles of inhibitor that hold a hydrate temperature down, computed in Rust.

    **The feed crosses in moles**, unlike every other hydrate model here: the secant adds an
    absolute amount to the inhibitor's entry, so a normalised feed would reproduce the
    equation and not the path.
    """
    spec = _models_gen.model("eos.hydrate_inhibitor_concentration")
    result = _core.hydrate_inhibitor_concentration(
        list(components),
        # A unit-bearing vector crosses as quantities; the kernel takes SI moles.
        [to_si(value, "mol", "moles") for value in moles],
        inhibitor,
        input_to_si(spec, "T_target", T_target),
        input_to_si(spec, "P", P),
        eos,
        hydrate_model,
    )
    return HydrateInhibitorConcentrationResult(
        inhibitor_moles=result.inhibitor_moles,
        weight_fraction=result.weight_fraction,
        hydrate_temperature=from_si(
            result.hydrate_temperature.magnitude_si, result.hydrate_temperature.unit
        ),
        iterations=result.iterations,
        residual=result.residual,
        warnings=_warnings(result.warnings),
    )


def hydrate_formation_pressure(
    components: Sequence[str],
    T: Q,
    z: Sequence[float],
    eos: str = "srk",
    hydrate_model: str = "pvtsim",
) -> HydrateFormationPressureResult:
    """The hydrate formation pressure of a fluid, computed in Rust.

    **The component names cross unresolved**, and the Rust side resolves them: the hydrate's
    guest tables are keyed by name, so a mixture built from constants alone could not carry
    them.
    """
    spec = _models_gen.model("eos.hydrate_formation_pressure")
    result = _core.hydrate_formation_pressure(
        list(components), input_to_si(spec, "T", T), list(z), eos, hydrate_model
    )
    return HydrateFormationPressureResult(
        pressure=from_si(result.pressure.magnitude_si, result.pressure.unit),
        structure=HydrateStructure(result.structure),
        iterations=result.iterations,
        residual=result.residual,
        warnings=_warnings(result.warnings),
    )


def tbp_fraction_properties(molar_mass: Q, density: Q) -> TbpFractionPropertiesResult:
    """A TBP pseudo-component's properties, computed in Rust."""
    result = _core.tbp_fraction_properties(
        float(molar_mass.to("kg/mol").magnitude), float(density.to("kg/m**3").magnitude)
    )
    return TbpFractionPropertiesResult(
        tc=from_si(result.tc.magnitude_si, result.tc.unit),
        pc=from_si(result.pc.magnitude_si, result.pc.unit),
        boiling_temperature=from_si(
            result.boiling_temperature.magnitude_si, result.boiling_temperature.unit
        ),
        acentric_factor=result.acentric_factor,
        attraction_exponent=result.attraction_exponent,
        warnings=_warnings(result.warnings),
    )


def wax_solid_fugacity(
    molar_mass: Q,
    tc: Q,
    pc: Q,
    omega: float,
    heat_of_fusion: Q,
    triple_point_temperature: Q,
    T: Q,
    P: Q,
    eos: str = "srk",
) -> WaxSolidFugacityResult:
    """A wax cut's solid fugacity coefficient, computed in Rust."""
    spec = _spec_for("eos.wax_solid_fugacity")
    result = _core.wax_solid_fugacity(
        input_to_si(spec, "molar_mass", molar_mass),
        input_to_si(spec, "tc", tc),
        input_to_si(spec, "pc", pc),
        omega,
        input_to_si(spec, "heat_of_fusion", heat_of_fusion),
        input_to_si(spec, "triple_point_temperature", triple_point_temperature),
        input_to_si(spec, "T", T),
        input_to_si(spec, "P", P),
        eos,
    )
    return WaxSolidFugacityResult(
        fugacity_coefficient=result.fugacity_coefficient,
        warnings=_warnings(result.warnings),
    )


def tp_multiflash_wax(
    components: Sequence[str], T: Q, P: Q, z: Sequence[float], eos: str = "srk"
) -> TpMultiflashWaxResult:
    """The wax fraction of a feed, computed in Rust.

    **The component names cross unresolved**, and the Rust side resolves them: the wax flag
    and the melt data are the databank's own columns.
    """
    spec = _models_gen.model("eos.tp_multiflash_wax")
    result = _core.tp_multiflash_wax(
        list(components),
        input_to_si(spec, "T", T),
        input_to_si(spec, "P", P),
        list(z),
        eos,
    )
    return TpMultiflashWaxResult(
        wax_fraction=result.wax_fraction,
        phase_count=result.phase_count,
        beta=tuple(result.beta),
        x=tuple(tuple(row) for row in result.x),
        iterations=result.iterations,
        residual=result.residual,
        converged=result.converged,
        warnings=_warnings(result.warnings),
    )


def tp_solid_flash(
    components: Sequence[str], solid: str, T: Q, P: Q, z: Sequence[float], eos: str = "srk"
) -> TpSolidFlashResult:
    """The solid fraction of a feed, computed in Rust.

    **The component names cross unresolved**, and the Rust side resolves them: the melt data
    and the density correlations are the databank's own columns.
    """
    spec = _models_gen.model("eos.tp_solid_flash")
    result = _core.tp_solid_flash(
        list(components),
        solid,
        input_to_si(spec, "T", T),
        input_to_si(spec, "P", P),
        list(z),
        eos,
    )
    return TpSolidFlashResult(
        solid_fraction=result.solid_fraction,
        phase_count=result.phase_count,
        beta=tuple(result.beta),
        x=tuple(tuple(row) for row in result.x),
        solid_fugacity_coefficient=result.solid_fugacity_coefficient,
        iterations=result.iterations,
        residual=result.residual,
        converged=result.converged,
        warnings=_warnings(result.warnings),
    )


def scale_saturation_ratio(
    salt: str,
    x1: float,
    x2: float,
    x_water: float,
    gamma1: float,
    gamma2: float,
    water_activity: float,
    T: Q,
    P: Q,
    h3o_molality: Q | None = None,
) -> ScaleSaturationRatioResult:
    """One salt's saturation ratio, computed in Rust."""
    spec = _spec_for("eos.scale_saturation_ratio")
    result = _core.scale_saturation_ratio(
        salt,
        x1,
        x2,
        x_water,
        gamma1,
        gamma2,
        water_activity,
        input_to_si(spec, "T", T),
        input_to_si(spec, "P", P),
        # By keyword: the boundary's own signature puts the defaulted one last, which is
        # the spec's order and not the position it has there.
        h3o_molality=(
            None if h3o_molality is None else input_to_si(spec, "h3o_molality", h3o_molality)
        ),
    )
    return ScaleSaturationRatioResult(
        saturation_ratio=result.saturation_ratio,
        ion_activity_product=result.ion_activity_product,
        solubility_product=result.solubility_product,
        warnings=_warnings(result.warnings),
    )


def salt_precipitation(
    components: Sequence[str], salt: str, T: Q, P: Q, z: Sequence[float]
) -> SaltPrecipitationResult:
    """The solid one mineral takes from a brine, computed in Rust.

    **The component names cross unresolved**, and the Rust side resolves them: the brine's
    coefficients come from the Pitzer phase over the same names.
    """
    spec = _models_gen.model("eos.salt_precipitation")
    result = _core.salt_precipitation(
        list(components),
        list(z),
        salt,
        input_to_si(spec, "T", T),
        input_to_si(spec, "P", P),
    )
    return SaltPrecipitationResult(
        precipitated_moles=result.precipitated_moles,
        initial_saturation_ratio=result.initial_saturation_ratio,
        final_saturation_ratio=result.final_saturation_ratio,
        iterations=result.iterations,
        extent_of_maximum=result.extent_of_maximum,
        warnings=_warnings(result.warnings),
    )


def solid_fugacity(
    heat_of_fusion: Q,
    triple_point_temperature: Q,
    delta_cp_sl: Q,
    delta_solid_volume: Q,
    tc: Q,
    pc: Q,
    omega: float,
    T: Q,
    P: Q,
    eos: str = "srk",
) -> SolidFugacityResult:
    """A pure solid's fugacity coefficient, computed in Rust."""
    spec = _spec_for("eos.solid_fugacity")
    result = _core.solid_fugacity(
        input_to_si(spec, "heat_of_fusion", heat_of_fusion),
        input_to_si(spec, "triple_point_temperature", triple_point_temperature),
        input_to_si(spec, "delta_cp_sl", delta_cp_sl),
        input_to_si(spec, "delta_solid_volume", delta_solid_volume),
        input_to_si(spec, "tc", tc),
        input_to_si(spec, "pc", pc),
        omega,
        input_to_si(spec, "T", T),
        input_to_si(spec, "P", P),
        eos,
    )
    return SolidFugacityResult(
        fugacity_coefficient=result.fugacity_coefficient,
        warnings=_warnings(result.warnings),
    )


def freezing_point(
    components: Sequence[str], z: Sequence[float], solid: str, P: Q
) -> FreezingPointResult:
    """A fluid's freezing-point temperature, computed in Rust.

    **The component names cross unresolved**, and the Rust side resolves them: the tabulated
    solid reads the databank's melting point, heat of fusion and density correlations, and the
    para-hydrogen route reads the reference equation on that side.
    """
    spec = _models_gen.model("eos.freezing_point")
    result = _core.freezing_point(list(components), list(z), solid, input_to_si(spec, "P", P))
    return FreezingPointResult(
        temperature=from_si(result.temperature.magnitude_si, result.temperature.unit),
        component=result.component,
        iterations=result.iterations,
        residual=result.residual,
        warnings=_warnings(result.warnings),
    )


def hybrid_eos_ge_flash(
    components: Sequence[str], cubic: str, T: Q, P: Q, moles: Sequence[float]
) -> HybridEosGeFlashResult:
    """The fixed-role gas-oil-brine flash, computed in Rust.

    The names cross **unresolved**: the two EoS roles' constants, the seeding's classes and
    the brine's ion mask all resolve on the Rust side from the same databank, so the two
    languages cannot disagree about which substance is which.
    """
    spec = _models_gen.model("eos.hybrid_eos_ge_flash")
    result = _core.hybrid_eos_ge_flash(
        list(components),
        cubic,
        input_to_si(spec, "T", T),
        input_to_si(spec, "P", P),
        # The unit-carrying vectors cross as SI magnitudes element by element, for the same
        # reason the reference converts them: a model's declared units are the boundary's.
        [_si(spec, "moles", value) for value in moles],
    )
    return HybridEosGeFlashResult(
        beta=tuple(result.beta),
        x=tuple(tuple(row) for row in result.x),
        ln_phi=tuple(tuple(row) for row in result.ln_phi),
        iterations=result.iterations,
        residual=result.residual,
        max_material_balance_residual=result.max_material_balance_residual,
        max_log_fugacity_residual=result.max_log_fugacity_residual,
        min_t_over_tc=result.min_t_over_tc,
        warnings=_warnings(result.warnings),
    )


def iapws_henry_law(gas: str, T: Q) -> IapwsHenryLawResult:
    """The Henry constant of a gas in water, computed in Rust.

    ``gas`` crosses as the table's own spelling, because the boundary a component name
    crosses is the table's lookup and not this one's: `ch4` is a row and `methane` is a
    name for it, and which of the two a caller has is the caller's business.
    """
    spec = _spec_for("eos.iapws_henry_law")
    result = _core.iapws_henry_law(gas, input_to_si(spec, "T", T))
    return IapwsHenryLawResult(
        henry=from_si(result.henry.magnitude_si, result.henry.unit),
        ln_henry=result.ln_henry,
        d_ln_henry_d_t=result.d_ln_henry_d_t,
        status=_HenryStatus(result.status),
        rms_log_residual=result.rms_log_residual,
        warnings=_warnings(result.warnings),
    )


def hydrogen_phase(
    T: Q, P: Q, hydrogen_type: str = "normal", compressed_phase: str = "vapour"
) -> HydrogenPhaseResult:
    """The Leachman hydrogen phase state, computed in Rust."""
    spec = _models_gen.model("eos.hydrogen_phase")
    result = _core.hydrogen_phase(
        input_to_si(spec, "T", T), input_to_si(spec, "P", P), hydrogen_type, compressed_phase
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
    components: list[_core.ComponentArguments] = [
        {
            "name": name,
            "tc": _si_of(stated.parameters.get("Tc")),
            "pc": _si_of(stated.parameters.get("Pc")),
            "omega": _plain(stated.parameters.get("omega")),
            "ion": stated.is_ion,
            "ionic_charge": _plain(stated.parameters.get("ionic_charge")),
            # Metres, which is the card's unit and the merge's input. The crossing to the
            # databank's ångström happens in one place, on the Rust side, so there is no
            # second conversion here to disagree with it.
            "deshmukh_mather_diameter": _si_of(stated.parameters.get("deshmukh_mather_diameter")),
            "dielectric": _dielectric(stated.parameters),
        }
        for name, stated in sorted(card.components.items())
    ]
    kij = [(a, b, value) for (a, b), value in sorted(card.kij.items())]
    return _core.overlay(components, kij)


def _si_of(value: Q | None) -> float | None:
    """A dimensioned card value as an SI magnitude, or ``None`` if the card omits it."""
    return None if value is None else float(value.to_base_units().magnitude)


def _dielectric(parameters: Mapping[str, Q]) -> list[float] | None:
    """The five coefficients the card states, each in its own unit, or ``None``."""
    names = ("dielectric_1", "dielectric_2", "dielectric_3", "dielectric_4", "dielectric_5")
    if not all(name in parameters for name in names):
        return None
    return [
        float(parameters["dielectric_1"].to("dimensionless").magnitude),
        float(parameters["dielectric_2"].to("K").magnitude),
        float(parameters["dielectric_3"].to("1/K").magnitude),
        float(parameters["dielectric_4"].to("1/K**2").magnitude),
        float(parameters["dielectric_5"].to("1/K**3").magnitude),
    ]


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


def tp_flash_saft(components: Sequence[str], T: Q, P: Q, z: Sequence[float]) -> TpFlashSaftResult:
    """The SAFT-VR-Mie flash, computed in Rust.

    The component names cross **unresolved**, so the Rust side resolves the fluid itself -
    the Mie set and the cubic constants the Wilson seed is built from.
    """
    spec = _models_gen.model("eos.tp_flash_saft")
    result = _core.tp_flash_saft(
        list(components),
        input_to_si(spec, "T", T),
        input_to_si(spec, "P", P),
        list(z),
    )
    return TpFlashSaftResult(
        beta=result.beta,
        x=tuple(result.x),
        y=tuple(result.y),
        k=tuple(result.k),
        ln_phi_liquid=tuple(result.ln_phi_liquid),
        ln_phi_vapour=tuple(result.ln_phi_vapour),
        z_liquid=result.z_liquid,
        z_vapour=result.z_vapour,
        phase=_Phase(result.phase),
        iterations=result.iterations,
        residual=result.residual,
        warnings=_warnings(result.warnings),
    )


def saft_vr_mie_phase(
    components: Sequence[str], T: Q, P: Q, z: Sequence[float], compressed_phase: str
) -> SaftVrMiePhaseResult:
    """The SAFT-VR-Mie phase state, computed in Rust.

    The component names cross **unresolved**, so the Rust side resolves the fluid itself -
    the five Mie columns, whose absence the table spells as zeros in three of them.
    """
    spec = _models_gen.model("eos.saft_vr_mie_phase")
    result = _core.saft_vr_mie_phase(
        list(components),
        input_to_si(spec, "T", T),
        input_to_si(spec, "P", P),
        list(z),
        compressed_phase,
    )
    return SaftVrMiePhaseResult(
        z_factor=result.z_factor,
        ln_phi=tuple(result.ln_phi),
        v=from_si(result.v.magnitude_si, result.v.unit),
        h_res=from_si(result.h_res.magnitude_si, result.h_res.unit),
        s_res=from_si(result.s_res.magnitude_si, result.s_res.unit),
        warnings=_warnings(result.warnings),
    )


def pcsaft_rahmat_phase(
    components: Sequence[str], T: Q, P: Q, z: Sequence[float], compressed_phase: str
) -> PcsaftRahmatPhaseResult:
    """The PC-SAFT phase state, computed in Rust.

    The component names cross **unresolved**, as they do for both CPA twins, so the Rust
    side resolves the fluid itself - the `mSAFT`/`sigmaSAFT`/`epsikSAFT` set, whose absence
    the table spells as zeros in all three columns, and the `KIJPCSAFT` column.
    """
    spec = _models_gen.model("eos.pcsaft_rahmat_phase")
    result = _core.pcsaft_rahmat_phase(
        list(components),
        input_to_si(spec, "T", T),
        input_to_si(spec, "P", P),
        list(z),
        compressed_phase,
    )
    return PcsaftRahmatPhaseResult(
        z_factor=result.z_factor,
        ln_phi=tuple(result.ln_phi),
        v=from_si(result.v.magnitude_si, result.v.unit),
        h_res=from_si(result.h_res.magnitude_si, result.h_res.unit),
        s_res=from_si(result.s_res.magnitude_si, result.s_res.unit),
        warnings=_warnings(result.warnings),
    )


def pr_cpa_phase(
    components: Sequence[str], T: Q, P: Q, z: Sequence[float], compressed_phase: str
) -> PrCpaPhaseResult:
    """The PR-CPA phase state, computed in Rust.

    The component names cross **unresolved**, as they do for the SRK twin, so the Rust side
    resolves the fluid itself and reads the PR family's fitted set and its `cpakij_PR`
    column - the first shipped model to read either.
    """
    spec = _models_gen.model("eos.pr_cpa_phase")
    result = _core.pr_cpa_phase(
        list(components),
        input_to_si(spec, "T", T),
        input_to_si(spec, "P", P),
        list(z),
        compressed_phase,
    )
    return PrCpaPhaseResult(
        z_factor=result.z_factor,
        ln_phi=tuple(result.ln_phi),
        h_res=from_si(result.h_res.magnitude_si, result.h_res.unit),
        s_res=from_si(result.s_res.magnitude_si, result.s_res.unit),
        warnings=_warnings(result.warnings),
    )


def umr_cpa_phase(
    components: Sequence[str], T: Q, P: Q, z: Sequence[float], compressed_phase: str
) -> UmrCpaPhaseResult:
    """The UMR-CPA phase state, computed in Rust.

    The component names cross **unresolved**, as they do for the other CPA models, so the
    Rust side resolves the UMR-CPA parameter set, the `UMRCPA_MC1..5` coefficients and the
    UNIFAC group decomposition itself - which is what makes the two-kernel comparison cover
    the resolution as well as the arithmetic.
    """
    spec = _models_gen.model("eos.umr_cpa_phase")
    result = _core.umr_cpa_phase(
        list(components),
        input_to_si(spec, "T", T),
        input_to_si(spec, "P", P),
        list(z),
        compressed_phase,
    )
    return UmrCpaPhaseResult(
        z_factor=result.z_factor,
        ln_phi=tuple(result.ln_phi),
        h_res=from_si(result.h_res.magnitude_si, result.h_res.unit),
        s_res=from_si(result.s_res.magnitude_si, result.s_res.unit),
        warnings=_warnings(result.warnings),
    )


def soreide_whitson_phase(
    components: Sequence[str],
    T: Q,
    P: Q,
    x: Sequence[float],
    salinity: Q,
    compressed_phase: str,
) -> SoreideWhitsonPhaseResult:
    """The Soreide-Whitson phase state, computed in Rust.

    The component names cross **unresolved**, and the Rust side resolves each one's role in
    the aqueous correlation and reads its row of `KIJWhitsonSoriede` itself. Both are name
    lookups on NeqSim's side too, so a caller could not state them as numbers without
    reimplementing the model.
    """
    spec = _models_gen.model("eos.soreide_whitson_phase")
    result = _core.soreide_whitson_phase(
        list(components),
        input_to_si(spec, "T", T),
        input_to_si(spec, "P", P),
        list(x),
        input_to_si(spec, "salinity", salinity),
        compressed_phase,
    )
    return SoreideWhitsonPhaseResult(
        z_factor=result.z_factor,
        ln_phi=tuple(result.ln_phi),
        warnings=_warnings(result.warnings),
    )


def furst_electrolyte_phase(
    components: Sequence[str],
    T: Q,
    P: Q,
    x: Sequence[float],
    compressed_phase: str,
) -> FurstElectrolytePhaseResult:
    """The Fürst electrolyte phase state, computed in Rust.

    The component names cross **unresolved**, and the Rust side resolves each one's
    dielectric coefficients, fitted covolume, Schwartzentruber parameters and short-range
    pair table. The salt is a component and not a scalar, so an ion crosses like any name.
    """
    spec = _models_gen.model("eos.furst_electrolyte_phase")
    result = _core.furst_electrolyte_phase(
        list(components),
        input_to_si(spec, "T", T),
        input_to_si(spec, "P", P),
        list(x),
        compressed_phase,
    )
    return FurstElectrolytePhaseResult(
        z_factor=result.z_factor,
        ln_phi=tuple(result.ln_phi),
        warnings=_warnings(result.warnings),
    )


def furst_electrolyte_mod2004_phase(
    components: Sequence[str],
    T: Q,
    P: Q,
    x: Sequence[float],
    compressed_phase: str,
) -> FurstElectrolyteMod2004PhaseResult:
    """The 2004 revision of the Fürst electrolyte phase state, computed in Rust.

    The same kernels as `furst_electrolyte_phase` in Python here, resolved to the Rust
    function that carries the difference.
    """
    spec = _models_gen.model("eos.furst_electrolyte_mod2004_phase")
    result = _core.furst_electrolyte_mod2004_phase(
        list(components),
        input_to_si(spec, "T", T),
        input_to_si(spec, "P", P),
        list(x),
        compressed_phase,
    )
    return FurstElectrolyteMod2004PhaseResult(
        z_factor=result.z_factor,
        ln_phi=tuple(result.ln_phi),
        warnings=_warnings(result.warnings),
    )


def srk_cpa_phase(
    components: Sequence[str], T: Q, P: Q, z: Sequence[float], compressed_phase: str
) -> SrkCpaPhaseResult:
    """The SRK-CPA phase state, computed in Rust.

    The component names cross **unresolved**, and the Rust side looks them up in its own
    databank. That is the `eos.eos_cg_phase` precedent, and for this model it is what makes
    the cross-implementation comparison cover the *resolution* as well as the arithmetic.
    """
    spec = _models_gen.model("eos.srk_cpa_phase")
    result = _core.srk_cpa_phase(
        list(components),
        input_to_si(spec, "T", T),
        input_to_si(spec, "P", P),
        list(z),
        compressed_phase,
    )
    return SrkCpaPhaseResult(
        z_factor=result.z_factor,
        ln_phi=tuple(result.ln_phi),
        h_res=from_si(result.h_res.magnitude_si, result.h_res.unit),
        s_res=from_si(result.s_res.magnitude_si, result.s_res.unit),
        warnings=_warnings(result.warnings),
    )


def reactive_phase_equilibrium(
    components: Sequence[str],
    source: str,
    phase: str,
    moles: Sequence[Q],
    phase_charge: Q,
    phase_moles: Q,
    whole_system: bool,
    log_activity: Sequence[float],
    T: Q,
    max_iterations: float,
    tolerance: float,
    concentration_basis: str = "mole_fraction",
    seed: str = "none",
) -> ReactivePhaseEquilibriumResult:
    """The phase's reactive equilibrium, computed in Rust.

    The component names cross **unresolved**, as the reference potentials' do, and the
    reaction name does for the same reason: the Rust side reads the element table, the
    component table and the reaction table itself, so no coefficient reaches Python and
    the two languages cannot disagree about which row answered.
    """
    spec = _models_gen.model("reactions.reactive_phase_equilibrium")
    result = _core.reactive_phase_equilibrium(
        list(components),
        source,
        phase,
        [_si(spec, "moles", value) for value in moles],
        input_to_si(spec, "phase_charge", phase_charge),
        input_to_si(spec, "phase_moles", phase_moles),
        whole_system,
        [_si(spec, "log_activity", value) for value in log_activity],
        input_to_si(spec, "T", T),
        int(max_iterations),
        float(tolerance),
        concentration_basis,
        seed,
    )
    return ReactivePhaseEquilibriumResult(
        skipped=result.skipped,
        a_matrix=tuple(tuple(row) for row in result.a_matrix),
        b=tuple(from_si(value.magnitude_si, value.unit) for value in result.b),
        chem_ref=tuple(from_si(value.magnitude_si, value.unit) for value in result.chem_ref),
        moles=tuple(from_si(value.magnitude_si, value.unit) for value in result.moles),
        iterations=result.iterations,
        error=result.error,
        converged=result.converged,
        refinements=result.refinements,
        certified=result.certified,
        max_reaction_log_residual=result.max_reaction_log_residual,
        net_charge_moles=from_si(
            result.net_charge_moles.magnitude_si, result.net_charge_moles.unit
        ),
        max_element_residual=from_si(
            result.max_element_residual.magnitude_si, result.max_element_residual.unit
        ),
        seed_applied=result.seed_applied,
        seed_moles=tuple(from_si(value.magnitude_si, value.unit) for value in result.seed_moles),
        warnings=_warnings(result.warnings),
    )


def mixer(
    components: Sequence[str],
    feed_n: Sequence[Q],
    feed_z: Sequence[Sequence[float]],
    feed_p: Sequence[Q],
    feed_t: Sequence[Q],
    outlet_pressure: Q | None = None,
) -> MixerResult:
    """`process.mixer`, computed in Rust.

    The feeds cross as vectors, one entry per feed, because a `many` port is one port: the
    outlet pressure is the only scalar here and it is optional, which is why the Rust side
    takes an `Option` and the case may omit it.
    """
    spec = _models_gen.model("process.mixer")
    result = _core.mixer(
        list(components),
        [input_to_si(spec, "feed_n", v) for v in feed_n],
        [[_si(spec, "feed_z", v) for v in row] for row in feed_z],
        [input_to_si(spec, "feed_p", v) for v in feed_p],
        [input_to_si(spec, "feed_t", v) for v in feed_t],
        None if outlet_pressure is None else input_to_si(spec, "outlet_pressure", outlet_pressure),
    )
    return MixerResult(
        product_n=from_si(result.product_n.magnitude_si, result.product_n.unit),
        product_z=tuple(result.product_z),
        product_p=from_si(result.product_p.magnitude_si, result.product_p.unit),
        product_t=from_si(result.product_t.magnitude_si, result.product_t.unit),
        product_h=from_si(result.product_h.magnitude_si, result.product_h.unit),
        warnings=_warnings(result.warnings),
    )


def separator(
    components: Sequence[str],
    feed_n: Q,
    feed_z: Sequence[float],
    feed_p: Q,
    feed_t: Q,
    pressure_drop: Q,
    gas_in_liquid: float,
    heat_input: Q | None = None,
) -> SeparatorResult:
    """`process.separator`, computed in Rust.

    The component names cross unresolved, as the other unit operations' do: the Rust side
    resolves the mixture itself, so the two languages cannot disagree about which databank
    row answered.
    """
    spec = _models_gen.model("process.separator")
    result = _core.separator(
        list(components),
        input_to_si(spec, "feed_n", feed_n),
        [_si(spec, "feed_z", v) for v in feed_z],
        input_to_si(spec, "feed_p", feed_p),
        input_to_si(spec, "feed_t", feed_t),
        input_to_si(spec, "pressure_drop", pressure_drop),
        _si(spec, "gas_in_liquid", gas_in_liquid),
        None if heat_input is None else input_to_si(spec, "heat_input", heat_input),
    )
    return SeparatorResult(
        vapour_n=from_si(result.vapour_n.magnitude_si, result.vapour_n.unit),
        vapour_z=tuple(result.vapour_z),
        vapour_p=from_si(result.vapour_p.magnitude_si, result.vapour_p.unit),
        vapour_t=from_si(result.vapour_t.magnitude_si, result.vapour_t.unit),
        vapour_h=from_si(result.vapour_h.magnitude_si, result.vapour_h.unit),
        liquid_n=from_si(result.liquid_n.magnitude_si, result.liquid_n.unit),
        liquid_z=tuple(result.liquid_z),
        liquid_p=from_si(result.liquid_p.magnitude_si, result.liquid_p.unit),
        liquid_t=from_si(result.liquid_t.magnitude_si, result.liquid_t.unit),
        liquid_h=from_si(result.liquid_h.magnitude_si, result.liquid_h.unit),
        warnings=_warnings(result.warnings),
    )


def shortcut_distillation_column(
    components: Sequence[str],
    feed_n: Q,
    feed_z: Sequence[float],
    feed_p: Q,
    feed_t: Q,
    light_key: str,
    heavy_key: str,
    light_key_recovery_distillate: float,
    heavy_key_recovery_bottoms: float,
    reflux_ratio_multiplier: float,
    condenser_pressure: Q | None = None,
    reboiler_pressure: Q | None = None,
) -> ShortcutDistillationColumnResult:
    """`process.shortcut_distillation_column`, computed in Rust.

    The component names cross unresolved, as the other unit operations' do: the Rust side
    resolves the mixture itself, so the two languages cannot disagree about which databank
    row answered. The two key names cross as strings, because they name rows of that same
    resolution rather than values of it.
    """
    spec = _models_gen.model("process.shortcut_distillation_column")
    result = _core.shortcut_distillation_column(
        list(components),
        input_to_si(spec, "feed_n", feed_n),
        [_si(spec, "feed_z", v) for v in feed_z],
        input_to_si(spec, "feed_p", feed_p),
        input_to_si(spec, "feed_t", feed_t),
        light_key,
        heavy_key,
        _si(spec, "light_key_recovery_distillate", light_key_recovery_distillate),
        _si(spec, "heavy_key_recovery_bottoms", heavy_key_recovery_bottoms),
        _si(spec, "reflux_ratio_multiplier", reflux_ratio_multiplier),
        None
        if condenser_pressure is None
        else input_to_si(spec, "condenser_pressure", condenser_pressure),
        None
        if reboiler_pressure is None
        else input_to_si(spec, "reboiler_pressure", reboiler_pressure),
    )
    return ShortcutDistillationColumnResult(
        distillate_n=from_si(result.distillate_n.magnitude_si, result.distillate_n.unit),
        distillate_z=tuple(result.distillate_z),
        distillate_p=from_si(result.distillate_p.magnitude_si, result.distillate_p.unit),
        distillate_t=from_si(result.distillate_t.magnitude_si, result.distillate_t.unit),
        distillate_h=from_si(result.distillate_h.magnitude_si, result.distillate_h.unit),
        bottoms_n=from_si(result.bottoms_n.magnitude_si, result.bottoms_n.unit),
        bottoms_z=tuple(result.bottoms_z),
        bottoms_p=from_si(result.bottoms_p.magnitude_si, result.bottoms_p.unit),
        bottoms_t=from_si(result.bottoms_t.magnitude_si, result.bottoms_t.unit),
        bottoms_h=from_si(result.bottoms_h.magnitude_si, result.bottoms_h.unit),
        minimum_stages=result.minimum_stages,
        minimum_reflux_ratio=result.minimum_reflux_ratio,
        actual_stages=result.actual_stages,
        actual_reflux_ratio=result.actual_reflux_ratio,
        feed_tray_number=result.feed_tray_number,
        condenser_duty=from_si(result.condenser_duty.magnitude_si, result.condenser_duty.unit),
        reboiler_duty=from_si(result.reboiler_duty.magnitude_si, result.reboiler_duty.unit),
        relative_volatility=result.relative_volatility,
        warnings=_warnings(result.warnings),
    )


def distillation_column(
    components: Sequence[str],
    feed_n: Q,
    feed_z: Sequence[float],
    feed_p: Q,
    feed_t: Q,
    number_of_stages: int,
    feed_stage: int,
    has_reboiler: bool,
    has_condenser: bool,
    top_pressure: Q,
    bottom_pressure: Q,
    temperature_tolerance: float = 1.0e-6,
    max_iterations: int = 200,
    reboiler_temperature: Q | None = None,
    condenser_temperature: Q | None = None,
    murphree_efficiency: float | None = None,
    solver_type: str | None = None,
    top_specification_type: str | None = None,
    top_specification_target: float | None = None,
    top_specification_component: str | None = None,
    bottom_specification_type: str | None = None,
    bottom_specification_target: float | None = None,
    bottom_specification_component: str | None = None,
) -> DistillationColumnResult:
    """`process.distillation_column`, computed in Rust.

    The component names cross unresolved, as the other unit operations' do. The four
    unported parameters cross as `None` or their value, and the Rust side refuses them: a
    refusal must be the kernel's, so that the Python reference and the extension refuse for the
    same stated reason rather than each in its own words.
    """
    spec = _models_gen.model("process.distillation_column")
    result = _core.distillation_column(
        list(components),
        input_to_si(spec, "feed_n", feed_n),
        [_si(spec, "feed_z", v) for v in feed_z],
        input_to_si(spec, "feed_p", feed_p),
        input_to_si(spec, "feed_t", feed_t),
        int(_si(spec, "number_of_stages", number_of_stages)),
        int(_si(spec, "feed_stage", feed_stage)),
        has_reboiler,
        has_condenser,
        input_to_si(spec, "top_pressure", top_pressure),
        input_to_si(spec, "bottom_pressure", bottom_pressure),
        _si(spec, "temperature_tolerance", temperature_tolerance),
        int(_si(spec, "max_iterations", max_iterations)),
        # **Absent means no pin**, which is what `setReboilerTemperature` not having been
        # called does - and what makes a duty specification reachable at all. **The optional
        # inputs cross last**, in the order the model's own signature declares them.
        None
        if reboiler_temperature is None
        else input_to_si(spec, "reboiler_temperature", reboiler_temperature),
        None
        if condenser_temperature is None
        else input_to_si(spec, "condenser_temperature", condenser_temperature),
        murphree_efficiency,
        solver_type,
        top_specification_type,
        top_specification_target,
        top_specification_component,
        bottom_specification_type,
        bottom_specification_target,
        bottom_specification_component,
    )
    return DistillationColumnResult(
        tray_temperature=tuple(from_si(q.magnitude_si, q.unit) for q in result.tray_temperature),
        tray_pressure=tuple(from_si(q.magnitude_si, q.unit) for q in result.tray_pressure),
        tray_gas_n=tuple(from_si(q.magnitude_si, q.unit) for q in result.tray_gas_n),
        tray_liquid_n=tuple(from_si(q.magnitude_si, q.unit) for q in result.tray_liquid_n),
        distillate_n=from_si(result.distillate_n.magnitude_si, result.distillate_n.unit),
        distillate_z=tuple(result.distillate_z),
        distillate_p=from_si(result.distillate_p.magnitude_si, result.distillate_p.unit),
        distillate_t=from_si(result.distillate_t.magnitude_si, result.distillate_t.unit),
        distillate_h=from_si(result.distillate_h.magnitude_si, result.distillate_h.unit),
        bottoms_n=from_si(result.bottoms_n.magnitude_si, result.bottoms_n.unit),
        bottoms_z=tuple(result.bottoms_z),
        bottoms_p=from_si(result.bottoms_p.magnitude_si, result.bottoms_p.unit),
        bottoms_t=from_si(result.bottoms_t.magnitude_si, result.bottoms_t.unit),
        bottoms_h=from_si(result.bottoms_h.magnitude_si, result.bottoms_h.unit),
        condenser_duty=from_si(result.condenser_duty.magnitude_si, result.condenser_duty.unit),
        reboiler_duty=from_si(result.reboiler_duty.magnitude_si, result.reboiler_duty.unit),
        iterations=int(result.iterations),
        temperature_residual=float(result.temperature_residual),
        mass_residual=float(result.mass_residual),
        energy_residual=float(result.energy_residual),
        warnings=_warnings(result.warnings),
    )


def throttling_valve(
    components: Sequence[str],
    inlet_n: Q,
    inlet_z: Sequence[float],
    inlet_p: Q,
    inlet_t: Q,
    outlet_pressure: Q,
) -> ThrottlingValveResult:
    """`process.throttling_valve`, computed in Rust.

    The component names cross unresolved, as the other unit operations' do: the Rust side
    resolves the mixture itself, so the two languages cannot disagree about which databank
    row answered.
    """
    spec = _models_gen.model("process.throttling_valve")
    result = _core.throttling_valve(
        list(components),
        input_to_si(spec, "inlet_n", inlet_n),
        [_si(spec, "inlet_z", v) for v in inlet_z],
        input_to_si(spec, "inlet_p", inlet_p),
        input_to_si(spec, "inlet_t", inlet_t),
        input_to_si(spec, "outlet_pressure", outlet_pressure),
    )
    return ThrottlingValveResult(
        outlet_n=from_si(result.outlet_n.magnitude_si, result.outlet_n.unit),
        outlet_z=tuple(result.outlet_z),
        outlet_p=from_si(result.outlet_p.magnitude_si, result.outlet_p.unit),
        outlet_t=from_si(result.outlet_t.magnitude_si, result.outlet_t.unit),
        outlet_h=from_si(result.outlet_h.magnitude_si, result.outlet_h.unit),
        warnings=_warnings(result.warnings),
    )


def heat_exchanger(
    hot_components: Sequence[str],
    hot_in_n: Q,
    hot_in_z: Sequence[float],
    hot_in_p: Q,
    hot_in_t: Q,
    cold_components: Sequence[str],
    cold_in_n: Q,
    cold_in_z: Sequence[float],
    cold_in_p: Q,
    cold_in_t: Q,
    ua: Q | None = None,
    flow_arrangement: str = "counterflow",
    hot_outlet_temperature: Q | None = None,
    cold_outlet_temperature: Q | None = None,
) -> HeatExchangerResult:
    """`process.heat_exchanger`, computed in Rust.

    **Two component lists**, because this is the only unit operation whose ports do not
    share one fluid. Both cross unresolved, so the Rust side resolves each through the
    databank and the two languages cannot disagree about which row answered.
    """
    spec = _models_gen.model("process.heat_exchanger")
    result = _core.heat_exchanger(
        list(hot_components),
        list(cold_components),
        input_to_si(spec, "hot_in_n", hot_in_n),
        [_si(spec, "hot_in_z", v) for v in hot_in_z],
        input_to_si(spec, "hot_in_p", hot_in_p),
        input_to_si(spec, "hot_in_t", hot_in_t),
        input_to_si(spec, "cold_in_n", cold_in_n),
        [_si(spec, "cold_in_z", v) for v in cold_in_z],
        input_to_si(spec, "cold_in_p", cold_in_p),
        input_to_si(spec, "cold_in_t", cold_in_t),
        flow_arrangement,
        None if ua is None else input_to_si(spec, "ua", ua),
        None
        if hot_outlet_temperature is None
        else input_to_si(spec, "hot_outlet_temperature", hot_outlet_temperature),
        None
        if cold_outlet_temperature is None
        else input_to_si(spec, "cold_outlet_temperature", cold_outlet_temperature),
    )
    return HeatExchangerResult(
        hot_out_n=from_si(result.hot_out_n.magnitude_si, result.hot_out_n.unit),
        hot_out_z=tuple(result.hot_out_z),
        hot_out_p=from_si(result.hot_out_p.magnitude_si, result.hot_out_p.unit),
        hot_out_t=from_si(result.hot_out_t.magnitude_si, result.hot_out_t.unit),
        hot_out_h=from_si(result.hot_out_h.magnitude_si, result.hot_out_h.unit),
        cold_out_n=from_si(result.cold_out_n.magnitude_si, result.cold_out_n.unit),
        cold_out_z=tuple(result.cold_out_z),
        cold_out_p=from_si(result.cold_out_p.magnitude_si, result.cold_out_p.unit),
        cold_out_t=from_si(result.cold_out_t.magnitude_si, result.cold_out_t.unit),
        cold_out_h=from_si(result.cold_out_h.magnitude_si, result.cold_out_h.unit),
        warnings=_warnings(result.warnings),
    )


def splitter(
    components: Sequence[str],
    feed_n: Q,
    feed_z: Sequence[float],
    feed_p: Q,
    feed_t: Q,
    split_factors: Sequence[float],
) -> SplitterResult:
    """`process.splitter`, computed in Rust.

    The component names cross unresolved, as the pump's do: the Rust side resolves the
    mixture itself, so the two languages cannot disagree about which databank row answered.
    """
    spec = _models_gen.model("process.splitter")
    result = _core.splitter(
        list(components),
        input_to_si(spec, "feed_n", feed_n),
        [_si(spec, "feed_z", v) for v in feed_z],
        input_to_si(spec, "feed_p", feed_p),
        input_to_si(spec, "feed_t", feed_t),
        [_si(spec, "split_factors", v) for v in split_factors],
    )
    return SplitterResult(
        products_n=tuple(from_si(value.magnitude_si, value.unit) for value in result.products_n),
        products_z=tuple(tuple(row) for row in result.products_z),
        products_p=tuple(from_si(value.magnitude_si, value.unit) for value in result.products_p),
        products_t=tuple(from_si(value.magnitude_si, value.unit) for value in result.products_t),
        products_h=tuple(from_si(value.magnitude_si, value.unit) for value in result.products_h),
        warnings=_warnings(result.warnings),
    )


def absorption_column(
    gas_components: Sequence[str],
    gas_n: Q,
    gas_z: Sequence[float],
    gas_p: Q,
    gas_t: Q,
    solvent_components: Sequence[str],
    solvent_n: Q,
    solvent_z: Sequence[float],
    solvent_p: Q,
    solvent_t: Q,
    number_of_stages: int,
    top_pressure: Q,
    bottom_pressure: Q,
    temperature_tolerance: float,
    max_iterations: int,
    tray_temperatures: Sequence[float] | None = None,
    murphree_efficiency: float | None = None,
    component_murphree_efficiency: Sequence[float] | None = None,
    max_allowable_gas_load_factor: float | None = None,
    solver_type: str | None = None,
) -> AbsorptionColumnResult:
    """Solve a tray absorber; see :func:`azoth.process.reference.absorption_column`.

    The declared inputs cross in SI magnitudes, and **the two that are refused cross as they
    are**: a refusal must be the kernel's, so that the two implementations refuse for the same
    stated reason rather than each in its own words.
    """
    spec = _models_gen.model("process.absorption_column")
    result = _core.absorption_column(
        list(gas_components),
        list(solvent_components),
        input_to_si(spec, "gas_n", gas_n),
        [_si(spec, "gas_z", v) for v in gas_z],
        input_to_si(spec, "gas_p", gas_p),
        input_to_si(spec, "gas_t", gas_t),
        input_to_si(spec, "solvent_n", solvent_n),
        [_si(spec, "solvent_z", v) for v in solvent_z],
        input_to_si(spec, "solvent_p", solvent_p),
        input_to_si(spec, "solvent_t", solvent_t),
        int(_si(spec, "number_of_stages", number_of_stages)),
        input_to_si(spec, "top_pressure", top_pressure),
        input_to_si(spec, "bottom_pressure", bottom_pressure),
        _si(spec, "temperature_tolerance", temperature_tolerance),
        int(_si(spec, "max_iterations", max_iterations)),
        None if tray_temperatures is None else [float(v) for v in tray_temperatures],
        murphree_efficiency,
        None
        if component_murphree_efficiency is None
        else [float(v) for v in component_murphree_efficiency],
        max_allowable_gas_load_factor,
        solver_type,
    )
    return AbsorptionColumnResult(
        tray_temperature=tuple(from_si(q.magnitude_si, q.unit) for q in result.tray_temperature),
        tray_pressure=tuple(from_si(q.magnitude_si, q.unit) for q in result.tray_pressure),
        tray_gas_n=tuple(from_si(q.magnitude_si, q.unit) for q in result.tray_gas_n),
        tray_liquid_n=tuple(from_si(q.magnitude_si, q.unit) for q in result.tray_liquid_n),
        gas_out_n=from_si(result.gas_out_n.magnitude_si, result.gas_out_n.unit),
        gas_out_z=tuple(result.gas_out_z),
        gas_out_p=from_si(result.gas_out_p.magnitude_si, result.gas_out_p.unit),
        gas_out_t=from_si(result.gas_out_t.magnitude_si, result.gas_out_t.unit),
        gas_out_h=from_si(result.gas_out_h.magnitude_si, result.gas_out_h.unit),
        liquid_out_n=from_si(result.liquid_out_n.magnitude_si, result.liquid_out_n.unit),
        liquid_out_z=tuple(result.liquid_out_z),
        liquid_out_p=from_si(result.liquid_out_p.magnitude_si, result.liquid_out_p.unit),
        liquid_out_t=from_si(result.liquid_out_t.magnitude_si, result.liquid_out_t.unit),
        liquid_out_h=from_si(result.liquid_out_h.magnitude_si, result.liquid_out_h.unit),
        iterations=result.iterations,
        temperature_residual=result.temperature_residual,
        mass_residual=result.mass_residual,
        energy_residual=result.energy_residual,
        warnings=_warnings(result.warnings),
    )


def filter(
    components: Sequence[str],
    inlet_n: Q,
    inlet_z: Sequence[float],
    inlet_p: Q,
    inlet_t: Q,
    pressure_drop: Q,
) -> FilterResult:
    """`process.filter`, computed in Rust.

    The drop is required rather than optional, because the class's own default of `0.01`
    bar is a number the palette entry does not declare: a caller states the drop this
    machine has.
    """
    spec = _models_gen.model("process.filter")
    result = _core.filter(
        list(components),
        input_to_si(spec, "inlet_n", inlet_n),
        [_si(spec, "inlet_z", v) for v in inlet_z],
        input_to_si(spec, "inlet_p", inlet_p),
        input_to_si(spec, "inlet_t", inlet_t),
        input_to_si(spec, "pressure_drop", pressure_drop),
    )
    return FilterResult(
        outlet_n=from_si(result.outlet_n.magnitude_si, result.outlet_n.unit),
        outlet_z=tuple(result.outlet_z),
        outlet_p=from_si(result.outlet_p.magnitude_si, result.outlet_p.unit),
        outlet_t=from_si(result.outlet_t.magnitude_si, result.outlet_t.unit),
        outlet_h=from_si(result.outlet_h.magnitude_si, result.outlet_h.unit),
        applied_drop=from_si(result.applied_drop.magnitude_si, result.applied_drop.unit),
        warnings=_warnings(result.warnings),
    )


def compressor(
    components: Sequence[str],
    inlet_n: Q,
    inlet_z: Sequence[float],
    inlet_p: Q,
    inlet_t: Q,
    outlet_pressure: Q,
    isentropic_efficiency: float,
) -> CompressorResult:
    """`process.compressor`, computed in Rust.

    The component names cross unresolved, as every process model's do: the Rust side resolves
    the mixture itself, so the two languages cannot disagree about which databank row answered.
    """
    spec = _models_gen.model("process.compressor")
    result = _core.compressor(
        list(components),
        input_to_si(spec, "inlet_n", inlet_n),
        [_si(spec, "inlet_z", v) for v in inlet_z],
        input_to_si(spec, "inlet_p", inlet_p),
        input_to_si(spec, "inlet_t", inlet_t),
        input_to_si(spec, "outlet_pressure", outlet_pressure),
        float(isentropic_efficiency),
    )
    return CompressorResult(
        outlet_n=from_si(result.outlet_n.magnitude_si, result.outlet_n.unit),
        outlet_z=tuple(result.outlet_z),
        outlet_p=from_si(result.outlet_p.magnitude_si, result.outlet_p.unit),
        outlet_t=from_si(result.outlet_t.magnitude_si, result.outlet_t.unit),
        outlet_h=from_si(result.outlet_h.magnitude_si, result.outlet_h.unit),
        warnings=_warnings(result.warnings),
    )


def expander(
    components: Sequence[str],
    inlet_n: Q,
    inlet_z: Sequence[float],
    inlet_p: Q,
    inlet_t: Q,
    outlet_pressure: Q,
    isentropic_efficiency: float,
) -> ExpanderResult:
    """`process.expander`, computed in Rust.

    The same arguments as `process.compressor`, and the same route: what differs is which side
    of the division the efficiency sits on, and that is the Rust kernel's business.
    """
    spec = _models_gen.model("process.expander")
    result = _core.expander(
        list(components),
        input_to_si(spec, "inlet_n", inlet_n),
        [_si(spec, "inlet_z", v) for v in inlet_z],
        input_to_si(spec, "inlet_p", inlet_p),
        input_to_si(spec, "inlet_t", inlet_t),
        input_to_si(spec, "outlet_pressure", outlet_pressure),
        float(isentropic_efficiency),
    )
    return ExpanderResult(
        outlet_n=from_si(result.outlet_n.magnitude_si, result.outlet_n.unit),
        outlet_z=tuple(result.outlet_z),
        outlet_p=from_si(result.outlet_p.magnitude_si, result.outlet_p.unit),
        outlet_t=from_si(result.outlet_t.magnitude_si, result.outlet_t.unit),
        outlet_h=from_si(result.outlet_h.magnitude_si, result.outlet_h.unit),
        warnings=_warnings(result.warnings),
    )


def pipe(
    components: Sequence[str],
    inlet_n: Q,
    inlet_z: Sequence[float],
    inlet_p: Q,
    inlet_t: Q,
    length: Q,
    diameter: Q,
    roughness: Q,
) -> PipeResult:
    """`process.pipe`, computed in Rust.

    The three geometry arguments cross as quantities and are converted here, once, so the
    two backends cannot disagree about what a metre is. The component names cross unresolved,
    as every process model's do.
    """
    spec = _models_gen.model("process.pipe")
    result = _core.pipe(
        list(components),
        input_to_si(spec, "inlet_n", inlet_n),
        [_si(spec, "inlet_z", v) for v in inlet_z],
        input_to_si(spec, "inlet_p", inlet_p),
        input_to_si(spec, "inlet_t", inlet_t),
        input_to_si(spec, "length", length),
        input_to_si(spec, "diameter", diameter),
        input_to_si(spec, "roughness", roughness),
    )
    return PipeResult(
        outlet_n=from_si(result.outlet_n.magnitude_si, result.outlet_n.unit),
        outlet_z=tuple(result.outlet_z),
        outlet_p=from_si(result.outlet_p.magnitude_si, result.outlet_p.unit),
        outlet_t=from_si(result.outlet_t.magnitude_si, result.outlet_t.unit),
        outlet_h=from_si(result.outlet_h.magnitude_si, result.outlet_h.unit),
        pressure_drop=from_si(result.pressure_drop.magnitude_si, result.pressure_drop.unit),
        warnings=_warnings(result.warnings),
    )


def tank(
    components: Sequence[str],
    feed_n: Sequence[Q],
    feed_z: Sequence[Sequence[float]],
    feed_p: Sequence[Q],
    feed_t: Sequence[Q],
) -> TankResult:
    """`process.tank`, computed in Rust.

    The feeds cross as vectors and a matrix, as `process.mixer`'s do: a tank's inlet is a
    mixer, and its two outlets are the record's five fields under `gas` and `liquid`.
    """
    spec = _models_gen.model("process.tank")
    result = _core.tank(
        list(components),
        [input_to_si(spec, "feed_n", v) for v in feed_n],
        [[_si(spec, "feed_z", x) for x in row] for row in feed_z],
        [input_to_si(spec, "feed_p", v) for v in feed_p],
        [input_to_si(spec, "feed_t", v) for v in feed_t],
    )
    return TankResult(
        gas_n=from_si(result.gas_n.magnitude_si, result.gas_n.unit),
        gas_z=tuple(result.gas_z),
        gas_p=from_si(result.gas_p.magnitude_si, result.gas_p.unit),
        gas_t=from_si(result.gas_t.magnitude_si, result.gas_t.unit),
        gas_h=from_si(result.gas_h.magnitude_si, result.gas_h.unit),
        liquid_n=from_si(result.liquid_n.magnitude_si, result.liquid_n.unit),
        liquid_z=tuple(result.liquid_z),
        liquid_p=from_si(result.liquid_p.magnitude_si, result.liquid_p.unit),
        liquid_t=from_si(result.liquid_t.magnitude_si, result.liquid_t.unit),
        liquid_h=from_si(result.liquid_h.magnitude_si, result.liquid_h.unit),
        warnings=_warnings(result.warnings),
    )


def flare(
    components: Sequence[str],
    inlet_n: Q,
    inlet_z: Sequence[float],
    inlet_p: Q,
    inlet_t: Q,
) -> FlareResult:
    """`process.flare`, computed in Rust.

    The component names cross unresolved, as the process models' do: the Rust side resolves
    each against the standard's table for the duty and the element table for the emission, so
    the two languages cannot disagree about which row answered.
    """
    spec = _models_gen.model("process.flare")
    result = _core.flare(
        list(components),
        input_to_si(spec, "inlet_n", inlet_n),
        [_si(spec, "inlet_z", value) for value in inlet_z],
        input_to_si(spec, "inlet_p", inlet_p),
        input_to_si(spec, "inlet_t", inlet_t),
    )
    return FlareResult(
        product_n=from_si(result.product_n.magnitude_si, result.product_n.unit),
        product_z=tuple(result.product_z),
        product_p=from_si(result.product_p.magnitude_si, result.product_p.unit),
        product_t=from_si(result.product_t.magnitude_si, result.product_t.unit),
        product_h=from_si(result.product_h.magnitude_si, result.product_h.unit),
        heat_duty=from_si(result.heat_duty.magnitude_si, result.heat_duty.unit),
        co2_emission=from_si(result.co2_emission.magnitude_si, result.co2_emission.unit),
        warnings=_warnings(result.warnings),
    )


def iso6976(
    components: Sequence[str],
    z: Sequence[float],
    volumetric_reference_temperature: Q,
    energy_reference_temperature: Q,
) -> Iso6976Result:
    """`standards.iso6976`, computed in Rust.

    The component names cross unresolved, as the process models' do: the Rust side resolves
    each against the standard's own table, so the two languages cannot disagree about which
    row answered. The two reference temperatures cross as kelvins, which is the unit the
    spec declares.
    """
    spec = _models_gen.model("standards.iso6976")
    result = _core.iso6976(
        list(components),
        [_si(spec, "z", value) for value in z],
        input_to_si(spec, "volumetric_reference_temperature", volumetric_reference_temperature),
        input_to_si(spec, "energy_reference_temperature", energy_reference_temperature),
    )
    return Iso6976Result(
        molar_mass=from_si(result.molar_mass.magnitude_si, result.molar_mass.unit),
        compression_factor=result.compression_factor,
        relative_density=result.relative_density,
        density_ideal=from_si(result.density_ideal.magnitude_si, result.density_ideal.unit),
        density_real=from_si(result.density_real.magnitude_si, result.density_real.unit),
        superior_calorific_value=from_si(
            result.superior_calorific_value.magnitude_si, result.superior_calorific_value.unit
        ),
        inferior_calorific_value=from_si(
            result.inferior_calorific_value.magnitude_si, result.inferior_calorific_value.unit
        ),
        warnings=_warnings(result.warnings),
    )


def ejector(
    motive_components: Sequence[str],
    motive_n: Q,
    motive_z: Sequence[float],
    motive_p: Q,
    motive_t: Q,
    suction_components: Sequence[str],
    suction_n: Q,
    suction_z: Sequence[float],
    suction_p: Q,
    suction_t: Q,
    discharge_pressure: Q,
    motive_nozzle_efficiency: float,
    suction_nozzle_efficiency: float,
    mixing_efficiency: float,
    diffuser_efficiency: float,
) -> EjectorResult:
    """`process.ejector`, computed in Rust.

    **Two component lists, as `process.heat_exchanger`'s does**: a motive stream and a suction
    stream need not carry the same fluid, so neither list is the plain `components` the
    single-inlet unit operations declare.
    """
    spec = _models_gen.model("process.ejector")
    # **The two component lists cross first**, which is the transport shape `_dispatch` gives
    # a model with more than one `components` input: the fluid is named once per port and the
    # rest of that port's fields follow - the arrangement `process.heat_exchanger` set.
    result = _core.ejector(
        list(motive_components),
        list(suction_components),
        input_to_si(spec, "motive_n", motive_n),
        [_si(spec, "motive_z", v) for v in motive_z],
        input_to_si(spec, "motive_p", motive_p),
        input_to_si(spec, "motive_t", motive_t),
        input_to_si(spec, "suction_n", suction_n),
        [_si(spec, "suction_z", v) for v in suction_z],
        input_to_si(spec, "suction_p", suction_p),
        input_to_si(spec, "suction_t", suction_t),
        input_to_si(spec, "discharge_pressure", discharge_pressure),
        _si(spec, "motive_nozzle_efficiency", motive_nozzle_efficiency),
        _si(spec, "suction_nozzle_efficiency", suction_nozzle_efficiency),
        _si(spec, "mixing_efficiency", mixing_efficiency),
        _si(spec, "diffuser_efficiency", diffuser_efficiency),
    )
    return EjectorResult(
        outlet_n=from_si(result.outlet_n.magnitude_si, result.outlet_n.unit),
        outlet_z=tuple(result.outlet_z),
        outlet_p=from_si(result.outlet_p.magnitude_si, result.outlet_p.unit),
        outlet_t=from_si(result.outlet_t.magnitude_si, result.outlet_t.unit),
        outlet_h=from_si(result.outlet_h.magnitude_si, result.outlet_h.unit),
        warnings=_warnings(result.warnings),
    )


def three_phase_separator(
    components: Sequence[str],
    feed_n: Q,
    feed_z: Sequence[float],
    feed_p: Q,
    feed_t: Q,
    pressure_drop: Q,
    gas_in_aqueous: float,
    gas_in_oil: float,
    oil_in_aqueous: float,
    oil_in_gas: float,
    aqueous_in_gas: float,
    aqueous_in_oil: float,
    heat_input: Q | None = None,
) -> ThreePhaseSeparatorResult:
    """`process.three_phase_separator`, computed in Rust.

    Six fractions cross as bare numbers, as `process.separator`'s one does: they are
    dimensionless parameters, and the spec declares them without a unit.
    """
    spec = _models_gen.model("process.three_phase_separator")
    result = _core.three_phase_separator(
        list(components),
        input_to_si(spec, "feed_n", feed_n),
        [_si(spec, "feed_z", v) for v in feed_z],
        input_to_si(spec, "feed_p", feed_p),
        input_to_si(spec, "feed_t", feed_t),
        input_to_si(spec, "pressure_drop", pressure_drop),
        _si(spec, "gas_in_aqueous", gas_in_aqueous),
        _si(spec, "gas_in_oil", gas_in_oil),
        _si(spec, "oil_in_aqueous", oil_in_aqueous),
        _si(spec, "oil_in_gas", oil_in_gas),
        _si(spec, "aqueous_in_gas", aqueous_in_gas),
        _si(spec, "aqueous_in_oil", aqueous_in_oil),
        None if heat_input is None else input_to_si(spec, "heat_input", heat_input),
    )
    return ThreePhaseSeparatorResult(
        vapour_n=from_si(result.vapour_n.magnitude_si, result.vapour_n.unit),
        vapour_z=tuple(result.vapour_z),
        vapour_p=from_si(result.vapour_p.magnitude_si, result.vapour_p.unit),
        vapour_t=from_si(result.vapour_t.magnitude_si, result.vapour_t.unit),
        vapour_h=from_si(result.vapour_h.magnitude_si, result.vapour_h.unit),
        light_liquid_n=from_si(result.light_liquid_n.magnitude_si, result.light_liquid_n.unit),
        light_liquid_z=tuple(result.light_liquid_z),
        light_liquid_p=from_si(result.light_liquid_p.magnitude_si, result.light_liquid_p.unit),
        light_liquid_t=from_si(result.light_liquid_t.magnitude_si, result.light_liquid_t.unit),
        light_liquid_h=from_si(result.light_liquid_h.magnitude_si, result.light_liquid_h.unit),
        heavy_liquid_n=from_si(result.heavy_liquid_n.magnitude_si, result.heavy_liquid_n.unit),
        heavy_liquid_z=tuple(result.heavy_liquid_z),
        heavy_liquid_p=from_si(result.heavy_liquid_p.magnitude_si, result.heavy_liquid_p.unit),
        heavy_liquid_t=from_si(result.heavy_liquid_t.magnitude_si, result.heavy_liquid_t.unit),
        heavy_liquid_h=from_si(result.heavy_liquid_h.magnitude_si, result.heavy_liquid_h.unit),
        warnings=_warnings(result.warnings),
    )


def manifold(
    components: Sequence[str],
    feed_n: Sequence[Q],
    feed_z: Sequence[Sequence[float]],
    feed_p: Sequence[Q],
    feed_t: Sequence[Q],
    split_factors: Sequence[float],
) -> ManifoldResult:
    """`process.manifold`, computed in Rust.

    The feeds cross as vectors and a matrix, as `process.mixer`'s do, and the factors as a
    vector, as `process.splitter`'s do - because the manifold is those two composed.
    """
    spec = _models_gen.model("process.manifold")
    result = _core.manifold(
        list(components),
        [input_to_si(spec, "feed_n", v) for v in feed_n],
        [[_si(spec, "feed_z", x) for x in row] for row in feed_z],
        [input_to_si(spec, "feed_p", v) for v in feed_p],
        [input_to_si(spec, "feed_t", v) for v in feed_t],
        [float(f) for f in split_factors],
    )
    return ManifoldResult(
        products_n=tuple(from_si(q.magnitude_si, q.unit) for q in result.products_n),
        products_z=tuple(tuple(row) for row in result.products_z),
        products_p=tuple(from_si(q.magnitude_si, q.unit) for q in result.products_p),
        products_t=tuple(from_si(q.magnitude_si, q.unit) for q in result.products_t),
        products_h=tuple(from_si(q.magnitude_si, q.unit) for q in result.products_h),
        warnings=_warnings(result.warnings),
    )


def gas_scrubber(
    components: Sequence[str],
    feed_n: Q,
    feed_z: Sequence[float],
    feed_p: Q,
    feed_t: Q,
    pressure_drop: Q,
    gas_in_liquid: float,
    heat_input: Q | None = None,
) -> GasScrubberResult:
    """`process.gas_scrubber`, computed in Rust.

    The separator's arguments, because the class does not override `run`.
    """
    spec = _models_gen.model("process.gas_scrubber")
    result = _core.gas_scrubber(
        list(components),
        input_to_si(spec, "feed_n", feed_n),
        [_si(spec, "feed_z", v) for v in feed_z],
        input_to_si(spec, "feed_p", feed_p),
        input_to_si(spec, "feed_t", feed_t),
        input_to_si(spec, "pressure_drop", pressure_drop),
        float(gas_in_liquid),
        None if heat_input is None else input_to_si(spec, "heat_input", heat_input),
    )
    return GasScrubberResult(
        vapour_n=from_si(result.vapour_n.magnitude_si, result.vapour_n.unit),
        vapour_z=tuple(result.vapour_z),
        vapour_p=from_si(result.vapour_p.magnitude_si, result.vapour_p.unit),
        vapour_t=from_si(result.vapour_t.magnitude_si, result.vapour_t.unit),
        vapour_h=from_si(result.vapour_h.magnitude_si, result.vapour_h.unit),
        liquid_n=from_si(result.liquid_n.magnitude_si, result.liquid_n.unit),
        liquid_z=tuple(result.liquid_z),
        liquid_p=from_si(result.liquid_p.magnitude_si, result.liquid_p.unit),
        liquid_t=from_si(result.liquid_t.magnitude_si, result.liquid_t.unit),
        liquid_h=from_si(result.liquid_h.magnitude_si, result.liquid_h.unit),
        warnings=_warnings(result.warnings),
    )


def component_splitter(
    components: Sequence[str],
    feed_n: Q,
    feed_z: Sequence[float],
    feed_p: Q,
    feed_t: Q,
    split_factors: Sequence[float],
) -> ComponentSplitterResult:
    """`process.component_splitter`, computed in Rust."""
    spec = _models_gen.model("process.component_splitter")
    result = _core.component_splitter(
        list(components),
        input_to_si(spec, "feed_n", feed_n),
        [_si(spec, "feed_z", v) for v in feed_z],
        input_to_si(spec, "feed_p", feed_p),
        input_to_si(spec, "feed_t", feed_t),
        [float(f) for f in split_factors],
    )
    return ComponentSplitterResult(
        overhead_n=from_si(result.overhead_n.magnitude_si, result.overhead_n.unit),
        overhead_z=tuple(result.overhead_z),
        overhead_p=from_si(result.overhead_p.magnitude_si, result.overhead_p.unit),
        overhead_t=from_si(result.overhead_t.magnitude_si, result.overhead_t.unit),
        overhead_h=from_si(result.overhead_h.magnitude_si, result.overhead_h.unit),
        bottoms_n=from_si(result.bottoms_n.magnitude_si, result.bottoms_n.unit),
        bottoms_z=tuple(result.bottoms_z),
        bottoms_p=from_si(result.bottoms_p.magnitude_si, result.bottoms_p.unit),
        bottoms_t=from_si(result.bottoms_t.magnitude_si, result.bottoms_t.unit),
        bottoms_h=from_si(result.bottoms_h.magnitude_si, result.bottoms_h.unit),
        warnings=_warnings(result.warnings),
    )


def cooler(
    components: Sequence[str],
    inlet_n: Q,
    inlet_z: Sequence[float],
    inlet_p: Q,
    inlet_t: Q,
    outlet_temperature: Q | None = None,
    duty: Q | None = None,
    pressure_drop: Q | None = None,
) -> CoolerResult:
    """`process.cooler`, computed in Rust.

    The same three optional arguments as `process.heater`, and the same refusals: the
    arithmetic behind both is `Heater.run`.
    """
    spec = _models_gen.model("process.cooler")
    result = _core.cooler(
        list(components),
        input_to_si(spec, "inlet_n", inlet_n),
        [_si(spec, "inlet_z", v) for v in inlet_z],
        input_to_si(spec, "inlet_p", inlet_p),
        input_to_si(spec, "inlet_t", inlet_t),
        None
        if outlet_temperature is None
        else input_to_si(spec, "outlet_temperature", outlet_temperature),
        None if duty is None else input_to_si(spec, "duty", duty),
        None if pressure_drop is None else input_to_si(spec, "pressure_drop", pressure_drop),
    )
    return CoolerResult(
        outlet_n=from_si(result.outlet_n.magnitude_si, result.outlet_n.unit),
        outlet_z=tuple(result.outlet_z),
        outlet_p=from_si(result.outlet_p.magnitude_si, result.outlet_p.unit),
        outlet_t=from_si(result.outlet_t.magnitude_si, result.outlet_t.unit),
        outlet_h=from_si(result.outlet_h.magnitude_si, result.outlet_h.unit),
        outlet_duty=from_si(result.outlet_duty.magnitude_si, result.outlet_duty.unit),
        warnings=_warnings(result.warnings),
    )


def heater(
    components: Sequence[str],
    inlet_n: Q,
    inlet_z: Sequence[float],
    inlet_p: Q,
    inlet_t: Q,
    outlet_temperature: Q | None = None,
    duty: Q | None = None,
    pressure_drop: Q | None = None,
) -> HeaterResult:
    """`process.heater`, computed in Rust.

    **The three optional arguments cross as ``None`` and not as a sentinel**, because each
    one is a branch of `Heater.run` rather than a value: a temperature with no duty, a duty
    with no temperature, a pressure drop on either, or none of the three - which leaves the
    drop isothermal. The Rust side refuses a temperature and a duty together.
    """
    spec = _models_gen.model("process.heater")
    result = _core.heater(
        list(components),
        input_to_si(spec, "inlet_n", inlet_n),
        [_si(spec, "inlet_z", v) for v in inlet_z],
        input_to_si(spec, "inlet_p", inlet_p),
        input_to_si(spec, "inlet_t", inlet_t),
        None
        if outlet_temperature is None
        else input_to_si(spec, "outlet_temperature", outlet_temperature),
        None if duty is None else input_to_si(spec, "duty", duty),
        None if pressure_drop is None else input_to_si(spec, "pressure_drop", pressure_drop),
    )
    return HeaterResult(
        outlet_n=from_si(result.outlet_n.magnitude_si, result.outlet_n.unit),
        outlet_z=tuple(result.outlet_z),
        outlet_p=from_si(result.outlet_p.magnitude_si, result.outlet_p.unit),
        outlet_t=from_si(result.outlet_t.magnitude_si, result.outlet_t.unit),
        outlet_h=from_si(result.outlet_h.magnitude_si, result.outlet_h.unit),
        outlet_duty=from_si(result.outlet_duty.magnitude_si, result.outlet_duty.unit),
        warnings=_warnings(result.warnings),
    )


def pump(
    components: Sequence[str],
    inlet_n: Q,
    inlet_z: Sequence[float],
    inlet_p: Q,
    inlet_t: Q,
    outlet_pressure: Q,
    isentropic_efficiency: float,
) -> PumpResult:
    """`process.pump`, computed in Rust.

    The component names cross unresolved, as the reaction ids' do: the Rust side resolves
    the mixture itself, so the two languages cannot disagree about which databank row
    answered.
    """
    spec = _models_gen.model("process.pump")
    result = _core.pump(
        list(components),
        input_to_si(spec, "inlet_n", inlet_n),
        [_si(spec, "inlet_z", v) for v in inlet_z],
        input_to_si(spec, "inlet_p", inlet_p),
        input_to_si(spec, "inlet_t", inlet_t),
        input_to_si(spec, "outlet_pressure", outlet_pressure),
        float(isentropic_efficiency),
    )
    return PumpResult(
        outlet_n=from_si(result.outlet_n.magnitude_si, result.outlet_n.unit),
        outlet_z=tuple(result.outlet_z),
        outlet_p=from_si(result.outlet_p.magnitude_si, result.outlet_p.unit),
        outlet_t=from_si(result.outlet_t.magnitude_si, result.outlet_t.unit),
        outlet_h=from_si(result.outlet_h.magnitude_si, result.outlet_h.unit),
        warnings=_warnings(result.warnings),
    )


def reactive_hybrid_eos_ge_flash(
    components: Sequence[str], cubic: str, T: Q, P: Q, moles: Sequence[float]
) -> ReactiveHybridEosGeFlashResult:
    """The coupled reactive hybrid flash, computed in Rust.

    The names cross **unresolved**: the two EoS roles' constants, the seeding's classes and
    the brine's ion mask all resolve on the Rust side from the same databank, so the two
    languages cannot disagree about which substance is which.
    """
    spec = _models_gen.model("reactions.reactive_hybrid_eos_ge_flash")
    result = _core.reactive_hybrid_eos_ge_flash(
        list(components),
        cubic,
        input_to_si(spec, "T", T),
        input_to_si(spec, "P", P),
        # The unit-carrying vector crosses as SI magnitudes element by element, for the same
        # reason the reference converts it: a model's declared units are the boundary's.
        [_si(spec, "moles", value) for value in moles],
    )
    return ReactiveHybridEosGeFlashResult(
        beta=tuple(result.beta),
        x=tuple(tuple(row) for row in result.x),
        coupled_moles=tuple(
            from_si(value.magnitude_si, value.unit) for value in result.coupled_moles
        ),
        aqueous_moles=tuple(
            from_si(value.magnitude_si, value.unit) for value in result.aqueous_moles
        ),
        passes=result.passes,
        chemical_deviation=result.chemical_deviation,
        residual=result.residual,
        max_material_balance_residual=result.max_material_balance_residual,
        max_log_fugacity_residual=result.max_log_fugacity_residual,
        element_residual=result.element_residual,
        charge_residual=result.charge_residual,
        warnings=_warnings(result.warnings),
    )


def reactive_tp_flash(
    components: Sequence[str],
    T: Q,
    P: Q,
    moles: Sequence[Q],
    max_phases: float,
) -> ReactiveTpFlashResult:
    """The reactive flash, computed in Rust.

    The component names cross **unresolved**, as the reference potentials' and the phase
    operation's do: the Rust side reads the element table, the component table and the
    formation columns itself, so no coefficient reaches Python and the two languages cannot
    disagree about which row answered.
    """
    spec = _models_gen.model("reactions.reactive_tp_flash")
    result = _core.reactive_tp_flash(
        list(components),
        input_to_si(spec, "T", T),
        input_to_si(spec, "P", P),
        [_si(spec, "moles", value) for value in moles],
        float(max_phases),
    )
    return ReactiveTpFlashResult(
        phase_count=result.phase_count,
        phase_moles=tuple(
            tuple(from_si(value.magnitude_si, value.unit) for value in row)
            for row in result.phase_moles
        ),
        phase_fraction=tuple(result.phase_fraction),
        converged=result.converged,
        total_iterations=result.total_iterations,
        equilibrium_total_moles=from_si(
            result.equilibrium_total_moles.magnitude_si,
            result.equilibrium_total_moles.unit,
        ),
        gibbs_energy=result.gibbs_energy,
        residual=result.residual,
        element_residual=result.element_residual,
        warnings=_warnings(result.warnings),
    )


def reactive_ph_flash(
    components: Sequence[str],
    T: Q,
    P: Q,
    moles: Sequence[Q],
    enthalpy: Q,
    max_phases: float,
) -> ReactivePhFlashResult:
    """The reactive PH flash, computed in Rust.

    The component names cross **unresolved**, as they do for the TP flash: the Rust side
    reads the element table, the formation columns and the heat-capacity polynomial itself.
    """
    spec = _models_gen.model("reactions.reactive_ph_flash")
    result = _core.reactive_ph_flash(
        list(components),
        input_to_si(spec, "T", T),
        input_to_si(spec, "P", P),
        [_si(spec, "moles", value) for value in moles],
        input_to_si(spec, "enthalpy", enthalpy),
        float(max_phases),
    )
    return ReactivePhFlashResult(
        temperature=from_si(result.temperature.magnitude_si, result.temperature.unit),
        converged=result.converged,
        outer_iterations=result.outer_iterations,
        total_inner_iterations=result.total_inner_iterations,
        warnings=_warnings(result.warnings),
    )


def kinetic_rate_law(
    law: str,
    T: Q,
    reference_rate: float,
    activation_energy: Q,
    reference_temperature: Q,
) -> KineticRateLawResult:
    """The kinetic rate factor, computed in Rust.

    The selector is parsed on the Rust side, so the law this library carries is decided in
    one place; the three Arrhenius parameters cross as magnitudes whether or not the branch
    reads them, because the legacy one ignores all three.
    """
    spec = _models_gen.model("reactions.kinetic_rate_law")
    result = _core.kinetic_rate_law(
        law,
        input_to_si(spec, "T", T),
        float(reference_rate),
        input_to_si(spec, "activation_energy", activation_energy),
        input_to_si(spec, "reference_temperature", reference_temperature),
    )
    return KineticRateLawResult(
        rate_factor=result.rate_factor,
        warnings=_warnings(result.warnings),
    )


def kinetics(
    components: Sequence[str],
    reaction_components: Sequence[str],
    reaction_lengths: Sequence[float],
    reaction_coefficients: Sequence[float],
    rate_factors: Sequence[float],
    equilibrium_constants: Sequence[float],
    fractions: Sequence[float],
    molar_masses: Sequence[float],
    density: Q,
    inter_fractions: Sequence[float],
    inter_density: Q,
    diffusion: Sequence[float],
) -> KineticsResult:
    """The Krishna-Standart rate matrix, computed in Rust.

    Every species and every reaction name crosses **unresolved**, as a string: which
    species a phase carries and which side of a reaction they sit on are read on the Rust
    side, so the two languages cannot disagree about which row answered.
    """
    spec = _models_gen.model("reactions.kinetics")
    result = _core.kinetics(
        list(components),
        list(reaction_components),
        [float(value) for value in reaction_lengths],
        [float(value) for value in reaction_coefficients],
        [float(value) for value in rate_factors],
        [float(value) for value in equilibrium_constants],
        [float(value) for value in fractions],
        [_si(spec, "molar_masses", value) for value in molar_masses],
        _si(spec, "density", density),
        [float(value) for value in inter_fractions],
        _si(spec, "inter_density", inter_density),
        [_si(spec, "diffusion", value) for value in diffusion],
    )
    return KineticsResult(
        coefficient=tuple(result.coefficient),
        phi_infinite=tuple(result.phi_infinite),
        irreversible=tuple(result.irreversible),
        warnings=_warnings(result.warnings),
    )


def effective_diffusion(
    binary_diffusion: Sequence[Sequence[Q]],
    x: Sequence[float],
) -> EffectiveDiffusionResult:
    """The effective diffusion coefficients, computed in Rust.

    The matrix's entries carry a unit, so each one crosses as its SI magnitude; which pair
    coefficient is which is the caller's matrix and the Rust side adds nothing to it.
    """
    spec = _models_gen.model("eos.effective_diffusion")
    result = _core.effective_diffusion(
        [[_si(spec, "binary_diffusion", value) for value in row] for row in binary_diffusion],
        [float(value) for value in x],
    )
    return EffectiveDiffusionResult(
        effective_diffusion=tuple(
            from_si(value.magnitude_si, value.unit) for value in result.effective_diffusion
        ),
        warnings=_warnings(result.warnings),
    )


def ge_flash(
    components: Sequence[str],
    cubic: str,
    liquid_model: str,
    T: Q,
    P: Q,
    z: Sequence[float],
) -> GeFlashResult:
    """The generalised gamma-phi flash, computed in Rust.

    The names cross **unresolved**, as they do for `eos.ge_nrtl_flash`'s parameters: the
    Rust side resolves whichever liquid model is named against its own databank, so the two
    languages cannot disagree about which row answered.
    """
    result = _core.ge_flash(
        list(components),
        liquid_model,
        input_to_si(_models_gen.model("eos.ge_flash"), "T", T),
        input_to_si(_models_gen.model("eos.ge_flash"), "P", P),
        list(z),
        cubic,
    )
    return GeFlashResult(
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
