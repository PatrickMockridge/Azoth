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

# **The generated half.** Every id whose signature is the public wrapper's and whose boundary is
# the spec's own inputs is emitted by `tools/gen_python_bridge.py`; what is left in this file is
# the ids that take a `Mixture` or a record the spec does not name. `resolve` derives the id's
# function name and looks it up in this module's globals, which the star import fills.
from azoth._rust_bridge_gen import *  # noqa: F403

# The helpers the generated adapters and the hand-written ones both call. Imported rather than
# defined here because `azoth._rust_bridge_gen` imports them too, and a module cannot import the
# file that re-exports it. `from ... import _warnings` keeps the name reachable as
# `azoth._rust_bridge._warnings`, which is where `azoth.batch._core` reads it.
from azoth.core.bridge_support import _warnings as _warnings
from azoth.core.errors import PropertyUnavailableError
from azoth.core.result import (
    AbsorptionColumnResult,
    AqueousViscosityResult,
    BubblePressureResult,
    BubbleTemperatureResult,
    BwrsPhaseResult,
    CapillaryDewPointResult,
    ChemicalEquilibriumResult,
    CriticalPointResult,
    DewPressureResult,
    DewTemperatureResult,
    DistillationColumnResult,
    EffectiveDiffusionResult,
    EjectorResult,
    EquilibriumConstantResult,
    GeNrtlFlashResult,
    GeNrtlPhaseResult,
    GeUnifacPhaseResult,
    GeUniquacPhaseResult,
    GeVanLaarAcidPhaseResult,
    GeWilsonPhaseResult,
    HeatExchangerResult,
    HybridEosGeFlashResult,
    HydrateInhibitorConcentrationResult,
    HydrateInhibitorWtResult,
    HydrogenPhaseResult,
    KineticsResult,
    LiquidConductivityPolynomResult,
    ManifoldResult,
    MasonSaxenaConductivityResult,
    MixerResult,
    MolarEnthalpyEntropyResult,
    NrtlActivityCoefficientsResult,
    OrificeFlowResult,
    PackedColumnResult,
    PhFlashResult,
    PsFlashResult,
    PtFlashResult,
    PtPhaseEnvelopeResult,
    PuFlashResult,
    PureSaturationResult,
    PvfFlashResult,
    PvFlashResult,
    PvRefluxFlashResult,
    RateBasedPackedColumnResult,
    ReactiveHybridEosGeFlashResult,
    ReactivePhaseEquilibriumResult,
    ReactivePhFlashResult,
    ReactiveTpFlashResult,
    SaltPrecipitationResult,
    StabilityTestResult,
    StrippingColumnResult,
    TankResult,
    ThermalConductivityResult,
    ThFlashResult,
    TpMultiflashResult,
    TsFlashResult,
    TuFlashResult,
    TvFlashResult,
    TvFractionFlashResult,
    UnifacActivityCoefficientsResult,
    UnifacPsrkActivityCoefficientsResult,
    UnifacUmrpruActivityCoefficientsResult,
    UniquacActivityCoefficientsResult,
    VanLaarAcidActivityCoefficientsResult,
    VhFlashResult,
    ViscosityResult,
    VsFlashResult,
    VuFlashResult,
    VuFlashSingleCompResult,
    WilkeViscosityResult,
    WilsonActivityCoefficientsResult,
)
from azoth.core.result import Phase as _Phase
from azoth.core.result import StabilityVerdict as _StabilityVerdict
from azoth.core.result import TpMultiflashSeed as _TpMultiflashSeed
from azoth.core.units import Q, from_si, input_to_si, to_si


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


def liquid_conductivity_polynom(
    liquid_conductivity: Sequence[Sequence[float]],
    molar_mass: Sequence[Q],
    z: Sequence[float],
    T: Q,
) -> LiquidConductivityPolynomResult:
    """A liquid's thermal conductivity, computed in Rust."""
    spec = _models_gen.model("eos.liquid_conductivity_polynom")
    result = _core.liquid_conductivity_polynom(
        [[float(value) for value in row] for row in liquid_conductivity],
        [_si(spec, "molar_mass", value) for value in molar_mass],
        [float(value) for value in z],
        _si(spec, "T", T),
    )
    return LiquidConductivityPolynomResult(
        k=from_si(result.k.magnitude_si, result.k.unit),
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

    `vapour_fraction` crosses as `Option<f64>` and becomes `None`, not a sentinel. That is the
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
        vapour_fraction=result.vapour_fraction,
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
        dew_temperature=tuple(from_si(q.magnitude_si, q.unit) for q in raw.dew_temperature),
        dew_pressure=tuple(from_si(q.magnitude_si, q.unit) for q in raw.dew_pressure),
        bubble_temperature=tuple(from_si(q.magnitude_si, q.unit) for q in raw.bubble_temperature),
        bubble_pressure=tuple(from_si(q.magnitude_si, q.unit) for q in raw.bubble_pressure),
        cricondenbar_temperature=from_si(raw.cricondenbar_temperature.magnitude_si, "K"),
        cricondenbar_pressure=from_si(raw.cricondenbar_pressure.magnitude_si, "Pa"),
        cricondentherm_temperature=from_si(raw.cricondentherm_temperature.magnitude_si, "K"),
        cricondentherm_pressure=from_si(raw.cricondentherm_pressure.magnitude_si, "Pa"),
        critical_temperature=from_si(raw.critical_temperature.magnitude_si, "K"),
        critical_pressure=from_si(raw.critical_pressure.magnitude_si, "Pa"),
        iterations=raw.iterations,
        residual=from_si(raw.residual.magnitude_si, raw.residual.unit),
        warnings=_warnings(raw.warnings),
    )


def ph_flash(mixture: Any, ideal_gas: Any, P: Q, H: Q, z: Sequence[float]) -> PhFlashResult:
    """The pressure-enthalpy flash of a mixture, solved in Rust.

    The ideal-gas vectors cross as six parallel lists plus the two reference states,
    exactly as `molar_enthalpy_entropy` sends them. They are the datum the requested
    enthalpy is a difference from, and a coefficient set without a reference state is
    not a thermodynamic model - which is why they are required here even though this
    model uses only four of the six.

    `vapour_fraction` crosses as `Option<f64>` and becomes `None`, not a sentinel: a single-phase
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
        vapour_fraction=result.vapour_fraction,
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
        vapour_fraction=result.vapour_fraction,
        x=tuple(result.x),
        y=tuple(result.y),
        k=tuple(result.k),
        phase=_Phase(result.phase),
        z_liquid=result.z_liquid,
        z_vapour=result.z_vapour,
        iterations=result.iterations,
        residual=from_si(result.residual.magnitude_si, result.residual.unit),
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
        vapour_fraction=result.vapour_fraction,
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
    mixture: Any, P: Q, vapour_fraction: float, temperature: Q, z: Sequence[float]
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
        float(vapour_fraction),
        input_to_si(spec, "temperature", temperature),
        list(z),
        mixture.cubic.name,
        mixture.alpha,
        [list(c.alpha_params) for c in mixture.components],
    )
    return PvfFlashResult(
        T=from_si(result.T.magnitude_si, result.T.unit),
        vapour_fraction=result.vapour_fraction,
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
        vapour_fraction=result.vapour_fraction,
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
        vapour_fraction=result.vapour_fraction,
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
        phase_fractions=tuple(result.phase_fractions),
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
        vapour_fraction=result.vapour_fraction,
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
        vapour_fraction=result.vapour_fraction,
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
        vapour_fraction=result.vapour_fraction,
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
        vapour_fraction=result.vapour_fraction,
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
        vapour_fraction=result.vapour_fraction,
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
        vapour_fraction=result.vapour_fraction,
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
        vapour_fraction=result.vapour_fraction,
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
        vapour_fraction=result.vapour_fraction,
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
        vapour_fraction=result.vapour_fraction,
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

    `vapour_fraction` crosses as `Option<f64>` and becomes `None`, not a sentinel.
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
        vapour_fraction=result.vapour_fraction,
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
        inhibitor_moles=from_si(result.inhibitor_moles.magnitude_si, result.inhibitor_moles.unit),
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
        inhibitor_moles=from_si(result.inhibitor_moles.magnitude_si, result.inhibitor_moles.unit),
        weight_fraction=result.weight_fraction,
        hydrate_temperature=from_si(
            result.hydrate_temperature.magnitude_si, result.hydrate_temperature.unit
        ),
        iterations=result.iterations,
        residual=from_si(result.residual.magnitude_si, result.residual.unit),
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
        phase_fractions=tuple(result.phase_fractions),
        x=tuple(tuple(row) for row in result.x),
        ln_phi=tuple(tuple(row) for row in result.ln_phi),
        iterations=result.iterations,
        residual=result.residual,
        max_material_balance_residual=result.max_material_balance_residual,
        max_log_fugacity_residual=result.max_log_fugacity_residual,
        min_t_over_tc=result.min_t_over_tc,
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
    tray_murphree_efficiency: Sequence[float] | None = None,
    solver_type: str | None = None,
    top_specification_type: str | None = None,
    top_specification_target: float | None = None,
    top_specification_component: str | None = None,
    bottom_specification_type: str | None = None,
    bottom_specification_target: float | None = None,
    bottom_specification_component: str | None = None,
    reactive: bool | None = None,
    reactive_start_tray: int | None = None,
    reactive_end_tray: int | None = None,
    gas_side_draw_fractions: Sequence[float] | None = None,
    liquid_side_draw_fractions: Sequence[float] | None = None,
    pumparound_fractions: Sequence[float] | None = None,
    side_draw_flow_tray: int | None = None,
    side_draw_flow_phase: str | None = None,
    side_draw_flow_target: object | None = None,
    side_draw_flow_tolerance: float | None = None,
    side_draw_flow_max_iterations: int | None = None,
    pumparound_return_tray: int | None = None,
    pumparound_draw_tray: int | None = None,
    pumparound_draw_fraction: float | None = None,
    pumparound_temperature_drop: object | None = None,
    pumparound_tolerance: float | None = None,
    pumparound_max_iterations: int | None = None,
    column_diameter: object | None = None,
    max_allowable_fs_factor: float | None = None,
    internals_type: str | None = None,
    tray_spacing: object | None = None,
    weir_height: object | None = None,
    hole_diameter: object | None = None,
    hole_area_fraction: float | None = None,
    downcommer_area_fraction: float | None = None,
    design_flood_fraction: float | None = None,
    column_diameter_override: object | None = None,
    hydraulic_pressure_drop_coupling: bool | None = None,
    hydraulic_pressure_drop_internals_type: str | None = None,
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
        None if tray_murphree_efficiency is None else [float(v) for v in tray_murphree_efficiency],
        solver_type,
        top_specification_type,
        top_specification_target,
        top_specification_component,
        bottom_specification_type,
        bottom_specification_target,
        bottom_specification_component,
        reactive,
        None if reactive_start_tray is None else int(reactive_start_tray),
        None if reactive_end_tray is None else int(reactive_end_tray),
        None if gas_side_draw_fractions is None else [float(v) for v in gas_side_draw_fractions],
        None
        if liquid_side_draw_fractions is None
        else [float(v) for v in liquid_side_draw_fractions],
        None if pumparound_fractions is None else [float(v) for v in pumparound_fractions],
        None if side_draw_flow_tray is None else int(side_draw_flow_tray),
        side_draw_flow_phase,
        None
        if side_draw_flow_target is None
        else input_to_si(spec, "side_draw_flow_target", side_draw_flow_target),
        side_draw_flow_tolerance,
        None if side_draw_flow_max_iterations is None else int(side_draw_flow_max_iterations),
        None if pumparound_return_tray is None else int(pumparound_return_tray),
        None if pumparound_draw_tray is None else int(pumparound_draw_tray),
        pumparound_draw_fraction,
        None
        if pumparound_temperature_drop is None
        else input_to_si(spec, "pumparound_temperature_drop", pumparound_temperature_drop),
        pumparound_tolerance,
        None if pumparound_max_iterations is None else int(pumparound_max_iterations),
        # **The capacity limits and the internals tree follow the same rule**: each crosses as
        # `None` or its SI magnitude, so the Rust half owns every default.
        None if column_diameter is None else input_to_si(spec, "column_diameter", column_diameter),
        max_allowable_fs_factor,
        internals_type,
        None if tray_spacing is None else input_to_si(spec, "tray_spacing", tray_spacing),
        None if weir_height is None else input_to_si(spec, "weir_height", weir_height),
        None if hole_diameter is None else input_to_si(spec, "hole_diameter", hole_diameter),
        hole_area_fraction,
        downcommer_area_fraction,
        design_flood_fraction,
        None
        if column_diameter_override is None
        else input_to_si(spec, "column_diameter_override", column_diameter_override),
        hydraulic_pressure_drop_coupling,
        hydraulic_pressure_drop_internals_type,
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
        gas_side_draw_n=tuple(
            from_si(value.magnitude_si, value.unit) for value in result.gas_side_draw_n
        ),
        liquid_side_draw_n=tuple(
            from_si(value.magnitude_si, value.unit) for value in result.liquid_side_draw_n
        ),
        pumparound_n=tuple(
            from_si(value.magnitude_si, value.unit) for value in result.pumparound_n
        ),
        condenser_duty=from_si(result.condenser_duty.magnitude_si, result.condenser_duty.unit),
        reboiler_duty=from_si(result.reboiler_duty.magnitude_si, result.reboiler_duty.unit),
        iterations=int(result.iterations),
        temperature_residual=float(result.temperature_residual),
        mass_residual=float(result.mass_residual),
        energy_residual=float(result.energy_residual),
        fs_factor=float(result.fs_factor),
        fs_factor_utilization=float(result.fs_factor_utilization),
        fs_factor_within_design_limit=bool(result.fs_factor_within_design_limit),
        minimum_diameter_for_fs_limit=from_si(
            result.minimum_diameter_for_fs_limit.magnitude_si,
            result.minimum_diameter_for_fs_limit.unit,
        ),
        required_diameter=from_si(
            result.required_diameter.magnitude_si, result.required_diameter.unit
        ),
        controlling_tray_index=int(result.controlling_tray_index),
        internals_design_ok=bool(result.internals_design_ok),
        max_percent_flood=float(result.max_percent_flood),
        min_percent_flood=float(result.min_percent_flood),
        average_tray_efficiency=float(result.average_tray_efficiency),
        total_pressure_drop=from_si(
            result.total_pressure_drop.magnitude_si, result.total_pressure_drop.unit
        ),
        total_pressure_drop_mbar=float(result.total_pressure_drop_mbar),
        tray_percent_flood=tuple(float(value) for value in result.tray_percent_flood),
        tray_pressure_drop=tuple(
            from_si(value.magnitude_si, value.unit) for value in result.tray_pressure_drop
        ),
        tray_efficiency=tuple(float(value) for value in result.tray_efficiency),
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
        duty=from_si(result.duty.magnitude_si, result.duty.unit),
        ntu=result.ntu,
        effectiveness=result.effectiveness,
        # **A quantity the run did not reach is `None`**, which the transport has already said:
        # `from_si` is only reached where the Rust side published a magnitude, so a pinned-outlet
        # run arrives with five `None`s rather than with zeros the dataclass would have to guess at.
        c_min=(
            None if result.c_min is None else from_si(result.c_min.magnitude_si, result.c_min.unit)
        ),
        c_max=(
            None if result.c_max is None else from_si(result.c_max.magnitude_si, result.c_max.unit)
        ),
        capacity_ratio=result.capacity_ratio,
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
    reactive: bool | None = None,
    reactive_start_tray: int | None = None,
    reactive_end_tray: int | None = None,
    solver_type: str | None = None,
    tray_murphree_efficiency: Sequence[float] | None = None,
    gas_side_draw_fractions: Sequence[float] | None = None,
    liquid_side_draw_fractions: Sequence[float] | None = None,
    pumparound_fractions: Sequence[float] | None = None,
    side_draw_flow_tray: int | None = None,
    side_draw_flow_phase: str | None = None,
    side_draw_flow_target: object | None = None,
    side_draw_flow_tolerance: float | None = None,
    side_draw_flow_max_iterations: int | None = None,
    pumparound_return_tray: int | None = None,
    pumparound_draw_tray: int | None = None,
    pumparound_draw_fraction: float | None = None,
    pumparound_temperature_drop: object | None = None,
    pumparound_tolerance: float | None = None,
    pumparound_max_iterations: int | None = None,
    column_diameter: Q | None = None,
    max_allowable_fs_factor: float | None = None,
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
        reactive,
        reactive_start_tray,
        reactive_end_tray,
        solver_type,
        None if tray_murphree_efficiency is None else [float(v) for v in tray_murphree_efficiency],
        None if gas_side_draw_fractions is None else [float(v) for v in gas_side_draw_fractions],
        None
        if liquid_side_draw_fractions is None
        else [float(v) for v in liquid_side_draw_fractions],
        None if pumparound_fractions is None else [float(v) for v in pumparound_fractions],
        side_draw_flow_tray,
        side_draw_flow_phase,
        None
        if side_draw_flow_target is None
        else input_to_si(spec, "side_draw_flow_target", side_draw_flow_target),
        side_draw_flow_tolerance,
        side_draw_flow_max_iterations,
        pumparound_return_tray,
        pumparound_draw_tray,
        pumparound_draw_fraction,
        None
        if pumparound_temperature_drop is None
        else input_to_si(spec, "pumparound_temperature_drop", pumparound_temperature_drop),
        pumparound_tolerance,
        pumparound_max_iterations,
        None if column_diameter is None else input_to_si(spec, "column_diameter", column_diameter),
        max_allowable_fs_factor,
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
        temperature_residual=from_si(
            result.temperature_residual.magnitude_si, result.temperature_residual.unit
        ),
        mass_residual=result.mass_residual,
        energy_residual=result.energy_residual,
        fs_factor=float(result.fs_factor),
        fs_factor_utilization=float(result.fs_factor_utilization),
        fs_factor_within_design_limit=bool(result.fs_factor_within_design_limit),
        minimum_diameter_for_fs_limit=from_si(
            result.minimum_diameter_for_fs_limit.magnitude_si,
            result.minimum_diameter_for_fs_limit.unit,
        ),
        gas_load_factor=from_si(result.gas_load_factor.magnitude_si, result.gas_load_factor.unit),
        gas_load_factor_utilization=float(result.gas_load_factor_utilization),
        gas_load_factor_within_design_limit=bool(result.gas_load_factor_within_design_limit),
        minimum_diameter_for_gas_load_limit=from_si(
            result.minimum_diameter_for_gas_load_limit.magnitude_si,
            result.minimum_diameter_for_gas_load_limit.unit,
        ),
        warnings=_warnings(result.warnings),
    )


def packed_column(
    components: Sequence[str],
    feed_n: Q,
    feed_z: Sequence[float],
    feed_p: Q,
    feed_t: Q,
    packed_height: Q,
    feed_stage: int,
    has_reboiler: bool,
    has_condenser: bool,
    top_pressure: Q,
    bottom_pressure: Q,
    temperature_tolerance: float,
    max_iterations: int,
    reboiler_temperature: Q | None = None,
    condenser_temperature: Q | None = None,
    packing_type: str | None = None,
    structured_packing: bool | None = None,
    design_flood_fraction: float | None = None,
    packing_hydraulic_capacity_factor: float | None = None,
    column_diameter: Q | None = None,
    murphree_efficiency: float | None = None,
    reactive: bool | None = None,
    reactive_start_tray: int | None = None,
    reactive_end_tray: int | None = None,
    solver_type: str | None = None,
    top_specification_type: str | None = None,
    top_specification_target: float | None = None,
    top_specification_component: str | None = None,
    bottom_specification_type: str | None = None,
    bottom_specification_target: float | None = None,
    bottom_specification_component: str | None = None,
    tray_murphree_efficiency: Sequence[float] | None = None,
    gas_side_draw_fractions: Sequence[float] | None = None,
    liquid_side_draw_fractions: Sequence[float] | None = None,
    pumparound_fractions: Sequence[float] | None = None,
    side_draw_flow_tray: int | None = None,
    side_draw_flow_phase: str | None = None,
    side_draw_flow_target: object | None = None,
    side_draw_flow_tolerance: float | None = None,
    side_draw_flow_max_iterations: int | None = None,
    pumparound_return_tray: int | None = None,
    pumparound_draw_tray: int | None = None,
    pumparound_draw_fraction: float | None = None,
    pumparound_temperature_drop: object | None = None,
    pumparound_tolerance: float | None = None,
    pumparound_max_iterations: int | None = None,
    max_allowable_fs_factor: float | None = None,
) -> PackedColumnResult:
    """Solve a packed column; see :func:`azoth.process.reference.packed_column`.

    The base column's own solve at the stage count this packed height's constructor derives.
    """
    spec = _models_gen.model("process.packed_column")
    result = _core.packed_column(
        list(components),
        input_to_si(spec, "feed_n", feed_n),
        [_si(spec, "feed_z", v) for v in feed_z],
        input_to_si(spec, "feed_p", feed_p),
        input_to_si(spec, "feed_t", feed_t),
        input_to_si(spec, "packed_height", packed_height),
        int(_si(spec, "feed_stage", feed_stage)),
        has_reboiler,
        has_condenser,
        input_to_si(spec, "top_pressure", top_pressure),
        input_to_si(spec, "bottom_pressure", bottom_pressure),
        _si(spec, "temperature_tolerance", temperature_tolerance),
        int(_si(spec, "max_iterations", max_iterations)),
        None
        if reboiler_temperature is None
        else input_to_si(spec, "reboiler_temperature", reboiler_temperature),
        None
        if condenser_temperature is None
        else input_to_si(spec, "condenser_temperature", condenser_temperature),
        packing_type,
        structured_packing,
        design_flood_fraction,
        packing_hydraulic_capacity_factor,
        None if column_diameter is None else input_to_si(spec, "column_diameter", column_diameter),
        murphree_efficiency,
        reactive,
        reactive_start_tray,
        reactive_end_tray,
        solver_type,
        top_specification_type,
        top_specification_target,
        top_specification_component,
        bottom_specification_type,
        bottom_specification_target,
        bottom_specification_component,
        None if tray_murphree_efficiency is None else [float(v) for v in tray_murphree_efficiency],
        None if gas_side_draw_fractions is None else [float(v) for v in gas_side_draw_fractions],
        None
        if liquid_side_draw_fractions is None
        else [float(v) for v in liquid_side_draw_fractions],
        None if pumparound_fractions is None else [float(v) for v in pumparound_fractions],
        side_draw_flow_tray,
        side_draw_flow_phase,
        None
        if side_draw_flow_target is None
        else input_to_si(spec, "side_draw_flow_target", side_draw_flow_target),
        side_draw_flow_tolerance,
        side_draw_flow_max_iterations,
        pumparound_return_tray,
        pumparound_draw_tray,
        pumparound_draw_fraction,
        None
        if pumparound_temperature_drop is None
        else input_to_si(spec, "pumparound_temperature_drop", pumparound_temperature_drop),
        pumparound_tolerance,
        pumparound_max_iterations,
        max_allowable_fs_factor,
    )
    return PackedColumnResult(
        gas_side_draw_n=tuple(
            from_si(value.magnitude_si, value.unit) for value in result.gas_side_draw_n
        ),
        liquid_side_draw_n=tuple(
            from_si(value.magnitude_si, value.unit) for value in result.liquid_side_draw_n
        ),
        pumparound_n=tuple(
            from_si(value.magnitude_si, value.unit) for value in result.pumparound_n
        ),
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
        iterations=result.iterations,
        temperature_residual=result.temperature_residual,
        mass_residual=result.mass_residual,
        energy_residual=result.energy_residual,
        fs_factor=float(result.fs_factor),
        fs_factor_utilization=float(result.fs_factor_utilization),
        fs_factor_within_design_limit=bool(result.fs_factor_within_design_limit),
        minimum_diameter_for_fs_limit=from_si(
            result.minimum_diameter_for_fs_limit.magnitude_si,
            result.minimum_diameter_for_fs_limit.unit,
        ),
        hetp=from_si(result.hetp.magnitude_si, result.hetp.unit),
        theoretical_stages=result.theoretical_stages,
        percent_flood=result.percent_flood,
        flooding_velocity=from_si(
            result.flooding_velocity.magnitude_si, result.flooding_velocity.unit
        ),
        packing_pressure_drop=from_si(
            result.packing_pressure_drop.magnitude_si, result.packing_pressure_drop.unit
        ),
        hydraulics_ok=result.hydraulics_ok,
        internal_diameter=from_si(
            result.internal_diameter.magnitude_si, result.internal_diameter.unit
        ),
        warnings=_warnings(result.warnings),
    )


def stripping_column(
    stripping_gas_components: Sequence[str],
    rich_liquid_components: Sequence[str],
    stripping_gas_n: Q,
    stripping_gas_z: Sequence[float],
    stripping_gas_p: Q,
    stripping_gas_t: Q,
    rich_liquid_n: Q,
    rich_liquid_z: Sequence[float],
    rich_liquid_p: Q,
    rich_liquid_t: Q,
    number_of_stages: int,
    top_pressure: Q,
    bottom_pressure: Q,
    temperature_tolerance: float,
    max_iterations: int,
    tray_temperatures: Sequence[float] | None = None,
    murphree_efficiency: float | None = None,
    component_murphree_efficiency: Sequence[float] | None = None,
    max_allowable_gas_load_factor: float | None = None,
    reactive: bool | None = None,
    reactive_start_tray: int | None = None,
    reactive_end_tray: int | None = None,
    solver_type: str | None = None,
    tray_murphree_efficiency: Sequence[float] | None = None,
    gas_side_draw_fractions: Sequence[float] | None = None,
    liquid_side_draw_fractions: Sequence[float] | None = None,
    pumparound_fractions: Sequence[float] | None = None,
    side_draw_flow_tray: int | None = None,
    side_draw_flow_phase: str | None = None,
    side_draw_flow_target: object | None = None,
    side_draw_flow_tolerance: float | None = None,
    side_draw_flow_max_iterations: int | None = None,
    pumparound_return_tray: int | None = None,
    pumparound_draw_tray: int | None = None,
    pumparound_draw_fraction: float | None = None,
    pumparound_temperature_drop: object | None = None,
    pumparound_tolerance: float | None = None,
    pumparound_max_iterations: int | None = None,
    column_diameter: Q | None = None,
    max_allowable_fs_factor: float | None = None,
) -> StrippingColumnResult:
    """Solve a tray stripper; see :func:`azoth.process.reference.stripping_column`.

    The base's own call with this class's names for the two inlets and the two products.
    """
    spec = _models_gen.model("process.stripping_column")
    result = _core.stripping_column(
        list(stripping_gas_components),
        list(rich_liquid_components),
        input_to_si(spec, "stripping_gas_n", stripping_gas_n),
        [_si(spec, "stripping_gas_z", v) for v in stripping_gas_z],
        input_to_si(spec, "stripping_gas_p", stripping_gas_p),
        input_to_si(spec, "stripping_gas_t", stripping_gas_t),
        input_to_si(spec, "rich_liquid_n", rich_liquid_n),
        [_si(spec, "rich_liquid_z", v) for v in rich_liquid_z],
        input_to_si(spec, "rich_liquid_p", rich_liquid_p),
        input_to_si(spec, "rich_liquid_t", rich_liquid_t),
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
        reactive,
        reactive_start_tray,
        reactive_end_tray,
        solver_type,
        None if tray_murphree_efficiency is None else [float(v) for v in tray_murphree_efficiency],
        None if gas_side_draw_fractions is None else [float(v) for v in gas_side_draw_fractions],
        None
        if liquid_side_draw_fractions is None
        else [float(v) for v in liquid_side_draw_fractions],
        None if pumparound_fractions is None else [float(v) for v in pumparound_fractions],
        side_draw_flow_tray,
        side_draw_flow_phase,
        None
        if side_draw_flow_target is None
        else input_to_si(spec, "side_draw_flow_target", side_draw_flow_target),
        side_draw_flow_tolerance,
        side_draw_flow_max_iterations,
        pumparound_return_tray,
        pumparound_draw_tray,
        pumparound_draw_fraction,
        None
        if pumparound_temperature_drop is None
        else input_to_si(spec, "pumparound_temperature_drop", pumparound_temperature_drop),
        pumparound_tolerance,
        pumparound_max_iterations,
        None if column_diameter is None else input_to_si(spec, "column_diameter", column_diameter),
        max_allowable_fs_factor,
    )
    return StrippingColumnResult(
        tray_temperature=tuple(from_si(q.magnitude_si, q.unit) for q in result.tray_temperature),
        tray_pressure=tuple(from_si(q.magnitude_si, q.unit) for q in result.tray_pressure),
        tray_gas_n=tuple(from_si(q.magnitude_si, q.unit) for q in result.tray_gas_n),
        tray_liquid_n=tuple(from_si(q.magnitude_si, q.unit) for q in result.tray_liquid_n),
        overhead_gas_n=from_si(result.overhead_gas_n.magnitude_si, result.overhead_gas_n.unit),
        overhead_gas_z=tuple(result.overhead_gas_z),
        overhead_gas_p=from_si(result.overhead_gas_p.magnitude_si, result.overhead_gas_p.unit),
        overhead_gas_t=from_si(result.overhead_gas_t.magnitude_si, result.overhead_gas_t.unit),
        overhead_gas_h=from_si(result.overhead_gas_h.magnitude_si, result.overhead_gas_h.unit),
        lean_liquid_n=from_si(result.lean_liquid_n.magnitude_si, result.lean_liquid_n.unit),
        lean_liquid_z=tuple(result.lean_liquid_z),
        lean_liquid_p=from_si(result.lean_liquid_p.magnitude_si, result.lean_liquid_p.unit),
        lean_liquid_t=from_si(result.lean_liquid_t.magnitude_si, result.lean_liquid_t.unit),
        lean_liquid_h=from_si(result.lean_liquid_h.magnitude_si, result.lean_liquid_h.unit),
        iterations=result.iterations,
        temperature_residual=from_si(
            result.temperature_residual.magnitude_si, result.temperature_residual.unit
        ),
        mass_residual=result.mass_residual,
        energy_residual=result.energy_residual,
        fs_factor=float(result.fs_factor),
        fs_factor_utilization=float(result.fs_factor_utilization),
        fs_factor_within_design_limit=bool(result.fs_factor_within_design_limit),
        minimum_diameter_for_fs_limit=from_si(
            result.minimum_diameter_for_fs_limit.magnitude_si,
            result.minimum_diameter_for_fs_limit.unit,
        ),
        gas_load_factor=from_si(result.gas_load_factor.magnitude_si, result.gas_load_factor.unit),
        gas_load_factor_utilization=float(result.gas_load_factor_utilization),
        gas_load_factor_within_design_limit=bool(result.gas_load_factor_within_design_limit),
        minimum_diameter_for_gas_load_limit=from_si(
            result.minimum_diameter_for_gas_load_limit.magnitude_si,
            result.minimum_diameter_for_gas_load_limit.unit,
        ),
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


def _optional_si(quantity: object) -> Q | None:
    """A quantity the run may not have reached, as the dataclass spells it.

    The Rust transport writes `None` for a number a run did not reach - an ejector handed no
    flow, an exchanger's rating where one outlet was pinned - so the bridge's job here is to
    keep that absence rather than to invent a zero for it.
    """
    if quantity is None:
        return None
    magnitude_si = quantity.magnitude_si  # type: ignore[attr-defined]
    unit = quantity.unit  # type: ignore[attr-defined]
    return from_si(magnitude_si, unit)


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
        # **A number the run did not reach is `None`**, which the transport has already said: a
        # machine handed no flow reached no mixing pressure and no velocity, and `from_si` is only
        # reached where the Rust side published a magnitude.
        mixing_pressure=_optional_si(result.mixing_pressure),
        motive_nozzle_velocity=_optional_si(result.motive_nozzle_velocity),
        suction_nozzle_velocity=_optional_si(result.suction_nozzle_velocity),
        mixing_velocity=_optional_si(result.mixing_velocity),
        diffuser_velocity=_optional_si(result.diffuser_velocity),
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
        phase_fractions=tuple(result.phase_fractions),
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
    cubic: str | None = None,
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
        "srk" if cubic is None else cubic,
    )
    return ReactiveTpFlashResult(
        phase_count=result.phase_count,
        phase_type=tuple(result.phase_type),
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
    cubic: str | None = None,
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
        "srk" if cubic is None else cubic,
    )
    return ReactivePhFlashResult(
        temperature=from_si(result.temperature.magnitude_si, result.temperature.unit),
        converged=result.converged,
        outer_iterations=result.outer_iterations,
        total_inner_iterations=result.total_inner_iterations,
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


def rate_based_packed_column(
    gas_components: Sequence[str],
    gas_n: Q,
    gas_z: Sequence[float],
    gas_p: Q,
    gas_t: Q,
    liquid_components: Sequence[str],
    liquid_n: Q,
    liquid_z: Sequence[float],
    liquid_p: Q,
    liquid_t: Q,
    transfer_components: Sequence[str] | None = None,
    column_diameter: Q | None = None,
    packed_height: Q | None = None,
    number_of_segments: int | None = None,
    packing_type: str | None = None,
    max_iterations: int | None = None,
    convergence_tolerance: Q | None = None,
    mass_transfer_correction: float | None = None,
    mass_transfer_correlation: str | None = None,
    film_model: str | None = None,
    heat_transfer_model: str | None = None,
    segment_solver: str | None = None,
    column_solver: str | None = None,
) -> RateBasedPackedColumnResult:
    """Solve a rate-based packed column; see the reference twin for the arithmetic."""
    spec = _models_gen.model("process.rate_based_packed_column")

    def si(name: str, value: Q | None) -> float | None:
        return None if value is None else input_to_si(spec, name, value)

    result = _core.rate_based_packed_column(
        list(gas_components),
        list(liquid_components),
        [] if transfer_components is None else list(transfer_components),
        input_to_si(spec, "gas_n", gas_n),
        [_si(spec, "gas_z", v) for v in gas_z],
        input_to_si(spec, "gas_p", gas_p),
        input_to_si(spec, "gas_t", gas_t),
        input_to_si(spec, "liquid_n", liquid_n),
        [_si(spec, "liquid_z", v) for v in liquid_z],
        input_to_si(spec, "liquid_p", liquid_p),
        input_to_si(spec, "liquid_t", liquid_t),
        si("column_diameter", column_diameter),
        si("packed_height", packed_height),
        None if number_of_segments is None else int(number_of_segments),
        packing_type,
        None if max_iterations is None else int(max_iterations),
        si("convergence_tolerance", convergence_tolerance),
        mass_transfer_correction,
        mass_transfer_correlation,
        film_model,
        heat_transfer_model,
        segment_solver,
        column_solver,
    )

    def q(value: Any) -> Q:
        return from_si(value.magnitude_si, value.unit)

    def qs(values: Any) -> tuple[Q, ...]:
        return tuple(from_si(value.magnitude_si, value.unit) for value in values)

    return RateBasedPackedColumnResult(
        gas_out_n=q(result.gas_out_n),
        gas_out_z=tuple(result.gas_out_z),
        gas_out_p=q(result.gas_out_p),
        gas_out_t=q(result.gas_out_t),
        gas_out_h=q(result.gas_out_h),
        liquid_out_n=q(result.liquid_out_n),
        liquid_out_z=tuple(result.liquid_out_z),
        liquid_out_p=q(result.liquid_out_p),
        liquid_out_t=q(result.liquid_out_t),
        liquid_out_h=q(result.liquid_out_h),
        iterations=result.iterations,
        convergence_residual=q(result.convergence_residual),
        converged=result.converged,
        total_absolute_molar_transfer=q(result.total_absolute_molar_transfer),
        component_transfer_totals=qs(result.component_transfer_totals),
        transfer_components=tuple(result.transfer_components),
        segment_height_from_bottom=qs(result.segment_height_from_bottom),
        segment_gas_temperature=qs(result.segment_gas_temperature),
        segment_liquid_temperature=qs(result.segment_liquid_temperature),
        segment_gas_pressure=qs(result.segment_gas_pressure),
        segment_liquid_pressure=qs(result.segment_liquid_pressure),
        segment_gas_molar_flow=qs(result.segment_gas_molar_flow),
        segment_liquid_molar_flow=qs(result.segment_liquid_molar_flow),
        segment_gas_density=qs(result.segment_gas_density),
        segment_liquid_density=qs(result.segment_liquid_density),
        segment_gas_viscosity=qs(result.segment_gas_viscosity),
        segment_liquid_viscosity=qs(result.segment_liquid_viscosity),
        segment_gas_diffusivity=qs(result.segment_gas_diffusivity),
        segment_liquid_diffusivity=qs(result.segment_liquid_diffusivity),
        segment_wetted_area=tuple(result.segment_wetted_area),
        segment_k_ga=tuple(result.segment_k_ga),
        segment_k_la=tuple(result.segment_k_la),
        segment_gas_heat_transfer_coefficient=tuple(result.segment_gas_heat_transfer_coefficient),
        segment_liquid_heat_transfer_coefficient=tuple(
            result.segment_liquid_heat_transfer_coefficient
        ),
        segment_overall_heat_transfer_coefficient=tuple(
            result.segment_overall_heat_transfer_coefficient
        ),
        segment_interface_temperature=qs(result.segment_interface_temperature),
        segment_heat_transfer_rate=qs(result.segment_heat_transfer_rate),
        segment_pressure_drop_per_meter=qs(result.segment_pressure_drop_per_meter),
        segment_percent_flood=tuple(result.segment_percent_flood),
        segment_net_molar_transfer=qs(result.segment_net_molar_transfer),
        segment_enthalpy_balance_residual=qs(result.segment_enthalpy_balance_residual),
        warnings=_warnings(result.warnings),
    )
