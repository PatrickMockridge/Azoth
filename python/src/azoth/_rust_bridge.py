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

# **The generated half.** Every id whose signature is the public wrapper's and whose boundary is
# the spec's own inputs is emitted by `tools/gen_python_bridge.py`; what is left in this file is
# the ids that take a `Mixture` or a record the spec does not name. `resolve` derives the id's
# function name and looks it up in this module's globals, which the star import fills.
from azoth._rust_bridge_gen import *  # noqa: F403

# The helpers the generated adapters and the hand-written ones both call. Imported rather than
# defined here because `azoth._rust_bridge_gen` imports them too, and a module cannot import the
# file that re-exports it. `from ... import _warnings` keeps the name reachable as
# `azoth._rust_bridge._warnings`, which is where `azoth.batch._core` reads it.
from azoth.core.bridge_support import _si as _si
from azoth.core.bridge_support import _warnings as _warnings
from azoth.core.errors import PropertyUnavailableError
from azoth.core.result import (
    AqueousViscosityResult,
    BubblePressureResult,
    BubbleTemperatureResult,
    BwrsPhaseResult,
    CapillaryDewPointResult,
    CriticalPointResult,
    DewPressureResult,
    DewTemperatureResult,
    GeNrtlFlashResult,
    GeNrtlPhaseResult,
    GeUnifacPhaseResult,
    GeUniquacPhaseResult,
    GeVanLaarAcidPhaseResult,
    GeWilsonPhaseResult,
    HydrogenPhaseResult,
    MolarEnthalpyEntropyResult,
    NrtlActivityCoefficientsResult,
    OrificeFlowResult,
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
    StabilityTestResult,
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
