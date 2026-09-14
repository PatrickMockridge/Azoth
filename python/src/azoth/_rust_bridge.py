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
from azoth.core.result import (
    BubblePressureResult,
    ChokedFlowAreaResult,
    ColebrookResult,
    CompressorResult,
    ConductionPlaneWallResult,
    ControlValveCvResult,
    CriticalPointResult,
    DarcyWeisbachResult,
    DewPressureResult,
    ExpanderResult,
    FlowRegime,
    HaalandResult,
    HeaterResult,
    IdealGasCpResult,
    KComponent,
    KFactorsResult,
    MixerResult,
    MolarEnthalpyEntropyResult,
    OrificeFlowResult,
    PhFlashResult,
    PrAlphaAbResult,
    PrDepartureResult,
    PrKappaResult,
    PrMassDensityResult,
    PrMolarVolumeResult,
    PrsvKappaResult,
    PrZFactorResult,
    PsFlashResult,
    PtFlashResult,
    PumpPowerResult,
    PumpResult,
    PureSaturationResult,
    RachfordRiceBinaryResult,
    ReynoldsNumberResult,
    RootStructure,
    SeparatorResult,
    SplitterResult,
    StabilityTestResult,
    SwameeJainResult,
    ThrottlingValveResult,
    Vdw1fMixBinaryResult,
)
from azoth.core.result import Phase as _Phase
from azoth.core.result import StabilityVerdict as _StabilityVerdict
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
    )
    return StabilityTestResult(
        verdict=_StabilityVerdict(result.verdict),
        tm=tuple(result.tm),
        w=tuple(tuple(row) for row in result.w),
        iterations=tuple(result.iterations),
        min_t_over_tc=result.min_t_over_tc,
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
    )
    return _boundary_result(raw, BubblePressureResult, liquid_first=True)  # type: ignore[no-any-return]


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
    )
    return _boundary_result(raw, DewPressureResult, liquid_first=False)  # type: ignore[no-any-return]


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


def separator(
    mixture: Any,
    ideal_gas: Any,
    T: Q,
    P: Q,
    n: Q,
    z: Sequence[float],
    pressure_drop: Q,
    heat_duty: Q,
) -> SeparatorResult:
    """One feed split into a gas and a liquid at a single state, run in Rust.

    The first bridge function for a unit operation, and the shape every one after it
    takes: the mixture crosses as the three per-component vectors plus a flattened
    `kij`, exactly as `pt_flash` sends them, because it *is* the same object.

    `n` crosses as `mol/s` and is the one input here that arrives as a bare `f64`
    magnitude with no `uom` type behind it on the Rust side - the units vocabulary has
    the name and the conversion is the identity, which is what makes that legitimate
    rather than a quantity smuggled across outside the vocabulary.
    """
    spec = _models_gen.model("process.separator")
    result = _core.separator(
        *_model_arguments(mixture, ideal_gas),
        input_to_si(spec, "T", T),
        input_to_si(spec, "P", P),
        input_to_si(spec, "n", n),
        list(z),
        input_to_si(spec, "pressure_drop", pressure_drop),
        input_to_si(spec, "heat_duty", heat_duty),
    )
    return SeparatorResult(
        T=from_si(result.T.magnitude_si, result.T.unit),
        P=from_si(result.P.magnitude_si, result.P.unit),
        beta=result.beta,
        gas_flow=result.gas_flow,
        gas_z=tuple(result.gas_z),
        liquid_flow=result.liquid_flow,
        liquid_z=tuple(result.liquid_z),
        phase=_Phase(result.phase),
        iterations=result.iterations,
        warnings=_warnings(result.warnings),
    )


def _mixture_arguments(
    mixture: Any,
) -> tuple[list[float], list[float], list[float], list[float]]:
    """The four vectors every model's mixture crosses as.

    One function rather than four lines repeated in eight unit operations. It is the same
    argument the crate's `build_mixture` makes on the Rust side: the component order is
    one thing to agree about, so it should be decided in one place.

    A fixed-length tuple rather than a list, because the caller splats it into a function
    of known arity: a list would leave the argument count unchecked at the call site.
    """
    return (
        [c.Tc.to_base_units().magnitude for c in mixture.components],
        [c.Pc.to_base_units().magnitude for c in mixture.components],
        [c.omega for c in mixture.components],
        mixture.flattened_kij(),
    )


def _ideal_gas_arguments(
    ideal_gas: Any,
) -> tuple[
    list[float],
    list[float],
    list[float],
    list[float],
    list[float],
]:
    """The five coefficients an ideal-gas model crosses as.

    There is no reference state among them: NeqSim's ideal-gas integrals start from a
    fixed temperature and pressure, so the coefficients are the whole of the model.
    """
    return (
        list(ideal_gas.cp_a),
        list(ideal_gas.cp_b),
        list(ideal_gas.cp_c),
        list(ideal_gas.cp_d),
        list(ideal_gas.cp_e),
    )


def _model_arguments(
    mixture: Any, ideal_gas: Any
) -> tuple[
    list[float],
    list[float],
    list[float],
    list[float],
    list[float],
    list[float],
    list[float],
    list[float],
    list[float],
]:
    """The leading arguments every unit operation takes, as one tuple: the mixture and
    the ideal-gas model.

    The concatenation of the two above, and not two separate unpacks at the call site:
    a call written `f(*a, *b, ...)` with both `a` and `b` of known length is one mypy
    cannot check the arity of - it aligns the first tuple and gives up on the second,
    reporting a mismatch that is not there. One tuple is checked correctly, so the
    argument count at the call site is verified rather than assumed.
    """
    return _mixture_arguments(mixture) + _ideal_gas_arguments(ideal_gas)


def mixer(
    mixture: Any,
    ideal_gas: Any,
    T: Sequence[Q],
    P: Sequence[Q],
    n: Sequence[Q],
    z: Sequence[Sequence[float]],
) -> MixerResult:
    """Several feeds blended into one, run in Rust.

    The first bridge function whose arguments are *vectors*: one entry per inlet, and a
    composition per inlet. `z` crosses flattened row-major, exactly as `kij` does and for
    the same reason - the boundary carries one list rather than a list of lists, and the
    Rust side derives the stream count from how many temperatures it was given.
    """
    spec = _models_gen.model("process.mixer")
    result = _core.mixer(
        *_model_arguments(mixture, ideal_gas),
        [input_to_si(spec, "T", value) for value in T],
        [input_to_si(spec, "P", value) for value in P],
        [input_to_si(spec, "n", value) for value in n],
        [float(value) for row in z for value in row],
    )
    return MixerResult(
        T=from_si(result.T.magnitude_si, result.T.unit),
        P=from_si(result.P.magnitude_si, result.P.unit),
        flow=result.flow,
        z_out=tuple(result.z_out),
        beta=result.beta,
        phase=_Phase(result.phase),
        iterations=result.iterations,
        warnings=_warnings(result.warnings),
    )


def splitter(
    mixture: Any,
    T: Q,
    P: Q,
    n: Q,
    z: Sequence[float],
    fractions: Sequence[float],
) -> SplitterResult:
    """One feed divided into branches, run in Rust.

    The only bridge function that sends no ideal-gas block, because the only unit
    operation that takes none - a splitter does no energy balance.
    """
    spec = _models_gen.model("process.splitter")
    result = _core.splitter(
        *_mixture_arguments(mixture),
        input_to_si(spec, "T", T),
        input_to_si(spec, "P", P),
        input_to_si(spec, "n", n),
        list(z),
        list(fractions),
    )
    return SplitterResult(
        T=from_si(result.T.magnitude_si, result.T.unit),
        P=from_si(result.P.magnitude_si, result.P.unit),
        phase=_Phase(result.phase),
        beta=result.beta,
        flows=tuple(result.flows),
        iterations=result.iterations,
        warnings=_warnings(result.warnings),
    )


def throttling_valve(
    mixture: Any,
    ideal_gas: Any,
    T: Q,
    P: Q,
    z: Sequence[float],
    pressure_drop: Q,
) -> ThrottlingValveResult:
    """A stream's pressure dropped at constant enthalpy, run in Rust.

    No flow rate crosses, because the model takes none: an isenthalpic flash is a molar
    property.
    """
    spec = _models_gen.model("process.throttling_valve")
    result = _core.throttling_valve(
        *_model_arguments(mixture, ideal_gas),
        input_to_si(spec, "T", T),
        input_to_si(spec, "P", P),
        list(z),
        input_to_si(spec, "pressure_drop", pressure_drop),
    )
    return ThrottlingValveResult(
        T=from_si(result.T.magnitude_si, result.T.unit),
        P=from_si(result.P.magnitude_si, result.P.unit),
        phase=_Phase(result.phase),
        beta=result.beta,
        iterations=result.iterations,
        warnings=_warnings(result.warnings),
    )


def heater(
    mixture: Any,
    ideal_gas: Any,
    T: Q,
    P: Q,
    n: Q,
    z: Sequence[float],
    pressure_drop: Q,
    heat_duty: Q,
) -> HeaterResult:
    """A duty applied at a fixed pressure, run in Rust. A negative duty is a cooler."""
    spec = _models_gen.model("process.heater")
    result = _core.heater(
        *_model_arguments(mixture, ideal_gas),
        input_to_si(spec, "T", T),
        input_to_si(spec, "P", P),
        input_to_si(spec, "n", n),
        list(z),
        input_to_si(spec, "pressure_drop", pressure_drop),
        input_to_si(spec, "heat_duty", heat_duty),
    )
    return HeaterResult(
        T=from_si(result.T.magnitude_si, result.T.unit),
        P=from_si(result.P.magnitude_si, result.P.unit),
        phase=_Phase(result.phase),
        beta=result.beta,
        iterations=result.iterations,
        warnings=_warnings(result.warnings),
    )


def compressor(
    mixture: Any,
    ideal_gas: Any,
    T: Q,
    P: Q,
    n: Q,
    z: Sequence[float],
    outlet_pressure: Q,
    efficiency: float,
) -> CompressorResult:
    """A pressure rise at a stated isentropic efficiency, run in Rust."""
    spec = _models_gen.model("process.compressor")
    result = _core.compressor(
        *_model_arguments(mixture, ideal_gas),
        input_to_si(spec, "T", T),
        input_to_si(spec, "P", P),
        input_to_si(spec, "n", n),
        list(z),
        input_to_si(spec, "outlet_pressure", outlet_pressure),
        efficiency,
    )
    return CompressorResult(
        T=from_si(result.T.magnitude_si, result.T.unit),
        P=from_si(result.P.magnitude_si, result.P.unit),
        power=result.power,
        beta=result.beta,
        phase=_Phase(result.phase),
        isentropic_temperature=from_si(
            result.isentropic_temperature.magnitude_si, result.isentropic_temperature.unit
        ),
        iterations=result.iterations,
        warnings=_warnings(result.warnings),
    )


def pump(
    mixture: Any,
    ideal_gas: Any,
    T: Q,
    P: Q,
    n: Q,
    z: Sequence[float],
    outlet_pressure: Q,
    efficiency: float,
) -> PumpResult:
    """A pressure rise in a liquid, run in Rust."""
    spec = _models_gen.model("process.pump")
    result = _core.pump(
        *_model_arguments(mixture, ideal_gas),
        input_to_si(spec, "T", T),
        input_to_si(spec, "P", P),
        input_to_si(spec, "n", n),
        list(z),
        input_to_si(spec, "outlet_pressure", outlet_pressure),
        efficiency,
    )
    return PumpResult(
        T=from_si(result.T.magnitude_si, result.T.unit),
        P=from_si(result.P.magnitude_si, result.P.unit),
        power=result.power,
        beta=result.beta,
        phase=_Phase(result.phase),
        isentropic_temperature=from_si(
            result.isentropic_temperature.magnitude_si, result.isentropic_temperature.unit
        ),
        iterations=result.iterations,
        warnings=_warnings(result.warnings),
    )


def expander(
    mixture: Any,
    ideal_gas: Any,
    T: Q,
    P: Q,
    n: Q,
    z: Sequence[float],
    outlet_pressure: Q,
    efficiency: float,
) -> ExpanderResult:
    """A pressure drop that produces work, run in Rust. `power` crosses negative."""
    spec = _models_gen.model("process.expander")
    result = _core.expander(
        *_model_arguments(mixture, ideal_gas),
        input_to_si(spec, "T", T),
        input_to_si(spec, "P", P),
        input_to_si(spec, "n", n),
        list(z),
        input_to_si(spec, "outlet_pressure", outlet_pressure),
        efficiency,
    )
    return ExpanderResult(
        T=from_si(result.T.magnitude_si, result.T.unit),
        P=from_si(result.P.magnitude_si, result.P.unit),
        power=result.power,
        beta=result.beta,
        phase=_Phase(result.phase),
        isentropic_temperature=from_si(
            result.isentropic_temperature.magnitude_si, result.isentropic_temperature.unit
        ),
        iterations=result.iterations,
        warnings=_warnings(result.warnings),
    )


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
