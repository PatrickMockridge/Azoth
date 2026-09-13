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

from azoth import _core
from azoth._registry_gen import spec as _spec_for
from azoth.core.result import (
    ChokedFlowAreaResult,
    ColebrookResult,
    ConductionPlaneWallResult,
    ControlValveCvResult,
    DarcyWeisbachResult,
    FlowRegime,
    HaalandResult,
    KComponent,
    KFactorsResult,
    OrificeFlowResult,
    PrAlphaAbResult,
    PrDepartureResult,
    PrKappaResult,
    PrsvKappaResult,
    PrZFactorResult,
    PumpPowerResult,
    RachfordRiceBinaryResult,
    ReynoldsNumberResult,
    RootStructure,
    SwameeJainResult,
    Vdw1fMixBinaryResult,
)
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

    Five dimensionless arguments and three dimensionless outputs, so there is
    nothing to convert and no unit string to keep in step.
    """
    result = _core.pr_departure(a_reduced, b_reduced, z, kappa, Tr)
    return PrDepartureResult(
        ln_phi=result.ln_phi,
        h_dep_rt=result.h_dep_rt,
        s_dep_r=result.s_dep_r,
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


#: Calc id -> the bridge function implementing it. Explicit rather than derived
#: from the function names, so a renamed id fails here at import rather than
#: resolving to the wrong calc.
_IMPLEMENTATIONS: dict[str, Callable[..., Any]] = {
    "hydraulics.reynolds_number": reynolds_number,
    "hydraulics.friction_factor_colebrook": friction_factor_colebrook,
    "hydraulics.friction_factor_swamee_jain": friction_factor_swamee_jain,
    "hydraulics.friction_factor_haaland": friction_factor_haaland,
    "hydraulics.crane_k_factors": crane_k_factors,
    "hydraulics.darcy_weisbach": darcy_weisbach,
    "hydraulics.pump_power": pump_power,
    "hydraulics.orifice_flow": orifice_flow,
    "hydraulics.control_valve_cv": control_valve_cv,
    "hydraulics.choked_flow_area": choked_flow_area,
    "thermal.conduction_plane_wall": conduction_plane_wall,
    "eos.pr_kappa": pr_kappa,
    "eos.pr_alpha_ab": pr_alpha_ab,
    "eos.pr_z_factor": pr_z_factor,
    "eos.prsv_kappa": prsv_kappa,
    "eos.pr_departure": pr_departure,
    "eos.vdw1f_mix_binary": vdw1f_mix_binary,
    "eos.rachford_rice_binary": rachford_rice_binary,
}


def resolve(calc_id: str) -> Callable[..., Any]:
    """The bridge function for a calc id.

    Raises:
        KeyError: if the calc has no Rust implementation. That is a programming
            error rather than a runtime condition - the reference and the
            extension are supposed to cover the same calcs, and a test asserts
            they do.
    """
    try:
        return _IMPLEMENTATIONS[calc_id]
    except KeyError:
        raise KeyError(
            f"no Rust implementation for {calc_id!r}; the extension covers "
            f"{sorted(_IMPLEMENTATIONS)}"
        ) from None
