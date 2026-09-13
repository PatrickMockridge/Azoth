"""Type stub for the compiled extension.

Hand-written, because the module is a Rust `cdylib` and there is nothing for a
type checker to read. Keeping it accurate is part of the contract: if the Rust
side gains a parameter, this file and the signature test both have to say so.

Nothing here is part of the public API. `azoth.hydraulics` is, and it returns
the result dataclasses in `azoth.core.result` rather than the transport types
below - see `azoth._rust_bridge`.
"""

from typing import final

@final
class DataFile:
    name: str
    path: str
    text: str

@final
class FittingRow:
    id: str
    family: str
    name: str
    n_ld: float
    f_t_basis: str
    citation: str
    verify_status: str
    source_ref: str | None
    source_locator: str | None

@final
class FluidRow:
    fluid: str
    temperature_c: float
    density_kg_m3: float
    dynamic_viscosity_pa_s: float
    citation: str
    verify_status: str

@final
class Qty:
    """A physical quantity: always SI magnitude, plus a canonical unit label."""

    magnitude_si: float
    unit: str
    def __init__(self, magnitude_si: float, unit: str) -> None: ...

@final
class Warning:
    """A warning, transported. `code` is the SCREAMING_SNAKE_CASE string."""

    code: str
    message: str
    field: str | None

@final
class BatchColumn:
    name: str
    unit: str
    values: list[float] | None
    labels: list[str | None] | None

@final
class BatchResult:
    columns: list[BatchColumn]
    warnings: list[list[Warning]]

@final
class KComponent:
    """One fitting's contribution to the total resistance coefficient."""

    fitting_id: str
    n_ld: float
    k: float

@final
class ReynoldsNumberResult:
    re: float
    regime: str
    warnings: list[Warning]

@final
class ColebrookResult:
    f: float
    iterations: int
    converged: bool
    residual: float
    warnings: list[Warning]

@final
class SwameeJainResult:
    f: float
    warnings: list[Warning]

@final
class HaalandResult:
    f: float
    warnings: list[Warning]

@final
class ConductionPlaneWallResult:
    q: Qty
    warnings: list[Warning]

@final
class PrKappaResult:
    # Dimensionless, so a bare float rather than a `Qty`: there is no unit for the
    # two implementations to disagree about, which is what writing the equation of
    # state in reduced variables buys at the boundary.
    kappa: float
    warnings: list[Warning]

@final
class PureSaturationResult:
    p_sat: Qty
    ln_phi: float
    iterations: int
    residual: float
    warnings: list[Warning]

@final
class PtFlashResult:
    beta: float | None
    x: list[float]
    y: list[float]
    k: list[float]
    ln_phi_liquid: list[float]
    ln_phi_vapour: list[float]
    z_liquid: float
    z_vapour: float
    min_t_over_tc: float
    phase: str
    iterations: int
    residual: float
    warnings: list[Warning]

@final
class MolarEnthalpyEntropyResult:
    h: Qty
    s: Qty
    h_ideal: Qty
    s_ideal: Qty
    h_departure: Qty
    s_departure: Qty
    psi_bar: float
    warnings: list[Warning]

@final
class IdealGasCpResult:
    cp_over_r: float
    cp: Qty
    warnings: list[Warning]

@final
class CriticalPointResult:
    tc: Qty
    pc: Qty
    vc: Qty
    z_c: float
    iterations: int
    residual: float
    warnings: list[Warning]

@final
class PhaseBoundaryResult:
    pressure: Qty
    incipient: list[float]
    k: list[float]
    z_liquid: float
    z_vapour: float
    min_t_over_tc: float
    iterations: int
    residual: float
    warnings: list[Warning]

@final
class StabilityTestResult:
    # The spec's spelling of the verdict, as `phase` is for the flash: the adapter
    # rebuilds the enum, so the transport carries a string.
    verdict: str
    # Two entries, always, in trial order - the vapour-like trial first.
    tm: list[float]
    w: list[list[float]]
    iterations: list[int]
    min_t_over_tc: float
    warnings: list[Warning]

@final
class PrMolarVolumeResult:
    v: Qty
    warnings: list[Warning]

@final
class PrMassDensityResult:
    rho: Qty
    warnings: list[Warning]

@final
class Vdw1fMixBinaryResult:
    a_mix: float
    b_mix: float
    warnings: list[Warning]

@final
class RachfordRiceBinaryResult:
    beta: float
    warnings: list[Warning]

@final
class PrDepartureResult:
    ln_phi: float
    h_dep_rt: float
    s_dep_r: float
    warnings: list[Warning]

@final
class PrsvKappaResult:
    kappa: float
    warnings: list[Warning]

@final
class PrAlphaAbResult:
    # All three dimensionless, so bare floats for the same reason `PrKappaResult`
    # carries one.
    alpha: float
    a_reduced: float
    b_reduced: float
    warnings: list[Warning]

@final
class PrZFactorResult:
    z_min: float
    z_max: float
    # The spec's spelling, rebuilt into `azoth.core.result.RootStructure` by the
    # bridge - the same arrangement `ReynoldsNumberResult.regime` uses.
    root_structure: str
    iterations: int
    converged: bool
    residual: float
    warnings: list[Warning]

@final
class PumpPowerResult:
    power: Qty
    warnings: list[Warning]

@final
class ChokedFlowAreaResult:
    a: Qty
    warnings: list[Warning]

@final
class ControlValveCvResult:
    q: Qty
    warnings: list[Warning]

@final
class OrificeFlowResult:
    q: Qty
    warnings: list[Warning]

@final
class KFactorsResult:
    k_total: float
    f_t: float
    components: list[KComponent]
    warnings: list[Warning]

@final
class DarcyWeisbachResult:
    dp: Qty
    f: float
    re: float | None
    regime: str | None
    warnings: list[Warning]

# --- calculations ---------------------------------------------------------
# All arguments and returns are SI magnitudes; unit handling happens once, in
# Python, before the call crosses this boundary. See crates/azoth-python.

def reynolds_number(rho: float, v: float, D: float, mu: float) -> ReynoldsNumberResult: ...
def friction_factor_colebrook(re: float, relative_roughness: float) -> ColebrookResult: ...
def friction_factor_swamee_jain(re: float, relative_roughness: float) -> SwameeJainResult: ...
def friction_factor_haaland(re: float, relative_roughness: float) -> HaalandResult: ...
def conduction_plane_wall(k: float, A: float, dT: float, L: float) -> ConductionPlaneWallResult: ...
def pr_kappa(omega: float) -> PrKappaResult: ...
def pr_alpha_ab(kappa: float, Tr: float, Pr: float) -> PrAlphaAbResult: ...
def pr_z_factor(a_reduced: float, b_reduced: float) -> PrZFactorResult: ...
def prsv_kappa(omega: float, Tr: float, kappa1: float) -> PrsvKappaResult: ...
def pr_departure(
    a_reduced: float, b_reduced: float, z: float, kappa: float, Tr: float
) -> PrDepartureResult: ...
def vdw1f_mix_binary(
    z1: float, a1: float, a2: float, b1: float, b2: float, k12: float
) -> Vdw1fMixBinaryResult: ...
def rachford_rice_binary(z1: float, K1: float, K2: float) -> RachfordRiceBinaryResult: ...
def pr_molar_volume(z: float, T: float, P: float) -> PrMolarVolumeResult: ...
def pr_mass_density(M: float, v: float) -> PrMassDensityResult: ...
def pure_saturation(Tc: float, Pc: float, omega: float, T: float) -> PureSaturationResult: ...
def molar_enthalpy_entropy(
    Tc: list[float],
    Pc: list[float],
    omega: list[float],
    kij: list[float],
    cp_a: list[float],
    cp_b: list[float],
    cp_c: list[float],
    cp_d: list[float],
    h_ref: list[float],
    s_ref: list[float],
    T_ref: float,
    P_ref: float,
    T: float,
    P: float,
    z: list[float],
    compressibility: float,
) -> MolarEnthalpyEntropyResult: ...
def ideal_gas_cp(a: float, b: float, c: float, d: float, T: float) -> IdealGasCpResult: ...
def critical_point(
    Tc: list[float],
    Pc: list[float],
    omega: list[float],
    kij: list[float],
    z: list[float],
) -> CriticalPointResult: ...
def bubble_pressure(
    Tc: list[float],
    Pc: list[float],
    omega: list[float],
    kij: list[float],
    T: float,
    held: list[float],
) -> PhaseBoundaryResult: ...
def dew_pressure(
    Tc: list[float],
    Pc: list[float],
    omega: list[float],
    kij: list[float],
    T: float,
    held: list[float],
) -> PhaseBoundaryResult: ...
def pt_flash(
    Tc: list[float],
    Pc: list[float],
    omega: list[float],
    kij: list[float],
    T: float,
    P: float,
    z: list[float],
) -> PtFlashResult: ...
def stability_test(
    Tc: list[float],
    Pc: list[float],
    omega: list[float],
    kij: list[float],
    T: float,
    P: float,
    z: list[float],
) -> StabilityTestResult: ...
def pump_power(rho: float, q: float, H: float, eta: float) -> PumpPowerResult: ...
def orifice_flow(d: float, dP: float, rho: float, Cd: float) -> OrificeFlowResult: ...
def control_valve_cv(Cv: float, dP: float, SG: float) -> ControlValveCvResult: ...
def choked_flow_area(m_dot: float, P0: float, rho0: float, k: float) -> ChokedFlowAreaResult: ...
def crane_k_factors(fittings: list[str], f_t: float) -> KFactorsResult: ...
def darcy_weisbach(
    f: float, L: float, D: float, rho: float, v: float, mu: float | None = None
) -> DarcyWeisbachResult: ...

# --- introspection --------------------------------------------------------

def batch_run(calc_id: str, inputs: dict[str, list[float]]) -> BatchResult: ...
def data_files() -> list[DataFile]: ...
def fittings_rows() -> list[FittingRow]: ...
def fluid_rows(name: str) -> list[FluidRow]: ...
def warning_codes() -> list[str]: ...
def unit_names() -> list[str]: ...
def solver_kinds() -> list[str]: ...
def model_ids() -> list[str]: ...
def model_schemes(model_id: str) -> list[str]: ...
def model_kind(model_id: str) -> str: ...
def result_fields(calc_id: str) -> list[str]: ...
def calc_ids() -> list[str]: ...
def version() -> str: ...

# --- exceptions -----------------------------------------------------------
# Re-exported from azoth.core.errors, so both backends raise the same class
# objects rather than two lookalike hierarchies.

class AzothError(Exception): ...
class InvalidInputError(AzothError, ValueError): ...
class OutOfRangeError(AzothError, ValueError): ...
class PropertyUnavailableError(AzothError, LookupError): ...
class SolverNotConvergedError(AzothError, RuntimeError): ...
class UnitMismatchError(AzothError, TypeError): ...
class UnknownFittingError(AzothError, LookupError): ...

class UnverifiedCalculationError(AzothError):
    # Declared but never raised, here or in the Rust core. Reserved for an
    # explicit opt-in strictness gate; a calc with an unconfirmed source does not
    # raise on that account, because warnings are not errors here.
    ...
