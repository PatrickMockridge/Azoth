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
class PumpPowerResult:
    power: Qty
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
def pump_power(rho: float, q: float, H: float, eta: float) -> PumpPowerResult: ...
def orifice_flow(d: float, dP: float, rho: float, Cd: float) -> OrificeFlowResult: ...
def crane_k_factors(fittings: list[str], f_t: float) -> KFactorsResult: ...
def darcy_weisbach(
    f: float, L: float, D: float, rho: float, v: float, mu: float | None = None
) -> DarcyWeisbachResult: ...

# --- introspection --------------------------------------------------------

def warning_codes() -> list[str]: ...
def unit_names() -> list[str]: ...
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
    # explicit opt-in strictness gate; a calc with an unconfirmed source warns
    # (UNVERIFIED_SOURCE) rather than raising, because warnings are not errors.
    ...
