"""``eos.vu_flash`` - the pressure and temperature a mixture reaches at a given volume and
internal energy.

Spec: ``specs/models/eos/vu_flash.toml``

A closed vessel at fixed volume and internal energy: both state variables are solved for,
by the decoupled 2x2 Newton in :mod:`azoth.eos.reference._flash_property`.
"""

from __future__ import annotations

from azoth import _models_gen
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import VuFlashResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning
from azoth.eos.mixture import Mixture
from azoth.eos.reference._flash_property import distinct_warnings, solve_pressure_temperature
from azoth.eos.reference.molar_enthalpy_entropy import IdealGasModel
from azoth.eos.reference.pr_molar_volume import MOLAR_GAS_CONSTANT

MODEL_ID = "eos.vu_flash"


def vu_flash(
    mixture: Mixture,
    ideal_gas: IdealGasModel,
    V: Q,
    U: Q,
    z: list[float],
) -> VuFlashResult:
    """The pressure and temperature at which a mixture has a given volume and internal energy.

    Raises:
        OutOfRangeError: if ``V`` is not positive.
        SolverNotConvergedError: if the iteration reaches its cap.
    """
    spec = _models_gen.model(MODEL_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    v_si = input_to_si(spec, "V", V)
    u_si = input_to_si(spec, "U", U)
    apply_checks(checks.on_input, {"V": v_si, "U": u_si}.get, warnings)

    start_t = float(spec["algorithm"].get("initial_temperature", 300.0))
    start_p = MOLAR_GAS_CONSTANT * start_t / v_si
    solved = solve_pressure_temperature(
        mixture, ideal_gas, v_si, u_si, list(z), spec["algorithm"], start_p, start_t
    )
    warnings.extend(solved["warnings"])
    flash = solved["state"]["flash"]

    return VuFlashResult(
        P=from_si(solved["P"], "Pa"),
        T=from_si(solved["T"], "K"),
        beta=flash.beta,
        x=tuple(flash.x),
        y=tuple(flash.y),
        k=tuple(flash.k),
        phase=flash.phase,
        z_liquid=flash.z_liquid,
        z_vapour=flash.z_vapour,
        iterations=solved["iterations"],
        residual=solved["residual"],
        warnings=distinct_warnings(warnings),
    )
