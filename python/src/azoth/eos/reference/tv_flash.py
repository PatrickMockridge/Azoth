"""``eos.tv_flash`` - the pressure a mixture reaches at a given temperature and volume.

Spec: ``specs/models/eos/tv_flash.toml``

An outer quasi-Newton on pressure, over the molar volume assembled from the phase split
at a trial pressure (:mod:`azoth.eos.reference.pt_flash`). The volume falls monotonically
with pressure at a fixed temperature. The loop lives in
:mod:`azoth.eos.reference._flash_property`, shared with the other flash-breadth models.
"""

from __future__ import annotations

from azoth import _models_gen
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import TvFlashResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning
from azoth.eos.mixture import Mixture
from azoth.eos.reference._flash_property import distinct_warnings, solve_pressure
from azoth.eos.reference.molar_enthalpy_entropy import IdealGasModel
from azoth.eos.reference.pr_molar_volume import MOLAR_GAS_CONSTANT

MODEL_ID = "eos.tv_flash"


def tv_flash(
    mixture: Mixture,
    ideal_gas: IdealGasModel,
    T: Q,
    V: Q,
    z: list[float],
) -> TvFlashResult:
    """The pressure at which a mixture has a given molar volume at a temperature.

    Raises:
        OutOfRangeError: if ``T`` or ``V`` is not positive.
        SolverNotConvergedError: if the iteration reaches its cap.
    """
    spec = _models_gen.model(MODEL_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    t_si = input_to_si(spec, "T", T)
    v_si = input_to_si(spec, "V", V)
    apply_checks(checks.on_input, {"T": t_si, "V": v_si}.get, warnings)

    start = MOLAR_GAS_CONSTANT * t_si / v_si
    solved = solve_pressure(mixture, ideal_gas, t_si, v_si, list(z), "v", spec["algorithm"], start)
    warnings.extend(solved["warnings"])
    flash = solved["state"]["flash"]

    apply_checks(
        checks.derived,
        lambda name: solved["P"] if name == "P" else None,
        warnings,
    )

    return TvFlashResult(
        P=from_si(solved["P"], "Pa"),
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
