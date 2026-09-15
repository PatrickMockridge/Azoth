"""``eos.pv_flash`` - the temperature a mixture reaches at a given pressure and volume.

Spec: ``specs/models/eos/pv_flash.toml``

An outer quasi-Newton on temperature, over the molar volume assembled from the phase
split at a trial temperature (:mod:`azoth.eos.reference.pt_flash`). The volume rises
monotonically with temperature at a fixed pressure. The loop lives in
:mod:`azoth.eos.reference._flash_property`.
"""

from __future__ import annotations

from azoth import _models_gen
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import PvFlashResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning
from azoth.eos.mixture import Mixture
from azoth.eos.reference._flash_property import distinct_warnings, solve_temperature
from azoth.eos.reference.molar_enthalpy_entropy import IdealGasModel

MODEL_ID = "eos.pv_flash"


def pv_flash(
    mixture: Mixture,
    ideal_gas: IdealGasModel,
    P: Q,
    V: Q,
    z: list[float],
) -> PvFlashResult:
    """The temperature at which a mixture has a given molar volume at a pressure.

    Raises:
        OutOfRangeError: if ``P`` or ``V`` is not positive.
        SolverNotConvergedError: if the iteration reaches its cap.
    """
    spec = _models_gen.model(MODEL_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    p_si = input_to_si(spec, "P", P)
    v_si = input_to_si(spec, "V", V)
    apply_checks(checks.on_input, {"P": p_si, "V": v_si}.get, warnings)

    solved = solve_temperature(mixture, ideal_gas, p_si, v_si, list(z), "v", spec["algorithm"])
    warnings.extend(solved["warnings"])
    flash = solved["state"]["flash"]

    apply_checks(
        checks.derived,
        lambda name: solved["T"] if name == "T" else None,
        warnings,
    )

    return PvFlashResult(
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
