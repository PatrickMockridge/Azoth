"""``eos.pu_flash`` - the temperature a mixture reaches at a given property.

Spec: ``specs/models/eos/pu_flash.toml``

A thin procedure model over the flash-property solver, mirroring ``eos.tv_flash`` /
``eos.pv_flash`` with the property exchanged.
"""

from __future__ import annotations

from azoth import _models_gen
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import PuFlashResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning
from azoth.eos.mixture import Mixture
from azoth.eos.reference._flash_property import distinct_warnings, solve_temperature
from azoth.eos.reference.molar_enthalpy_entropy import IdealGasModel

MODEL_ID = "eos.pu_flash"


def pu_flash(
    mixture: Mixture,
    ideal_gas: IdealGasModel,
    P: Q,
    U: Q,
    z: list[float],
) -> PuFlashResult:
    """The temperature at which a mixture has a given property.

    Raises:
        OutOfRangeError: if a state input is not positive.
        SolverNotConvergedError: if the iteration reaches its cap.
    """
    spec = _models_gen.model(MODEL_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    p_si = input_to_si(spec, "P", P)
    u_si = input_to_si(spec, "U", U)
    apply_checks(checks.on_input, {"P": p_si, "U": u_si}.get, warnings)

    solved = solve_temperature(mixture, ideal_gas, p_si, u_si, list(z), "u", spec["algorithm"])
    warnings.extend(solved["warnings"])
    flash = solved["state"]["flash"]

    apply_checks(
        checks.derived,
        lambda name: solved["T"] if name == "T" else None,
        warnings,
    )

    return PuFlashResult(
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
