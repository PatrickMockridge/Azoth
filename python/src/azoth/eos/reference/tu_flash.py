"""``eos.tu_flash`` - the pressure a mixture reaches at a given property.

Spec: ``specs/models/eos/tu_flash.toml``

A thin procedure model over the flash-property solver, mirroring ``eos.tv_flash`` /
``eos.pv_flash`` with the property exchanged.
"""

from __future__ import annotations

from azoth import _models_gen
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import TuFlashResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning
from azoth.eos.mixture import Mixture
from azoth.eos.reference._flash_property import distinct_warnings, solve_pressure
from azoth.eos.reference.molar_enthalpy_entropy import IdealGasModel

MODEL_ID = "eos.tu_flash"


def tu_flash(
    mixture: Mixture,
    ideal_gas: IdealGasModel,
    T: Q,
    U: Q,
    z: list[float],
) -> TuFlashResult:
    """The pressure at which a mixture has a given property.

    Raises:
        OutOfRangeError: if a state input is not positive.
        SolverNotConvergedError: if the iteration reaches its cap.
    """
    spec = _models_gen.model(MODEL_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    t_si = input_to_si(spec, "T", T)
    u_si = input_to_si(spec, "U", U)
    apply_checks(checks.on_input, {"T": t_si, "U": u_si}.get, warnings)

    solved = solve_pressure(mixture, ideal_gas, t_si, u_si, list(z), "u", spec["algorithm"], 1.0e5)
    warnings.extend(solved["warnings"])
    flash = solved["state"]["flash"]

    apply_checks(
        checks.derived,
        lambda name: solved["P"] if name == "P" else None,
        warnings,
    )

    return TuFlashResult(
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
