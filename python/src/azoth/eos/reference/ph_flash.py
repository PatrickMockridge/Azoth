"""``eos.ph_flash`` - the temperature a mixture reaches at a given pressure and enthalpy.

Spec: ``specs/models/eos/ph_flash.toml``

An outer bisection on temperature, over the enthalpy assembled from two things this
library already has: the phase split at a trial temperature, from
:mod:`azoth.eos.reference.pt_flash`, and each phase's enthalpy, from
``eos.molar_enthalpy_entropy``. The enthalpy rises monotonically with temperature at a
fixed pressure, which is what makes the bisection well posed.

The bracket, the loop and the phase branch live in
:mod:`azoth.eos.reference._flash_property`, shared with ``eos.ps_flash``.
"""

from __future__ import annotations

from typing import Any

from azoth import _models_gen
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import PhFlashResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning
from azoth.eos.mixture import Mixture
from azoth.eos.reference._flash_property import Property, distinct_warnings, property_at
from azoth.eos.reference._flash_property import solve_temperature as _solve
from azoth.eos.reference.molar_enthalpy_entropy import IdealGasModel

MODEL_ID = "eos.ph_flash"

#: Which property this model inverts. The one thing it does not share with `ps_flash`.
WHICH: Property = "h"


def enthalpy_at(
    mixture: Mixture,
    ideal_gas: IdealGasModel,
    t_si: float,
    p_si: float,
    z: list[float],
) -> tuple[float, dict[str, Any]]:
    """The molar enthalpy of a mixture at a temperature and pressure, and its split.

    The composition the flash settles on is the equilibrium one, so this is the enthalpy
    of the *feed* at that state - which is what makes it comparable with a duty a caller
    supplied.
    """
    return property_at(mixture, ideal_gas, t_si, p_si, z, WHICH)


def solve_temperature(
    mixture: Mixture,
    ideal_gas: IdealGasModel,
    p_si: float,
    h_si: float,
    z: list[float],
    algorithm: dict[str, Any],
) -> dict[str, Any]:
    """Invert ``H(T, P) = h_si`` for ``T`` by bracketing, then bisecting."""
    return _solve(mixture, ideal_gas, p_si, h_si, z, WHICH, algorithm)


def ph_flash(
    mixture: Mixture,
    ideal_gas: IdealGasModel,
    P: Q,
    H: Q,
    z: list[float],
) -> PhFlashResult:
    """The temperature at which a mixture has a given molar enthalpy at a pressure.

    ``H`` is a *difference* from the datum ``ideal_gas`` carries, not an absolute
    quantity: two calls with different reference values are not comparable, and their
    difference is a plausible number rather than an error.

    Raises:
        InvalidInputError: if the composition or any ideal-gas vector is the wrong
            length, or if ``z`` is not a composition.
        OutOfRangeError: if ``P`` is not positive, or a range check on the answer fails.
        SolverNotConvergedError: if no temperature on the bracket covers the requested
            enthalpy, or if the bisection reaches its cap.
    """
    spec = _models_gen.model(MODEL_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    p_si = input_to_si(spec, "P", P)
    h_si = input_to_si(spec, "H", H)
    apply_checks(checks.on_input, {"P": p_si, "H": h_si}.get, warnings)

    solved = solve_temperature(mixture, ideal_gas, p_si, h_si, list(z), spec["algorithm"])
    warnings.extend(solved["warnings"])
    flash = solved["state"]["flash"]

    apply_checks(
        checks.derived,
        lambda name: solved["T"] if name == "T" else None,
        warnings,
    )

    return PhFlashResult(
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
