"""``eos.ps_flash`` - the temperature a mixture reaches at a given pressure and entropy.

Spec: ``specs/models/eos/ps_flash.toml``

An outer bisection on temperature, for an **isentropic** unit operation - a compressor,
an expander, a turbine or a nozzle, assumed ideal - whose outlet temperature is not
known because entropy is conserved and temperature is not. The bracket, the loop, the
phase branch and the warning handling live in
:mod:`azoth.eos.reference._flash_property`, shared with ``eos.ph_flash``.

Each phase's entropy is taken at *that phase's own composition*, so the weighted sum
carries the entropy of mixing; summing the phases at the feed composition instead would
conserve entropy across a phase change, which is wrong.
"""

from __future__ import annotations

from typing import Any

from azoth import _models_gen
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import PsFlashResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning
from azoth.eos.mixture import Mixture
from azoth.eos.reference._flash_property import Property, distinct_warnings, property_at
from azoth.eos.reference._flash_property import solve_temperature as _solve
from azoth.eos.reference.molar_enthalpy_entropy import IdealGasModel

MODEL_ID = "eos.ps_flash"

#: Which property this model inverts. The one thing it does not share with `ph_flash`.
WHICH: Property = "s"


def entropy_at(
    mixture: Mixture,
    ideal_gas: IdealGasModel,
    t_si: float,
    p_si: float,
    z: list[float],
) -> tuple[float, dict[str, Any]]:
    """The molar entropy of a mixture at a temperature and pressure, and its split.

    The composition the flash settles on is the equilibrium one, so this is the entropy
    of the *feed* at that state - which is what makes it comparable with the entropy a
    caller says the stream arrived with.
    """
    return property_at(mixture, ideal_gas, t_si, p_si, z, WHICH)


def solve_temperature(
    mixture: Mixture,
    ideal_gas: IdealGasModel,
    p_si: float,
    s_si: float,
    z: list[float],
    algorithm: dict[str, Any],
) -> dict[str, Any]:
    """Invert ``S(T, P) = s_si`` for ``T`` by bracketing, then bisecting."""
    return _solve(mixture, ideal_gas, p_si, s_si, z, WHICH, algorithm)


def ps_flash(
    mixture: Mixture,
    ideal_gas: IdealGasModel,
    P: Q,
    S: Q,
    z: list[float],
) -> PsFlashResult:
    """The temperature at which a mixture has a given molar entropy at a pressure.

    ``S`` is a *difference* from the datum ``ideal_gas`` carries, not an absolute
    quantity - the same caveat ``eos.ph_flash`` and ``eos.molar_enthalpy_entropy``
    carry, and for the same reason.

    Raises:
        InvalidInputError: if the composition or any ideal-gas vector is the wrong
            length, or if ``z`` is not a composition.
        OutOfRangeError: if ``P`` is not positive, or a range check on the answer fails.
        SolverNotConvergedError: if no temperature on the bracket covers the requested
            entropy, or if the bisection reaches its cap.
    """
    spec = _models_gen.model(MODEL_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    p_si = input_to_si(spec, "P", P)
    s_si = input_to_si(spec, "S", S)
    apply_checks(checks.on_input, {"P": p_si, "S": s_si}.get, warnings)

    solved = solve_temperature(mixture, ideal_gas, p_si, s_si, list(z), spec["algorithm"])
    warnings.extend(solved["warnings"])
    flash = solved["state"]["flash"]

    apply_checks(
        checks.derived,
        lambda name: solved["T"] if name == "T" else None,
        warnings,
    )

    return PsFlashResult(
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
