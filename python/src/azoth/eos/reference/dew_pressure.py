"""``eos.dew_pressure`` - the pressure at which a vapour first condenses.

Spec: ``specs/models/eos/dew_pressure.toml``

The companion of ``eos.bubble_pressure`` and the same iteration with the phases
exchanged. The loop lives in :mod:`azoth.eos.reference._phase_boundary` so that the
guard against the trivial solution is written once; what is here is the spec lookup,
the range checks and the result.
"""

from __future__ import annotations

from azoth import _models_gen
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import DewPressureResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning
from azoth.eos.mixture import Mixture
from azoth.eos.reference._phase_boundary import LIQUID, phase_boundary_pressure

MODEL_ID = "eos.dew_pressure"


def dew_pressure(mixture: Mixture, T: Q, y: list[float]) -> DewPressureResult:
    """The pressure at which a vapour of composition ``y`` first condenses.

    ``y`` is the vapour's composition and is taken as given: this model does not ask
    whether that vapour is stable, only where its dew point is.

    Raises:
        InvalidInputError: if the mixture has one component, or if ``y`` is the wrong
            length, has a negative entry, or does not sum to one.
        OutOfRangeError: if ``T`` is not positive, or if the mixture has no dew point
            at this temperature - the mixture is at or above its critical condition.
            Reported on ``min_t_over_tc``.
        SolverNotConvergedError: if the iteration hits its cap.
    """
    spec = _models_gen.model(MODEL_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    temperature = input_to_si(spec, "T", T)
    apply_checks(checks.on_input, {"T": temperature}.get, warnings)

    min_t_over_tc = min(temperature / c.Tc.to_base_units().magnitude for c in mixture.components)
    apply_checks(
        checks.derived,
        lambda name: min_t_over_tc if name == "min_t_over_tc" else None,
        warnings,
    )

    boundary = phase_boundary_pressure(mixture, temperature, list(y), LIQUID, spec["algorithm"])
    warnings.extend(boundary["warnings"])

    return DewPressureResult(
        pressure=from_si(boundary["pressure"], "Pa"),
        incipient=tuple(boundary["incipient"]),
        k=tuple(boundary["k"]),
        z_liquid=boundary["z_incipient"],
        z_vapour=boundary["z_held"],
        min_t_over_tc=min_t_over_tc,
        iterations=boundary["iterations"],
        residual=boundary["residual"],
        warnings=tuple(warnings),
    )
