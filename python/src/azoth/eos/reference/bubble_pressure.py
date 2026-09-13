"""``eos.bubble_pressure`` - the pressure at which a liquid first gives off vapour.

Spec: ``specs/models/eos/bubble_pressure.yaml``

The iteration lives in :mod:`azoth.eos.reference._phase_boundary`, because
``eos.dew_pressure`` runs the same loop with the phases exchanged. What is here is
the spec lookup, the range checks and the result.
"""

from __future__ import annotations

from azoth import _models_gen
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import BubblePressureResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning
from azoth.eos.mixture import Mixture
from azoth.eos.reference._phase_boundary import VAPOUR, phase_boundary_pressure

MODEL_ID = "eos.bubble_pressure"


def bubble_pressure(mixture: Mixture, T: Q, x: list[float]) -> BubblePressureResult:
    """The pressure at which a liquid of composition ``x`` first gives off vapour.

    ``x`` is the liquid's composition and is taken as given: this model does not ask
    whether that liquid is stable, only where its bubble point is.

    Raises:
        InvalidInputError: if the mixture has one component, or if ``x`` is the wrong
            length, has a negative entry, or does not sum to one.
        OutOfRangeError: if ``T`` is not positive, or if the mixture has no bubble
            point at this temperature - the mixture is at or above its critical
            condition. Reported on ``min_t_over_tc``.
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

    boundary = phase_boundary_pressure(mixture, temperature, list(x), VAPOUR, spec["algorithm"])
    warnings.extend(boundary["warnings"])

    return BubblePressureResult(
        pressure=from_si(boundary["pressure"], "Pa"),
        incipient=tuple(boundary["incipient"]),
        k=tuple(boundary["k"]),
        z_liquid=boundary["z_held"],
        z_vapour=boundary["z_incipient"],
        min_t_over_tc=min_t_over_tc,
        iterations=boundary["iterations"],
        residual=boundary["residual"],
        warnings=tuple(warnings),
    )
