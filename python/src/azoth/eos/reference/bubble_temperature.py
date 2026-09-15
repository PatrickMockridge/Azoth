"""``eos.bubble_temperature`` - the temperature at which a liquid first gives off vapour.

Spec: ``specs/models/eos/bubble_temperature.toml``

The iteration lives in :mod:`azoth.eos.reference._phase_boundary_temperature`, because
``eos.dew_temperature`` runs the same loop with the phases exchanged. What is here is
the spec lookup, the range checks and the result.
"""

from __future__ import annotations

from azoth import _models_gen
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import BubbleTemperatureResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning
from azoth.eos.mixture import Mixture
from azoth.eos.reference._phase_boundary import VAPOUR
from azoth.eos.reference._phase_boundary_temperature import phase_boundary_temperature

MODEL_ID = "eos.bubble_temperature"


def bubble_temperature(mixture: Mixture, P: Q, x: list[float]) -> BubbleTemperatureResult:
    """The temperature at which a liquid of composition ``x`` first gives off vapour.

    ``x`` is the liquid's composition and is taken as given: this model does not ask
    whether that liquid is stable, only where its bubble point is.

    Raises:
        InvalidInputError: if the mixture has one component, or if ``x`` is the wrong
            length, has a negative entry, or does not sum to one.
        OutOfRangeError: if ``P`` is not positive, or if the mixture has no bubble
            point at this pressure. Reported on ``min_t_over_tc``.
        SolverNotConvergedError: if the iteration hits its cap.
    """
    spec = _models_gen.model(MODEL_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    pressure = input_to_si(spec, "P", P)
    apply_checks(checks.on_input, {"P": pressure}.get, warnings)

    boundary = phase_boundary_temperature(mixture, pressure, list(x), VAPOUR, spec["algorithm"])
    warnings.extend(boundary["warnings"])

    min_t_over_tc = min(
        boundary["temperature"] / c.Tc.to_base_units().magnitude for c in mixture.components
    )
    apply_checks(
        checks.derived,
        lambda name: min_t_over_tc if name == "min_t_over_tc" else None,
        warnings,
    )

    return BubbleTemperatureResult(
        temperature=from_si(boundary["temperature"], "K"),
        incipient=tuple(boundary["incipient"]),
        k=tuple(boundary["k"]),
        z_liquid=boundary["z_held"],
        z_vapour=boundary["z_incipient"],
        min_t_over_tc=min_t_over_tc,
        iterations=boundary["iterations"],
        residual=boundary["residual"],
        warnings=tuple(warnings),
    )
