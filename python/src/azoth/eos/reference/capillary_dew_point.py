"""``eos.capillary_dew_point`` - the temperature at which a vapour condenses inside a pore.

Spec: ``specs/models/eos/capillary_dew_point.toml``. NeqSim's ``CapillaryDewPointFlash``,
reached through ``ThermodynamicOperations.capillaryDewPointTemperatureFlash(r[, theta])``.

The same boundary :mod:`azoth.eos.reference.dew_temperature` finds, drawn on a curved interface
instead of a flat one. Young-Laplace puts the liquid inside a pore at a pressure above the
vapour outside it by ``2 sigma cos(theta) / r``, and the Kelvin equation turns that into a
shift on every K-value:

    K_cap,i = K_i exp(-Vm_L dP_cap / (R T))

so the incipient liquid is *stabilised* - the dew point moves **up**, not down, which is the
direction that matters for condensation in tight rock.

**The surface tension is the caller's, not the mixture's.** NeqSim reads it off the phase's
interphase properties and falls back to a constant ``0.005 N/m`` when that returns nothing;
this takes it as an argument, because a surface tension is a fitted quantity with its own
provenance and reading one inside a model whose other inputs are all stated would hide which
was used. That is also why the numbers differ from NeqSim's by about 5%: the shift is linear in
the molar volume, and NeqSim takes that volume from its phase's ``getMolarVolume`` where this
takes it from the cubic's own lower root. ``validation/neqsim/CapillaryDew.java`` measures both.
"""

from __future__ import annotations

import math

from azoth import _models_gen
from azoth.core.errors import InvalidInputError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import CapillaryDewPointResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning
from azoth.eos.mixture import Mixture
from azoth.eos.reference._phase_boundary import LIQUID
from azoth.eos.reference._phase_boundary_temperature import phase_boundary_temperature

MODEL_ID = "eos.capillary_dew_point"


def capillary_dew_point(
    mixture: Mixture,
    P: Q,
    y: list[float],
    pore_radius: float,
    contact_angle: float,
    surface_tension: float,
) -> CapillaryDewPointResult:
    """The dew-point temperature of a vapour held in a pore of a stated radius.

    ``y`` is the vapour's composition and is taken as given: this model does not ask whether
    that vapour is stable, only where its dew point is. ``pore_radius`` is in metres,
    ``contact_angle`` in radians - zero is a perfectly wetting liquid - and ``surface_tension``
    in newtons per metre.

    Raises:
        InvalidInputError: if the contact angle is not finite, or ``y`` is the wrong length,
            has a negative entry, or does not sum to one.
        OutOfRangeError: if ``P`` is not positive, ``pore_radius`` is not positive,
            ``surface_tension`` is negative, or the mixture has no dew point at this pressure.
        SolverNotConvergedError: if the iteration hits its cap.
    """
    spec = _models_gen.model(MODEL_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    pressure = input_to_si(spec, "P", P)
    radius = input_to_si(spec, "pore_radius", pore_radius)
    sigma = input_to_si(spec, "surface_tension", surface_tension)
    # **The angle is a bare number, and it is the only scalar input here that is.** The spec
    # declares it `dimensionless` because the vocabulary carries no angle; every *other*
    # dimensionless scalar in this library is a bare number too - a composition, a mole
    # fraction - because a dimensionless quantity is the number. `input_to_si` would demand a
    # `Q(x, "")` for it, which is a wrapper around nothing. A caller who writes the angle as a
    # quantity is taken at its magnitude rather than refused, because `Q(0.5, "rad")` is what a
    # reader reaches for and it is not wrong.
    angle = (
        contact_angle.to_base_units().magnitude
        if hasattr(contact_angle, "to_base_units")
        else float(contact_angle)
    )
    apply_checks(
        checks.on_input,
        {"P": pressure, "pore_radius": radius, "surface_tension": sigma}.get,
        warnings,
    )
    # The radius and the tension are range-checked above, from the spec's `valid_range` blocks.
    # The angle has no range - every real angle is admissible, and the two outside the first
    # quadrant are the same interface by symmetry - so only a non-finite one is refused.
    if not math.isfinite(angle):
        raise InvalidInputError("contact_angle", "the contact angle is not a finite number")

    capillary_pressure = 2.0 * sigma * math.cos(angle) / radius
    boundary = phase_boundary_temperature(
        mixture, pressure, list(y), LIQUID, spec["algorithm"], capillary_pressure
    )
    warnings.extend(boundary["warnings"])

    min_t_over_tc = min(
        boundary["temperature"] / c.Tc.to_base_units().magnitude for c in mixture.components
    )
    apply_checks(
        checks.derived,
        lambda name: min_t_over_tc if name == "min_t_over_tc" else None,
        warnings,
    )

    return CapillaryDewPointResult(
        temperature=from_si(boundary["temperature"], "K"),
        incipient=tuple(boundary["incipient"]),
        k=tuple(boundary["k"]),
        z_liquid=boundary["z_incipient"],
        z_vapour=boundary["z_held"],
        capillary_pressure=from_si(capillary_pressure, "Pa"),
        min_t_over_tc=min_t_over_tc,
        iterations=boundary["iterations"],
        residual=boundary["residual"],
        warnings=tuple(warnings),
    )
