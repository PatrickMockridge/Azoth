"""``eos.sound_speed`` - the cubic's isentropic speed of sound.

Spec: ``specs/models/eos/sound_speed.toml``.

.. code-block:: text

    c^2 = -gamma (v^2 / M) (dP/dv)_T        gamma = cp/cv
    cp - cv = -T (dP/dT)^2_v / (dP/dv)_T

NeqSim's ``PhaseEos.getSoundSpeed`` in SI: its ``bs = -(1/V)(Cv/Cp)/(dPdVTn)`` and
``sqrt(V/Mw/bs)`` reduce to the expression above once the class's 1e5 volume unit and the
1e-10 inside ``dPdVTn`` cancel. The two derivatives are closed forms of the same residual the
departure functions evaluate:

.. code-block:: text

    P = R T/(v - b) - a/((v + d1 b)(v + d2 b))
    (dP/dv)_T = -R T/(v - b)^2 + a (2v + (d1 + d2) b) / f^2
    (dP/dT)_v = R/(v - b) - (da/dT)/f
    da/dT = a psi_bar/T

``psi_bar`` is ``dln a/dln T``, which the phase state already carries for the departure heat
capacity. The system value is the class's beta-weighted mean, not a mixture speed.
"""

from __future__ import annotations

import math
from typing import Any

from azoth.core.errors import InvalidInputError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import SoundSpeedResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning
from azoth.eos import components as _components
from azoth.eos.reference._heat_capacity import phase_capacity
from azoth.eos.reference._mixture_state import phase_state, reduced_parameters
from azoth.eos.reference.pr_molar_volume import MOLAR_GAS_CONSTANT
from azoth.eos.reference.pt_flash import pt_flash

MODEL_ID = "eos.sound_speed"


def sound_speed(
    components: list[str],
    T: Q,
    P: Q,
    z: list[float],
) -> SoundSpeedResult:
    """The speed of sound in a mixture at a state.

    Args:
        components: the substances, by name, resolved against the databank.
        T: the state's absolute temperature.
        P: the state's absolute pressure.
        z: the mixture's mole fractions, checked rather than renormalised.

    Returns:
        The phase's speed, the class's beta-weighted system value, and what the flash found.

    Raises:
        InvalidInputError: if the composition's length disagrees, or the state gives no
            admissible derivative.
        OutOfRangeError: for a non-positive temperature or pressure.

    See :func:`azoth.eos.sound_speed`.
    """
    from azoth import _models_gen

    spec = _models_gen.model(MODEL_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    t_si = input_to_si(spec, "T", T)
    p_si = input_to_si(spec, "P", P)
    apply_checks(checks.on_input, {"T": t_si, "P": p_si}.get, warnings)

    mixture, ideal_gas = _components.mixture_of(components, eos="pr")
    if len(z) != len(components):
        raise InvalidInputError(
            "z",
            f"a mixture of {len(components)} components needs {len(components)} mole "
            f"fractions, and {z!r} has {len(z)}",
        )

    reduced = reduced_parameters(mixture, t_si, p_si)
    kij = mixture.effective_kij(reduced.t_kelvin)
    flash = pt_flash(mixture, from_si(t_si, "K"), from_si(p_si, "Pa"), list(z))

    if flash.phase == "two_phase":
        beta = flash.vapour_fraction
        if beta is None:
            raise InvalidInputError(
                "state",
                "the flash split the feed and reported no vapour fraction, so the system "
                "speed has no weight for either phase",
            )
        vapour = _phase_speed(
            mixture, ideal_gas, reduced, kij, t_si, p_si, list(flash.y), flash.z_vapour
        )
        liquid = _phase_speed(
            mixture, ideal_gas, reduced, kij, t_si, p_si, list(flash.x), flash.z_liquid
        )
        phase_speed = vapour
        system_speed = beta * vapour + (1.0 - beta) * liquid
    else:
        # **The feed's own root, at the feed's composition.** `pt_flash` on a superheated
        # stream reports the *incipient* phase's root; the class reads the one phase it holds.
        liquid = flash.phase != "all_vapour"
        state = phase_state(reduced, kij, list(z), liquid=liquid)
        phase_speed = _phase_speed(mixture, ideal_gas, reduced, kij, t_si, p_si, list(z), state.z)
        system_speed = phase_speed

    return SoundSpeedResult(
        sound_speed=from_si(phase_speed, "m/s"),
        system_sound_speed=from_si(system_speed, "m/s"),
        phase=flash.phase,
        warnings=tuple(warnings),
    )


def _phase_speed(
    mixture: Any,
    ideal_gas: Any,
    reduced: Any,
    kij: Any,
    t: float,
    p: float,
    composition: list[float],
    root: float,
) -> float:
    """One phase's isentropic speed, m/s."""
    capacity = phase_capacity(mixture, ideal_gas, reduced, kij, t, p, composition, root)
    molar_mass = sum(
        component.molar_mass.to("kg/mol").magnitude * fraction
        for component, fraction in zip(mixture.components, composition, strict=True)
    )
    v = root * MOLAR_GAS_CONSTANT * t / p
    gamma = capacity.cp / capacity.cv
    return v * math.sqrt(-gamma * capacity.d_p_d_v / molar_mass)
