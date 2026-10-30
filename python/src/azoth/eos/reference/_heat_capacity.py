"""One phase's heat capacities at constant pressure and volume.

The Python twin of ``crates/azoth-eos/src/heat_capacity.rs``. ``cp - cv = -T (dP/dT)^2 /
(dP/dV)``, the cubic's own relation, with ``cp`` the molar departure capacity the
enthalpy/entropy reference already carries. ``eos.sound_speed`` needs ``gamma = cp/cv`` for
its speed; ``process.leak_blowdown`` needs the same capacities, weighted over the phases,
for the class's ``getGamma()``.
"""

from __future__ import annotations

import math
from typing import Any, NamedTuple

from azoth.core.errors import InvalidInputError
from azoth.core.units import from_si
from azoth.eos.reference._mixture_state import phase_state_at
from azoth.eos.reference.molar_enthalpy_entropy import molar_enthalpy_entropy
from azoth.eos.reference.pr_molar_volume import MOLAR_GAS_CONSTANT

#: Peng-Robinson's two geometry deltas, the mixture the model is stated on.
DELTA1 = 1.0 + math.sqrt(2.0)
DELTA2 = 1.0 - math.sqrt(2.0)


class PhaseCapacity(NamedTuple):
    """One phase's capacities and the pressure derivative they rest on."""

    cp: float
    cv: float
    d_p_d_v: float


def phase_capacity(
    mixture: Any,
    ideal_gas: Any,
    reduced: Any,
    kij: Any,
    t: float,
    p: float,
    composition: list[float],
    root: float,
) -> PhaseCapacity:
    """One phase's ``cp`` and ``cv`` from the cubic's closed-form derivatives.

    Raises:
        InvalidInputError: when the state gives no admissible derivative, or the capacities
            are not finite and positive - a sound speed and a blowdown both need a stable
            phase.
    """
    constant = MOLAR_GAS_CONSTANT
    state = phase_state_at(reduced, kij, composition, root)
    v = root * constant * t / p

    # `A` is dimensionless (`aP/(RT)^2`), so the dimensional `a` and `b` are these.
    a = state.a_mix * (constant * t) * (constant * t) / p
    b = state.b_mix * constant * t / p
    d_a_d_t = a * state.psi_bar / t
    f = (v + DELTA1 * b) * (v + DELTA2 * b)
    d_f = 2.0 * v + (DELTA1 + DELTA2) * b
    d_p_d_v = -constant * t / ((v - b) * (v - b)) + a * d_f / (f * f)
    d_p_d_t = constant / (v - b) - d_a_d_t / f
    if not d_p_d_v < 0.0 or not math.isfinite(d_p_d_v) or not math.isfinite(d_p_d_t):
        raise InvalidInputError(
            "state",
            f"the cubic's pressure-volume derivative is {d_p_d_v} Pa m3/mol at {t:.3} K and "
            f"{p:.3} Pa, and a sound speed needs a stable phase",
        )

    cp = float(
        molar_enthalpy_entropy(
            mixture, ideal_gas, from_si(t, "K"), from_si(p, "Pa"), composition, root
        )
        .cp.to_base_units()
        .magnitude
    )
    cv = cp + t * d_p_d_t * d_p_d_t / d_p_d_v
    if not cv > 0.0 or not math.isfinite(cp):
        raise InvalidInputError(
            "state", f"the phase's heat capacities are cp={cp}, cv={cv}, and both must be positive"
        )
    return PhaseCapacity(cp=cp, cv=cv, d_p_d_v=d_p_d_v)
