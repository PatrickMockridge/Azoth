"""The enthalpy of an outlet that is *a phase*.

The Python twin of ``azoth_process::Stream::from_side``, and it exists for the same
reason: **re-flashing an outlet's composition is a different question from the state the
vessel produced.**

``Stream::from_pt`` - and :func:`azoth.eos.reference.ph_flash.enthalpy_at`, which is what
this namespace used for its outlets until the tank found otherwise - answers with the state
a composition settles on *on its own*. A phase of a split is on the root the split put it
on, and where the two differ the difference is not small. Measured on a water-bearing gas
leaving a tank, the outlet's composition re-flashes to two phases at a vapour fraction of
``0.767`` and an enthalpy of ``-10043.91`` J/mol, where the phase's own root gives
``441.77`` - and NeqSim's ``setThermoSystemFromPhase`` reports ``441.83``.

Where the composition *is* stable on its own the two agree, which is why this was
invisible until a feed whose phase separates when its parent does not.

**Which one an outlet is, is the class's decision and not this module's.** A vessel that
re-runs an outlet as a stream - ``Separator.run`` does, for the liquid, under
``gasInLiquid != 0.0`` and for the vapour under ``oilInGas``/``waterInGas`` - re-flashes
it, and that branch is the caller's to take. This is the other branch.
"""

from __future__ import annotations

from typing import Any

from azoth.core.units import from_si
from azoth.eos.reference._mixture_state import phase_state, reduced_parameters
from azoth.eos.reference.molar_enthalpy_entropy import molar_enthalpy_entropy

__all__ = ["phase_enthalpy"]


def phase_enthalpy(
    mixture: Any,
    ideal_gas: Any,
    t_si: float,
    p_si: float,
    composition: list[float],
    *,
    liquid: bool,
) -> float:
    """One phase's molar enthalpy at a state, in J/mol, on its own side of the cubic.

    Args:
        mixture: the components and their interaction parameters.
        ideal_gas: the heat-capacity coefficients and the datum.
        t_si: the state's temperature in kelvin.
        p_si: the state's pressure in pascals.
        composition: the *phase's* mole fractions.
        liquid: which side of the cubic to root on - the smallest for the liquid, the
            largest for the vapour, `eos.pr_z_factor`'s rule.

    Returns:
        The molar enthalpy in J/mol, on the datum ``ideal_gas`` carries.
    """
    reduced = reduced_parameters(mixture, t_si, p_si)
    root = phase_state(reduced, mixture.kij, list(composition), liquid=liquid).z
    state = molar_enthalpy_entropy(
        mixture,
        ideal_gas,
        from_si(t_si, "K"),
        from_si(p_si, "Pa"),
        list(composition),
        root,
    )
    return float(state.h.to_base_units().magnitude)
