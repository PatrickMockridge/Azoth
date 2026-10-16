"""``eos.calc_ionic_composition`` - molalities and activity coefficients of one phase.

Spec: ``specs/models/eos/calc_ionic_composition.toml``. A *direct* model: the class solves
nothing, it reads a phase. What makes it more than arithmetic is the third column.

This is the pure-Python reference: a second, independent expression of the same physics as
the Rust kernel, held to it by ``test_cross_impl.py``.

The activity coefficient is ``exp(ln phi_i - ln phi_i^reference)``, and the reference is a
**second phase of the same fluid at a different composition** - pure for a ``solvent``,
``1e-10`` mol beside ``10`` mol of water for anything else. So an ion's coefficient is
against the ion at infinite dilution in water and not against water in the same phase.
"""

from __future__ import annotations

import math
from collections.abc import Sequence

from azoth import _models_gen
from azoth.core.errors import InvalidInputError, PropertyUnavailableError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import CalcIonicCompositionResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning
from azoth.eos.components import SOLVENT, electrolyte_cpa_statoil_mixture_of, entry
from azoth.eos.reference._mixture_state import phase_state, reduced_parameters

MODEL_ID = "eos.calc_ionic_composition"

#: The solute's amount in the reference binary, the class's own.
_DILUTE_MOLES = 1.0e-10

#: The solvent's amount beside it, the class's own.
_SOLVENT_MOLES = 10.0


def calc_ionic_composition(
    components: Sequence[str],
    T: Q,
    P: Q,
    x: Sequence[float],
) -> CalcIonicCompositionResult:
    """Molalities and activity coefficients for each component of one phase.

    Args:
        components: the substances, by name, in the order every output vector is indexed by.
            ``water`` must be one of them.
        T: absolute temperature.
        P: absolute pressure.
        x: the phase's mole fractions, checked rather than renormalised.

    Returns:
        Each component's molality in mol/kg and mass ratio in kg/kg of water, and its
        activity coefficient at its own reference state.

    Raises:
        InvalidInputError: if ``x`` is not a composition, a name is not in the databank, or no
            component is named ``water``.
        PropertyUnavailableError: if the databank states no molar mass for a component.
        OutOfRangeError: if ``T`` or ``P`` is not positive, or no volume root exists.

    See :func:`azoth.eos.calc_ionic_composition`.
    """
    spec = _models_gen.model(MODEL_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    t_si = input_to_si(spec, "T", T)
    p_si = input_to_si(spec, "P", P)
    apply_checks(checks.on_input, {"T": t_si, "P": p_si}.get, warnings)

    n = len(components)
    if len(x) != n:
        raise InvalidInputError(
            "x", f"{n} components but {len(x)} mole fractions; the two must match"
        )
    total = sum(x)
    if abs(total - 1.0) > 1.0e-9:
        raise InvalidInputError(
            "x",
            f"the mole fractions sum to {total}, not to one. Renormalising them here "
            "would make a composition error invisible in every number downstream, so "
            "it is refused instead",
        )

    # **The solvent is found by name**, as the class finds it: it asks the phase for its
    # `water` component rather than for one whose reference state is `solvent`.
    water = next((i for i, name in enumerate(components) if name.strip().lower() == "water"), None)
    if water is None:
        raise InvalidInputError(
            "components",
            "the two concentration columns are per kilogram of water, and no component is "
            "named `water`",
        )

    entries = [entry(name) for name in components]
    if entries[water].molar_mass is None:
        raise PropertyUnavailableError(
            components[water], "a molar mass", "the concentration columns are per kilogram of it"
        )
    water_mass = entries[water].molar_mass.to("kg/mol").magnitude

    fluid, _ = electrolyte_cpa_statoil_mixture_of(list(components))
    reduced = reduced_parameters(fluid, t_si, p_si)
    warnings.extend(reduced.warnings)
    kij = fluid.effective_kij(t_si)
    state = phase_state(reduced, kij, list(x), liquid=True)

    molality: list[float] = []
    mass_ratio: list[float] = []
    activity: list[float] = []
    for i in range(n):
        # The class's `mmol/kgSolvent` column is this times 1000 and its `mg/kgSolvent` is
        # this times the molar mass times `1e6`; the vocabulary states the quantity itself.
        per_kg = x[i] / x[water] / water_mass
        molality.append(per_kg)
        if entries[i].molar_mass is None:
            raise PropertyUnavailableError(
                components[i], "a molar mass", "the mass column is the molar one times it"
            )
        mass_ratio.append(per_kg * entries[i].molar_mass.to("kg/mol").magnitude)
        reference = (
            _pure_log_fugacity(components[i], t_si, p_si)
            if entries[i].reference_state == SOLVENT
            else _infinite_dilute_log_fugacity(components[i], components[water], t_si, p_si)
        )
        activity.append(math.exp(state.ln_phi[i] - reference))

    return CalcIonicCompositionResult(
        molality=tuple(from_si(value, "mol/kg") for value in molality),
        mass_ratio=tuple(mass_ratio),
        activity_coefficient=tuple(activity),
        warnings=tuple(warnings),
    )


def _pure_log_fugacity(name: str, t_si: float, p_si: float) -> float:
    """``ln phi`` of one component in a phase of itself, at the same ``T`` and ``P``."""
    fluid, _ = electrolyte_cpa_statoil_mixture_of([name])
    reduced = reduced_parameters(fluid, t_si, p_si)
    state = phase_state(reduced, fluid.effective_kij(t_si), [1.0], liquid=True)
    return state.ln_phi[0]


def _infinite_dilute_log_fugacity(name: str, solvent: str, t_si: float, p_si: float) -> float:
    """``ln phi`` of one component at infinite dilution in the solvent, at the same state.

    The binary is ``1e-10`` mol of the component beside ``10`` mol of the solvent, so the
    composition is that ratio and not a limit: the phase is solved there.
    """
    total = _DILUTE_MOLES + _SOLVENT_MOLES
    x = [_DILUTE_MOLES / total, _SOLVENT_MOLES / total]
    fluid, _ = electrolyte_cpa_statoil_mixture_of([name, solvent])
    reduced = reduced_parameters(fluid, t_si, p_si)
    state = phase_state(reduced, fluid.effective_kij(t_si), x, liquid=True)
    return state.ln_phi[0]


__all__ = ["calc_ionic_composition"]
