"""Chemical reaction equilibrium.

The first namespace whose subject is not a fluid's own state. A reaction reads a
substance's elemental composition and a fitted equilibrium constant, and the two
together are what decide whether a fluid can react and how far.

* :func:`equilibrium_constant` - one reaction's ``K``, its derivative and its heat

# Which implementation answers

As in the other namespaces, each function dispatches to the Rust extension when it is
built and to :mod:`azoth.reactions.reference` otherwise. Both are always reachable -
see :func:`azoth.backends` and :func:`azoth.use_backend`.
"""

from __future__ import annotations

from azoth._dispatch import resolve
from azoth.core.result import (
    ChemicalEquilibriumResult,
    EquilibriumConstantResult,
    ReactivePhaseEquilibriumResult,
    ReferencePotentialsResult,
)
from azoth.core.units import Q

__all__ = [
    "chemical_equilibrium",
    "equilibrium_constant",
    "reactive_phase_equilibrium",
    "reference_potentials",
]

_CHEMICAL_EQUILIBRIUM = "reactions.chemical_equilibrium"
_EQUILIBRIUM_CONSTANT = "reactions.equilibrium_constant"
_REACTIVE_PHASE_EQUILIBRIUM = "reactions.reactive_phase_equilibrium"
_REFERENCE_POTENTIALS = "reactions.reference_potentials"


def equilibrium_constant(source: str, reaction: str, T: Q) -> EquilibriumConstantResult:
    """One reaction's equilibrium constant, its temperature derivative and its heat.

    ``source`` is required rather than defaulting to ``"standard"``: the three reaction
    tables are different standard states whose answers part company away from 298 K, so
    a default here would be a silent choice of convention. ``reaction`` names the row,
    and ``T`` is an absolute temperature.

    Raises:
        InvalidInputError: if ``source`` is not one of ``"standard"``, ``"pitzer"`` or
            ``"kent-eisenberg"``.
        PropertyUnavailableError: if the source carries no row by that name.

    See :func:`azoth.reactions.reference.equilibrium_constant`.
    """
    return resolve(_EQUILIBRIUM_CONSTANT)(source=source, reaction=reaction, T=T)  # type: ignore[no-any-return]


def reference_potentials(components: list[str], source: str, T: Q) -> ReferencePotentialsResult:
    """The standard-state reference potentials of a fluid's reactive components.

    ``components`` is the reactive set by name, **in the order the potentials are wanted
    back** - the column indices the basis is chosen over are those positions. ``source``
    is required rather than defaulting, for the reason :func:`equilibrium_constant` gives.

    Raises:
        InvalidInputError: if ``source`` is unknown, or the propagation cannot reach a
            component, or a rank test meets a non-integral coefficient.

    See :func:`azoth.reactions.reference.reference_potentials`.
    """
    return resolve(_REFERENCE_POTENTIALS)(components=components, source=source, T=T)  # type: ignore[no-any-return]


def chemical_equilibrium(
    a_matrix: list[list[float]],
    b: list[float],
    whole_system: bool,
    moles: list[float],
    chem_ref: list[float],
    log_activity: list[float],
    T: Q,
    max_iterations: float,
    tolerance: float,
) -> ChemicalEquilibriumResult:
    """The reactive equilibrium composition of a phase, by the Smith-Missen method.

    ``a_matrix`` is the element matrix with the electroneutrality row last, ``b`` the
    element amounts it conserves, ``whole_system`` whether the phase is the only one (which
    is NeqSim's ``getNumberOfPhases() == 1`` and decides whether the conservation coupling
    is corrected), ``moles`` the starting composition, ``chem_ref`` the
    reduced standard-state potentials and ``log_activity`` the ``ln(gamma)``.

    Returns a result whether or not the solve converged; ``converged`` says which.

    Raises:
        InvalidInputError: where the shapes disagree.

    See :func:`azoth.reactions.reference.chemical_equilibrium`.
    """
    return resolve(_CHEMICAL_EQUILIBRIUM)(  # type: ignore[no-any-return]
        a_matrix=a_matrix,
        b=b,
        whole_system=whole_system,
        moles=moles,
        chem_ref=chem_ref,
        log_activity=log_activity,
        T=T,
        max_iterations=max_iterations,
        tolerance=tolerance,
    )


def reactive_phase_equilibrium(
    components: list[str],
    source: str,
    phase: str,
    moles: list[Q],
    phase_charge: Q,
    phase_moles: Q,
    whole_system: bool,
    log_activity: list[float],
    T: Q,
    max_iterations: float,
    tolerance: float,
) -> ReactivePhaseEquilibriumResult:
    """The reactive equilibrium composition of one phase, from its own state.

    ``components`` and ``moles`` are the phase's reactive substances and their amounts.
    ``phase`` is its type name as NeqSim spells it: ``aqueous``, ``liquid`` or ``oil``
    takes the solve and **anything else skips it**, which is how NeqSim's ``-1`` reaches
    a caller here. ``phase_charge`` and ``phase_moles`` are the whole phase's own, which
    is what the electroneutrality row's correction is read from, and ``whole_system``
    whether this phase is the only one.

    The matrix, the element amounts and the reference potentials come back either way.
    ``skipped`` says whether a solve ran, and it is what separates a skip from a solve
    that ran and did not converge.

    Raises:
        InvalidInputError: where the shapes disagree, or a component has no element row
            or no charge.

    See :func:`azoth.reactions.reference.reactive_phase_equilibrium`.
    """
    return resolve(_REACTIVE_PHASE_EQUILIBRIUM)(  # type: ignore[no-any-return]
        components=components,
        source=source,
        phase=phase,
        moles=moles,
        phase_charge=phase_charge,
        phase_moles=phase_moles,
        whole_system=whole_system,
        log_activity=log_activity,
        T=T,
        max_iterations=max_iterations,
        tolerance=tolerance,
    )
