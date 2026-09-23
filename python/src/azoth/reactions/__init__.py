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

from collections.abc import Sequence

from azoth._dispatch import resolve
from azoth.core.result import (
    ChemicalEquilibriumResult,
    EquilibriumConstantResult,
    ReactivePhaseEquilibriumResult,
    ReactiveTpFlashResult,
    ReferencePotentialsResult,
)
from azoth.core.units import Q

__all__ = [
    "chemical_equilibrium",
    "equilibrium_constant",
    "reactive_phase_equilibrium",
    "reactive_tp_flash",
    "reference_potentials",
]

_CHEMICAL_EQUILIBRIUM = "reactions.chemical_equilibrium"
_EQUILIBRIUM_CONSTANT = "reactions.equilibrium_constant"
_REACTIVE_PHASE_EQUILIBRIUM = "reactions.reactive_phase_equilibrium"
_REACTIVE_TP_FLASH = "reactions.reactive_tp_flash"
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
    concentration_basis: str,
    solvent_weight: Q,
    solvent_mask: Sequence[float],
    phase_moles: Q,
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
        concentration_basis=concentration_basis,
        solvent_weight=solvent_weight,
        solvent_mask=solvent_mask,
        phase_moles=phase_moles,
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
    seed: str = "none",
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
        seed=seed,
        phase_moles=phase_moles,
        whole_system=whole_system,
        log_activity=log_activity,
        T=T,
        max_iterations=max_iterations,
        tolerance=tolerance,
    )


def reactive_tp_flash(
    components: Sequence[str],
    T: Q,
    P: Q,
    moles: Sequence[Q],
    max_phases: float,
) -> ReactiveTpFlashResult:
    """Simultaneous chemical and phase equilibrium at fixed temperature and pressure.

    ``moles`` is the **overall component amounts and not a composition**, which is what
    the driver's ``getOverallMoles`` reads and what its frozen element inventory is built
    from. ``max_phases`` is its effective phase ceiling; a ceiling of one skips the VLE
    initialisation and collapses the phase list, so it changes the answer rather than the
    format.

    **The composition is the answer and the split is not.** Where two phases converge to
    the same composition the Gibbs energy is flat along the direction that trades moles
    between them, and the split a run reports is decided by its path. The moles summed
    over the returned rows are what the element balance and the equilibrium fix.

    ``gibbs_energy`` is ``computeGibbsEnergy``, which weighs each phase by the fraction
    its list carries - so on the class's own two-phase states it is twice one phase's
    worth, and on the single-phase branch it is ``0.0``.

    Raises:
        InvalidInputError: for a charged component, a shape disagreement, or a component
            the databank does not carry.

    See :func:`azoth.reactions.reference.reactive_tp_flash`.
    """
    return resolve(_REACTIVE_TP_FLASH)(  # type: ignore[no-any-return]
        components=list(components),
        T=T,
        P=P,
        moles=list(moles),
        max_phases=max_phases,
    )
