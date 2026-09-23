"""``reactions.reactive_hybrid_eos_ge_flash`` - the reactive coupling on the hybrid flash.

The pure-Python implementation, the twin of
``crates/azoth-reactions/src/reactive_hybrid_eos_ge_flash.rs``. ``eos.hybrid_eos_ge_flash``
solves one gas-oil-brine split at a *fixed* species inventory and NeqSim wraps it in a second
loop when ``system.isChemicalSystem()``: each pass re-equilibrates the brine, projects the
reaction delta onto the stoichiometric conservation null space, writes the adjusted inventory
back, and solves the fractions again.

**The chemistry's input is the phase's activity vector**, which ``solveChemEq`` reads once
before its solve. :func:`log_activities` recovers it from ``eos.pitzer_phase``, and what that
recovery does and does not claim is the spec's business.
"""

from __future__ import annotations

import math
from collections.abc import Sequence

from azoth.core.errors import InvalidInputError, SolverNotConvergedError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import ReactiveHybridEosGeFlashResult
from azoth.core.units import Q, input_to_si, quantity
from azoth.core.warnings import Warning
from azoth.eos import components as _components
from azoth.eos.reference import hybrid_eos_ge_flash as _hybrid
from azoth.eos.reference.pitzer_phase import pitzer_phase
from azoth.reactions.reference import _linalg, _tables
from azoth.reactions.reference.reactive_phase_equilibrium import (
    _element_matrix,
    reactive_phase_equilibrium,
)
from azoth.reactions.reference.reference_potentials import side_is_present

__all__ = ["log_activities", "reactive_components", "reactive_hybrid_eos_ge_flash"]

#: The reaction source every fluid of this model runs, from
#: ``SystemPitzer.getChemicalReactionDataSource()``. The hybrid route is reachable only from an
#: EoS/GE system, so the source is a property of the model and not an input.
SOURCE = "pitzer"

#: The coupled loop's cap and floor, ``MAXIMUM_REACTIVE_ITERATIONS`` and
#: ``MINIMUM_REACTIVE_ITERATIONS``.
MAXIMUM_REACTIVE_PASSES = 100
MINIMUM_REACTIVE_PASSES = 3

#: The two stopping rules: the brine's composition and the fraction solve's own residual.
REACTIVE_COMPOSITION_TOLERANCE = 1.0e-10
HYBRID_SOLVER_TOLERANCE = 1.0e-10

#: The floor the coupled inventory is raised to, and the two scales on a reaction delta.
MINIMUM_COUPLED_MOLES = 1.0e-45
REACTION_DELTA_CONSERVATION_TOLERANCE_MOLES = 1.0e-8
NEGATIVE_INVENTORY_TOLERANCE_MOLES = 1.0e-9

#: ``ChemicalEquilibrium``'s own defaults, which ``solveChemEq``'s caller passes through.
CHEMISTRY_MAX_ITERATIONS = 100
CHEMISTRY_TOLERANCE = 1.0e-8

#: The class's two-component reference phase: the solvent at ten mol, the solute at ``1e-10``.
REFERENCE_SOLVENT_MOLES = 10.0
REFERENCE_SOLUTE_MOLES = 1.0e-10

#: The aqueous role's index, in ``[gas, oil, aqueous]`` order.
AQUEOUS_ROLE = 2


def _spec() -> dict[str, object]:
    """The generated spec, imported here so the module has no import-time cycle."""
    from azoth._models_gen import model

    return model("reactions.reactive_hybrid_eos_ge_flash")


def _si(spec: dict[str, object], name: str, value: float | Q) -> float:
    """A declared input as its SI magnitude, quantity or not."""
    from azoth.reactions.reference.chemical_equilibrium import _si as shared

    return shared(spec, name, value)


def reactive_hybrid_eos_ge_flash(
    components: list[str],
    cubic: str,
    T: Q,
    P: Q,
    moles: list[float],
) -> ReactiveHybridEosGeFlashResult:
    """The coupled reactive flash: chemistry and a fixed gas-oil-brine topology at once.

    Args:
        components: the substances the feed is made of, by name.
        cubic: the equation of state **both** EoS roles are.
        T: absolute temperature.
        P: absolute pressure.
        moles: the feed's mole numbers. The inventory the passes hand each other is not this
            vector: a reaction changes the number of species moles, and the answer reports
            the adjusted one.

    Returns:
        The roles' fractions and compositions, the coupled inventory, the brine's species
        amounts, and the numbers the loop stopped on.

    Raises:
        InvalidInputError: if a name is unknown, a component has no formula or no charge row,
            or the chemistry drives an amount past its tolerance.
        OutOfRangeError: if ``T`` or ``P`` is not positive.
        PropertyUnavailableError: if a phase model cannot answer, or the source's evidence
            gate refuses an active reaction row.
        SolverNotConvergedError: if the coupled loop reaches its cap without certifying.

    Example:
        >>> import azoth
        >>> q = azoth.ureg.Quantity
        >>> result = reactive_hybrid_eos_ge_flash(
        ...     ["methane", "CO2", "water", "Ca++", "Cl-", "HCO3-", "CO3--", "OH-", "H3O+"],
        ...     "srk",
        ...     q(313.15, "K"),
        ...     q(50.0e5, "Pa"),
        ...     [5.0, 0.05, 55.5, 6.0e-4, 2.0e-4, 1.0e-3, 1.0e-10, 1.0e-10, 1.0e-10],
        ... )
        >>> round(result.beta[0], 6)
        0.083462
    """
    spec = _spec()
    warnings: list[Warning] = []

    t_si = input_to_si(spec, "T", T)
    p_si = input_to_si(spec, "P", P)
    apply_checks(checks_for(spec).on_input, {"T": t_si, "P": p_si}.get, warnings)

    if len(components) != len(moles):
        raise InvalidInputError(
            "moles", f"{len(moles)} entry(ies) against {len(components)} component(s)"
        )

    fraction_algorithm = _fraction_algorithm()
    outcome = solve_coupled(
        components,
        cubic,
        t_si,
        p_si,
        # A model's declared units are the boundary's: the vector carries a unit, so it crosses
        # as SI magnitudes element by element and no arithmetic is done on a quantity.
        [_si(spec, "moles", value) for value in moles],
        tolerance=float(fraction_algorithm["tolerance"]),
        cap=int(fraction_algorithm["max_iterations"]),
    )

    return ReactiveHybridEosGeFlashResult(
        beta=tuple(outcome["beta"]),
        x=tuple(tuple(row) for row in outcome["x"]),
        coupled_moles=tuple(quantity(value, "mol") for value in outcome["coupled_moles"]),
        aqueous_moles=tuple(quantity(value, "mol") for value in outcome["aqueous_moles"]),
        passes=outcome["passes"],
        chemical_deviation=outcome["chemical_deviation"],
        residual=outcome["residual"],
        max_material_balance_residual=outcome["max_material_balance_residual"],
        max_log_fugacity_residual=outcome["max_log_fugacity_residual"],
        element_residual=outcome["element_residual"],
        charge_residual=outcome["charge_residual"],
        warnings=tuple(warnings),
    )


def _fraction_algorithm() -> dict[str, object]:
    """The fraction solve's own algorithm, which is the underlying model's.

    This loop re-enters that solve rather than running a second one, so it takes its tolerance
    and cap from ``eos.hybrid_eos_ge_flash``'s spec instead of declaring a second copy.
    """
    from azoth._models_gen import model

    return model("eos.hybrid_eos_ge_flash")["algorithm"]


def solve_coupled(
    components: list[str],
    cubic: str,
    t_si: float,
    p_si: float,
    moles: list[float],
    *,
    tolerance: float,
    cap: int,
    chemistry_max_iterations: int = CHEMISTRY_MAX_ITERATIONS,
    chemistry_tolerance: float = CHEMISTRY_TOLERANCE,
) -> dict[str, object]:
    """The coupled loop itself, over names resolved inside it.

    Kept separate from the public function so the two kernels' difference is the arithmetic
    and not the boundary.
    """
    if not components or sum(moles) <= 0.0 or not math.isfinite(sum(moles)):
        raise InvalidInputError("moles", "the feed sums to nothing that is a mole count")

    entries = _hybrid._resolved(components)
    ion = [_hybrid._is_ion(entry) for entry in entries]
    mixture, _ = _components.mixture_of_with_ions(components, eos=cubic)
    reduced = _hybrid.reduced_parameters(mixture, t_si, p_si)

    reactive = reactive_components(components)
    if not reactive:
        raise InvalidInputError(
            "components", "no surviving reaction of the pitzer source names any of these"
        )
    a_matrix = _element_matrix(reactive)
    basis = _linalg.row_space_basis(a_matrix)
    reactive_at = [components.index(name) for name in reactive]

    total = sum(moles)
    z = [value / total for value in moles]
    phases = _hybrid._seed(entries, z, ion)

    # `run`'s first solve, outside the loop: chemistry would otherwise meet the role seeds.
    _hybrid.solve_fixed_topology_from(
        reduced, mixture.kij, phases, moles, ion, components, t_si, p_si, tolerance, cap
    )

    coupled = [float(value) for value in moles]
    passes = 0
    coupled_converged = False
    chemical_deviation = math.inf
    residual = math.nan
    for index in range(MAXIMUM_REACTIVE_PASSES):
        passes = index + 1
        chemical_deviation = _chemical_step(
            components,
            reactive,
            reactive_at,
            phases,
            coupled,
            a_matrix,
            basis,
            t_si,
            p_si,
            initialise=index == 0,
            chemistry_max_iterations=chemistry_max_iterations,
            chemistry_tolerance=chemistry_tolerance,
        )
        # The adjusted inventory is what the next fraction solve balances, which is what
        # `synchronizeHybridEosGeOverallComposition` does on NeqSim's own system.
        _, residual, gradient_norm = _hybrid.solve_fixed_topology_from(
            reduced, mixture.kij, phases, coupled, ion, components, t_si, p_si, tolerance, cap
        )
        residual = max(residual, gradient_norm)

        coupled_converged = (
            passes >= MINIMUM_REACTIVE_PASSES
            and chemical_deviation <= REACTIVE_COMPOSITION_TOLERANCE
            and math.isfinite(residual)
            and residual <= HYBRID_SOLVER_TOLERANCE
        )
        if coupled_converged:
            break

    if not coupled_converged:
        raise SolverNotConvergedError(
            passes, max(residual, chemical_deviation), HYBRID_SOLVER_TOLERANCE
        )

    element_residual, charge_residual = _conservation_residuals(
        components, a_matrix, reactive_at, phases, coupled
    )
    coupled_total = sum(coupled)

    return {
        "beta": [phase.fraction for phase in phases],
        "x": [list(phase.composition) for phase in phases],
        "coupled_moles": coupled,
        "aqueous_moles": [
            phases[AQUEOUS_ROLE].composition[index] * phases[AQUEOUS_ROLE].fraction * coupled_total
            for index in reactive_at
        ],
        "passes": passes,
        "chemical_deviation": chemical_deviation,
        "residual": residual,
        "max_material_balance_residual": _material_balance(phases, coupled),
        "max_log_fugacity_residual": _log_fugacity_residual(
            mixture, reduced, phases, components, coupled, t_si, p_si, ion
        ),
        "element_residual": element_residual,
        "charge_residual": charge_residual,
    }


def _chemical_step(
    components: list[str],
    reactive: list[str],
    reactive_at: list[int],
    phases: list[object],
    coupled: list[float],
    a_matrix: list[list[float]],
    basis: list[list[float]],
    t_si: float,
    p_si: float,
    *,
    initialise: bool,
    chemistry_max_iterations: int,
    chemistry_tolerance: float,
) -> float:
    """One pass's chemistry: solve the brine, project the delta, write the inventory back.

    Returns the composition deviation the convergence test reads, taken after the chemistry
    wrote the brine back and before the fractions are solved again.
    """
    aqueous = phases[AQUEOUS_ROLE]
    total = sum(coupled)
    phase_moles = aqueous.fraction * total

    phase_amounts = [max(fraction * phase_moles, 0.0) for fraction in aqueous.composition]
    old_fractions = [aqueous.composition[index] for index in reactive_at]
    old_moles = [phase_amounts[index] for index in reactive_at]

    # The phase's own charge over *every* component it holds, which is the quantity
    # `calcBVector` subtracts the reactive set's share from.
    phase_charge = 0.0
    for index, name in enumerate(components):
        charge = _tables.ionic_charge(name)
        if charge is not None:
            phase_charge += charge * phase_amounts[index]

    activity = log_activities(components, t_si, p_si, aqueous.composition)
    log_activity = [activity[index] for index in reactive_at]

    solved = reactive_phase_equilibrium(
        reactive,
        SOURCE,
        "aqueous",
        [quantity(value, "mol") for value in old_moles],
        quantity(phase_charge, "mol"),
        quantity(phase_moles, "mol"),
        False,
        log_activity,
        quantity(t_si, "K"),
        chemistry_max_iterations,
        chemistry_tolerance,
        "solute_molality",
        "linear_programming" if initialise else "none",
    )

    # `updateMoles`'s answer, back in SI moles: the result carries quantities and the loop
    # works in magnitudes.
    solved_moles = [float(value.to("mol").magnitude) for value in solved.moles]
    raw = [new - old for new, old in zip(solved_moles, old_moles, strict=True)]
    delta = _conservative_delta(a_matrix, basis, coupled, reactive_at, raw)
    _apply_delta(coupled, reactive_at, delta)

    # `updateMoles` writes the answer into the phase and the deviation is read from there, so
    # the phase's own total moves with the reactive components and the others do not.
    for position, index in enumerate(reactive_at):
        phase_amounts[index] = max(old_moles[position] + raw[position], MINIMUM_COUPLED_MOLES)
    new_total = sum(phase_amounts)
    if new_total <= 0.0:
        return 0.0
    return sum(
        abs(old_fractions[position] - phase_amounts[index] / new_total)
        for position, index in enumerate(reactive_at)
    )


def _conservative_delta(
    a_matrix: list[list[float]],
    basis: list[list[float]],
    coupled: list[float],
    reactive_at: list[int],
    raw: list[float],
) -> list[float]:
    """The conservative part of a reaction delta, with the two shortcuts in front of it.

    A delta that already satisfies the element balance to
    :data:`REACTION_DELTA_CONSERVATION_TOLERANCE_MOLES` and leaves no inventory negative is
    used as it stands; otherwise it is projected onto the conservation null space. **Measured:
    the shortcut is the branch the captured fluids take on every pass** - the worst raw
    residual is ``2.0e-10`` - so the reactive oracle does not exercise the projection.

    The last step scales the delta whole when an inventory would go negative, which is
    ``feasibleStep``: scaling together keeps the delta a reaction direction, where clipping one
    component would not.
    """
    worst = max(
        abs(sum(coefficient * delta for coefficient, delta in zip(row, raw, strict=True)))
        for row in a_matrix
    )
    non_negative = all(
        coupled[index] + delta >= -NEGATIVE_INVENTORY_TOLERANCE_MOLES
        for delta, index in zip(raw, reactive_at, strict=True)
    )

    if worst <= REACTION_DELTA_CONSERVATION_TOLERANCE_MOLES and non_negative:
        delta = list(raw)
    else:
        delta = _linalg.project_onto_null_space(basis, raw)

    step = 1.0
    for value, index in zip(delta, reactive_at, strict=True):
        if value >= 0.0:
            continue
        available = max(coupled[index] - MINIMUM_COUPLED_MOLES, 0.0)
        step = min(step, available / -value)
    if step < 1.0:
        delta = [value * max(step, 0.0) for value in delta]
    return delta


def _apply_delta(coupled: list[float], reactive_at: list[int], delta: list[float]) -> None:
    """Add a reaction delta to the coupled inventory, floored as `updateCoupledOverallComposition`
    floors it and refusing a negative amount past its own tolerance."""
    for position, index in enumerate(reactive_at):
        updated = coupled[index] + delta[position]
        if not math.isfinite(updated) or updated < -NEGATIVE_INVENTORY_TOLERANCE_MOLES:
            raise InvalidInputError(
                "moles",
                f"aqueous chemical equilibrium produced an invalid overall amount "
                f"({updated}) for component {index}",
            )
        coupled[index] = max(updated, MINIMUM_COUPLED_MOLES)


def log_activities(
    components: list[str],
    t_si: float,
    p_si: float,
    x: Sequence[float],
) -> tuple[float, ...]:
    """``ChemicalEquilibrium.calcRefPot``'s vector, recovered from ``eos.pitzer_phase``.

    ``log_activity_i = ln phi_i - ln phi_i^reference``, where the reference is the pure
    component for a solvent and a **two-component reference phase** - the solute at
    :data:`REFERENCE_SOLUTE_MOLES`, water at :data:`REFERENCE_SOLVENT_MOLES` - for everything
    else. ``pitzer_phase`` reports both halves of it: ``ln_gamma`` is the live arm's
    ``ln phi`` up to the branch's own factors, and ``gamma_inf`` is the infinite-dilution
    coefficient the ionic arm divides by.

    **Measured against the captured vector**, on two passes and every component: water matches
    to ``3e-16`` as ``ln_gamma`` alone, every ion to ``3e-15`` as ``ln_gamma - ln gamma_inf``,
    and a neutral that is not water takes the same rule plus the reference phase's own
    ``ln(m/x)`` shift, ``-ln x_water + ln x_water^reference``.

    What this does not claim: the neutral arm's ``gamma_inf`` is reported as one where
    ``pitzer_phase`` does not divide by it, so a strongly non-ideal neutral keeps its full
    ``ln gamma`` where NeqSim's reference phase would divide it out. No captured fluid
    exercises that.

    Raises:
        InvalidInputError: if the fluid carries no water, which is the component the reference
            is taken against.
    """
    activity = pitzer_phase(components, quantity(t_si, "K"), quantity(p_si, "Pa"), list(x))
    water = next((index for index, name in enumerate(components) if name.lower() == "water"), None)
    if water is None:
        raise InvalidInputError(
            "components", "the activity vector is taken against water and this fluid carries none"
        )

    shift = -math.log(x[water]) + math.log(
        REFERENCE_SOLVENT_MOLES / (REFERENCE_SOLVENT_MOLES + REFERENCE_SOLUTE_MOLES)
    )
    out: list[float] = []
    for index, name in enumerate(components):
        value = activity.ln_gamma[index] - math.log(activity.gamma_inf[index])
        charge = _tables.ionic_charge(name)
        if not (charge is not None and charge != 0.0) and index != water:
            value += shift
        out.append(value)
    return tuple(out)


def reactive_components(components: list[str]) -> list[str]:
    """The components a fluid drives out of :data:`SOURCE`'s table, **in the fluid's order**.

    NeqSim reaches this through ``readReactions`` -> ``removeJunkReactions(componentNames)`` ->
    ``getAllComponents()``: the reactions the fluid can support are kept, and the set is the
    union of *their* components. A component the fluid carries but no surviving reaction names
    is not reactive and must not be solved for - which is how ``Ca++`` and ``Cl-``, both of
    which the element table covers and neither of which the pitzer source's reactions name for
    a carbonate brine, stay out of its chemistry.

    **The order is this library's own.** NeqSim's is a ``HashSet`` iteration, a property of
    Java's string hash and not of the chemistry; every vector that depends on this order moves
    with it, so the order is stated rather than reproduced.
    """
    present = {name: index for index, name in enumerate(components)}
    named: list[str] = []
    for row in _tables.reactions(SOURCE):
        if not row.use_reaction:
            continue
        coefficients = _tables.stoichiometry(row.name)
        if not side_is_present(coefficients, present, negative=True) and not side_is_present(
            coefficients, present, negative=False
        ):
            continue
        for component, _ in coefficients:
            if component in present and component not in named:
                named.append(component)
    return [name for name in components if name in named]


def _conservation_residuals(
    components: list[str],
    a_matrix: list[list[float]],
    reactive_at: list[int],
    phases: list[object],
    coupled: list[float],
) -> tuple[float, float]:
    """The split's element and charge residuals against the coupled inventory it was solved at.

    An element row is taken over the reactive components alone, and the charge row over every
    component's charge, which is how ``SystemHybridEosGeFlashTest`` builds the quantity.
    """
    summed = _phase_sums(phases, components, sum(coupled))
    rows = len(a_matrix)
    worst = 0.0
    for row, coefficients in enumerate(a_matrix):
        if row + 1 == rows:
            continue
        difference = sum(
            coefficient * (coupled[reactive_at[column]] - summed[reactive_at[column]])
            for column, coefficient in enumerate(coefficients)
        )
        worst = max(worst, abs(difference))

    charge = 0.0
    for index, name in enumerate(components):
        value = _tables.ionic_charge(name)
        if value is not None:
            charge += value * (coupled[index] - summed[index])
    return worst, abs(charge)


def _phase_sums(phases: list[object], components: list[str], total: float) -> list[float]:
    """Each component's amount summed over the roles, `sum_k beta_k x_ik n`."""
    return [
        sum(phase.fraction * phase.composition[index] for phase in phases) * total
        for index in range(len(components))
    ]


def _material_balance(phases: list[object], coupled: list[float]) -> float:
    """The worst `|z_i - sum_k beta_k x_ik|` over the coupled inventory."""
    total = sum(coupled)
    z = [value / total for value in coupled]
    return max(
        abs(z[index] - sum(phase.fraction * phase.composition[index] for phase in phases))
        for index in range(len(coupled))
    )


def _log_fugacity_residual(
    mixture: object,
    reduced: object,
    phases: list[object],
    components: list[str],
    coupled: list[float],
    t_si: float,
    p_si: float,
    ion: list[bool],
) -> float:
    """The worst spread of ``ln(x_i phi_i P)`` over the roles that are there.

    The underlying model's own acceptance quantity, over the coupled inventory rather than the
    feed - the roles were solved at the adjusted one.
    """
    ln_phi = [
        [math.log(value) for value in row]
        for row in _hybrid._coefficients(
            reduced, mixture.kij, phases, components, t_si, p_si, AQUEOUS_ROLE
        )
    ]
    total = sum(coupled)
    _, fugacity = _hybrid._acceptance(
        [phase.fraction for phase in phases],
        [list(phase.composition) for phase in phases],
        ln_phi,
        [value / total for value in coupled],
        ion,
        p_si,
    )
    return float(fugacity)
