"""``reactions.reactive_phase_equilibrium`` - the solve as a phase operation.

Spec: ``specs/models/reactions/reactive_phase_equilibrium.toml``

The Python twin of ``crates/azoth-reactions/src/reactive_phase_equilibrium.rs``, written
to mirror it line for line.

# What is built and what is skipped

NeqSim's ``ChemicalReactionOperations`` chooses the phase, builds the element matrix and
the element amounts from it, solves, and writes the composition back. The matrix, the
amounts and the reference potentials are built **whether or not the phase is reactive**,
because that is what NeqSim's own constructor does; only the solve is conditional.

# The charge row is a correction and not a zero

```text
inert = sum(z_i n_i over the phase) - (A n)[charge]
b[charge] = |inert| <= 1e-10 max(1, n_phase) ? 0 : -inert
```

The spectator ions are not in ``components``, which is why the phase's own charge and
mole count are inputs. The row is zero only when every ion in the phase is in the
reaction set; the captured bicarbonate brine is the fluid where it is not.
"""

from __future__ import annotations

import math
from collections.abc import Sequence

from azoth.core.errors import InvalidInputError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import ReactivePhaseEquilibriumResult
from azoth.core.units import Q, from_si, input_to_si, quantity
from azoth.core.warnings import Warning
from azoth.eos import components as _components
from azoth.reactions.reference import _lp_seed, _tables
from azoth.reactions.reference.chemical_equilibrium import MIN_MOLES, chemical_equilibrium
from azoth.reactions.reference.equilibrium_constant import (
    GAS_CONSTANT,
    equilibrium_constant,
)
from azoth.reactions.reference.reference_potentials import reference_potentials

#: The charge, in moles of elementary charge, below which the phase is treated as
#: neutral, from ``calcBVector``'s ``1e-10 * Math.max(1.0, phase moles)``.
CHARGE_NOISE_MOLE_FRACTION = 1e-10

#: The floor ``updateMoles`` raises every mole number to before writing it back, from
#: ``Math.max(newMoles[i], 1e-45)``.
MIN_WRITTEN_MOLES = 1e-45

#: The reaction log residual `solveChemEq` requires, from `REACTION_LOG_RESIDUAL_TOLERANCE`.
REACTION_LOG_RESIDUAL_TOLERANCE = 2.0e-6

#: The net charge it requires, in moles of elementary charge.
REACTIVE_PHASE_CHARGE_TOLERANCE_MOLES = 1.0e-8

#: The element-balance residual it requires, in mol.
ELEMENT_BALANCE_RESIDUAL_TOLERANCE_MOLES = 1.0e-8

#: The three phase type names NeqSim solves reactions in.
REACTIVE_PHASE_LABELS = ("aqueous", "liquid", "oil")


def is_reactive_phase(label: str) -> bool:
    """Whether a phase of this type takes the reaction solve.

    The comparison is case-insensitive, as ``equalsIgnoreCase`` is. A name that is none of
    the three - ``gas`` above all - is the skip.
    """
    return label.strip().lower() in REACTIVE_PHASE_LABELS


def reactive_phase_index(labels: list[str]) -> int | None:
    """The first phase of ``labels`` that takes the reaction solve, or ``None``.

    **``aqueous`` is looked for over every phase before ``liquid`` or ``oil`` is looked
    for over any**, which is two passes and not one: a liquid phase listed ahead of an
    aqueous one does not win. Within each pass the first match in the caller's order is
    the answer, so the order is part of the input.

    ``None`` is the skip, and the skip is not a failure: it is what NeqSim returns as
    ``-1``, and a solve that failed is a different result.
    """
    folded = [label.strip().lower() for label in labels]
    if "aqueous" in folded:
        return folded.index("aqueous")
    for index, label in enumerate(folded):
        if label in ("liquid", "oil"):
            return index
    return None


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
    concentration_basis: str = "mole_fraction",
    seed: str = "none",
) -> ReactivePhaseEquilibriumResult:
    """The reactive equilibrium composition of one phase.

    Args:
        components: the phase's reactive substances - not every substance in it.
        source: which of the three reaction tables the reference potentials come from.
        phase: the phase's type name as NeqSim spells it. ``aqueous``, ``liquid`` or
            ``oil`` takes the solve; anything else skips it.
        moles: their amounts in the phase, ``A n``'s operand and the starting
            composition wherever the linear-program estimate is absent.
        phase_charge: the phase's own ``sum(z_i n_i)``, over **every** component it
            holds, in moles of elementary charge.
        phase_moles: the phase's total moles, read for the charge-row tolerance.
        whole_system: whether this phase is the only one, which is NeqSim's
            ``getNumberOfPhases() == 1`` and the condition its solve corrects the
            conservation coupling under.
        log_activity: each component's ``ln(gamma)``, held fixed for the solve.
        T: absolute temperature.
        max_iterations: the solve's pass cap.
        tolerance: the error the solve is trying to reach.

    Returns:
        The composition with ``updateMoles``'s floor applied, or the caller's own where
        the solve was skipped, beside the matrix, the element amounts and the potentials
        that were built either way. **A skip is an answer**: ``skipped`` is what separates
        it from a solve that ran and did not converge.

    Raises:
        InvalidInputError: where the shapes disagree, or a component has no element row or
            no charge.
        OutOfRangeError: for a temperature, tolerance or cap the spec bounds.

    Example:
        >>> import azoth
        >>> q = azoth.ureg.Quantity
        >>> r = reactive_phase_equilibrium(
        ...     ["CO2", "water"],
        ...     "standard",
        ...     "gas",
        ...     [q(0.01, "mol"), q(9.0, "mol")],
        ...     q(0.0, "mol"),
        ...     q(9.01, "mol"),
        ...     [0.0, 0.0],
        ...     q(298.15, "K"),
        ...     100.0,
        ...     1e-8,
        ... )
        >>> r.skipped
        True
    """
    spec = _spec()
    checks = checks_for(spec)
    warnings: list[Warning] = []
    temperature = input_to_si(spec, "T", T)

    apply_checks(
        checks.on_input,
        {
            "T": temperature,
            "tolerance": tolerance,
            "max_iterations": max_iterations,
        }.get,
        warnings,
    )

    if len(components) != len(moles) or len(components) != len(log_activity):
        raise InvalidInputError(
            "components",
            f"{len(components)} component(s), {len(moles)} mole(s) and "
            f"{len(log_activity)} activity coefficient(s)",
        )

    # **Every declared input goes through `_si`, element by element.** The vectors carry
    # units, so the case hands them over as pint quantities and arithmetic on one against
    # a bare float raises rather than converts. Everything below is an SI magnitude;
    # `from_si` puts the units back on the way out.
    magnitudes = [_si(spec, "moles", value) for value in moles]
    activity = [_si(spec, "log_activity", value) for value in log_activity]

    a_matrix = _element_matrix(components)
    b = _element_amounts(a_matrix, magnitudes)

    charge_row = len(b) - 1
    inert_charge = input_to_si(spec, "phase_charge", phase_charge) - b[charge_row]
    noise = CHARGE_NOISE_MOLE_FRACTION * max(input_to_si(spec, "phase_moles", phase_moles), 1.0)
    b[charge_row] = 0.0 if abs(inert_charge) <= noise else -inert_charge

    potentials = reference_potentials(components, source, T)
    warnings.extend(potentials.warnings)
    chem_ref = potentials.potentials

    if not is_reactive_phase(phase):
        return ReactivePhaseEquilibriumResult(
            skipped=True,
            a_matrix=tuple(tuple(row) for row in a_matrix),
            b=tuple(from_si(value, "mol") for value in b),
            chem_ref=chem_ref,
            moles=tuple(from_si(value, "mol") for value in magnitudes),
            iterations=0,
            error=0.0,
            converged=False,
            # Nothing was solved, so nothing is certified: NeqSim's net-charge reader is NaN
            # without a reactive phase and its two other residuals have no phase to be
            # computed on.
            refinements=0,
            certified=False,
            max_reaction_log_residual=math.nan,
            net_charge_moles=from_si(math.nan, "mol"),
            max_element_residual=from_si(math.nan, "mol"),
            # Nothing was solved, so no estimate was applied: the caller's own composition
            # comes back as both the seed and the answer.
            seed_applied=False,
            seed_moles=tuple(from_si(value, "mol") for value in magnitudes),
            warnings=tuple(warnings),
        )

    # The solver takes the potentials reduced, which is the form its own iteration adds a
    # logarithm to. NeqSim divides by `R T` at the same boundary.
    reduced = [
        float(potential.to("J/mol").magnitude) / (GAS_CONSTANT * temperature)
        for potential in chem_ref
    ]
    # NeqSim's seed: the estimate is written back into the phase through `updateMoles`
    # before the solve starts, so it is the starting composition where one exists.
    estimate = None if seed == "none" else _lp_seed.initial_estimate(a_matrix, b, reduced)
    seed_moles = list(magnitudes) if estimate is None else estimate
    start = (
        list(magnitudes)
        if estimate is None
        else [max(value, MIN_WRITTEN_MOLES) for value in seed_moles]
    )

    # **The basis is the caller's, and its other two facts are derived here.** NeqSim reads
    # all three off its system; this operation is handed vectors, so the mask is the
    # components whose reference state is `solvent` and the solvent's mass is `sum(n_j M_j)`
    # over them - both properties of the moles and the component databank, and neither
    # something a caller should have to restate.
    mask = _solvent_mask(components)
    solvent_weight = _solvent_weight(
        components, [value.to("mol").magnitude for value in moles], mask
    )
    solved = chemical_equilibrium(
        a_matrix,
        b,
        whole_system,
        start,
        reduced,
        activity,
        T,
        max_iterations,
        tolerance,
        concentration_basis,
        from_si(solvent_weight, "kg"),
        mask,
        phase_moles,
    )
    warnings.extend(solved.warnings)

    # NeqSim's certificate, on the composition the solve left. **`converged` is this and not
    # the solver's own flag**: `solveChemEq` returns false wherever one of the three
    # residuals is over its tolerance, and on all five captured fluids it is.
    certificate = _certify(
        components,
        source,
        temperature,
        a_matrix,
        b,
        solved.moles,
        magnitudes,
        activity,
        input_to_si(spec, "phase_charge", phase_charge),
        input_to_si(spec, "phase_moles", phase_moles),
    )
    max_reaction_log_residual, net_charge_moles, max_element_residual = certificate
    certified = (
        solved.converged
        and max_reaction_log_residual <= REACTION_LOG_RESIDUAL_TOLERANCE
        and abs(net_charge_moles) <= REACTIVE_PHASE_CHARGE_TOLERANCE_MOLES
        and max_element_residual <= ELEMENT_BALANCE_RESIDUAL_TOLERANCE_MOLES
    )

    return ReactivePhaseEquilibriumResult(
        skipped=False,
        seed_applied=estimate is not None,
        seed_moles=tuple(from_si(value, "mol") for value in seed_moles),
        a_matrix=tuple(tuple(row) for row in a_matrix),
        b=tuple(from_si(value, "mol") for value in b),
        chem_ref=chem_ref,
        # `updateMoles`'s floor, applied to the value the phase is left holding.
        moles=tuple(
            from_si(max(float(value.to("mol").magnitude), MIN_WRITTEN_MOLES), "mol")
            for value in solved.moles
        ),
        iterations=solved.iterations,
        error=solved.error,
        converged=certified,
        refinements=1,
        certified=certified,
        max_reaction_log_residual=max_reaction_log_residual,
        net_charge_moles=from_si(net_charge_moles, "mol"),
        max_element_residual=from_si(max_element_residual, "mol"),
        warnings=tuple(warnings),
    )


def _si(spec: dict[str, object], name: str, value: float | Q) -> float:
    """A declared input as its SI magnitude, quantity or not.

    **The boundary is mixed by design.** A vector whose spec declares a unit arrives as a
    pint quantity and has to be converted; a dimensionless one arrives as the bare number
    it is, and ``input_to_si`` refuses that rather than passing it through. The rule is
    the solve's own and is imported rather than written twice.
    """
    from azoth.reactions.reference.chemical_equilibrium import _si as shared

    return shared(spec, name, value)


def _certify(
    components: list[str],
    source: str,
    temperature: float,
    a_matrix: list[list[float]],
    b: list[float],
    answer: Sequence[Q],
    input_moles: list[float],
    log_activity: list[float],
    phase_charge: float,
    phase_moles: float,
) -> tuple[float, float, float]:
    """NeqSim's three-residual certificate, evaluated at the answer.

    The three are `solveChemEq`'s own gate (``ChemicalReactionOperations.java:652-654``), and
    one of them needs the reaction set back: ``max |ln Q - ln K|`` over the reactions this
    fluid can run, with ``Q`` from the same activity term the solve used. The other two are
    ``A n - b`` over the element rows and the phase's net charge.

    **A species a surviving reaction names but the caller did not supply is refused**, not
    skipped: NeqSim reads it off the live phase, which holds every product the reaction
    machinery added, and a shorter list here would report a smaller residual than the one
    that exists.
    """
    answer_magnitudes = [float(value.to("mol").magnitude) for value in answer]
    worst_reaction = 0.0
    for reaction in _tables.reactions(source):
        if not reaction.use_reaction:
            continue
        species = _tables.stoichiometry(reaction.name)
        if not species:
            continue
        # `removeJunkReactions`: a reaction is kept only when every *reactant* the fluid could
        # hold is one it holds. The reactants are the negative coefficients.
        if not all(name in components for name, coefficient in species if coefficient < 0.0):
            continue
        ln_k = equilibrium_constant(source, reaction.name, quantity(temperature, "K")).ln_k
        quotient = 0.0
        for name, coefficient in species:
            if name not in components:
                raise InvalidInputError(
                    "components",
                    f"the reaction `{reaction.name}` names {name}, which this fluid does not "
                    f"carry, so its log residual cannot be evaluated here",
                )
            index = components.index(name)
            x = max(answer_magnitudes[index] / phase_moles, MIN_MOLES)
            quotient += coefficient * (math.log(x) + log_activity[index])
        worst_reaction = max(worst_reaction, abs(quotient - ln_k))

    # The charge row is the phase's net charge, and the correction `b` carries is applied to
    # the reactive set only - so the phase's own charge moves by what the set gained.
    charge_row = len(a_matrix) - 1
    net_charge = phase_charge + sum(
        a_matrix[charge_row][i] * (answer_magnitudes[i] - input_moles[i])
        for i in range(len(answer_magnitudes))
    )

    worst_element = 0.0
    for element in range(charge_row):
        amount = sum(
            a_matrix[element][i] * answer_magnitudes[i] for i in range(len(answer_magnitudes))
        )
        worst_element = max(worst_element, abs(amount - b[element]))

    return worst_reaction, net_charge, worst_element


def _element_matrix(components: list[str]) -> list[list[float]]:
    """The element matrix of a reaction set.

    Rows are the elements the substances carry, **sorted by symbol**, with the
    electroneutrality row last. NeqSim collects them into a ``HashSet`` and iterates it,
    so its row order is Java's hash order; the two matrices are the same up to a row
    permutation.
    """
    elements: list[str] = []
    composition: list[tuple[tuple[str, float], ...]] = []
    charge: list[float] = []

    for component in components:
        counts = _tables.element_composition(component)
        if counts is None:
            raise InvalidInputError(
                "components",
                f"`{component}` has no row in the element table, so it has no formula to "
                f"put in the balance. NeqSim's `Element` is empty for it and the row is "
                f"built from whatever components do have one; two fluids can therefore "
                f"carry different matrices under the same name",
            )
        for element, _ in counts:
            if element not in elements:
                elements.append(element)
        composition.append(counts)
        z = _tables.ionic_charge(component)
        if z is None:
            raise InvalidInputError(
                "components",
                f"`{component}` has no row in the component databank, so its charge is "
                f"unknown and cannot be defaulted to zero",
            )
        charge.append(z)

    elements.sort()
    matrix = [[0.0] * len(components) for _ in range(len(elements) + 1)]
    for column, counts in enumerate(composition):
        for element, count in counts:
            matrix[elements.index(element)][column] = count
    for column, z in enumerate(charge):
        matrix[len(elements)][column] = z
    return matrix


def _element_amounts(a_matrix: list[list[float]], moles: list[float]) -> list[float]:
    """``A n``: the element amounts the composition carries."""
    return [
        sum(coefficient * amount for coefficient, amount in zip(row, moles, strict=True))
        for row in a_matrix
    ]


def _spec() -> dict[str, object]:
    """The generated spec, imported here so the module has no import-time cycle."""
    from azoth._models_gen import model

    return model("reactions.reactive_phase_equilibrium")


def _solvent_mask(components: list[str]) -> tuple[float, ...]:
    """Which components keep the mole-fraction form on the solute-molality basis.

    The reference state is the component databank's, for the reason the ionic charge is: a
    caller that stated the mask could state a different one from the table's, and the branch
    it selects is a property of the substance.
    """
    return tuple(
        1.0 if _components.entry(name).reference_state == _components.SOLVENT else 0.0
        for name in components
    )


def _solvent_weight(components: list[str], moles: list[float], mask: tuple[float, ...]) -> float:
    """``sum(n_j M_j)`` over the masked components, in kg."""
    weight = 0.0
    for name, amount, marked in zip(components, moles, mask, strict=True):
        if marked < 0.5:
            continue
        record = _components.entry(name)
        if record.molar_mass is not None:
            weight += max(amount, 0.0) * record.molar_mass.to("kg/mol").magnitude
    return weight
