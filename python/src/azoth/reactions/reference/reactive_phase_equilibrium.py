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

from azoth.core.errors import InvalidInputError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import ReactivePhaseEquilibriumResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning
from azoth.reactions.reference import _tables
from azoth.reactions.reference.chemical_equilibrium import chemical_equilibrium
from azoth.reactions.reference.equilibrium_constant import GAS_CONSTANT
from azoth.reactions.reference.reference_potentials import reference_potentials

#: The charge, in moles of elementary charge, below which the phase is treated as
#: neutral, from ``calcBVector``'s ``1e-10 * Math.max(1.0, phase moles)``.
CHARGE_NOISE_MOLE_FRACTION = 1e-10

#: The floor ``updateMoles`` raises every mole number to before writing it back, from
#: ``Math.max(newMoles[i], 1e-45)``.
MIN_WRITTEN_MOLES = 1e-45

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
) -> ReactivePhaseEquilibriumResult:
    """The reactive equilibrium composition of one phase.

    Args:
        components: the phase's reactive substances - not every substance in it.
        source: which of the three reaction tables the reference potentials come from.
        phase: the phase's type name as NeqSim spells it. ``aqueous``, ``liquid`` or
            ``oil`` takes the solve; anything else skips it.
        moles: their amounts in the phase, and the starting composition.
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
            warnings=tuple(warnings),
        )

    # The solver takes the potentials reduced, which is the form its own iteration adds a
    # logarithm to. NeqSim divides by `R T` at the same boundary.
    reduced = [
        float(potential.to("J/mol").magnitude) / (GAS_CONSTANT * temperature)
        for potential in chem_ref
    ]
    solved = chemical_equilibrium(
        a_matrix,
        b,
        whole_system,
        magnitudes,
        reduced,
        activity,
        T,
        max_iterations,
        tolerance,
    )
    warnings.extend(solved.warnings)

    return ReactivePhaseEquilibriumResult(
        skipped=False,
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
        converged=solved.converged,
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
