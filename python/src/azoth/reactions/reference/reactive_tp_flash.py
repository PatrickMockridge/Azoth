"""``reactions.reactive_tp_flash`` - the reactive flash as a model.

The pure-Python implementation, the twin of ``crates/azoth-reactions/src/
reactive_tp_flash.rs``. Everything the driver does is in
:mod:`azoth.reactions.reference._reactive_flash`, which takes its phase model and its
reactive solve as closures; this module is what wires those closures to the databank and
a cubic, which is why the reaction namespace reaches into ``azoth.eos`` here and nowhere
else.

**SRK, and every phase takes the vapour root** - and the second is a measurement rather
than a shortcut: NeqSim clones phase 0, a ``gas`` phase, for the phase its VLE
initialisation adds, so both phases of a split are gas-typed there. Rooting the water-rich
phase liquid collapses the 300 K split to a single phase; the vapour root reproduces it.

**Nothing here names the phases.** NeqSim's phase *types* come from ``init(1)``'s
assignment and are its system's bookkeeping, so the rows come back in the driver's order
with their mole numbers and fractions and nothing said about what they are.
"""

from __future__ import annotations

from azoth.core.errors import InvalidInputError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import ReactiveTpFlashResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning
from azoth.eos.components import entry, from_names
from azoth.eos.reference._mixture_state import phase_state, reduced_parameters
from azoth.reactions.reference import _reactive_flash as driver
from azoth.reactions.reference import _tables
from azoth.reactions.reference._rand_solver import (
    PhaseFeed,
    ThermoData,
    solve_single_phase,
    standard_potentials,
)
from azoth.reactions.reference._reactive_stability import CriticalConstants

__all__ = ["reactive_tp_flash"]

#: The unit NeqSim's ``ln(P/P_ref)`` term is taken against, and the one its
#: ``getPressure()`` reports: the potentials are reduced at bara and the cubic at pascals.
BARA = 1.0e5


def _spec() -> dict[str, object]:
    """The generated spec, imported here so the module has no import-time cycle."""
    from azoth._models_gen import model

    return model("reactions.reactive_tp_flash")


def _si(spec: dict[str, object], name: str, value: float | Q) -> float:
    """A declared input as its SI magnitude, quantity or not."""
    from azoth.reactions.reference.chemical_equilibrium import _si as shared

    return shared(spec, name, value)


def _element_matrix(names: list[str]) -> list[list[float]]:
    """The element matrix, with the rows in **first-encountered** order.

    ``FormulaMatrix.build`` walks the components and appends an element the first time it
    sees it, so the WGS fluid's rows are ``C O H`` and not a sorted ``C H O``. The matrix
    is internal - it is what ``b`` is built from and what the solve is handed - so the
    order is only the two kernels' agreement with each other.
    """
    element_names: list[str] = []
    per_component: list[tuple[tuple[str, float], ...]] = []
    for name in names:
        found = _tables.element_composition(name)
        if found is None:
            raise InvalidInputError(
                "components",
                f"`{name}` has no element row, so it cannot take part in an element "
                f"balance. NeqSim builds the row from whatever components do have one; "
                f"this refuses rather than answering a smaller balance",
            )
        per_component.append(found)
        for element, _count in found:
            if element not in element_names:
                element_names.append(element)

    matrix = [[0.0] * len(names) for _ in element_names]
    for column, counts in enumerate(per_component):
        for element, count in counts:
            matrix[element_names.index(element)][column] = count
    return matrix


def _formation(names: list[str]) -> list[ThermoData]:
    """The formation columns and the heat-capacity polynomial the potentials are built on.

    The three columns come from the reaction databank and the polynomial from the
    component row, which is where ``computeG0`` reads them in NeqSim.
    """
    out: list[ThermoData] = []
    for name in names:
        formation = _tables.formation_properties(name)
        if formation is None:
            raise InvalidInputError(
                "components",
                f"`{name}` has no formation row, so its standard potential is unknown. "
                f"NeqSim's `getIdealGasEnthalpyOfFormation` returns zero there and the "
                f"potential comes out as `ln(P)` alone; this refuses rather than answering "
                f"that",
            )
        coefficients = entry(name).cp
        if coefficients is None:
            raise InvalidInputError(
                "components",
                f"`{name}` has no heat-capacity polynomial, so its potential cannot be "
                f"corrected from 298.15 K to the state's temperature",
            )
        out.append(
            ThermoData(
                enthalpy_of_formation=formation.enthalpy_of_formation,
                absolute_entropy=formation.absolute_entropy,
                gibbs_energy_of_formation=formation.gibbs_energy_of_formation,
                cp=tuple(coefficients),  # type: ignore[arg-type]
            )
        )
    return out


def reactive_tp_flash(
    components: list[str],
    T: Q,
    P: Q,
    moles: list[Q],
    max_phases: float,
    cubic: str | None = None,
) -> ReactiveTpFlashResult:
    """Simultaneous chemical and phase equilibrium at fixed temperature and pressure.

    Args:
        components: the substances the fluid is made of, by name. **A charged one is
            refused**: the RAND solve's ionic branch is not ported.
        T: absolute temperature.
        P: absolute pressure.
        moles: the overall component amounts, which is what the driver's
            ``getOverallMoles`` reads and what its frozen element inventory is built
            from. Not a composition and not renormalised.
        max_phases: the driver's effective phase ceiling. Two is the captured states' own,
            and **one changes the answer rather than the format**: a ceiling of one skips
            the VLE initialisation and collapses the phase list.
        cubic: the cubic the fluid is flashed with, ``"srk"`` or ``"pr"``. Omitted
            means ``"srk"``, which is what the class flashes and what every captured
            state is; ``"pr"`` is the cubic the process layer's streams carry.

    Returns:
        The phases, in the driver's own order and with no type promised, beside the
        driver's own numbers. **The composition is the answer and the split is not**:
        where two phases converge to the same composition the Gibbs energy is flat along
        the direction that trades moles between them.

    Raises:
        InvalidInputError: for a charged component, a shape disagreement, or a component
            the databank does not carry.
        OutOfRangeError: for a temperature, pressure or ceiling the spec bounds.

    Example:
        >>> import azoth
        >>> q = azoth.ureg.Quantity
        >>> r = reactive_tp_flash(
        ...     ["CO", "water", "CO2", "hydrogen"],
        ...     q(600.0, "K"),
        ...     q(1.0, "bar"),
        ...     [q(0.25, "mol")] * 4,
        ...     2.0,
        ... )
        >>> round(r.phase_count, 0)
        2
    """
    spec = _spec()
    checks = checks_for(spec)
    warnings: list[Warning] = []
    temperature = input_to_si(spec, "T", T)
    pressure = input_to_si(spec, "P", P)
    ceiling = int(max_phases)

    apply_checks(
        checks.on_input,
        {"T": temperature, "P": pressure, "max_phases": float(ceiling)}.get,
        warnings,
    )

    if not components:
        raise InvalidInputError("components", "a fluid with no component has no state")
    if len(components) != len(moles):
        raise InvalidInputError(
            "moles",
            f"{len(moles)} entry(ies) against {len(components)} component(s)",
        )
    if ceiling < 1:
        raise InvalidInputError("max_phases", f"a phase ceiling of {ceiling} is not a ceiling")

    amounts = [_si(spec, "moles", value) for value in moles]
    total_moles = sum(amounts)
    if total_moles <= 0.0:
        raise InvalidInputError("moles", "the feed holds no moles")

    # The ionic branch is refused rather than approximated, and the charge comes from the
    # component databank rather than being defaulted to zero.
    for name in components:
        charge = _tables.ionic_charge(name)
        if charge is None:
            raise InvalidInputError(
                "components",
                f"`{name}` has no row in the component databank, so its charge is unknown "
                f"and cannot be defaulted to zero",
            )
        if charge != 0.0:
            raise InvalidInputError(
                "components",
                f"`{name}` is charged, and the RAND solve's ionic branch is not ported",
            )

    fluid = from_names(list(components), eos="srk" if cubic is None else cubic)
    reduced = reduced_parameters(fluid, temperature, pressure)
    constants = [
        CriticalConstants(
            tc=component.Tc.to("K").magnitude,
            pc=component.Pc.to("bar").magnitude,
            omega=component.omega,
        )
        for component in fluid.components
    ]

    matrix = _element_matrix(list(components))
    b = [sum(row[i] * amounts[i] for i in range(len(amounts))) for row in matrix]

    # NeqSim's `getPressure()` is in bara, which is the unit its `ln(P/P_ref)` term is taken
    # against; the model's pressure is in pascals and is converted here and nowhere else.
    g0 = standard_potentials(_formation(list(components)), temperature, pressure / BARA)
    charges = [0.0] * len(components)

    def ln_phi_at(x: list[float]) -> list[float]:
        """The vapour root at a composition, which is the root every phase takes."""
        return list(phase_state(reduced, fluid.kij, x, liquid=False).ln_phi)

    def phase_ln_phi(_index: int, x: list[float]) -> list[float]:
        return ln_phi_at(x)

    def ce(x: list[float]) -> list[float]:
        """The single-phase reactive solve the analysis brings its reference to.

        It answers with its **input** where the solve did not converge, which is what both
        of the class's callers do with a failed solve. The element inventory is the
        driver's frozen one; NeqSim's stability step recomputes it from a clone of phase 0
        instead, and the two agree on the first call and differ on a recheck.
        """
        solved = solve_single_phase(matrix, g0, b, x, ln_phi_at)
        if not solved.converged:
            return list(x)
        total = sum(solved.moles)
        return [value / total for value in solved.moles]

    outcome = driver.run(
        driver.DriverState(
            feed_moles=amounts,
            a_matrix=matrix,
            g0=g0,
            b=b,
            total_moles=total_moles,
            constants=constants,
            charges=charges,
            temperature=temperature,
            pressure=pressure / BARA,
            max_phases=ceiling,
            phases=[
                PhaseFeed(fractions=[value / total_moles for value in amounts], beta=1.0),
                PhaseFeed(fractions=[value / total_moles for value in amounts], beta=1.0),
            ],
        ),
        phase_ln_phi,
        ln_phi_at,
        ce,
    )

    solution = outcome.solution
    # The moles per phase are the solve's own `n[j][i]`, and the `NR = 0` branch has none:
    # its phases carry a composition at `beta = 1.0` and the amounts follow from it.
    rows: list[tuple[Q, ...]] = []
    for index, phase in enumerate(outcome.phases):
        held = (
            solution.phase_moles[index]
            if solution is not None and index < len(solution.phase_moles)
            else [x * total_moles * phase.beta for x in phase.fractions]
        )
        rows.append(tuple(from_si(value, "mol") for value in held))

    # **The label is the driver's own index, where it built a pair, and nothing where it did
    # not.** The `NR = 0` delegation and the VLE initialisation index their phases by vapour and
    # liquid; a reacting solve's phases come out of the loop as two copies of the feed, so naming
    # one of them the vapour would be picking a row.
    phase_type = (
        tuple(
            "vapour" if index == outcome.gas_index else "liquid"
            for index in range(len(outcome.phases))
        )
        if outcome.gas_index is not None
        else ()
    )

    return ReactiveTpFlashResult(
        phase_count=len(outcome.phases),
        phase_type=phase_type,
        phase_moles=tuple(rows),
        phase_fraction=tuple(phase.beta for phase in outcome.phases),
        converged=outcome.converged,
        total_iterations=outcome.total_iterations,
        equilibrium_total_moles=from_si(outcome.equilibrium_total_moles, "mol"),
        gibbs_energy=outcome.gibbs_energy,
        residual=solution.final_residual if solution is not None else 0.0,
        element_residual=solution.element_residual if solution is not None else 0.0,
        warnings=tuple(warnings),
    )
