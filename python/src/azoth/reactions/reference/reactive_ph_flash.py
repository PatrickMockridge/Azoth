"""``reactions.reactive_ph_flash`` - the reactive PH flash as a model.

The pure-Python implementation, the twin of ``crates/azoth-reactions/src/
reactive_ph_flash.rs``. The loop is
:func:`azoth.reactions.reference._reactive_ph_flash.reactive_ph_flash`, and everything it
walks over is a piece the namespace already carries: the reactive TP flash for the inner step,
and ``azoth.eos``' molar enthalpy and heat capacity for the state.

**SRK, and the vapour root on every phase**, for the reason the TP model gives.

**The specification is thermochemical and the caller states it.** NeqSim's
``getEnthalpy()`` is a sensible enthalpy that excludes the formation enthalpies; the class's
constructor adds the inventory ``sum n_i dHf_i`` at whatever composition the system holds when
it is built - which its own test has flashed first. A model has no such system to read, so a
caller states the thermochemical number, and a caller holding a sensible enthalpy adds the
inventory at their own composition before calling.
"""

from __future__ import annotations

from azoth.core.errors import InvalidInputError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import ReactivePhFlashResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning
from azoth.eos import IdealGasModel, molar_enthalpy_entropy
from azoth.eos.components import entry, from_names
from azoth.reactions.reference import _tables
from azoth.reactions.reference._reactive_ph_flash import PhState
from azoth.reactions.reference._reactive_ph_flash import reactive_ph_flash as loop
from azoth.reactions.reference.reactive_tp_flash import reactive_tp_flash as tp_flash

__all__ = ["reactive_ph_flash"]

#: The unit NeqSim's ``getPressure()`` reports and its potentials are reduced at.
BARA = 1.0e5


def _spec() -> dict[str, object]:
    """The generated spec, imported here so the module has no import-time cycle."""
    from azoth._models_gen import model

    return model("reactions.reactive_ph_flash")


def _si(spec: dict[str, object], name: str, value: float | Q) -> float:
    """A declared input as its SI magnitude, quantity or not."""
    from azoth.reactions.reference.chemical_equilibrium import _si as shared

    return shared(spec, name, value)


def reactive_ph_flash(
    components: list[str],
    T: Q,
    P: Q,
    moles: list[Q],
    enthalpy: Q,
    max_phases: float,
) -> ReactivePhFlashResult:
    """The temperature at which a reactive fluid's enthalpy matches a specification.

    Args:
        components: the substances the fluid is made of, by name. **A charged one is
            refused**, because the inner flash's ionic branch is not ported.
        T: the temperature the search starts from. It is not a bound and not a guess at the
            answer: the loop is a secant, so a different start is a different path.
        P: absolute pressure, held fixed while the temperature moves.
        moles: the overall component amounts, the same input ``reactive_tp_flash`` takes.
        enthalpy: **the thermochemical specification** - the fluid's sensible enthalpy plus
            the formation inventory. A sensible enthalpy without the inventory is a different
            number and finds a different temperature.
        max_phases: the inner flash's phase ceiling.

    Returns:
        The temperature, whether the loop converged, and the passes it cost. **The counts are
        path quantities**: the same temperature is reached in a different number of steps by a
        different implementation.

    Raises:
        InvalidInputError: for a charged component, a shape disagreement, or a component the
            databank does not carry.
        OutOfRangeError: for a temperature, pressure or ceiling the spec bounds.

    Example:
        >>> import azoth
        >>> q = azoth.ureg.Quantity
        >>> r = reactive_ph_flash(
        ...     ["CO", "water", "CO2", "hydrogen"],
        ...     q(500.0, "K"),
        ...     q(1.0, "bar"),
        ...     [q(0.25, "mol")] * 4,
        ...     q(-182008.71, "J"),
        ...     2.0,
        ... )
        >>> round(float(r.temperature.magnitude), 3)
        600.0
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
            "moles", f"{len(moles)} entry(ies) against {len(components)} component(s)"
        )

    # The amounts are the inner flash's to read - it converts them itself - but the shape is
    # checked here, where the caller's list is, rather than one call deeper.
    if not moles:
        raise InvalidInputError("moles", "the feed holds no component")
    specified = input_to_si(spec, "enthalpy", enthalpy)

    fluid = from_names(list(components), eos="srk")
    coefficients = []
    formation = []
    for name in components:
        cp = entry(name).cp
        if cp is None:
            raise InvalidInputError(
                "components",
                f"`{name}` has no heat-capacity polynomial, so its enthalpy cannot be "
                f"evaluated at a trial temperature",
            )
        coefficients.append(cp)
        properties = _tables.formation_properties(name)
        if properties is None:
            raise InvalidInputError(
                "components",
                f"`{name}` has no formation row, so its thermochemical enthalpy is unknown",
            )
        formation.append(properties.enthalpy_of_formation)

    ideal = IdealGasModel(
        cp_a=tuple(c[0] for c in coefficients),
        cp_b=tuple(c[1] for c in coefficients),
        cp_c=tuple(c[2] for c in coefficients),
        cp_d=tuple(c[3] for c in coefficients),
        cp_e=tuple(c[4] for c in coefficients),
    )

    def inner(trial: float) -> PhState:
        """The inner flash at a trial temperature: its passes, its thermochemical enthalpy
        and its heat capacity."""
        outcome = tp_flash(list(components), _quantity(trial, "K"), P, moles, float(ceiling))
        sensible = 0.0
        heat_capacity = 0.0
        inventory = 0.0
        for row in outcome.phase_moles:
            held = sum(float(value.to("mol").magnitude) for value in row)
            if held <= 0.0:
                continue
            composition = [float(value.to("mol").magnitude) / held for value in row]
            z = _compressibility(fluid, trial, pressure, composition)
            state = molar_enthalpy_entropy(
                fluid, ideal, _quantity(trial, "K"), _quantity(pressure, "Pa"), composition, z
            )
            sensible += held * float(state.h.to("J/mol").magnitude)
            heat_capacity += held * float(state.cp.to("J/(mol*K)").magnitude)
            for index, value in enumerate(row):
                inventory += float(value.to("mol").magnitude) * formation[index]
        return PhState(
            iterations=outcome.total_iterations,
            thermochemical_enthalpy=sensible + inventory,
            cp=heat_capacity,
        )

    outcome = loop(temperature, specified, inner)
    return ReactivePhFlashResult(
        temperature=from_si(outcome.temperature, "K"),
        converged=outcome.converged,
        outer_iterations=outcome.outer_iterations,
        total_inner_iterations=outcome.total_inner_iterations,
        warnings=tuple(warnings),
    )


def _quantity(magnitude: float, unit: str) -> Q:
    """A pint quantity, built without importing the registry into this module's namespace."""
    import azoth

    return azoth.ureg.Quantity(magnitude, unit)


def _compressibility(fluid: object, temperature: float, pressure: float, z: list[float]) -> float:
    """The cubic's vapour root at a composition, which is the root every phase takes."""
    from azoth.eos.reference._mixture_state import phase_state, reduced_parameters

    reduced = reduced_parameters(fluid, temperature, pressure)  # type: ignore[arg-type]
    return float(phase_state(reduced, fluid.kij, z, liquid=False).z)  # type: ignore[attr-defined]
