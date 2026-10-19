"""``reactions.reactive_co2_brine_phase_equilibrium`` - the phase/reaction coupling.

Spec: ``specs/models/reactions/reactive_co2_brine_phase_equilibrium.toml``

The Python twin of ``crates/azoth-reactions/src/reactive_co2_brine_phase_equilibrium.rs``,
written to mirror it statement for statement.

# The loop, and why it is a loop

``ReactiveCO2BrinePhaseEquilibrium.run`` alternates two solves. A phase selection holds the
species inventory fixed and confines the ions to the aqueous phase; chemistry is then solved
on a *copy* of that aqueous phase, and only its conservative species changes are written back
to the full inventory. Feeding a generic coupled multiphase iteration with stale overall
species amounts is what the class exists to avoid. A result requires molecular fugacity
equality (the phase solve's own acceptance), reaction equilibrium, element conservation and
electroneutrality at once.

# The residual is computed on the flashed state

`getMaximumAbsoluteReactionLogResidual` is **lazy**: it walks the reaction list against the
system's current reactive phase every time it is read. So the accept test is evaluated on the
phase split the loop has just produced, not on the isolated chemistry solve's own residual -
which is why the accepted state's residual is `2e-7` and not the solver's `1e-12`.

# The activity coefficient is the unsymmetric one

`ChemicalReaction.getLogReactionActivity` adds
`getLogActivityCoefficient(k, water)`: the phase's `ln phi_k` less a reference. For a
solvent the reference is `getLogPureComponentFugacity` - the component alone; for a solute it
is `getLogInfiniteDiluteFugacity(k, water)` - a *binary* phase of component `k` at `1e-11`
against water. Both are NeqSim's `Phase.initRefPhases`, and omitting them moves the
water-gas-shift speciation by a third: measured, `H3O+` comes out `1.0e-4` against the
capture's `1.5e-4`.
"""

from __future__ import annotations

import math
from collections.abc import Sequence

from azoth.core.errors import InvalidInputError, SolverNotConvergedError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import (
    Co2BrinePhaseEquilibriumResult,
    ReactiveCo2BrinePhaseEquilibriumResult,
)
from azoth.core.units import Q, from_si, input_to_si, quantity
from azoth.core.warnings import Warning
from azoth.eos import components as _components
from azoth.eos.reference.co2_brine_phase_equilibrium import co2_brine_phase_equilibrium
from azoth.eos.reference.electrolyte_cpa_statoil_phase import electrolyte_cpa_statoil_phase
from azoth.reactions.reference import _lp_seed, _tables
from azoth.reactions.reference.chemical_equilibrium import chemical_equilibrium
from azoth.reactions.reference.equilibrium_constant import GAS_CONSTANT
from azoth.reactions.reference.reactive_phase_equilibrium import (
    _certify,
    _element_amounts,
    _element_matrix,
    _solvent_mask,
    _solvent_weight,
)
from azoth.reactions.reference.reference_potentials import reference_potentials

CALC_ID = "reactions.reactive_co2_brine_phase_equilibrium"

#: The class's `MAXIMUM_ITERATIONS`, and the two gates its accept test applies.
MAXIMUM_ITERATIONS = 30
REACTION_TOLERANCE = 2.0e-6
CHARGE_TOLERANCE_MOLES = 1.0e-8

#: The element-conservation tolerance `verifyElements` uses: `max(1e-9, |expected| * 1e-10)`.
ELEMENT_TOLERANCE_FLOOR = 1.0e-9
ELEMENT_TOLERANCE_RELATIVE = 1.0e-10

#: The chemistry solve's own parameters, from `solveChemEq`'s defaults.
CHEMISTRY_MAX_ITERATIONS = 100
CHEMISTRY_TOLERANCE = 1.0e-8

#: The finite-difference the reference phase dilutes a solute to: `1e-10 mol` against
#: `10 mol` of water, which is `Phase.initRefPhases`' own binary.
REFERENCE_SOLUTE_MOLES = 1.0e-10
REFERENCE_WATER_MOLES = 10.0


def _reference_log_fugacity(name: str, temperature: Q, pressure: Q) -> float:
    """``getLogPureComponentFugacity`` or ``getLogInfiniteDiluteFugacity`` for one name."""
    entry = _components.entry(name)
    if entry.reference_state == _components.SOLVENT:
        state = electrolyte_cpa_statoil_phase([name], temperature, pressure, [1.0], "liquid")
        return state.ln_phi[0]
    fraction = REFERENCE_SOLUTE_MOLES / (REFERENCE_SOLUTE_MOLES + REFERENCE_WATER_MOLES)
    state = electrolyte_cpa_statoil_phase(
        [name, "water"], temperature, pressure, [fraction, 1.0 - fraction], "liquid"
    )
    return state.ln_phi[0]


def _log_activity(
    components: Sequence[str],
    actual: Sequence[float],
    temperature: Q,
    pressure: Q,
) -> list[float]:
    """`ln gamma_k` per component: the phase's own less its unsymmetric reference."""
    return [
        float(actual[index]) - _reference_log_fugacity(name, temperature, pressure)
        for index, name in enumerate(components)
    ]


def reactive_co2_brine_phase_equilibrium(
    components: list[str],
    moles: list[Q],
    T: Q,
    P: Q,
) -> ReactiveCo2BrinePhaseEquilibriumResult:
    """The simultaneous phase and chemical equilibrium of a water-rich CO2/brine feed.

    Args:
        components: the fluid's substances, by name: CO2, water and its ions.
        moles: each component's overall amount.
        T: absolute temperature.
        P: absolute pressure.

    Returns:
        The phase inventory, the accepted root, the residual the answer was accepted on, the
        loop's pass count and the minimum tangent-plane distance.

    Raises:
        InvalidInputError: where the shapes disagree, a component carries no element row, or a
            species inventory the loop produces is not one.
        SolverNotConvergedError: if the coupling reaches its pass cap, a chemistry solve does
            not converge, or the element balance moves.
        Whatever the phase selection raises for a feed outside its scope.

    See :func:`azoth.reactions.reactive_co2_brine_phase_equilibrium`.
    """
    spec = _spec()
    checks = checks_for(spec)
    warnings: list[Warning] = []

    temperature = input_to_si(spec, "T", T)
    pressure = input_to_si(spec, "P", P)
    apply_checks(checks.on_input, {"T": temperature, "P": pressure}.get, warnings)

    feed = [input_to_si(spec, "moles", value) for value in moles]
    if len(components) != len(feed):
        raise InvalidInputError(
            "moles",
            f"{len(components)} component(s) and {len(feed)} mole number(s); the two must match",
        )
    if any(not math.isfinite(value) or value < 0.0 for value in feed):
        raise InvalidInputError("moles", "every mole number must be finite and not negative")

    a_matrix = _element_matrix(list(components))
    expected = _element_amounts(a_matrix, feed)
    mask = _solvent_mask(list(components))
    charges = [_charge(name) for name in components]

    work = list(feed)
    flash = None
    reaction_residual = math.nan
    for iteration in range(1, MAXIMUM_ITERATIONS + 1):
        flash = co2_brine_phase_equilibrium(
            list(components),
            quantity(temperature, "K"),
            quantity(pressure, "Pa"),
            [quantity(value, "mol") for value in work],
        )
        aqueous = 1 if flash.phase_count == 2 else 0
        beta = flash.phase_fractions[aqueous]
        total = sum(work)
        before = [value * beta * total for value in flash.x[aqueous]]
        b = _element_amounts(a_matrix, before)
        b[-1] = 0.0
        log_activity = _log_activity(
            components, flash.ln_phi[aqueous], quantity(temperature, "K"), quantity(pressure, "Pa")
        )
        phase_charge = sum(z * n for z, n in zip(charges, before, strict=True))
        reaction_residual, _, _ = _certify(
            list(components),
            "standard",
            temperature,
            a_matrix,
            b,
            [quantity(value, "mol") for value in before],
            before,
            log_activity,
            phase_charge,
            sum(before),
        )

        if math.isfinite(reaction_residual) and reaction_residual <= REACTION_TOLERANCE:
            charge = sum(z * n for z, n in zip(charges, work, strict=True))
            if not math.isfinite(charge) or abs(charge) > CHARGE_TOLERANCE_MOLES:
                raise SolverNotConvergedError(iteration, charge, CHARGE_TOLERANCE_MOLES)
            return _result(
                flash,
                iteration,
                reaction_residual,
                charge,
                _element_balance(a_matrix, expected, work),
                warnings,
            )

        # **The chemistry runs on an isolated copy of the aqueous phase.** Only its
        # conservative species changes are written back; the rest of the inventory keeps the
        # amounts the phase selection left it.
        potentials = reference_potentials(list(components), "standard", quantity(temperature, "K"))
        warnings.extend(potentials.warnings)
        reduced = [
            float(potential.to("J/mol").magnitude) / (GAS_CONSTANT * temperature)
            for potential in potentials.potentials
        ]
        estimate = _lp_seed.initial_estimate(a_matrix, b, reduced)
        start = (
            list(before)
            if estimate is None
            else [max(value, _lp_seed.MIN_MOLES) for value in estimate]
        )
        solved = chemical_equilibrium(
            a_matrix,
            b,
            True,
            start,
            reduced,
            log_activity,
            quantity(temperature, "K"),
            CHEMISTRY_MAX_ITERATIONS,
            CHEMISTRY_TOLERANCE,
            "mole_fraction",
            quantity(_solvent_weight(list(components), before, mask), "kg"),
            mask,
            quantity(sum(before), "mol"),
        )
        warnings.extend(solved.warnings)
        # **A non-converged chemistry solve contributes its last iterate.** NeqSim's
        # `solveChemEq` throws here; this solver fails from the floored post-flash
        # composition where NeqSim's converges, and the loop's own gates - the element
        # check below and the flashed-state residual above - are what decide the answer.
        # The spec states the divergence.
        after = [float(value.to("mol").magnitude) for value in solved.moles]
        updated: list[float] = []
        for index in range(len(work)):
            value = work[index] + (after[index] - before[index])
            if not math.isfinite(value) or value < -1.0e-12:
                raise InvalidInputError(
                    "moles",
                    f"the coupling produced {value} mol of component {index}, which is not an "
                    f"inventory",
                )
            updated.append(max(0.0, value))
        residual = _element_balance(a_matrix, expected, updated)
        if not math.isfinite(residual):
            raise SolverNotConvergedError(iteration, residual, ELEMENT_TOLERANCE_FLOOR)
        actual_elements = _element_amounts(a_matrix, updated)
        for expected_value, actual in zip(expected, actual_elements, strict=True):
            tolerance = max(
                ELEMENT_TOLERANCE_FLOOR, abs(expected_value) * ELEMENT_TOLERANCE_RELATIVE
            )
            if abs(actual - expected_value) > tolerance:
                raise SolverNotConvergedError(iteration, actual - expected_value, tolerance)
        work = updated

    raise SolverNotConvergedError(MAXIMUM_ITERATIONS, reaction_residual, REACTION_TOLERANCE)


def _charge(name: str) -> float:
    """The component's ionic charge, refused rather than defaulted when absent."""
    value = _tables.ionic_charge(name)
    if value is None:
        raise InvalidInputError(
            "components", f"`{name}` has no row in the component databank, so its charge is unknown"
        )
    return value


def _element_balance(
    a_matrix: list[list[float]], expected: list[float], moles: list[float]
) -> float:
    """The worst `|A n - b|` over the element rows, excluding the charge row."""
    actual = _element_amounts(a_matrix, moles)
    return max(
        (abs(actual[index] - expected[index]) for index in range(len(expected) - 1)),
        default=0.0,
    )


def _result(
    flash: Co2BrinePhaseEquilibriumResult,
    iterations: int,
    reaction_residual: float,
    charge: float,
    element_residual: float,
    warnings: list[Warning],
) -> ReactiveCo2BrinePhaseEquilibriumResult:
    """The accepted state as the record."""
    phase_fractions = tuple(float(value) for value in flash.phase_fractions)
    x = tuple(tuple(float(value) for value in row) for row in flash.x)
    phase_type: tuple[str, ...]
    if len(phase_fractions) == 2:
        phase_type = (str(flash.accepted_root), "aqueous")
    else:
        phase_type = ("aqueous",)
    return ReactiveCo2BrinePhaseEquilibriumResult(
        phase_count=len(phase_fractions),
        phase_type=phase_type,
        phase_fractions=phase_fractions,
        phase_composition=x,
        accepted_root=str(flash.accepted_root),
        minimum_trial_distance=float(flash.minimum_trial_distance),
        iterations=iterations,
        reaction_log_residual=reaction_residual,
        charge_moles=from_si(charge, "mol"),
        element_balance_residual=from_si(element_residual, "mol"),
        warnings=tuple(warnings),
    )


def _spec() -> dict[str, object]:
    """The generated spec, imported here so the module has no import-time cycle."""
    from azoth._models_gen import model

    return model(CALC_ID)
