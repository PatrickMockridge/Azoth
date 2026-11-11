"""``process.separator_capacity_assessment`` - the separator's detached capacity report.

Spec: ``specs/models/process/separator_capacity_assessment.toml``

The Python twin of ``crates/azoth-process/src/models/separator_capacity_assessment.rs``
over ``separator.rs``'s ``_route`` and ``_sizing_states``.

# Fail closed, not refuse

``SeparatorCapacityAssessment.from`` answers ``UNASSESSED`` for a blank provenance, an
orientation that is neither ``horizontal`` nor ``vertical``, a non-positive or non-finite
diameter or design K, no gas phase, a design liquid level outside ``[0, 1)``, and non-finite
evidence. Those are statuses, not exceptions: the snapshot's purpose is to report that the
evidence is incomplete. ``process.separator`` refuses several of the same inputs at the
boundary, and the two behaviours are deliberate - one computes, the other reports.

# One arithmetic, two ids

The flash and the sizing family are :mod:`azoth.process.reference.separator`'s own
functions, called rather than re-derived, so this id cannot disagree with the capacity
surface on a field they share. Its checks run in the class's order, so a state failing two of
them carries the earlier diagnostic.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from azoth.core.errors import InvalidInputError
from azoth.core.result import SeparatorCapacityAssessmentResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.process.reference.separator import _route, _sizing_states

MODEL_ID = "process.separator_capacity_assessment"

#: The class's fixed diagnostics, verbatim.
PROVENANCE_REQUIRED = "Design-basis provenance is required"
BAD_ORIENTATION = "Separator orientation must be horizontal or vertical"
GEOMETRY_REQUIRED = "Finite geometry, design K-factor, and a calculated gas phase are required"
BAD_LEVEL = "Design liquid level fraction must be in [0, 1)"
BAD_EVIDENCE = "Separator capacity calculation returned non-finite evidence"
FEASIBLE = "Gas capacity is within the evidenced design limit"
INFEASIBLE = "Gas capacity exceeds the evidenced design limit"


@dataclass(frozen=True, slots=True)
class AssessmentStates:
    """The snapshot, in the order the class builds it.

    Every numeric field is ``None`` on ``unassessed``: the class holds ``NaN`` there, and the
    port does not publish one.
    """

    status: str
    diagnostic: str
    design_basis_source: str
    governing_capacity_basis: str
    gas_area: float | None
    gas_density: float | None
    liquid_density: float | None
    maximum_gas_velocity: float | None
    current_gas_flow: float | None
    mechanical_design_maximum_gas_flow: float | None
    souders_brown_maximum_gas_flow: float | None
    maximum_gas_flow: float | None
    utilization: float | None


def separator_capacity_assessment(
    components: list[str],
    feed_n: Q,
    feed_z: list[float],
    feed_p: Q,
    feed_t: Q,
    pressure_drop: Q,
    gas_in_liquid: float,
    internal_diameter: Q,
    design_gas_load_factor: Q,
    orientation: str,
    design_liquid_level_fraction: float,
    design_basis_provenance: str,
    heat_input: Q | None = None,
    mechanical_design_max_gas_flow: Q | None = None,
) -> SeparatorCapacityAssessmentResult:
    """Assess a separator's gas capacity and report whether the evidence is complete.

    Args:
        components: the fluid's substances, by name.
        feed_n: the inlet's molar flow.
        feed_z: the inlet composition.
        feed_p: the inlet pressure.
        feed_t: the inlet temperature.
        pressure_drop: the pressure drop across the vessel.
        gas_in_liquid: the fraction of the vapour's moles carried into the liquid outlet.
        internal_diameter: the vessel's internal diameter.
        design_gas_load_factor: the Souders-Brown ``K`` the vessel is designed to.
        orientation: ``"horizontal"`` or ``"vertical"``; anything else is ``unassessed``.
        design_liquid_level_fraction: the horizontal gas area's derating, in ``[0, 1)``.
        design_basis_provenance: where the geometry and design ``K`` come from; blank is
            ``unassessed``.
        heat_input: heat added before the flash, or ``None`` for the feed's temperature.
        mechanical_design_max_gas_flow: the mechanical design's gas limit; absent or zero
            means no mechanical basis.

    Returns:
        The snapshot, with the numeric evidence absent on ``unassessed``.

    Example:
        >>> import azoth
        >>> q = azoth.ureg.Quantity
        >>> r = separator_capacity_assessment(
        ...     ["methane", "n-butane"], q(1.0, "mol/s"), [0.7, 0.3], q(20.0, "bar"),
        ...     q(300.0, "K"), q(0.0, "Pa"), 0.0, q(1.0, "m"), q(0.11, "m/s"),
        ...     "vertical", 0.8, "design basis A",
        ... )
        >>> r.status
        'feasible'
    """
    spec = _spec()
    n = input_to_si(spec, "feed_n", feed_n)
    # Converted for the unit check alone: the route takes the quantity itself.
    input_to_si(spec, "feed_t", feed_t)
    p = input_to_si(spec, "feed_p", feed_p)
    drop = input_to_si(spec, "pressure_drop", pressure_drop)
    diameter = input_to_si(spec, "internal_diameter", internal_diameter)
    load_factor = input_to_si(spec, "design_gas_load_factor", design_gas_load_factor)
    level = float(design_liquid_level_fraction)
    provenance = str(design_basis_provenance).strip()
    kind = str(orientation)
    mechanical = (
        None
        if mechanical_design_max_gas_flow is None
        else input_to_si(spec, "mechanical_design_max_gas_flow", mechanical_design_max_gas_flow)
    )
    duty = None if heat_input is None else input_to_si(spec, "heat_input", heat_input)

    # **The run precedes the assessment**, as it does in the class: a feed that cannot flash
    # is an error, not an `unassessed` snapshot. The checks below then run in the class's order,
    # so a state failing two carries the earlier diagnostic.
    states = _route(components, n, feed_z, p, feed_t, drop, gas_in_liquid, duty)
    outlet_t = states.temperature.to("K").magnitude

    if not provenance:
        return _result(_unassessed(PROVENANCE_REQUIRED))
    if kind not in ("horizontal", "vertical"):
        return _result(_unassessed(BAD_ORIENTATION))
    if not (math.isfinite(diameter) and diameter > 0.0):
        return _result(_unassessed(GEOMETRY_REQUIRED))
    if not (math.isfinite(load_factor) and load_factor > 0.0):
        return _result(_unassessed(GEOMETRY_REQUIRED))
    # The class checks the gas phase before the level, so a liquid-only feed with a level of
    # one answers the geometry sentence rather than the level's.
    if not _has_gas(components, outlet_t, states.p_out, list(feed_z)):
        return _result(_unassessed(GEOMETRY_REQUIRED))
    if not math.isfinite(level) or level < 0.0 or level >= 1.0:
        return _result(_unassessed(BAD_LEVEL))

    try:
        sizing = _sizing_states(
            components,
            outlet_t,
            states.p_out,
            list(feed_z),
            states.n_vapour,
            list(states.liquid_z),
            diameter,
            load_factor,
            kind,
            level,
        )
    except InvalidInputError:
        return _result(_unassessed(BAD_EVIDENCE))
    if sizing is None:  # pragma: no cover - the probe above already refused this state
        return _result(_unassessed(GEOMETRY_REQUIRED))

    if mechanical is not None and not math.isfinite(mechanical):
        raise InvalidInputError(
            "mechanical_design_max_gas_flow",
            f"a stated mechanical gas limit is a finite flow, and {mechanical} is not",
        )
    mechanical_flow = 0.0 if mechanical is None else mechanical
    mechanical_positive = math.isfinite(mechanical_flow) and mechanical_flow > 0.0
    source = "mechanicalDesign" if mechanical_positive else "equipment"
    governs = mechanical_positive and mechanical_flow <= sizing.max_flow
    maximum_flow = mechanical_flow if governs else sizing.max_flow
    governing = "mechanicalDesignGasFlow" if governs else "soudersBrown"
    if not (math.isfinite(maximum_flow) and maximum_flow > 0.0):
        return _result(_unassessed(BAD_EVIDENCE))
    utilization = sizing.volumetric_flow / maximum_flow
    if not math.isfinite(utilization) or utilization < 0.0:
        return _result(_unassessed(BAD_EVIDENCE))

    return _result(
        AssessmentStates(
            status="feasible" if utilization <= 1.0 else "infeasible",
            diagnostic=FEASIBLE if utilization <= 1.0 else INFEASIBLE,
            design_basis_source=source,
            governing_capacity_basis=governing,
            gas_area=sizing.gas_area,
            gas_density=sizing.gas_density,
            liquid_density=sizing.liquid_density,
            maximum_gas_velocity=sizing.max_velocity,
            current_gas_flow=sizing.volumetric_flow,
            mechanical_design_maximum_gas_flow=mechanical_flow,
            souders_brown_maximum_gas_flow=sizing.max_flow,
            maximum_gas_flow=maximum_flow,
            utilization=utilization,
        )
    )


def _has_gas(components: list[str], t: float, p: float, feed_z: list[float]) -> bool:
    """Whether the flashed system carries a gas phase - the class's ``hasPhaseType("gas")``."""
    from azoth.eos.reference.hydrate_inhibitor_wt import GAS
    from azoth.process.reference.capacity import _labelled_phases

    _mixture, labelled, _shares = _labelled_phases(components, t, p, feed_z, from_feed=True)
    return any(entry[0] == GAS for entry in labelled)


def _result(states: AssessmentStates) -> SeparatorCapacityAssessmentResult:
    """The snapshot as the result dataclass, with its dimensioned fields as quantities."""

    def quantity(value: float | None, unit: str) -> Q | None:
        return None if value is None else from_si(value, unit)

    return SeparatorCapacityAssessmentResult(
        status=states.status,
        diagnostic=states.diagnostic,
        design_basis_source=states.design_basis_source,
        governing_capacity_basis=states.governing_capacity_basis,
        gas_area=quantity(states.gas_area, "m**2"),
        gas_density=quantity(states.gas_density, "kg/m**3"),
        liquid_density=quantity(states.liquid_density, "kg/m**3"),
        maximum_gas_velocity=quantity(states.maximum_gas_velocity, "m/s"),
        current_gas_flow=quantity(states.current_gas_flow, "m**3/s"),
        mechanical_design_maximum_gas_flow=quantity(
            states.mechanical_design_maximum_gas_flow, "m**3/s"
        ),
        souders_brown_maximum_gas_flow=quantity(states.souders_brown_maximum_gas_flow, "m**3/s"),
        maximum_gas_flow=quantity(states.maximum_gas_flow, "m**3/s"),
        utilization=states.utilization,
        warnings=(),
    )


def _unassessed(diagnostic: str) -> AssessmentStates:
    """The class's ``unavailable`` snapshot: the status, the sentence, and no numbers."""
    return AssessmentStates(
        status="unassessed",
        diagnostic=diagnostic,
        design_basis_source="",
        governing_capacity_basis="",
        gas_area=None,
        gas_density=None,
        liquid_density=None,
        maximum_gas_velocity=None,
        current_gas_flow=None,
        mechanical_design_maximum_gas_flow=None,
        souders_brown_maximum_gas_flow=None,
        maximum_gas_flow=None,
        utilization=None,
    )


def _spec() -> dict[str, object]:
    """The generated spec, imported here so the module has no import-time cycle."""
    from azoth._models_gen import model

    return model(MODEL_ID)
