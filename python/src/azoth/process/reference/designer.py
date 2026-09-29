"""``ColumnInternalsDesigner.calculateTrayed``: the trayed internals tree, as a report.

The Python twin of ``crates/azoth-process/src/column/designer.rs``. Every step mirrors the Rust
line for line, and the two are compared case by case by ``python/tests/test_cross_impl.py``.

**The class the tray calculator is driven *by*.** It walks the column's trays, picks the
**controlling** one by the largest vapour mass flow, sizes a diameter off that tray, then runs one
calculator per tray at the diameter it resolved and sums them. Everything the class publishes
comes out of that second loop.

**Three of the class's own lines are reproduced rather than corrected**, each measured in
``validation/neqsim/captures/process_internals_designer.tsv``:

**The controlling tray's calculator is discarded.** It exists only to size the diameter from, and
it is built **without a relative volatility** while every tray in the summation loop carries one.

**The extraction's fallbacks are the class's.** ``getTrayProperties`` answers
``1.0``/``800.0``/``0.001``/``0.02``/``2.0`` when the tray's fluid is absent or single-phase, and
on a single-phase tray it takes **phase 0** as the vapour density and substitutes ``600.0`` for
the liquid's where that phase is a gas. **Which trays those are is read off the tray's own two
traffics and not off a rebuilt stream's flash**: the reboiler's vapour stream and its liquid
stream each flash to one phase here, so a kind comparison would call a two-phase tray
single-phase and substitute ``600.0`` for a liquid whose density is its own.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, Any

from azoth.core.units import quantity
from azoth.hydraulics.reference.tray_hydraulics import size_column_diameter, tray_hydraulics

if TYPE_CHECKING:
    # **The solved column's own record, imported for its type alone.** The column imports this
    # module's report, so a runtime import here would be a cycle - and every function below
    # duck-types the record rather than constructing one.
    from azoth.process.reference.distillation_column import _States

#: The tray spacing and the weir height the class's own constructor carries, one per geometry
#: input. They are `ColumnInternalsDesigner`'s initialisers and not correlated values.
DEFAULT_TRAY_SPACING_M = 0.6
DEFAULT_WEIR_HEIGHT_M = 0.05

#: The relative volatility the class's own default stands at, and the one the *sizing* branch runs
#: at because its calculator is never given the tray's.
DEFAULT_RELATIVE_VOLATILITY = 2.0
#: The liquid density a **single gas phase** substitutes for a missing liquid.
SINGLE_PHASE_GAS_LIQUID_DENSITY = 600.0
#: The liquid density the extraction leaves in place on a single non-gas phase.
SINGLE_PHASE_LIQUID_DENSITY_KG_PER_M3 = 800.0
#: The liquid viscosity the extraction falls back to where the phase answers nothing.
DEFAULT_LIQUID_VISCOSITY_PA_S = 0.001
#: The surface tension the extraction falls back to when the interface answers nothing, which on
#: every captured tray it does: the interphase route reads ``0.0``.
DESIGNER_SURFACE_TENSION_N_PER_M = 0.02
#: The average the class answers when its loop walked no tray at all.
EMPTY_AVERAGE_TRAY_EFFICIENCY = 0.65
#: The relative-volatility spread's own ceiling.
MAX_RELATIVE_VOLATILITY = 20.0
#: The column diameter that means "size it", ``ColumnInternalsDesigner``'s own ``-1.0``.
UNSIZED_COLUMN_DIAMETER_M = -1.0


@dataclass(frozen=True, slots=True)
class DesignerGeometry:
    """The designer's own geometry: its setters, and the class's constructor defaults."""

    #: ``sieve``, ``valve`` or ``bubble-cap``. The class's own constructor default is ``sieve``.
    internals_type: str = "sieve"
    #: The tray spacing, m. The class's own default is ``0.6``.
    tray_spacing_m: float = DEFAULT_TRAY_SPACING_M
    #: The weir height, m. The class's own default is ``0.05``.
    weir_height_m: float = DEFAULT_WEIR_HEIGHT_M
    #: The hole diameter, m. The class's field is millimetres.
    hole_diameter_m: float = 12.7e-3
    #: The hole area over the active area. The class's own default is ``0.1``.
    hole_area_fraction: float = 0.1
    #: The downcomer area over the total. The class's own default is ``0.1``.
    downcommer_area_fraction: float = 0.1
    #: The fraction of flood a tray is sized to. The class's own default is ``0.8``.
    design_flood_fraction: float = 0.8
    #: A stated diameter, which **replaces the sizing branch entirely** at or below zero.
    column_diameter_override_m: float = UNSIZED_COLUMN_DIAMETER_M


@dataclass(frozen=True, slots=True)
class DesignerTrayResult:
    """One tray's own results, which the designer keeps as a list of calculators."""

    #: The tray's load, in per cent of flood.
    percent_flood: float
    #: The tray's own total pressure drop, Pa.
    total_pressure_drop_pa: float
    #: The tray's efficiency.
    tray_efficiency: float
    #: Whether this tray is inside its own design window.
    design_ok: bool


@dataclass(frozen=True, slots=True)
class DesignerReport:
    """``ColumnInternalsDesigner.calculateTrayed``'s whole answer."""

    #: The diameter the sizing resolved, or the stated override where one was.
    required_diameter_m: float
    #: Which tray the diameter was sized from: the largest vapour mass flow, the first of equals.
    controlling_tray_index: int
    #: Every tray's own verdict, and-ed together.
    design_ok: bool
    #: The largest per-tray load, in per cent of flood.
    max_percent_flood: float
    #: The smallest load above zero, which is not the smallest: an exact zero is skipped.
    min_percent_flood: float
    #: The mean of the per-tray efficiencies, or the class's own ``0.65`` with no tray walked.
    average_tray_efficiency: float
    #: The trays' pressure drops, summed, Pa.
    total_pressure_drop_pa: float
    #: One entry per tray, in the column's own order.
    trays: tuple[DesignerTrayResult, ...]


@dataclass(frozen=True, slots=True)
class TrayReadings:
    """What ``getTrayFlows`` and ``getTrayProperties`` read off one tray."""

    vapor_mass_flow: float
    liquid_mass_flow: float
    vapor_density: float
    liquid_density: float
    liquid_viscosity: float
    surface_tension: float
    relative_volatility: float


def _readings(states: _States, index: int, components: list[str]) -> TrayReadings:
    """The class's own two extraction helpers, at one tray."""
    vapour = _phase(components, states, index, "gas")
    liquid = _phase(components, states, index, "liquid")

    # **The tray's own two traffics decide**, and they are the class's ``getNumberOfPhases()``:
    # a tray with no vapour or no liquid has one phase, and one with both has two.
    single_phase = states.tray_gas_n[index] <= 0.0 or states.tray_liquid_n[index] <= 0.0
    if single_phase:
        density = vapour["density"]
        substituted = (
            SINGLE_PHASE_GAS_LIQUID_DENSITY
            if vapour["kind"] == "gas"
            else SINGLE_PHASE_LIQUID_DENSITY_KG_PER_M3
        )
        vapor_density, liquid_density = density, substituted
        relative_volatility = DEFAULT_RELATIVE_VOLATILITY
    else:
        vapor_density, liquid_density = vapour["density"], liquid["density"]
        relative_volatility = _relative_volatility(vapour["z"], liquid["z"])

    viscosity = float(liquid["transport"].mu.to("Pa*s").magnitude)
    if viscosity <= 0.0:
        viscosity = DEFAULT_LIQUID_VISCOSITY_PA_S

    return TrayReadings(
        vapor_mass_flow=states.tray_gas_n[index] * vapour["molar_mass"],
        liquid_mass_flow=states.tray_liquid_n[index] * liquid["molar_mass"],
        vapor_density=vapor_density,
        liquid_density=liquid_density,
        liquid_viscosity=viscosity,
        # **The interface route answers nothing on these states** - every tray of the capture
        # prints `0.0` - so the extraction's own fallback is what stands in.
        surface_tension=DESIGNER_SURFACE_TENSION_N_PER_M,
        relative_volatility=relative_volatility,
    )


def _phase(components: list[str], states: _States, index: int, pick: str) -> dict[str, Any]:
    """The phase one tray's stream **carries**, and not what a re-flash of it answers.

    A pinned end is saturated by construction, so a re-flash of either side lands on the knife
    edge - the segment model's own readers, stated for that case and measured on the capture's
    reboiler.
    """
    from azoth.process.reference.rate_based_packed_column import _liquid_view, _vapour_view

    view = _vapour_view if pick == "gas" else _liquid_view
    return view(
        components,
        states.tray_temperature[index],
        states.tray_pressure[index],
        list(states.tray_gas_z[index] if pick == "gas" else states.tray_liquid_z[index]),
        states.tray_gas_n[index] if pick == "gas" else states.tray_liquid_n[index],
    )


def _relative_volatility(vapor_z: Any, liquid_z: Any) -> float:
    """``getTrayProperties``'s own spread: the largest ``y_i/x_i`` over the smallest, capped."""
    max_ratio = 0.0
    min_ratio = 1.0e10
    for y, x in zip(vapor_z, liquid_z, strict=True):
        if x > 1.0e-10 and y > 1.0e-10:
            ratio = y / x
            if ratio > max_ratio:
                max_ratio = ratio
            if ratio < min_ratio:
                min_ratio = ratio
    relative = max_ratio / min_ratio if min_ratio > 0.0 else DEFAULT_RELATIVE_VOLATILITY
    return min(relative, MAX_RELATIVE_VOLATILITY)


def _tray_outcome(
    geometry: DesignerGeometry,
    tray: TrayReadings,
    diameter_m: float,
    relative_volatility: float,
) -> Any:
    """One tray's calculator, at the designer's geometry and the diameter it resolved."""
    return tray_hydraulics(
        geometry.internals_type,
        quantity(diameter_m, "m"),
        quantity(geometry.tray_spacing_m, "m"),
        quantity(geometry.weir_height_m, "m"),
        # **The designer has no weir-length field at all**, so the calculator's own `-1.0` - the
        # rule that derives `0.73 D` - is what runs.
        quantity(-1.0, "m"),
        geometry.downcommer_area_fraction,
        quantity(geometry.hole_diameter_m, "m"),
        geometry.hole_area_fraction,
        geometry.design_flood_fraction,
        quantity(tray.vapor_mass_flow, "kg/s"),
        quantity(tray.liquid_mass_flow, "kg/s"),
        quantity(tray.vapor_density, "kg/m**3"),
        quantity(tray.liquid_density, "kg/m**3"),
        quantity(tray.liquid_viscosity, "Pa*s"),
        quantity(tray.surface_tension, "N/m"),
        relative_volatility,
    )


def designer_report(
    states: _States,
    components: list[str],
    geometry: DesignerGeometry,
) -> DesignerReport:
    """``ColumnInternalsDesigner.calculateTrayed``, over a solved column's trays.

    Raises:
        InvalidInputError: where the column carries no tray at all, which the class logs and
            returns from - here it is a refusal, because a report of nothing is not an answer.
    """
    from azoth.core.errors import InvalidInputError

    count = len(states.tray_temperature)
    if count == 0:
        raise InvalidInputError(
            "trays",
            "the column has no trays, so its internals have nothing to be designed against",
        )

    # ---- The controlling tray: the largest vapour mass flow, the first where two are equal.
    readings = [_readings(states, index, components) for index in range(count)]
    controlling = 0
    largest = float("-inf")
    for index, tray in enumerate(readings):
        if tray.vapor_mass_flow > largest:
            largest = tray.vapor_mass_flow
            controlling = index

    # ---- The diameter, which is the class's sizing branch or the stated override.
    if geometry.column_diameter_override_m > 0.0:
        required_diameter = geometry.column_diameter_override_m
    else:
        controller = readings[controlling]
        # **Without the tray's relative volatility**, which is the class's own omission: this
        # calculator is built before the loop and never given one.
        required_diameter = float(
            size_column_diameter(
                geometry.internals_type,
                quantity(1.0, "m"),
                quantity(geometry.tray_spacing_m, "m"),
                quantity(geometry.weir_height_m, "m"),
                quantity(-1.0, "m"),
                geometry.downcommer_area_fraction,
                quantity(geometry.hole_diameter_m, "m"),
                geometry.hole_area_fraction,
                geometry.design_flood_fraction,
                quantity(controller.vapor_mass_flow, "kg/s"),
                quantity(controller.liquid_mass_flow, "kg/s"),
                quantity(controller.vapor_density, "kg/m**3"),
                quantity(controller.liquid_density, "kg/m**3"),
                quantity(controller.liquid_viscosity, "Pa*s"),
                quantity(controller.surface_tension, "N/m"),
                DEFAULT_RELATIVE_VOLATILITY,
            )
            .to("m")
            .magnitude
        )

    # ---- The summation loop, at the diameter just resolved.
    results: list[DesignerTrayResult] = []
    total_pressure_drop = 0.0
    max_percent_flood = 0.0
    min_percent_flood = 100.0
    design_ok = True
    efficiency_sum = 0.0
    for tray in readings:
        out = _tray_outcome(geometry, tray, required_diameter, tray.relative_volatility)
        drop = float(out.total_tray_pressure_drop.to("Pa").magnitude)
        total_pressure_drop += drop
        if out.percent_flood > max_percent_flood:
            max_percent_flood = out.percent_flood
        if 0.0 < out.percent_flood < min_percent_flood:
            min_percent_flood = out.percent_flood
        if not out.design_ok:
            design_ok = False
        efficiency_sum += out.tray_efficiency
        results.append(
            DesignerTrayResult(
                percent_flood=out.percent_flood,
                total_pressure_drop_pa=drop,
                tray_efficiency=out.tray_efficiency,
                design_ok=out.design_ok,
            )
        )

    return DesignerReport(
        required_diameter_m=required_diameter,
        controlling_tray_index=controlling,
        design_ok=design_ok,
        max_percent_flood=max_percent_flood,
        min_percent_flood=min_percent_flood,
        average_tray_efficiency=(
            efficiency_sum / len(results) if results else EMPTY_AVERAGE_TRAY_EFFICIENCY
        ),
        total_pressure_drop_pa=total_pressure_drop,
        trays=tuple(results),
    )


__all__ = [
    "DesignerGeometry",
    "DesignerReport",
    "DesignerTrayResult",
    "designer_report",
]
