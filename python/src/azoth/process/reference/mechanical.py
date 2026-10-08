"""``DistillationColumnMechanicalDesign.calcDesign``: the vessel around the internals tree.

The Python twin of ``crates/azoth-process/src/column/mechanical.rs``. Every step mirrors the Rust
line for line, and **the two halves are held to the same capture**: the Rust module's own ``report``
test and ``python/tests/test_mechanical_report.py`` compare the same rows of
``validation/neqsim/captures/process_column_mechanical_design.tsv``.

**Nothing calls this yet, and neither does anything call the Rust module.** ``calcDesign`` is
ported and not wired - what no model declares is its inputs - so this is the other half of a pair
that waits for ``process.distillation_column`` to publish the vessel.

**The class runs two sizings and publishes the second one's diameter beside the first one's
numbers**, which is the shape ``validation/neqsim/captures/process_column_mechanical_design.tsv``
was written to measure. The first pass is Souders-Brown on tray 0's two outlets; at that diameter
it also fixes the flooding factor and takes a weir loading from ``0.7 D``. Then
``calculateContactorCapacity`` builds the internals designer, hands it
:func:`resolve_rating_diameter` as an override, and **replaces the diameter with the designer's
answer** - while the three quantities above keep the values they were computed at.

**And the designer is never allowed to size.** The rating diameter answers the stated override
where there is one and the first pass's own diameter otherwise, and the first pass is never below
the class's ``0.5`` floor - so ``getRequiredDiameter()`` always returns the override, the sizing
branch is dead on this path, and the class's own ``contactorDesignFloodFraction`` of ``0.70`` is
never read. Both are stated here rather than ported as knobs.

**One quantity of ``calcDesign`` is computed and then thrown away**: ``totalPressureDrop`` is
assigned the tray estimate at the class's own line 376 and unconditionally overwritten by the
designer's sum at 512, so this computes only the second.

**The two tray-0 density reads are upstream defect #4140**, and this is where the port states the
difference rather than reproducing it: the class takes ``getFluid().getDensity()`` without ever
asking the fluid to initialise, so on an absorber's *liquid* outlet it reads ``0.0`` - and the weir
loading dividing by it is ``Infinity`` while the tray pressure drop drops its liquid-head term.
azoth always has a density for a ``(T, P, z)``, so this port publishes the number the capture
holds after one ``initProperties()``. **The difference is the defect and not a correction of it** -
the capture carries both sides on every row, and this is held to the repaired one.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import TYPE_CHECKING

from azoth.core.errors import InvalidInputError
from azoth.hydraulics.reference.tray_hydraulics import round_to_standard_diameter
from azoth.process.reference.designer import DesignerGeometry, designer_report
from azoth.process.reference.rate_based_packed_column import _liquid_view, _vapour_view

if TYPE_CHECKING:
    from azoth.process.reference.distillation_column import _States

#: ``trayEfficiency``'s own field initialiser, and the efficiency ``calcDesign`` divides the tray
#: count by.
DEFAULT_TRAY_EFFICIENCY = 0.65
#: ``maxFloodingFactor``'s own initialiser - the fraction of flood the design velocity is.
DEFAULT_MAX_FLOODING_FACTOR = 0.85
#: ``trayType``'s own initialiser. **It answers twice**: the Souders-Brown ``kFactor`` below and,
#: through ``resolveInternalsType``, the designer's internals type.
DEFAULT_TRAY_TYPE = "sieve"
#: ``contactorInternalsType``'s own initialiser.
DEFAULT_CONTACTOR_INTERNALS_TYPE = "auto"
#: ``contactorDesignFloodFraction``'s own initialiser. **Carried and not declared**, because the
#: designer it reaches is always given an override and so never runs the branch that reads it.
DEFAULT_CONTACTOR_DESIGN_FLOOD_FRACTION = 0.7
#: ``materialGrade``'s own initialiser.
DEFAULT_MATERIAL_GRADE = "SA-516-70"
#: ``MechanicalDesign.maxOperationPressure``'s own field initialiser, bara.
#:
#: **The equipment path overwrites it**: ``DistillationColumn`` sets it from
#: ``getEconomicDesignPressure()`` - the larger of the two endpoint pressures - before it calls
#: ``calcDesign``, and the capture's rows construct the design directly and so never take that
#: path. This declares it, so a case can state either.
DEFAULT_MAX_OPERATION_PRESSURE_BARA = 100.0
#: ``MechanicalDesign.tensileStrength``'s own initialiser, MPa. Not declared: no path this
#: module's caller runs reaches its setter.
TENSILE_STRENGTH_MPA = 483.0
#: ``MechanicalDesign.jointEfficiency``'s own initialiser, fully radiographed.
JOINT_EFFICIENCY = 1.0
#: ``MechanicalDesign.corrosionAllowance``'s own initialiser, mm.
CORROSION_ALLOWANCE_MM = 0.0
#: ``getMaxOperationPressure() * 1.1``, the class's own design margin.
DESIGN_PRESSURE_MARGIN = 1.1
#: ``getTensileStrength() * 0.4``, the allowable stress in ``calculateRequiredWallThickness``.
ALLOWABLE_STRESS_FRACTION = 0.4
#: The floor ``calculateRequiredWallThickness`` takes its answer at, mm.
MINIMUM_WALL_THICKNESS_MM = 6.0
#: The Souders-Brown ``kFactor`` per tray type, m/s, **compared case-sensitively** as the class's
#: own ``"valve".equals(trayType)`` does.
SIEVE_K_FACTOR = 0.1
#: See :data:`SIEVE_K_FACTOR`.
VALVE_K_FACTOR = 0.12
#: See :data:`SIEVE_K_FACTOR`.
BUBBLE_CAP_K_FACTOR = 0.08
#: ``dryTrayDp``, the class's own constant for a sieve tray, mbar/tray. It also stands in for the
#: whole tray pressure drop wherever #4140 zeroes the liquid density.
DRY_TRAY_PRESSURE_DROP_MBAR = 5.0
#: The weir length as a fraction of the diameter: ``weirLength = columnDiameter * 0.7``.
WEIR_LENGTH_FRACTION = 0.7
#: The vapour disengagement section, m.
TOP_SECTION_M = 1.0
#: The liquid holdup section, m.
BOTTOM_SECTION_M = 2.0
#: One head, m - taken twice.
HEAD_HEIGHT_M = 0.5
#: The diameter the class falls back to where the vapour volume flow is not positive, before
#: rounding it through the table: ``roundToStandardDiameter(0.5)``, which is ``0.5``.
DEGENERATE_DIAMETER_M = 0.5
#: The class's own ``500.0`` in ``calculateRequiredWallThickness``: ``D/2`` in millimetres.
THICKNESS_DIAMETER_MM_PER_M = 500.0
#: ``(int)``'s saturation, which the class's own tray-count cast takes its ceiling at.
I32_MAX = 2147483647


@dataclass(frozen=True, slots=True)
class MechanicalGeometry:
    """``DistillationColumnMechanicalDesign``'s own geometry and design statements.

    **Six of these are the class's fields and the rest are the internals tree's**, which the class
    passes straight through to the designer it builds - so a caller states each of them once. The
    designer's own ``designFloodFraction`` is the exception: the class overwrites it with its own
    :data:`DEFAULT_CONTACTOR_DESIGN_FLOOD_FRACTION`, which is why it is not here.
    """

    #: ``trayType``: the Souders-Brown ``kFactor`` and, through ``resolveInternalsType``, the
    #: designer's internals type.
    tray_type: str = DEFAULT_TRAY_TYPE
    #: ``contactorInternalsType``. ``auto`` resolves to :attr:`tray_type` on this equipment.
    contactor_internals_type: str = DEFAULT_CONTACTOR_INTERNALS_TYPE
    #: ``trayEfficiency``, which the tray count is divided by.
    tray_efficiency: float = DEFAULT_TRAY_EFFICIENCY
    #: ``maxFloodingFactor``, the fraction of flood the design velocity is.
    max_flooding_factor: float = DEFAULT_MAX_FLOODING_FACTOR
    #: ``materialGrade``, a vessel statement the class carries and reports.
    material_grade: str = DEFAULT_MATERIAL_GRADE
    #: ``getMaxOperationPressure()``, bara - the wall thickness's design basis.
    max_operation_pressure_bara: float = DEFAULT_MAX_OPERATION_PRESSURE_BARA
    #: The tray spacing, m, passed to the designer.
    tray_spacing_m: float = 0.6
    #: The weir height, m: the tray pressure drop's liquid head and the designer's own field.
    weir_height_m: float = 0.05
    #: The sieve hole diameter, m.
    hole_diameter_m: float = 12.7e-3
    #: The hole area over the active area.
    hole_area_fraction: float = 0.1
    #: The downcomer area over the total, which the sizing's total area is divided by.
    downcommer_area_fraction: float = 0.1
    #: ``columnDiameterOverride``, m: a stated diameter, or ``-1`` for the class's own sizing.
    column_diameter_override_m: float = -1.0


@dataclass(frozen=True, slots=True)
class MechanicalReport:
    """``calcDesign``'s own report, minus the capacity result, the demister and the costs."""

    #: ``getColumnDiameter()``: the diameter the designer answered, which is
    #: :func:`resolve_rating_diameter`'s own value - its first pass's rounded diameter where none
    #: was stated.
    vessel_diameter_m: float
    #: ``getColumnHeight()``: ``actualTrays * spacing + 1 + 2 + 2 * 0.5``, m.
    vessel_height_m: float
    #: ``getColumnWallThickness()``, mm - taken at the *final* diameter.
    vessel_wall_thickness_mm: float
    #: ``getActualTrays()``: ``ceil(trays / efficiency)``, the class's own integer cast saturated
    #: at ``i32::MAX``.
    actual_trays: int
    #: ``getFloodingFactor()``: the actual velocity over the flooding velocity, **at the first
    #: pass's diameter** - the class's own order of statements, and not at the one it publishes.
    flooding_factor: float
    #: ``getWeirLoading()``: the liquid volume flow over ``0.7 D``, m3/hr per m, at the same first
    #: pass's diameter.
    weir_loading: float
    #: ``getTrayPressureDrop()``, mbar/tray: the class's own ``5.0`` constant plus
    #: ``weirHeight * liquidDensity * 9.81 / 100``.
    tray_pressure_drop_mbar: float
    #: ``getTotalPressureDrop()``, **bar**: the designer's summed pressure drop over ``1e5``, and
    #: not the tray estimate the class computes and discards.
    total_pressure_drop_bar: float
    #: ``getReboilerDuty()``, kW.
    reboiler_duty_kw: float
    #: ``getCondenserDuty()``, kW - **the absolute value**, which is the class's own.
    condenser_duty_kw: float
    #: ``getMaterialGrade()``, carried and reported.
    material_grade: str


def mechanical_design(
    states: _States,
    components: list[str],
    geometry: MechanicalGeometry,
) -> MechanicalReport:
    """``calcDesign``, over a solved column.

    Raises:
        InvalidInputError: where the column carries no tray at all, where ``tray_efficiency`` is
            not positive and finite, where ``contactor_internals_type`` is one this port does not
            carry, and whatever a tray's flash, its two streams, their densities or the designer
            raise.
    """
    count = len(states.tray_temperature)
    if count == 0:
        raise InvalidInputError(
            "trays",
            "`DistillationColumnMechanicalDesign.calcDesign` returns early on a column with no "
            "tray and leaves every quantity at its own initialiser - a vessel of zero diameter "
            "and zero height, which is not a design of anything",
        )
    if not math.isfinite(geometry.tray_efficiency) or geometry.tray_efficiency <= 0.0:
        raise InvalidInputError(
            "tray_efficiency",
            f"a tray efficiency of {geometry.tray_efficiency} is not positive and finite, and "
            f"`calcDesign` divides the tray count by it with no guard - `(int) Math.ceil` "
            f"answers {I32_MAX} there, which is a saturation of its own cast rather than a tray "
            f"count",
        )
    active_type = resolve_internals_type(geometry)

    # ---- The tray count, and what the efficiency makes of it.
    actual_trays = int(min(max(math.ceil(count / geometry.tray_efficiency), 0.0), float(I32_MAX)))

    # ---- Tray 0's two outlets, which is the whole of the class's input surface. **The density
    # is the *phase's*, not the stream's system average**, and that is a measurement rather than
    # a preference: on the capture's row the label rule splits the tray's vapour outlet into two
    # phases while its flash reports a vapour fraction of `1.0`, so the system density is
    # `84.32487` where the class reads `44.17833075904439` - the vapour's own.
    gas_n, gas_z, gas_t, gas_p = _tray_zero(states, "gas")
    liquid_n, liquid_z, liquid_t, liquid_p = _tray_zero(states, "liquid")
    vapour = _vapour_view(components, gas_t, gas_p, gas_z, gas_n)
    liquid = _liquid_view(components, liquid_t, liquid_p, liquid_z, liquid_n)
    vapor_molar_flow = gas_n * 3600.0
    liquid_molar_flow = liquid_n * 3600.0
    # The class's own `getMolarMass() * 1000`, which cancels against the `/ 1000` below and is
    # kept so this reads against the source.
    vapor_mw = float(vapour["molar_mass"]) * 1000.0
    liquid_mw = float(liquid["molar_mass"]) * 1000.0
    vapor_density = float(vapour["density"])
    liquid_density = float(liquid["density"])
    vapor_mass_flow = vapor_molar_flow * vapor_mw / 1000.0
    liquid_mass_flow = liquid_molar_flow * liquid_mw / 1000.0
    vapor_volume_flow = vapor_mass_flow / vapor_density
    liquid_volume_flow = liquid_mass_flow / liquid_density

    # ---- The Souders-Brown pass, whose diameter the second sizing then throws away.
    k_factor_value = k_factor(geometry.tray_type)
    u_flood_raw = k_factor_value * math.sqrt((liquid_density - vapor_density) / vapor_density)
    u_flood = (
        k_factor_value
        if math.isnan(u_flood_raw) or math.isinf(u_flood_raw) or u_flood_raw <= 0.0
        else u_flood_raw
    )
    u_design = u_flood * geometry.max_flooding_factor
    vapor_volume_flow_m3s = vapor_volume_flow / 3600.0
    if (
        vapor_volume_flow_m3s <= 0.0
        or math.isnan(vapor_volume_flow_m3s)
        or math.isinf(vapor_volume_flow_m3s)
        or u_design <= 0.0
    ):
        first_pass = round_to_standard_diameter(DEGENERATE_DIAMETER_M)
    else:
        required_vapor_area = vapor_volume_flow_m3s / u_design
        total_area = required_vapor_area / (1.0 - geometry.downcommer_area_fraction)
        first_pass = round_to_standard_diameter(math.sqrt(4.0 * total_area / math.pi))

    # ---- The three quantities read *at* that diameter and published beside the one below.
    actual_area = math.pi * (first_pass / 2.0) ** 2
    actual_vapor_area = actual_area * (1.0 - geometry.downcommer_area_fraction)
    flooding_factor = vapor_volume_flow_m3s / actual_vapor_area / u_flood
    weir_length = first_pass * WEIR_LENGTH_FRACTION
    weir_loading = liquid_volume_flow / weir_length

    tray_pressure_drop_mbar = (
        DRY_TRAY_PRESSURE_DROP_MBAR + geometry.weir_height_m * liquid_density * 9.81 / 100.0
    )

    # ---- The internals tree, given the rating diameter as an override so it never sizes.
    rating_diameter = resolve_rating_diameter(geometry, first_pass)
    internals = designer_report(
        states,
        components,
        DesignerGeometry(
            internals_type=active_type,
            tray_spacing_m=geometry.tray_spacing_m,
            weir_height_m=geometry.weir_height_m,
            hole_diameter_m=geometry.hole_diameter_m,
            hole_area_fraction=geometry.hole_area_fraction,
            downcommer_area_fraction=geometry.downcommer_area_fraction,
            design_flood_fraction=DEFAULT_CONTACTOR_DESIGN_FLOOD_FRACTION,
            column_diameter_override_m=rating_diameter,
        ),
    )
    # **The designer's answer replaces the diameter**, which the class guards with `> 0.0`.
    vessel_diameter = (
        internals.required_diameter_m if internals.required_diameter_m > 0.0 else first_pass
    )

    return MechanicalReport(
        vessel_diameter_m=vessel_diameter,
        vessel_height_m=(
            actual_trays * geometry.tray_spacing_m
            + TOP_SECTION_M
            + BOTTOM_SECTION_M
            + 2.0 * HEAD_HEIGHT_M
        ),
        # Taken at the final diameter, because the class recomputes it there.
        vessel_wall_thickness_mm=required_wall_thickness(
            vessel_diameter, geometry.max_operation_pressure_bara
        ),
        actual_trays=actual_trays,
        flooding_factor=flooding_factor,
        weir_loading=weir_loading,
        tray_pressure_drop_mbar=tray_pressure_drop_mbar,
        total_pressure_drop_bar=internals.total_pressure_drop_pa / 1.0e5,
        reboiler_duty_kw=states.reboiler_duty / 1000.0,
        condenser_duty_kw=abs(states.condenser_duty) / 1000.0,
        material_grade=geometry.material_grade,
    )


def _tray_zero(states: _States, pick: str) -> tuple[float, list[float], float, float]:
    """Tray 0's outlet for one side, **with the class's own fallback applied**.

    ``calcDesign`` reads the tray's traffic and, where it is not positive, takes the column's own
    product instead. Its guard is the *flow* and not the stream - the product is read whenever
    the tray's traffic is not positive, even where that is also zero.
    """
    gas = pick == "gas"
    tray_n = states.tray_gas_n[0] if gas else states.tray_liquid_n[0]
    if tray_n > 0.0:
        return (
            tray_n,
            list(states.tray_gas_z[0] if gas else states.tray_liquid_z[0]),
            states.tray_temperature[0],
            states.tray_pressure[0],
        )
    if gas:
        return (
            states.distillate_n,
            list(states.distillate_z),
            states.distillate_t,
            states.distillate_p,
        )
    return states.bottoms_n, list(states.bottoms_z), states.bottoms_t, states.bottoms_p


def resolve_internals_type(geometry: MechanicalGeometry) -> str:
    """``resolveInternalsType``: the class's own ``contactorInternalsType`` unless it is ``auto``,
    in which case the equipment's own - which on a ``DistillationColumn`` is ``trayType``.

    Raises:
        InvalidInputError: for ``packed``. The class carries that branch, and reaching it here
            would need a packed preset, a bed height and a capacity factor and would put the
            height on the class's ``packedHeight + 4`` path instead of the trayed one - a second
            machine, refused by name rather than half-carried.
    """
    if geometry.contactor_internals_type.lower() == "auto":
        return geometry.tray_type
    if geometry.contactor_internals_type.lower() == "packed":
        # The *id* refuses this value from its own declaration, and `distillation_column` is
        # where that refusal is written - a refusal is a claim about a declaration, and this
        # module carries none. This guard is for a caller that reaches the vessel directly.
        raise InvalidInputError(
            "contactor_internals_type",
            "a packed contactor has no trayed design to build, so this vessel cannot be sized "
            "from it - `process.distillation_column` refuses the value from its declaration",
        )
    return geometry.contactor_internals_type


def resolve_rating_diameter(geometry: MechanicalGeometry, first_pass: float) -> float:
    """``resolveRatingDiameter``: the stated override where there is one, and otherwise the
    diameter the Souders-Brown pass resolved.

    **The class's packed arm is not reachable here**, and its ``column instanceof PackedColumn``
    test could not be reached on this equipment even where a packed contactor type were carried -
    so the two conditions that remain are these.
    """
    if geometry.column_diameter_override_m > 0.0:
        return geometry.column_diameter_override_m
    return first_pass


def required_wall_thickness(diameter_m: float, max_operation_pressure_bara: float) -> float:
    """``calculateRequiredWallThickness``, whose design pressure is ``bara``, whose stress is
    ``MPa`` and whose answer is millimetres - the class's own unit mixing, reproduced rather than
    corrected.
    """
    design_pressure = max_operation_pressure_bara * DESIGN_PRESSURE_MARGIN
    allowable_stress = TENSILE_STRENGTH_MPA * ALLOWABLE_STRESS_FRACTION
    denominator = allowable_stress * JOINT_EFFICIENCY - 0.6 * design_pressure
    calculated = (
        design_pressure * diameter_m * THICKNESS_DIAMETER_MM_PER_M / denominator
        if denominator > 0.0
        else 0.0
    )
    return max(calculated + CORROSION_ALLOWANCE_MM, MINIMUM_WALL_THICKNESS_MM)


def k_factor(tray_type: str) -> float:
    """The class's three Souders-Brown constants, compared with ``String.equals`` and so
    case-sensitively: an unrecognised type is the sieve tray's ``0.1``, which is what its own
    constructor default is.
    """
    if tray_type == "valve":
        return VALVE_K_FACTOR
    if tray_type == "bubble-cap":
        return BUBBLE_CAP_K_FACTOR
    return SIEVE_K_FACTOR
