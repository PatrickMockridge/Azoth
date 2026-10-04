"""``hydraulics.tray_hydraulics`` - a tray's flooding, weeping, entrainment, pressure drop and
efficiency.

Spec: ``specs/calcs/hydraulics/tray_hydraulics.toml``.

The Python twin of ``crates/azoth-hydraulics/src/tray_hydraulics.rs``, and both mirror
``TrayHydraulicsCalculator``. **Five correlations, the class's own**: Fair's capacity factor and
its entrainment fit, Sinnott's weeping minimum, O'Connell's efficiency, and a three-term tray
pressure drop.

# The tray type is read at four sites and its fallback differs at each

The hole area is ``hole_area_fraction`` / ``0.14`` / ``0.12``, the flooding factor ``1.0`` /
``1.05`` / ``0.85``, the orifice coefficient ``0.73`` / ``0.80`` / ``0.65``, and only the weeping
check returns for ``bubble-cap`` alone. So an unrecognised type is *not* uniformly a bubble-cap,
which is why the type is compared at each site rather than resolved once.
"""

from __future__ import annotations

import math
from itertools import pairwise

from azoth._registry_gen import spec as _spec_for
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import TrayHydraulicsResult
from azoth.core.units import Q, from_si, input_to_si, quantity
from azoth.core.warnings import Warning

CALC_ID = "hydraulics.tray_hydraulics"

#: The gravitational acceleration the class writes inline, in m/s**2.
G = 9.81

#: The diameter ``sizeColumnDiameter`` sizes *from*: the class writes ``1.0`` into
#: ``columnDiameter``, re-derives the areas and the flooding velocity there, and sizes from that.
TRIAL_DIAMETER_M = 1.0
#: The floor on the vapour's density in that sizing, ``max(vaporDensity, 0.01)``.
VAPOR_DENSITY_FLOOR = 0.01
#: ``roundToStandardDiameter``'s table, m - the class's own thirty-one sizes.
STANDARD_DIAMETERS_M = (
    0.5,
    0.6,
    0.7,
    0.8,
    0.9,
    1.0,
    1.1,
    1.2,
    1.4,
    1.5,
    1.6,
    1.8,
    2.0,
    2.2,
    2.4,
    2.6,
    2.8,
    3.0,
    3.2,
    3.4,
    3.6,
    3.8,
    4.0,
    4.5,
    5.0,
    5.5,
    6.0,
    7.0,
    8.0,
    9.0,
    10.0,
)

#: The surface tension the flooding velocity's correction is referenced to, in N/m.
REFERENCE_SURFACE_TENSION = 0.020

#: The Fair capacity factor's correction, both exponents the class's own.
FLV_CORRECTION_SCALE = 1.463
FLV_CORRECTION_POWER = 0.842
#: The flow parameter the tabulated factors are stated at, which normalises them.
FLV_REFERENCE = 0.1
MIN_FLV, MAX_FLV = 0.01, 2.0

#: The flow parameter's bands and Fair's ``(a, b)`` pairs.
ENTRAINMENT_LOW_FLV, ENTRAINMENT_HIGH_FLV = 0.02, 0.2
ENTRAINMENT_LOW, ENTRAINMENT_MID, ENTRAINMENT_HIGH = (0.085, 4.2), (0.05, 3.8), (0.03, 3.5)
ENTRAINMENT_LIMIT = 0.1

#: The orifice discharge coefficients, one per tray type.
ORIFICE_SIEVE, ORIFICE_VALVE, ORIFICE_OTHER = 0.73, 0.80, 0.65
#: The flooding velocity's tray-type factors.
FLOOD_FACTOR_VALVE, FLOOD_FACTOR_BUBBLE_CAP = 1.05, 0.85
#: The hole area's fallback fractions.
HOLE_FRACTION_VALVE, HOLE_FRACTION_OTHER = 0.14, 0.12
#: The weir length's auto-derived fraction of the diameter.
WEIR_LENGTH_FRACTION = 0.73

#: Sinnott's weeping-chart fit and its two references.
WEEPING_CHART_BASE, WEEPING_CHART_SLOPE = 7.0, 0.27
WEEPING_HOLE_SLOPE, WEEPING_REFERENCE_HOLE_MM = 0.90, 25.4
WEEPING_NUMERATOR_FLOOR = 0.5

#: O'Connell's efficiency constants and its two clamps.
OCONNELL_BASE, OCONNELL_SLOPE = 51.0, 32.5
OCONNELL_ALPHA_MU_MIN, OCONNELL_ALPHA_MU_MAX = 0.1, 10.0
OCONNELL_EFFICIENCY_MIN, OCONNELL_EFFICIENCY_MAX = 10.0, 100.0

#: The design window the flood is judged inside, in per cent.
DESIGN_FLOOD_MIN, DESIGN_FLOOD_MAX = 50.0, 85.0
#: The downcomer backup's allowance, as a fraction of ``spacing + weir height``.
BACKUP_LIMIT_FRACTION = 0.5
#: The apron gap and its two constants.
APRON_GAP_FRACTION, APRON_LOSS_COEFFICIENT, APRON_AREA_FLOOR_M2 = 0.025, 166.0, 0.001

#: The Fair capacity factors at ``FLV = 0.1``, tabulated against tray spacing.
CAPACITY_FACTOR_TABLE = (
    (0.15, 0.025),
    (0.23, 0.040),
    (0.30, 0.050),
    (0.46, 0.070),
    (0.61, 0.085),
    (0.91, 0.105),
)


def round_to_standard_diameter(diameter: float) -> float:
    """``roundToStandardDiameter``: the first standard size at or above ``diameter``, m.

    **Past the table's end the class rounds to the nearest half metre** rather than refusing, so
    a column wider than ``10`` m answers ``ceil(2 d) / 2`` - its own arithmetic, not a policy.
    """
    for standard in STANDARD_DIAMETERS_M:
        if standard >= diameter:
            return standard
    return math.ceil(diameter * 2.0) / 2.0


def size_column_diameter(
    tray_type: str,
    column_diameter: Q,
    tray_spacing: Q,
    weir_height: Q,
    weir_length: Q,
    downcommer_area_fraction: float,
    hole_diameter: Q,
    hole_area_fraction: float,
    design_flood_fraction: float,
    vapor_mass_flow: Q,
    liquid_mass_flow: Q,
    vapor_density: Q,
    liquid_density: Q,
    liquid_viscosity: Q,
    surface_tension: Q,
    relative_volatility: float,
) -> Q:
    """``TrayHydraulicsCalculator.sizeColumnDiameter``: the diameter a tray sizes to.

    **The trial diameter is ``1.0`` m, and it is not a guess at the answer.** The class writes
    ``1.0`` into ``columnDiameter``, re-derives the areas and the **flooding velocity** there, and
    sizes from *that* velocity - so the flooding velocity this reads is the one at a ``1.0`` m
    tray and not at the diameter it resolves. The rest is the class's own arithmetic: the vapour's
    volumetric flow at ``max(rho_v, 0.01)``, the design velocity at ``designFloodFraction``, the
    area from the two, the **net** area by dividing out ``1 - downcommerAreaFraction``, and the
    round up to the standard table.

    **A design velocity that is not positive answers ``1.0`` and leaves the trial diameter
    behind**, which is the class's own early return rather than a refusal.

    Args:
        tray_type: the tray type, which the flooding velocity's factor and the orifice read.
        column_diameter: the stated diameter, which **the trial overwrites before anything is
            sized** - carried so that these are the calculation's own sixteen inputs and read by
            nothing.
        tray_spacing: the tray spacing.
        weir_height: the weir height.
        weir_length: at or below zero derives ``0.73 D``.
        downcommer_area_fraction: the downcomer fraction the net area is taken over.
        hole_diameter: the hole diameter, in millimetres.
        hole_area_fraction: the hole area over the active area.
        design_flood_fraction: the fraction of flood the design velocity is.
        vapor_mass_flow: the vapour's mass flow.
        liquid_mass_flow: the liquid's mass flow.
        vapor_density: the vapour's mass density.
        liquid_density: the liquid's mass density.
        liquid_viscosity: the liquid's dynamic viscosity.
        surface_tension: the interfacial surface tension.
        relative_volatility: the relative volatility of the key components.

    Returns:
        The standard diameter the tray sizes to, m.

    See :func:`azoth.hydraulics.reference.tray_hydraulics.tray_hydraulics`.
    """
    # **The stated diameter is overwritten and never read**: the class writes `1.0` into
    # `columnDiameter` before it sizes, so a caller's own diameter changes nothing here.
    _ = column_diameter
    # **The trial run is the whole calculation at that diameter**, which is what the class does
    # too: it re-derives the areas and the flooding velocity there and sizes from *that*
    # velocity rather than from the one it resolves.
    trial = tray_hydraulics(
        tray_type,
        quantity(TRIAL_DIAMETER_M, "m"),
        tray_spacing,
        weir_height,
        weir_length,
        downcommer_area_fraction,
        hole_diameter,
        hole_area_fraction,
        design_flood_fraction,
        vapor_mass_flow,
        liquid_mass_flow,
        vapor_density,
        liquid_density,
        liquid_viscosity,
        surface_tension,
        relative_volatility,
    )
    # **The result carries its velocities as quantities**, which is the shape the dataclass
    # declares for them, so this reads the value in m/s rather than the object.
    flooding_velocity = float(trial.flooding_velocity.to("m/s").magnitude)
    density = float(vapor_density.to("kg/m**3").magnitude)
    volumetric_flow = float(vapor_mass_flow.to("kg/s").magnitude) / max(
        density, VAPOR_DENSITY_FLOOR
    )
    design_velocity = flooding_velocity * float(design_flood_fraction)
    if not design_velocity > 0.0:
        return from_si(TRIAL_DIAMETER_M, "m")
    area = volumetric_flow / design_velocity
    net_area = area / (1.0 - float(downcommer_area_fraction))
    return from_si(round_to_standard_diameter(math.sqrt(4.0 * net_area / math.pi)), "m")


def capacity_factor(spacing: float, flv: float) -> float:
    """``getCapacityFactor``: the tabulated Fair factors interpolated, then corrected.

    The tabulation is the class's own, with its linear extrapolation past the last row; the
    correction is ``exp(-1.463 FLV**0.842)`` normalised at ``FLV = 0.1``, and **the clamp bites**,
    so a flow parameter below ``0.01`` or above ``2`` is read at the bound.
    """
    if spacing <= CAPACITY_FACTOR_TABLE[0][0]:
        k_at_flv_01 = CAPACITY_FACTOR_TABLE[0][1]
    elif spacing > CAPACITY_FACTOR_TABLE[-1][0]:
        k_at_flv_01 = CAPACITY_FACTOR_TABLE[-1][1] + (spacing - CAPACITY_FACTOR_TABLE[-1][0]) * 0.02
    else:
        k_at_flv_01 = CAPACITY_FACTOR_TABLE[-1][1]
        for (low_spacing, low_k), (high_spacing, high_k) in pairwise(CAPACITY_FACTOR_TABLE):
            if spacing <= high_spacing:
                k_at_flv_01 = low_k + (spacing - low_spacing) / (high_spacing - low_spacing) * (
                    high_k - low_k
                )
                break
    flv_clamped = min(max(flv, MIN_FLV), MAX_FLV)
    correction = math.exp(-FLV_CORRECTION_SCALE * math.pow(flv_clamped, FLV_CORRECTION_POWER))
    normalizer = math.exp(-FLV_CORRECTION_SCALE * math.pow(FLV_REFERENCE, FLV_CORRECTION_POWER))
    return k_at_flv_01 * correction / normalizer


def tray_hydraulics(
    tray_type: str,
    column_diameter: Q,
    tray_spacing: Q,
    weir_height: Q,
    weir_length: Q,
    downcommer_area_fraction: float,
    hole_diameter: Q,
    hole_area_fraction: float,
    design_flood_fraction: float,
    vapor_mass_flow: Q,
    liquid_mass_flow: Q,
    vapor_density: Q,
    liquid_density: Q,
    liquid_viscosity: Q,
    surface_tension: Q,
    relative_volatility: float,
) -> TrayHydraulicsResult:
    """A tray's hydraulics at a stated geometry and load.

    Args:
        tray_type: ``sieve``, ``valve`` or ``bubble-cap``, which four branches read.
        column_diameter: the column's internal diameter, which the four areas come from.
        tray_spacing: the tray spacing, which the capacity factor and the backup allowance read.
        weir_height: the weir height, which the liquid head and the weeping chart take in mm.
        weir_length: **at or below zero derives** ``0.73 D``, and the derived value is what every
            crest term reads.
        downcommer_area_fraction: the downcomer area over the total, which fixes the active area.
        hole_diameter: the hole diameter, **in millimetres**.
        hole_area_fraction: the hole area over the active area, read for a ``sieve`` tray alone.
        design_flood_fraction: stored and read by nothing, which is the class's own state.
        vapor_mass_flow: the vapour's mass flow.
        liquid_mass_flow: the liquid's mass flow.
        vapor_density: the vapour's mass density.
        liquid_density: the liquid's mass density.
        liquid_viscosity: the liquid's dynamic viscosity.
        surface_tension: the interfacial surface tension.
        relative_volatility: the relative volatility of the key components.

    Returns:
        The four verdicts, the three pressure-drop terms, the four areas and the class's rates.

    Raises:
        OutOfRangeError: for any input outside the range its own division requires.

    Example:
        >>> import azoth
        >>> q = azoth.ureg.Quantity
        >>> r = tray_hydraulics(
        ...     "sieve", q(1.0, "m"), q(0.6, "m"), q(0.05, "m"), q(-1.0, "m"), 0.1,
        ...     q(12.7, "mm"), 0.1, 0.8, q(1.0, "kg/s"), q(3.0, "kg/s"), q(2.0, "kg/m**3"),
        ...     q(800.0, "kg/m**3"), q(1.0e-3, "Pa*s"), q(0.02, "N/m"), 2.0)
        >>> round(r.percent_flood, 10)
        45.9273367045
    """
    spec = _spec_for(CALC_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    # **A declared input arrives in one of two shapes**: a unit-carrying one is a quantity and a
    # dimensionless one is the bare number it is. `input_to_si` refuses the second, so the three
    # dimensionless inputs of this calculation go straight into the checks dict - the shape
    # `hydraulics.packing_hydraulics` gives its capacity factor.
    values = {
        "column_diameter": input_to_si(spec, "column_diameter", column_diameter),
        "tray_spacing": input_to_si(spec, "tray_spacing", tray_spacing),
        "weir_height": input_to_si(spec, "weir_height", weir_height),
        "weir_length": input_to_si(spec, "weir_length", weir_length),
        "hole_diameter": input_to_si(spec, "hole_diameter", hole_diameter),
        "vapor_mass_flow": input_to_si(spec, "vapor_mass_flow", vapor_mass_flow),
        "liquid_mass_flow": input_to_si(spec, "liquid_mass_flow", liquid_mass_flow),
        "vapor_density": input_to_si(spec, "vapor_density", vapor_density),
        "liquid_density": input_to_si(spec, "liquid_density", liquid_density),
        "liquid_viscosity": input_to_si(spec, "liquid_viscosity", liquid_viscosity),
        "surface_tension": input_to_si(spec, "surface_tension", surface_tension),
    }
    # **`design_flood_fraction` is stored and read by nothing**, which is the class's own state:
    # `designFloodFraction` is a field no branch of `calculate` touches. It is named here so the
    # parameter is not dropped, and there is nothing to do with it after that.
    _ = design_flood_fraction
    downcommer_fraction = float(downcommer_area_fraction)
    applied_hole_fraction = float(hole_area_fraction)
    volatility = float(relative_volatility)
    apply_checks(
        checks.on_input,
        {
            **values,
            "downcommer_area_fraction": downcommer_fraction,
            "hole_area_fraction": applied_hole_fraction,
            "relative_volatility": volatility,
        }.get,
        warnings,
    )

    diameter = values["column_diameter"]
    spacing = values["tray_spacing"]
    weir_height_m = values["weir_height"]
    vapor_mass_flow_si = values["vapor_mass_flow"]
    liquid_mass_flow_si = values["liquid_mass_flow"]
    vapor_density_si = values["vapor_density"]
    liquid_density_si = values["liquid_density"]
    liquid_viscosity_si = values["liquid_viscosity"]
    surface_tension_si = values["surface_tension"]
    # **The spec declares millimetres and the calcs work in SI**, so the two correlations that
    # read the hole diameter - Sinnott's minimum and the residual head - convert back here.
    hole_diameter_mm = values["hole_diameter"] * 1000.0
    sieve = tray_type.lower() == "sieve"
    valve = tray_type.lower() == "valve"
    bubble_cap = tray_type.lower() == "bubble-cap"

    # ---- `calculateAreas`. `activeAreaFraction` is derived here and not an input.
    total_area = math.pi / 4.0 * diameter * diameter
    downcommer_area = total_area * downcommer_fraction
    active_area_fraction = 1.0 - 2.0 * downcommer_fraction
    active_area = total_area * active_area_fraction
    if sieve:
        hole_area = active_area * applied_hole_fraction
    elif valve:
        hole_area = active_area * HOLE_FRACTION_VALVE
    else:
        hole_area = active_area * HOLE_FRACTION_OTHER
    # **The derived weir length is what every crest term reads**, so a stated one replaces it
    # everywhere below rather than being reported back.
    calculated_weir_length = (
        WEIR_LENGTH_FRACTION * diameter if values["weir_length"] <= 0.0 else values["weir_length"]
    )

    # ---- `calculateFloodingVelocity`.
    flv = 0.0
    if vapor_mass_flow_si > 0.0 and vapor_density_si > 0.0 and liquid_density_si > 0.0:
        flv = (liquid_mass_flow_si / vapor_mass_flow_si) * math.sqrt(
            vapor_density_si / liquid_density_si
        )
    k_base = capacity_factor(spacing, flv)
    sigma_corrected = k_base * math.pow(surface_tension_si / REFERENCE_SURFACE_TENSION, 0.2)
    if valve:
        tray_factor = FLOOD_FACTOR_VALVE
    elif bubble_cap:
        tray_factor = FLOOD_FACTOR_BUBBLE_CAP
    else:
        tray_factor = 1.0
    k_final = sigma_corrected * tray_factor
    delta_rho = max(liquid_density_si - vapor_density_si, 0.1)
    flooding_velocity = k_final * math.sqrt(delta_rho / vapor_density_si)

    # ---- `calculateActualVelocity`. The net area is the total less **one** downcomer.
    net_area = total_area - downcommer_area
    if net_area <= 0.0 or vapor_density_si <= 0.0:
        actual_vapor_velocity = 0.0
        fs_factor = 0.0
        percent_flood = 0.0
    else:
        vapor_volume_flow = vapor_mass_flow_si / vapor_density_si
        actual_vapor_velocity = vapor_volume_flow / net_area
        fs_factor = actual_vapor_velocity * math.sqrt(vapor_density_si)
        percent_flood = (
            (actual_vapor_velocity / flooding_velocity) * 100.0 if flooding_velocity > 0.0 else 0.0
        )

    # ---- `calculateWeepingCheck`.
    liquid_vol_flow = liquid_mass_flow_si / liquid_density_si

    def crest(weir: float) -> float:
        if weir > 0.0:
            return 750.0 * math.pow(liquid_vol_flow / weir, 2.0 / 3.0)
        return 0.0

    actual_hole_velocity = (
        (vapor_mass_flow_si / vapor_density_si) / hole_area
        if hole_area > 0.0 and vapor_density_si > 0.0
        else 0.0
    )
    # The class returns before the check for a bubble-cap tray, leaving `turndownRatio` at its
    # field initialiser - which is what the captured row's `0.0` is.
    if bubble_cap:
        minimum_vapor_velocity = 0.0
        weeping_ok = True
        turndown_written = None
    else:
        how_mm = crest(calculated_weir_length)
        hw_plus_how = weir_height_m * 1000.0 + how_mm
        kw = (
            WEEPING_CHART_BASE
            if hw_plus_how < 5.0
            else WEEPING_CHART_BASE + WEEPING_CHART_SLOPE * math.sqrt(hw_plus_how)
        )
        numerator = kw - WEEPING_HOLE_SLOPE * (WEEPING_REFERENCE_HOLE_MM - hole_diameter_mm)
        if numerator < 0.0:
            numerator = WEEPING_NUMERATOR_FLOOR
        u_min_hole = numerator / math.sqrt(max(vapor_density_si, 0.01))
        minimum_vapor_velocity = (
            u_min_hole * hole_area / (total_area - downcommer_area) if hole_area > 0.0 else 0.0
        )
        weeping_ok = actual_hole_velocity >= u_min_hole
        turndown_written = actual_hole_velocity / u_min_hole if u_min_hole > 0.0 else 1.0

    # ---- `calculateEntrainment`.
    if percent_flood / 100.0 <= 0.0:
        entrainment = 0.0
    else:
        flv_entrainment = 0.0
        if vapor_mass_flow_si > 0.0:
            flv_entrainment = (liquid_mass_flow_si / vapor_mass_flow_si) * math.sqrt(
                vapor_density_si / liquid_density_si
            )
        if flv_entrainment < ENTRAINMENT_LOW_FLV:
            a, b = ENTRAINMENT_LOW
        elif flv_entrainment < ENTRAINMENT_HIGH_FLV:
            a, b = ENTRAINMENT_MID
        else:
            a, b = ENTRAINMENT_HIGH
        entrainment = min(a * math.pow(percent_flood / 100.0, b), 1.0)
    entrainment_ok = entrainment < ENTRAINMENT_LIMIT

    # ---- `calculatePressureDrop`.
    if sieve:
        orifice_coefficient = ORIFICE_SIEVE
    elif valve:
        orifice_coefficient = ORIFICE_VALVE
    else:
        orifice_coefficient = ORIFICE_OTHER
    if orifice_coefficient > 0.0 and liquid_density_si > 0.0:
        dry_mm = (
            51.0
            * math.pow(actual_hole_velocity / orifice_coefficient, 2)
            * (vapor_density_si / liquid_density_si)
        )
    else:
        dry_mm = 0.0
    dry_tray_pressure_drop = dry_mm * liquid_density_si * G / 1000.0
    h_liquid_mm = weir_height_m * 1000.0 + crest(calculated_weir_length)
    liquid_head_pressure_drop = h_liquid_mm * liquid_density_si * G / 1000.0
    residual_mm = (
        6.0 * surface_tension_si * 1000.0 / (liquid_density_si * G * hole_diameter_mm / 1000.0)
        if hole_diameter_mm > 0.0 and liquid_density_si > 0.0
        else 0.0
    )
    residual_head_pressure_drop = residual_mm * liquid_density_si * G / 1000.0
    total_tray_pressure_drop = (
        dry_tray_pressure_drop + liquid_head_pressure_drop + residual_head_pressure_drop
    )

    # ---- `calculateDowncommerBackup`.
    ht_liquid = total_tray_pressure_drop / (liquid_density_si * G)
    area_apron = max(calculated_weir_length * APRON_GAP_FRACTION, APRON_AREA_FLOOR_M2)
    apron_flow = liquid_mass_flow_si / max(liquid_density_si, 1.0)
    hap_m = APRON_LOSS_COEFFICIENT * math.pow(apron_flow / area_apron, 2) / 1000.0
    liquid_in_downcomer = (weir_height_m * 1000.0 + crest(calculated_weir_length)) / 1000.0
    downcommer_backup = ht_liquid + liquid_in_downcomer + hap_m
    max_allowed = spacing + weir_height_m
    downcommer_backup_fraction = downcommer_backup / max_allowed if max_allowed > 0.0 else 0.0
    downcommer_backup_ok = downcommer_backup < BACKUP_LIMIT_FRACTION * max_allowed

    # ---- `calculateTrayEfficiency`. O'Connell takes the viscosity in centipoise.
    alpha_mu = min(
        max(volatility * liquid_viscosity_si * 1000.0, OCONNELL_ALPHA_MU_MIN),
        OCONNELL_ALPHA_MU_MAX,
    )
    eo_percent = min(
        max(OCONNELL_BASE - OCONNELL_SLOPE * math.log10(alpha_mu), OCONNELL_EFFICIENCY_MIN),
        OCONNELL_EFFICIENCY_MAX,
    )
    tray_efficiency = eo_percent / 100.0

    # ---- `calculateTurndownRatio`. **The same number as the weeping check's ratio** wherever
    # both are defined, so it only bites where `minimum_vapor_velocity` is zero - which is a
    # bubble-cap tray, and there the field keeps its initialiser.
    if turndown_written is None:
        turndown_ratio = 0.0
    elif actual_vapor_velocity > 0.0 and minimum_vapor_velocity > 0.0:
        turndown_ratio = actual_vapor_velocity / minimum_vapor_velocity
    else:
        turndown_ratio = turndown_written

    # ---- `assessDesign`. Five conditions.
    design_ok = (
        weeping_ok
        and entrainment_ok
        and downcommer_backup_ok
        and DESIGN_FLOOD_MIN <= percent_flood <= DESIGN_FLOOD_MAX
    )

    return TrayHydraulicsResult(
        flooding_velocity=from_si(flooding_velocity, "m/s"),
        actual_vapor_velocity=from_si(actual_vapor_velocity, "m/s"),
        percent_flood=percent_flood,
        minimum_vapor_velocity=from_si(minimum_vapor_velocity, "m/s"),
        fs_factor=fs_factor,
        weeping_ok=weeping_ok,
        entrainment=entrainment,
        entrainment_ok=entrainment_ok,
        downcommer_backup=from_si(downcommer_backup, "m"),
        downcommer_backup_fraction=downcommer_backup_fraction,
        downcommer_backup_ok=downcommer_backup_ok,
        total_tray_pressure_drop=from_si(total_tray_pressure_drop, "Pa"),
        total_tray_pressure_drop_mbar=total_tray_pressure_drop / 100.0,
        dry_tray_pressure_drop=from_si(dry_tray_pressure_drop, "Pa"),
        liquid_head_pressure_drop=from_si(liquid_head_pressure_drop, "Pa"),
        residual_head_pressure_drop=from_si(residual_head_pressure_drop, "Pa"),
        tray_efficiency=tray_efficiency,
        turndown_ratio=turndown_ratio,
        calculated_weir_length=from_si(calculated_weir_length, "m"),
        active_area=from_si(active_area, "m**2"),
        total_area=from_si(total_area, "m**2"),
        hole_area=from_si(hole_area, "m**2"),
        downcommer_area=from_si(downcommer_area, "m**2"),
        design_ok=design_ok,
        warnings=tuple(warnings),
    )
