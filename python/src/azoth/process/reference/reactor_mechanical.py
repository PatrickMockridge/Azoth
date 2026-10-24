"""``ReactorMechanicalDesign.calcDesign``: the vessel around a ``GibbsReactor``.

The Python twin of ``crates/azoth-process/src/reactor/mechanical.rs``. Every step mirrors the Rust
line for line, and **the two halves are held to the same capture**: the Rust module's own report
tests and ``python/tests/test_reactor_mechanical_report.py`` compare the same three rows of
``validation/neqsim/captures/process_reactor_mechanical_design.tsv``.

**The class reads the inlet stream and nothing else.** ``GibbsReactor`` is a ``TwoPortEquipment``,
so the design takes ``getInletStream()``'s pressure, temperature, volumetric flow and density, and
from those sizes a catalyst volume, a vessel, a shell, an Ergun bed drop and a weight. The
reactor's own solve is not on this path.

**The two stream reads are system quantities, not phase ones.** ``getFlowRate("m3/hr")`` and
``getThermoSystem().getDensity("kg/m3")`` are :func:`system_volumetric_flow` and
:func:`system_mass_density`'s own definitions, so this design follows the class onto a two-phase
feed where a single-phase cubic volume would refuse.

**Three branches cannot answer, and the capture measures it.** The L/D clamps are unreachable
because the diameter is built at the target ``4.0`` and the length is ``4.0 * D``; the head
formula divides by a larger denominator than the shell's at every positive design pressure, so
``max(head, shell)`` always reports the shell; and the catalyst and ``0.5`` m diameter floors
answer on every default-GHSV state. The capture's CSTR row states an LHSV that lifts the diameter
off its floor so the cube root is measured rather than clamped.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from azoth.core.errors import InvalidInputError
from azoth.process.reference.capacity import system_mass_density, system_volumetric_flow

#: ``reactorType``'s own initialiser. **Only ``CSTR`` is read**, as the class's own
#: ``"CSTR".equals(reactorType)`` - a case-sensitive comparison this reproduces.
DEFAULT_REACTOR_TYPE = "FIXED_BED"
#: ``catalystBulkDensity``'s own initialiser, kg/m**3.
DEFAULT_CATALYST_BULK_DENSITY = 800.0
#: ``catalystParticleDiameter``'s own initialiser, mm.
DEFAULT_CATALYST_PARTICLE_DIAMETER_MM = 3.0
#: ``bedVoidFraction``'s own initialiser.
DEFAULT_BED_VOID_FRACTION = 0.4
#: ``ghsv``'s own initialiser, 1/hr.
DEFAULT_GHSV_PER_HOUR = 5000.0
#: ``lhsv``'s own initialiser, 1/hr.
DEFAULT_LHSV_PER_HOUR = 1.0
#: ``numberOfBeds``'s own initialiser.
DEFAULT_NUMBER_OF_BEDS = 1
#: ``MechanicalDesign.maxOperationPressure``'s own initialiser, bara. The equipment path would
#: overwrite it from an economic design pressure this port does not carry.
DEFAULT_MAX_OPERATION_PRESSURE_BARA = 100.0
#: ``MechanicalDesign.maxOperationTemperature``'s own initialiser, K.
DEFAULT_MAX_OPERATION_TEMPERATURE_K = 100.0
#: ``MechanicalDesign.jointEfficiency``'s own initialiser, fully radiographed.
JOINT_EFFICIENCY = 1.0
#: ``MechanicalDesign.corrosionAllowance``'s own initialiser, mm.
CORROSION_ALLOWANCE_MM = 0.0
#: ``designPressureMargin``: the design pressure over the operating pressure.
DESIGN_PRESSURE_MARGIN = 1.1
#: ``designTemperatureMarginC``: degrees Celsius above the operating temperature.
DESIGN_TEMPERATURE_MARGIN_C = 25.0
#: ``bedFillFraction``: the catalyst volume over the vessel volume.
BED_FILL_FRACTION = 0.7
#: ``targetLD``: the length over the diameter the vessel is built at.
TARGET_LD_RATIO = 4.0
#: ``minLDRatio``, unreachable - see the module doc.
MIN_LD_RATIO = 2.0
#: ``maxLDRatio``, unreachable - see the module doc.
MAX_LD_RATIO = 8.0
#: The floor ``calcDesign`` takes the catalyst volume at, m**3.
MINIMUM_CATALYST_VOLUME_M3 = 0.01
#: The floor ``calcDesign`` takes the vessel diameter at, m.
MINIMUM_VESSEL_DIAMETER_M = 0.5
#: ``gasViscosity``, **hard-coded in the class** rather than read from the fluid, Pa*s.
GAS_VISCOSITY_PA_S = 1.5e-5
#: ``allowableStressMPa``, SA-516-70.
ALLOWABLE_STRESS_MPA = 137.9
#: The floor ``calcDesign`` takes the shell thickness at, mm.
MINIMUM_SHELL_THICKNESS_MM = 6.0
#: ``steelDensity``, kg/m**3.
STEEL_DENSITY_KG_PER_M3 = 7850.0
#: ``distributionPlateThickness``, mm.
DISTRIBUTION_PLATE_THICKNESS_MM = 12.0
#: The 2:1 ellipsoidal head's area factor over a hemisphere.
ELLIPSOIDAL_HEAD_AREA_FACTOR = 1.084
#: The nozzle weight as a fraction of the empty vessel's.
NOZZLE_WEIGHT_FRACTION = 0.1


@dataclass(frozen=True, slots=True)
class ReactorMechanicalGeometry:
    """``calcDesign``'s own configuration, every field the class's setters reach."""

    #: ``reactorType``. **Only ``CSTR`` is read**, case-sensitively.
    reactor_type: str = DEFAULT_REACTOR_TYPE
    #: ``catalystBulkDensity``, kg/m**3.
    catalyst_bulk_density: float = DEFAULT_CATALYST_BULK_DENSITY
    #: ``catalystParticleDiameter``, mm.
    catalyst_particle_diameter_mm: float = DEFAULT_CATALYST_PARTICLE_DIAMETER_MM
    #: ``bedVoidFraction``.
    bed_void_fraction: float = DEFAULT_BED_VOID_FRACTION
    #: ``ghsv``, 1/hr.
    ghsv_per_hour: float = DEFAULT_GHSV_PER_HOUR
    #: ``lhsv``, 1/hr.
    lhsv_per_hour: float = DEFAULT_LHSV_PER_HOUR
    #: ``numberOfBeds``.
    number_of_beds: int = DEFAULT_NUMBER_OF_BEDS
    #: ``getMaxOperationPressure()``, bara.
    max_operation_pressure_bara: float = DEFAULT_MAX_OPERATION_PRESSURE_BARA
    #: ``getMaxOperationTemperature()``, K.
    max_operation_temperature_k: float = DEFAULT_MAX_OPERATION_TEMPERATURE_K
    #: ``getJointEfficiency()``.
    joint_efficiency: float = JOINT_EFFICIENCY
    #: ``getCorrosionAllowance()``, mm.
    corrosion_allowance_mm: float = CORROSION_ALLOWANCE_MM


@dataclass(frozen=True, slots=True)
class ReactorMechanicalReport:
    """``calcDesign``'s own report."""

    #: ``getVesselDiameter()``, m.
    vessel_diameter_m: float
    #: ``getVesselLength()``, m.
    vessel_length_m: float
    #: ``getShellThickness()``, mm.
    shell_thickness_mm: float
    #: ``getHeadThickness()``, mm - **always the shell's**, see the module doc.
    head_thickness_mm: float
    #: ``getDesignPressureBara()``.
    design_pressure_bara: float
    #: ``getDesignTemperatureC()``, the class's own Celsius.
    design_temperature_c: float
    #: ``getCatalystVolume()``, m**3.
    catalyst_volume_m3: float
    #: ``getCatalystMass()``, kg.
    catalyst_mass_kg: float
    #: ``getBedPressureDrop()``, bar.
    bed_pressure_drop_bar: float
    #: ``getEmptyVesselWeight()``, kg.
    empty_vessel_weight_kg: float
    #: ``getTotalEquippedWeight()``, kg.
    total_equipped_weight_kg: float
    #: ``getReactorType()``, the class's own spelling.
    reactor_type: str
    #: ``getNumberOfBeds()``.
    number_of_beds: int


def reactor_mechanical_design(
    components: list[str],
    t: float,
    p: float,
    z: list[float],
    n: float,
    geometry: ReactorMechanicalGeometry,
) -> ReactorMechanicalReport:
    """``calcDesign``, over the reactor's inlet state.

    Args:
        components: the inlet's substances.
        t: the inlet's temperature, K.
        p: the inlet's pressure, Pa.
        z: the inlet's composition.
        n: the inlet's molar flow, mol/s.
        geometry: the class's configuration.

    Raises:
        InvalidInputError: where a configuration value makes the arithmetic undefined - a
            non-positive bulk density, particle diameter, joint efficiency or design pressure, a
            void fraction outside ``(0, 1)``, no bed at all, or a negative corrosion allowance -
            and whatever the inlet's flash, density or volumetric flow raise.
    """
    _validate(geometry)

    # ---- What it reads. Both are the *system* quantities the class asks for.
    feed_flow_m3s = system_volumetric_flow(components, t, p, z, n)
    gas_density = system_mass_density(components, t, p, z)

    # ---- The design conditions, each a max with the class's own operation limit.
    inlet_pressure_bara = p / 1.0e5
    operating_pressure_bara = max(inlet_pressure_bara, geometry.max_operation_pressure_bara)
    inlet_temperature_c = t - 273.15
    operation_limit_c = (
        geometry.max_operation_temperature_k - 273.15
        if geometry.max_operation_temperature_k > 0.0
        else 0.0
    )
    operating_temperature_c = max(inlet_temperature_c, operation_limit_c)
    design_pressure_bara = operating_pressure_bara * DESIGN_PRESSURE_MARGIN
    design_temperature_c = operating_temperature_c + DESIGN_TEMPERATURE_MARGIN_C

    # ---- The catalyst, from the space velocity that reactor type selects.
    feed_flow_m3hr = feed_flow_m3s * 3600.0
    if geometry.reactor_type == "CSTR":
        space_velocity = max(geometry.lhsv_per_hour, 0.01)
    else:
        space_velocity = max(geometry.ghsv_per_hour, 1.0)
    catalyst_volume = max(feed_flow_m3hr / space_velocity, MINIMUM_CATALYST_VOLUME_M3)
    catalyst_mass = catalyst_volume * geometry.catalyst_bulk_density

    # ---- The vessel, at the target L/D and then the clamps that cannot fire.
    total_volume = catalyst_volume / BED_FILL_FRACTION
    vessel_diameter = max(
        (4.0 * total_volume / (math.pi * TARGET_LD_RATIO)) ** (1.0 / 3.0),
        MINIMUM_VESSEL_DIAMETER_M,
    )
    vessel_length = TARGET_LD_RATIO * vessel_diameter
    ld_ratio = vessel_length / vessel_diameter
    if ld_ratio < MIN_LD_RATIO:
        vessel_length = vessel_diameter * MIN_LD_RATIO
    elif ld_ratio > MAX_LD_RATIO:
        vessel_diameter = vessel_length / MAX_LD_RATIO

    # ---- The Ergun bed drop, over the bed's own length inside the vessel.
    cross_section = math.pi * vessel_diameter * vessel_diameter / 4.0
    bed_length = catalyst_volume / cross_section
    velocity = feed_flow_m3s / cross_section
    particle_diameter_m = geometry.catalyst_particle_diameter_mm / 1000.0
    eps = geometry.bed_void_fraction
    viscous = (
        150.0
        * GAS_VISCOSITY_PA_S
        * (1.0 - eps) ** 2
        * velocity
        / (particle_diameter_m * particle_diameter_m * eps**3)
    )
    inertial = (
        1.75 * gas_density * (1.0 - eps) * velocity * velocity / (particle_diameter_m * eps**3)
    )
    bed_pressure_drop_bar = (viscous + inertial) * bed_length / 1.0e5

    # ---- The shell and the head, ASME VIII Div 1, in the class's own mm.
    p_mpa = design_pressure_bara * 0.1
    radius_m = vessel_diameter / 2.0
    e = geometry.joint_efficiency
    ca_mm = geometry.corrosion_allowance_mm
    shell_thickness_mm = max(
        (p_mpa * radius_m * 1000.0) / (ALLOWABLE_STRESS_MPA * e - 0.6 * p_mpa) + ca_mm,
        MINIMUM_SHELL_THICKNESS_MM,
    )
    head_thickness_mm = max(
        (p_mpa * vessel_diameter * 1000.0) / (2.0 * ALLOWABLE_STRESS_MPA * e - 0.2 * p_mpa) + ca_mm,
        shell_thickness_mm,
    )

    # ---- The weights.
    shell_area = math.pi * vessel_diameter * vessel_length
    head_area = 2.0 * math.pi * (vessel_diameter / 2.0) ** 2 * ELLIPSOIDAL_HEAD_AREA_FACTOR
    empty_vessel_weight = (
        (shell_area + head_area) * (shell_thickness_mm / 1000.0) * STEEL_DENSITY_KG_PER_M3
    )
    internals_weight = (
        geometry.number_of_beds
        * cross_section
        * (DISTRIBUTION_PLATE_THICKNESS_MM / 1000.0)
        * STEEL_DENSITY_KG_PER_M3
        * 2.0
    )
    nozzle_weight = empty_vessel_weight * NOZZLE_WEIGHT_FRACTION
    total_equipped_weight = empty_vessel_weight + catalyst_mass + internals_weight + nozzle_weight

    return ReactorMechanicalReport(
        vessel_diameter_m=vessel_diameter,
        vessel_length_m=vessel_length,
        shell_thickness_mm=shell_thickness_mm,
        head_thickness_mm=head_thickness_mm,
        design_pressure_bara=design_pressure_bara,
        design_temperature_c=design_temperature_c,
        catalyst_volume_m3=catalyst_volume,
        catalyst_mass_kg=catalyst_mass,
        bed_pressure_drop_bar=bed_pressure_drop_bar,
        empty_vessel_weight_kg=empty_vessel_weight,
        total_equipped_weight_kg=total_equipped_weight,
        reactor_type=geometry.reactor_type,
        number_of_beds=geometry.number_of_beds,
    )


def _positive(value: float, name: str) -> None:
    if not (math.isfinite(value) and value > 0.0):
        raise InvalidInputError(name, f"{name} is {value}; the design divides by it")


def _validate(geometry: ReactorMechanicalGeometry) -> None:
    """The eight values whose arithmetic is undefined rather than merely odd.

    **The two space velocities are absent deliberately**: the class floors them
    (``max(ghsv, 1.0)``, ``max(lhsv, 0.01)``), and the floor is the behaviour this reproduces.
    """
    _positive(geometry.catalyst_bulk_density, "catalyst_bulk_density")
    _positive(geometry.catalyst_particle_diameter_mm, "catalyst_particle_diameter")
    _positive(geometry.joint_efficiency, "joint_efficiency")
    _positive(geometry.max_operation_pressure_bara, "max_operation_pressure")
    _positive(geometry.max_operation_temperature_k, "max_operation_temperature")
    if not (0.0 < geometry.bed_void_fraction < 1.0):
        raise InvalidInputError(
            "bed_void_fraction",
            f"a bed void fraction of {geometry.bed_void_fraction} is not in (0, 1), and the "
            f"Ergun terms divide by its cube",
        )
    if geometry.number_of_beds == 0:
        raise InvalidInputError(
            "number_of_beds", "a reactor with no catalyst bed has no bed weight to carry"
        )
    if not (
        math.isfinite(geometry.corrosion_allowance_mm) and geometry.corrosion_allowance_mm >= 0.0
    ):
        raise InvalidInputError(
            "corrosion_allowance",
            f"a corrosion allowance of {geometry.corrosion_allowance_mm} mm is negative",
        )
