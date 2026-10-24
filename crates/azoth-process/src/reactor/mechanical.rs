//! `ReactorMechanicalDesign.calcDesign`: the vessel around a `GibbsReactor`.
//!
//! **The class reads the inlet stream and nothing else.** `GibbsReactor` is a
//! `TwoPortEquipment`, so the design takes `getInletStream()`'s pressure, temperature,
//! volumetric flow and density, and from those sizes a catalyst volume, a vessel, a shell, an
//! Ergun bed drop and a weight. The reactor's own solve is not on this path.
//!
//! **The two stream reads are system quantities, not phase ones.**
//! `getFlowRate("m3/hr")` and `getThermoSystem().getDensity("kg/m3")` are
//! [`system_volumetric_flow`] and [`system_mass_density`]'s own definitions, so this design
//! follows the class onto a two-phase feed where [`Stream::density`] would refuse - the cubic
//! volume of one phase is not what the class asks for.
//!
//! **Three branches cannot answer, and `process_reactor_mechanical_design.tsv` measures it.**
//! The L/D clamps are unreachable because the diameter is built at the target `4.0` and the
//! length is `4.0 * D`; the head formula divides by a larger denominator than the shell's at
//! every positive design pressure, so `max(head, shell)` always reports the shell; and the
//! catalyst floor and the `0.5` m diameter floor both answer on every default-GHSV state. The
//! capture's CSTR row states an LHSV that lifts the diameter off its floor so the cube root is
//! measured rather than clamped.
//!
//! **`super.calcDesign()` is not ported.** It calls `setCompanySpecificDesignStandards` and
//! `readDesignSpecifications`, which would overwrite the joint efficiency and the corrosion
//! allowance from a design-standard tree this workspace does not carry; the capture shows the
//! two at the base class's own `1.0` and `0.0` at the pin.

use std::f64::consts::PI;

use azoth_core::units::{Length, Mass, Volume, cubic_meters, kilograms, meters};
use azoth_core::{AzothError, Result};

use crate::segment::phase::{system_mass_density, system_volumetric_flow};
use crate::stream::Stream;

/// `reactorType`'s own initialiser. **Only `CSTR` is read**, as the class's own
/// `"CSTR".equals(reactorType)` - a case-sensitive comparison this reproduces.
pub const DEFAULT_REACTOR_TYPE: &str = "FIXED_BED";
/// `catalystBulkDensity`'s own initialiser, kg/m3.
pub const DEFAULT_CATALYST_BULK_DENSITY: f64 = 800.0;
/// `catalystParticleDiameter`'s own initialiser, mm.
pub const DEFAULT_CATALYST_PARTICLE_DIAMETER_MM: f64 = 3.0;
/// `bedVoidFraction`'s own initialiser.
pub const DEFAULT_BED_VOID_FRACTION: f64 = 0.4;
/// `ghsv`'s own initialiser, 1/hr.
pub const DEFAULT_GHSV_PER_HOUR: f64 = 5000.0;
/// `lhsv`'s own initialiser, 1/hr.
pub const DEFAULT_LHSV_PER_HOUR: f64 = 1.0;
/// `numberOfBeds`'s own initialiser.
pub const DEFAULT_NUMBER_OF_BEDS: usize = 1;
/// `MechanicalDesign.maxOperationPressure`'s own initialiser, bara. The equipment path would
/// overwrite it from an economic design pressure this port does not carry.
pub const DEFAULT_MAX_OPERATION_PRESSURE_BARA: f64 = 100.0;
/// `MechanicalDesign.maxOperationTemperature`'s own initialiser, K.
pub const DEFAULT_MAX_OPERATION_TEMPERATURE_K: f64 = 100.0;
/// `MechanicalDesign.jointEfficiency`'s own initialiser, fully radiographed.
pub const JOINT_EFFICIENCY: f64 = 1.0;
/// `MechanicalDesign.corrosionAllowance`'s own initialiser, mm.
pub const CORROSION_ALLOWANCE_MM: f64 = 0.0;
/// `designPressureMargin`: the design pressure over the operating pressure.
pub const DESIGN_PRESSURE_MARGIN: f64 = 1.1;
/// `designTemperatureMarginC`: degrees Celsius above the operating temperature.
pub const DESIGN_TEMPERATURE_MARGIN_C: f64 = 25.0;
/// `bedFillFraction`: the catalyst volume over the vessel volume.
pub const BED_FILL_FRACTION: f64 = 0.7;
/// `targetLD`: the length over the diameter the vessel is built at.
pub const TARGET_LD_RATIO: f64 = 4.0;
/// `minLDRatio`, unreachable - see the module doc.
pub const MIN_LD_RATIO: f64 = 2.0;
/// `maxLDRatio`, unreachable - see the module doc.
pub const MAX_LD_RATIO: f64 = 8.0;
/// The floor `calcDesign` takes the catalyst volume at, m3.
pub const MINIMUM_CATALYST_VOLUME_M3: f64 = 0.01;
/// The floor `calcDesign` takes the vessel diameter at, m.
pub const MINIMUM_VESSEL_DIAMETER_M: f64 = 0.5;
/// `gasViscosity`, **hard-coded in the class** rather than read from the fluid, Pa·s.
pub const GAS_VISCOSITY_PA_S: f64 = 1.5e-5;
/// `allowableStressMPa`, SA-516-70.
pub const ALLOWABLE_STRESS_MPA: f64 = 137.9;
/// The floor `calcDesign` takes the shell thickness at, mm.
pub const MINIMUM_SHELL_THICKNESS_MM: f64 = 6.0;
/// `steelDensity`, kg/m3.
pub const STEEL_DENSITY_KG_PER_M3: f64 = 7850.0;
/// `distributionPlateThickness`, mm.
pub const DISTRIBUTION_PLATE_THICKNESS_MM: f64 = 12.0;
/// The 2:1 ellipsoidal head's area factor over a hemisphere.
pub const ELLIPSOIDAL_HEAD_AREA_FACTOR: f64 = 1.084;
/// The nozzle weight as a fraction of the empty vessel's.
pub const NOZZLE_WEIGHT_FRACTION: f64 = 0.1;

/// `calcDesign`'s own configuration, every field the class's setters reach.
#[derive(Debug, Clone, PartialEq)]
pub struct ReactorMechanicalGeometry {
    /// `reactorType`. **Only `CSTR` is read**, case-sensitively.
    pub reactor_type: String,
    /// `catalystBulkDensity`, kg/m3.
    pub catalyst_bulk_density: f64,
    /// `catalystParticleDiameter`, mm.
    pub catalyst_particle_diameter_mm: f64,
    /// `bedVoidFraction`.
    pub bed_void_fraction: f64,
    /// `ghsv`, 1/hr.
    pub ghsv_per_hour: f64,
    /// `lhsv`, 1/hr.
    pub lhsv_per_hour: f64,
    /// `numberOfBeds`.
    pub number_of_beds: usize,
    /// `getMaxOperationPressure()`, bara.
    pub max_operation_pressure_bara: f64,
    /// `getMaxOperationTemperature()`, K.
    pub max_operation_temperature_k: f64,
    /// `getJointEfficiency()`.
    pub joint_efficiency: f64,
    /// `getCorrosionAllowance()`, mm.
    pub corrosion_allowance_mm: f64,
}

impl Default for ReactorMechanicalGeometry {
    fn default() -> Self {
        Self {
            reactor_type: DEFAULT_REACTOR_TYPE.to_string(),
            catalyst_bulk_density: DEFAULT_CATALYST_BULK_DENSITY,
            catalyst_particle_diameter_mm: DEFAULT_CATALYST_PARTICLE_DIAMETER_MM,
            bed_void_fraction: DEFAULT_BED_VOID_FRACTION,
            ghsv_per_hour: DEFAULT_GHSV_PER_HOUR,
            lhsv_per_hour: DEFAULT_LHSV_PER_HOUR,
            number_of_beds: DEFAULT_NUMBER_OF_BEDS,
            max_operation_pressure_bara: DEFAULT_MAX_OPERATION_PRESSURE_BARA,
            max_operation_temperature_k: DEFAULT_MAX_OPERATION_TEMPERATURE_K,
            joint_efficiency: JOINT_EFFICIENCY,
            corrosion_allowance_mm: CORROSION_ALLOWANCE_MM,
        }
    }
}

/// `calcDesign`'s own report.
#[derive(Debug, Clone, PartialEq)]
pub struct ReactorMechanicalReport {
    /// `getVesselDiameter()`, m.
    pub vessel_diameter: Length,
    /// `getVesselLength()`, m.
    pub vessel_length: Length,
    /// `getShellThickness()`, mm.
    pub shell_thickness_mm: f64,
    /// `getHeadThickness()`, mm - **always the shell's**, see the module doc.
    pub head_thickness_mm: f64,
    /// `getDesignPressureBara()`.
    pub design_pressure_bara: f64,
    /// `getDesignTemperatureC()`, the class's own Celsius.
    pub design_temperature_c: f64,
    /// `getCatalystVolume()`, m3.
    pub catalyst_volume: Volume,
    /// `getCatalystMass()`, kg.
    pub catalyst_mass: Mass,
    /// `getBedPressureDrop()`, bar.
    pub bed_pressure_drop_bar: f64,
    /// `getEmptyVesselWeight()`, kg.
    pub empty_vessel_weight: Mass,
    /// `getTotalEquippedWeight()`, kg.
    pub total_equipped_weight: Mass,
    /// `getReactorType()`, the class's own spelling.
    pub reactor_type: String,
    /// `getNumberOfBeds()`.
    pub number_of_beds: usize,
}

/// `ReactorMechanicalDesign.calcDesign`, over the reactor's inlet stream.
///
/// # Errors
/// * [`AzothError::InvalidInput`] where a configuration value makes the arithmetic undefined -
///   a non-positive bulk density, particle diameter, joint efficiency or design pressure, a
///   void fraction outside `(0, 1)`, no bed at all, or a negative corrosion allowance;
/// * whatever the inlet's flash, density or volumetric flow raise.
pub fn reactor_mechanical_design(
    feed: &Stream,
    geometry: &ReactorMechanicalGeometry,
) -> Result<ReactorMechanicalReport> {
    validate(geometry)?;

    // ---- What it reads. Both are the *system* quantities the class asks for.
    let feed_flow_m3s = system_volumetric_flow(feed)?;
    let gas_density = system_mass_density(feed)?;

    // ---- The design conditions, each a max with the class's own operation limit.
    let inlet_pressure_bara = feed.p.value / 1.0e5;
    let operating_pressure_bara = inlet_pressure_bara.max(geometry.max_operation_pressure_bara);
    let inlet_temperature_c = feed.t.value - 273.15;
    let operation_limit_c = if geometry.max_operation_temperature_k > 0.0 {
        geometry.max_operation_temperature_k - 273.15
    } else {
        0.0
    };
    let operating_temperature_c = inlet_temperature_c.max(operation_limit_c);
    let design_pressure_bara = operating_pressure_bara * DESIGN_PRESSURE_MARGIN;
    let design_temperature_c = operating_temperature_c + DESIGN_TEMPERATURE_MARGIN_C;

    // ---- The catalyst, from the space velocity that reactor type selects.
    let feed_flow_m3hr = feed_flow_m3s * 3600.0;
    let space_velocity = if geometry.reactor_type == "CSTR" {
        geometry.lhsv_per_hour.max(0.01)
    } else {
        geometry.ghsv_per_hour.max(1.0)
    };
    let catalyst_volume = (feed_flow_m3hr / space_velocity).max(MINIMUM_CATALYST_VOLUME_M3);
    let catalyst_mass = catalyst_volume * geometry.catalyst_bulk_density;

    // ---- The vessel, at the target L/D and then the clamps that cannot fire.
    let total_volume = catalyst_volume / BED_FILL_FRACTION;
    let mut vessel_diameter = (4.0 * total_volume / (PI * TARGET_LD_RATIO))
        .cbrt()
        .max(MINIMUM_VESSEL_DIAMETER_M);
    let mut vessel_length = TARGET_LD_RATIO * vessel_diameter;
    let ld_ratio = vessel_length / vessel_diameter;
    if ld_ratio < MIN_LD_RATIO {
        vessel_length = vessel_diameter * MIN_LD_RATIO;
    } else if ld_ratio > MAX_LD_RATIO {
        vessel_diameter = vessel_length / MAX_LD_RATIO;
    }

    // ---- The Ergun bed drop, over the bed's own length inside the vessel.
    let cross_section = PI * vessel_diameter * vessel_diameter / 4.0;
    let bed_length = catalyst_volume / cross_section;
    let velocity = feed_flow_m3s / cross_section;
    let particle_diameter_m = geometry.catalyst_particle_diameter_mm / 1000.0;
    let eps = geometry.bed_void_fraction;
    let viscous = 150.0 * GAS_VISCOSITY_PA_S * (1.0 - eps).powi(2) * velocity
        / (particle_diameter_m * particle_diameter_m * eps.powi(3));
    let inertial = 1.75 * gas_density * (1.0 - eps) * velocity * velocity
        / (particle_diameter_m * eps.powi(3));
    let bed_pressure_drop_bar = (viscous + inertial) * bed_length / 1.0e5;

    // ---- The shell and the head, ASME VIII Div 1, in the class's own mm.
    let p_mpa = design_pressure_bara * 0.1;
    let radius_m = vessel_diameter / 2.0;
    let e = geometry.joint_efficiency;
    let ca_mm = geometry.corrosion_allowance_mm;
    let shell_thickness_mm =
        ((p_mpa * radius_m * 1000.0) / (ALLOWABLE_STRESS_MPA * e - 0.6 * p_mpa) + ca_mm)
            .max(MINIMUM_SHELL_THICKNESS_MM);
    let head_thickness_mm = ((p_mpa * vessel_diameter * 1000.0)
        / (2.0 * ALLOWABLE_STRESS_MPA * e - 0.2 * p_mpa)
        + ca_mm)
        .max(shell_thickness_mm);

    // ---- The weights.
    let shell_area = PI * vessel_diameter * vessel_length;
    let head_area = 2.0 * PI * (vessel_diameter / 2.0).powi(2) * ELLIPSOIDAL_HEAD_AREA_FACTOR;
    let empty_vessel_weight =
        (shell_area + head_area) * (shell_thickness_mm / 1000.0) * STEEL_DENSITY_KG_PER_M3;
    let internals_weight = geometry.number_of_beds as f64
        * cross_section
        * (DISTRIBUTION_PLATE_THICKNESS_MM / 1000.0)
        * STEEL_DENSITY_KG_PER_M3
        * 2.0;
    let nozzle_weight = empty_vessel_weight * NOZZLE_WEIGHT_FRACTION;
    let total_equipped_weight =
        empty_vessel_weight + catalyst_mass + internals_weight + nozzle_weight;

    Ok(ReactorMechanicalReport {
        vessel_diameter: meters(vessel_diameter),
        vessel_length: meters(vessel_length),
        shell_thickness_mm,
        head_thickness_mm,
        design_pressure_bara,
        design_temperature_c,
        catalyst_volume: cubic_meters(catalyst_volume),
        catalyst_mass: kilograms(catalyst_mass),
        bed_pressure_drop_bar,
        empty_vessel_weight: kilograms(empty_vessel_weight),
        total_equipped_weight: kilograms(total_equipped_weight),
        reactor_type: geometry.reactor_type.clone(),
        number_of_beds: geometry.number_of_beds,
    })
}

/// The eight values whose arithmetic is undefined rather than merely odd.
///
/// **The two space velocities are absent deliberately**: the class floors them
/// (`max(ghsv, 1.0)`, `max(lhsv, 0.01)`), and the floor is the behaviour this reproduces.
fn validate(geometry: &ReactorMechanicalGeometry) -> Result<()> {
    let positive = |value: f64, name: &str| -> Result<()> {
        if value.is_finite() && value > 0.0 {
            Ok(())
        } else {
            Err(AzothError::invalid_input(
                name,
                format!("{name} is {value}; the design divides by it"),
            ))
        }
    };
    positive(geometry.catalyst_bulk_density, "catalyst_bulk_density")?;
    positive(
        geometry.catalyst_particle_diameter_mm,
        "catalyst_particle_diameter",
    )?;
    positive(geometry.joint_efficiency, "joint_efficiency")?;
    positive(
        geometry.max_operation_pressure_bara,
        "max_operation_pressure",
    )?;
    positive(
        geometry.max_operation_temperature_k,
        "max_operation_temperature",
    )?;
    if !(geometry.bed_void_fraction > 0.0 && geometry.bed_void_fraction < 1.0) {
        return Err(AzothError::invalid_input(
            "bed_void_fraction",
            format!(
                "a bed void fraction of {} is not in (0, 1), and the Ergun terms divide by its \
                 cube",
                geometry.bed_void_fraction
            ),
        ));
    }
    if geometry.number_of_beds == 0 {
        return Err(AzothError::invalid_input(
            "number_of_beds",
            "a reactor with no catalyst bed has no bed weight to carry",
        ));
    }
    if !(geometry.corrosion_allowance_mm.is_finite() && geometry.corrosion_allowance_mm >= 0.0) {
        return Err(AzothError::invalid_input(
            "corrosion_allowance",
            format!(
                "a corrosion allowance of {} mm is negative",
                geometry.corrosion_allowance_mm
            ),
        ));
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use azoth_core::units::{kelvins, pascals};

    /// A feed at a state, with the probe's `addComponent` amounts scaled to one mole per second.
    fn stream(names: &[&str], amounts: &[f64], temperature_k: f64, pressure_bara: f64) -> Stream {
        let total: f64 = amounts.iter().sum();
        let z: Vec<f64> = amounts.iter().map(|amount| amount / total).collect();
        Stream::from_pt(
            names.iter().map(|name| (*name).to_string()).collect(),
            z,
            1.0,
            pascals(pressure_bara * 1.0e5),
            kelvins(temperature_k),
        )
        .expect("the feed resolves")
    }

    /// **`1e-6`, because the volumetric-flow identity is not exact between the libraries.**
    /// azoth's `system_volumetric_flow` and NeqSim's `getFlowRate("m3/hr")` differ at about
    /// `2e-8` relative on these rows, and the catalyst volume and the bed drop carry it; every
    /// formula error this test exists to catch is orders larger.
    fn close(got: f64, want: f64, what: &str) {
        assert!(
            (got - want).abs() / want.abs() < 1e-6,
            "{what}: {got} against {want}"
        );
    }

    /// **`process_reactor_mechanical_design.tsv`'s `gibbs_fixed_bed_defaults` row**: every class
    /// default, on the ammonia feed, where the catalyst volume is floored at `0.01` m3.
    #[test]
    fn the_defaults_row_floors_the_catalyst_and_sizes_the_shell() {
        let feed = stream(
            &["hydrogen", "nitrogen", "ammonia"],
            &[1.5, 0.5, 0.0],
            450.0,
            300.0,
        );
        let report = reactor_mechanical_design(&feed, &ReactorMechanicalGeometry::default())
            .expect("the design computes");

        assert_eq!(report.reactor_type, "FIXED_BED");
        assert_eq!(report.number_of_beds, 1);
        close(report.vessel_diameter.value, 0.5, "vessel_diameter_m");
        close(report.vessel_length.value, 2.0, "vessel_length_m");
        close(
            report.shell_thickness_mm,
            69.856_054_191_363_24,
            "shell_thickness_mm",
        );
        close(
            report.head_thickness_mm,
            69.856_054_191_363_24,
            "head_thickness_mm",
        );
        close(report.design_pressure_bara, 330.0, "design_pressure_bara");
        close(
            report.design_temperature_c,
            201.850_000_000_000_02,
            "design_temperature_c",
        );
        close(report.catalyst_volume.value, 0.01, "catalyst_volume_m3");
        close(report.catalyst_mass.value, 8.0, "catalyst_mass_kg");
        close(
            report.bed_pressure_drop_bar,
            6.101_079_671_327_595e-7,
            "bed_pressure_drop_bar",
        );
        close(
            report.empty_vessel_weight.value,
            1956.188_578_713_105_2,
            "empty_vessel_weight_kg",
        );
        close(
            report.total_equipped_weight.value,
            2196.799_690_080_435_5,
            "total_equipped_weight_kg",
        );
    }

    /// **The `gibbs_cstr` row**: the LHSV branch, and the one row whose diameter is off the
    /// `0.5` m floor, so the cube root is measured rather than clamped.
    #[test]
    fn the_cstr_row_takes_the_lhsv_branch_and_leaves_the_diameter_floor() {
        let feed = stream(
            &["hydrogen", "nitrogen", "ammonia"],
            &[1.5, 0.5, 0.0],
            450.0,
            300.0,
        );
        let geometry = ReactorMechanicalGeometry {
            reactor_type: "CSTR".to_string(),
            lhsv_per_hour: 1.5,
            ..ReactorMechanicalGeometry::default()
        };
        let report = reactor_mechanical_design(&feed, &geometry).expect("the design computes");

        assert_eq!(report.reactor_type, "CSTR");
        close(
            report.vessel_diameter.value,
            0.538_551_169_195_642_5,
            "vessel_diameter_m",
        );
        close(
            report.vessel_length.value,
            2.154_204_676_782_57,
            "vessel_length_m",
        );
        close(
            report.shell_thickness_mm,
            75.242_119_320_305_68,
            "shell_thickness_mm",
        );
        close(
            report.catalyst_volume.value,
            0.343_501_653_505_556_8,
            "catalyst_volume_m3",
        );
        close(
            report.catalyst_mass.value,
            274.801_322_804_445_4,
            "catalyst_mass_kg",
        );
        close(
            report.bed_pressure_drop_bar,
            1.526_048_947_657_984_3e-5,
            "bed_pressure_drop_bar",
        );
        close(
            report.empty_vessel_weight.value,
            2444.452_627_150_544_2,
            "empty_vessel_weight_kg",
        );
        close(
            report.total_equipped_weight.value,
            3006.615_755_120_348,
            "total_equipped_weight_kg",
        );
    }

    /// **The `gibbs_high_pressure` row**: the feed is at 1 bara and the stated maximum operation
    /// pressure is 150, so `max` takes the limit and the GHSV division is above its floor.
    #[test]
    fn the_high_pressure_row_takes_the_max_branch() {
        let feed = stream(
            &["CO2", "hydrogen", "CO", "water"],
            &[1.0, 1.0, 0.0, 0.0],
            1000.0,
            1.0,
        );
        let geometry = ReactorMechanicalGeometry {
            max_operation_pressure_bara: 150.0,
            ..ReactorMechanicalGeometry::default()
        };
        let report = reactor_mechanical_design(&feed, &geometry).expect("the design computes");

        close(report.vessel_diameter.value, 0.5, "vessel_diameter_m");
        close(report.vessel_length.value, 2.0, "vessel_length_m");
        close(report.design_pressure_bara, 165.0, "design_pressure_bara");
        close(report.design_temperature_c, 751.85, "design_temperature_c");
        close(
            report.shell_thickness_mm,
            32.226_562_5,
            "shell_thickness_mm",
        );
        close(
            report.catalyst_volume.value,
            0.059_878_658_517_044_936,
            "catalyst_volume_m3",
        );
        close(
            report.catalyst_mass.value,
            47.902_926_813_635_95,
            "catalyst_mass_kg",
        );
        close(
            report.bed_pressure_drop_bar,
            0.002_644_322_348_722_190_6,
            "bed_pressure_drop_bar",
        );
        close(
            report.empty_vessel_weight.value,
            902.444_809_164_131_8,
            "empty_vessel_weight_kg",
        );
        close(
            report.total_equipped_weight.value,
            1077.584_470_390_200_7,
            "total_equipped_weight_kg",
        );
    }

    /// **The clamps that cannot fire**, held as a property rather than a number: the target
    /// ratio is inside the class's own bounds, so the built ratio is exactly `4.0`.
    #[test]
    fn the_ld_clamps_are_unreachable() {
        for catalyst in [0.01, 0.1, 1.0, 100.0] {
            let total = catalyst / BED_FILL_FRACTION;
            let diameter = (4.0 * total / (PI * TARGET_LD_RATIO))
                .cbrt()
                .max(MINIMUM_VESSEL_DIAMETER_M);
            let length = TARGET_LD_RATIO * diameter;
            let ratio = length / diameter;
            assert!(ratio >= MIN_LD_RATIO && ratio <= MAX_LD_RATIO);
        }
    }

    /// Every refusal, on the field it names.
    #[test]
    fn the_undefined_configurations_are_refused_by_name() {
        let feed = stream(&["methane"], &[1.0], 300.0, 10.0);
        let cases: [(ReactorMechanicalGeometry, &str); 9] = [
            (
                ReactorMechanicalGeometry {
                    catalyst_bulk_density: 0.0,
                    ..Default::default()
                },
                "catalyst_bulk_density",
            ),
            (
                ReactorMechanicalGeometry {
                    catalyst_particle_diameter_mm: 0.0,
                    ..Default::default()
                },
                "catalyst_particle_diameter",
            ),
            (
                ReactorMechanicalGeometry {
                    bed_void_fraction: 0.0,
                    ..Default::default()
                },
                "bed_void_fraction",
            ),
            (
                ReactorMechanicalGeometry {
                    bed_void_fraction: 1.0,
                    ..Default::default()
                },
                "bed_void_fraction",
            ),
            (
                ReactorMechanicalGeometry {
                    number_of_beds: 0,
                    ..Default::default()
                },
                "number_of_beds",
            ),
            (
                ReactorMechanicalGeometry {
                    joint_efficiency: 0.0,
                    ..Default::default()
                },
                "joint_efficiency",
            ),
            (
                ReactorMechanicalGeometry {
                    corrosion_allowance_mm: -1.0,
                    ..Default::default()
                },
                "corrosion_allowance",
            ),
            (
                ReactorMechanicalGeometry {
                    max_operation_pressure_bara: 0.0,
                    ..Default::default()
                },
                "max_operation_pressure",
            ),
            (
                ReactorMechanicalGeometry {
                    max_operation_temperature_k: 0.0,
                    ..Default::default()
                },
                "max_operation_temperature",
            ),
        ];
        for (geometry, field) in cases {
            let error = reactor_mechanical_design(&feed, &geometry)
                .expect_err("the configuration is refused");
            assert_eq!(error.field(), Some(field), "{field}");
        }
    }
}
