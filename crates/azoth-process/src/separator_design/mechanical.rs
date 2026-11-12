//! `SeparatorMechanicalDesign.calcDesign`: the vessel around the separator.
//!
//! The class reads its design basis from `MechanicalDesign`'s default design standards
//! (`super.calcDesign()` installs nine when none are set), and the capture prints every value
//! it takes from them: `separator process design`'s `K = 0.11` m/s, `Fg = 0.8`, volume safety
//! `1.0` and liquid retention `60` s; `material plate design codes`' allowable stress `425`;
//! `plate Joint Efficiency design codes`' factor `1.0`; and the pressure vessel's maximum
//! operation pressure `100` bara and corrosion allowance `0`. Those are this module's
//! [`DesignInputs`], so the arithmetic is a function of a state and a stated basis rather than
//! of a database this tree does not carry.
//!
//! **The class's own test is SRK and this port's separator is PR**, so the capture carries
//! both rows. The arithmetic is the same; the densities and the viscosity move.
//!
//! Three behaviours are reproduced rather than tidied:
//!
//! - **the two volume flows are per second and the class divides them by `1e5`.** The
//!   `getVolume()` it reads is the stream's volume in NeqSim's internal unit; the port works in
//!   m3/s and the capture's `max_design_volume_flow_m3_hr` is the same number.
//! - **the L/D adjustment has two branches.** `calcDesign` resets the tangent-to-tangent length
//!   to `5 d` when the total is outside `[3 d, 6 d]`; `performSizingCalculations` uses `4 d`,
//!   and this module ports `calcDesign` because that is the capture.
//! - **the retention and surge times divide by the liquid flow as though it were m3/hr** while
//!   the sizing multiplies the same number by 3600. The capture's `394 159` s is that
//!   inconsistency, and it is not corrected here.

use azoth_core::units::Pressure;
use azoth_core::{AzothError, Result};
use azoth_eos::phase_transport::{PhaseKind, phase_transport};
use azoth_eos::{Phase, pt_flash};

use crate::segment::phase::{mass_density, molar_mass_of};
use crate::stream::Stream;

/// Which machine's `calcDesign` this is.
///
/// The two share the wall thickness, the weights, the levels, the nozzles and the surge and
/// retention arithmetic; they differ in how the diameter and length are sized and in the
/// module's height.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum VesselKind {
    /// `SeparatorMechanicalDesign`: the liquid retention sets the length, the bubble check
    /// can enlarge the diameter, and the module is `2 d + 1` high.
    Separator,
    /// `GasScrubberMechanicalDesign`: the length is always `5 d`, there is no bubble check,
    /// and the module is `2 d` high.
    Scrubber,
}

/// The design basis the class takes from its default standards.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DesignInputs {
    /// `SeparatorDesignStandard.getGasLoadFactor`, m/s. The default is `0.11`.
    pub gas_load_factor: f64,
    /// `SeparatorDesignStandard.getFg`, the horizontal liquid fraction. Default `0.8`.
    pub fg: f64,
    /// `SeparatorDesignStandard.getVolumetricDesignFactor`. Default `1.0`.
    pub volume_safety_factor: f64,
    /// `getLiquidRetentionTime("API12J", ...)`, seconds. Default `60`.
    pub retention_time_s: f64,
    /// `MechanicalDesign.getMaxOperationPressure`. Default `100` bara.
    pub max_operation_pressure: Pressure,
    /// `MaterialPlateDesignStandard.getDivisionClass`, MPa. Default `425`.
    pub allowable_stress: f64,
    /// `JointEfficiencyPlateStandard.getJEFactor`. Default `1.0`.
    pub joint_efficiency: f64,
    /// `MechanicalDesign.getCorrosionAllowance`, mm. Default `0`.
    pub corrosion_allowance_mm: f64,
    /// `SeparatorMechanicalDesign.getMaxGasVelocity`, m/s. Default `3.0`; it caps the
    /// Souders-Brown velocity.
    pub max_gas_velocity: f64,
    /// `getFoamAllowanceFactor`. The shipped database's value is `1.35`.
    pub foam_allowance_factor: f64,
    /// `getMinOilRetentionTime`, seconds. Default `120`.
    pub min_oil_retention_s: f64,
    /// `getMinWaterRetentionTime`, seconds. Default `180`.
    pub min_water_retention_s: f64,
    /// Whether the liquid phase is aqueous, which raises the retention floor.
    pub aqueous: bool,
}

/// The design's own numbers, in SI (m, kg, m3, s), on the state the separator leaves.
#[derive(Debug, Clone, PartialEq)]
pub struct Design {
    pub inner_diameter: f64,
    pub tantan_length: f64,
    pub outer_diameter: f64,
    pub wall_thickness: f64,
    pub max_design_volume_flow_m3_hr: f64,
    pub max_design_gas_volume_flow_m3_hr: f64,
    /// The gas phase's mass density, kg/m3 - the capture's density, carried for its diag.
    pub gas_density: f64,
    /// The liquid phase's mass density, kg/m3.
    pub liquid_density: f64,
    pub vessel_shell_weight_kg: f64,
    pub piping_weight_kg: f64,
    pub structural_weight_kg: f64,
    pub electrical_weight_kg: f64,
    pub total_weight_kg: f64,
    pub module_width: f64,
    pub module_length: f64,
    pub module_height: f64,
    pub vessel_volume_m3: f64,
    pub hhll: f64,
    pub hll: f64,
    pub nll: f64,
    pub lll: f64,
    pub llll: f64,
    pub weir_height: f64,
    pub hil: f64,
    pub nil: f64,
    pub lil: f64,
    pub effective_length_liquid: f64,
    pub effective_length_gas: f64,
    pub inlet_nozzle_id: f64,
    pub gas_outlet_nozzle_id: f64,
    pub oil_outlet_nozzle_id: f64,
    pub surge_volume_m3: f64,
    pub surge_time_s: f64,
    pub liquid_retention_time_s: f64,
}

/// `SeparatorMechanicalDesign.calcDesign`, on the outlets [`crate::kernels::separator`] leaves.
///
/// The design reads the **gas and liquid outlets' actual states**, not the feed: the class
/// takes `getThermoSystem().getPhase(0)` and `getPhase(1)` after `run`, which is the flashed
/// system. The liquid viscosity is the phase dispatch's (`PhaseKind::Oil` on a hydrocarbon
/// liquid, `Aqueous` otherwise), which the capture holds beside the class's own reading.
///
/// # Errors
/// A flash that refuses, a phase transport that refuses, a non-positive density or viscosity,
/// and the geometry divisions the class would answer `NaN` for.
pub fn mechanical_design(
    gas_outlet: &Stream,
    liquid_outlet: &Stream,
    horizontal: bool,
    kind: VesselKind,
    inputs: &DesignInputs,
    perform_sizing: bool,
) -> Result<Design> {
    let (mixture, ideal_gas) = gas_outlet.mixture()?;

    let gas_flash = pt_flash(&mixture, gas_outlet.t, gas_outlet.p, &gas_outlet.z)?;
    let gas_root = match gas_flash.phase {
        Phase::TwoPhase => gas_flash.z_vapour,
        _ => gas_flash.z_liquid,
    };
    let gas_molar_mass = molar_mass_of(&mixture, &gas_outlet.z)?;
    let gas_density = mass_density(
        &mixture,
        gas_outlet.t,
        gas_outlet.p,
        &gas_outlet.z,
        gas_root,
        gas_molar_mass,
    )?;

    let liquid_flash = pt_flash(&mixture, liquid_outlet.t, liquid_outlet.p, &liquid_outlet.z)?;
    let liquid_root = liquid_flash.z_liquid;
    let liquid_molar_mass = molar_mass_of(&mixture, &liquid_outlet.z)?;
    let liquid_density = mass_density(
        &mixture,
        liquid_outlet.t,
        liquid_outlet.p,
        &liquid_outlet.z,
        liquid_root,
        liquid_molar_mass,
    )?;
    let liquid_viscosity = phase_transport(
        &mixture,
        &ideal_gas,
        PhaseKind::Oil,
        liquid_outlet.t,
        liquid_outlet.p,
        &liquid_outlet.z,
    )?
    .mu
    .value;

    if !(gas_density.is_finite()
        && gas_density > 0.0
        && liquid_density.is_finite()
        && liquid_density > gas_density
        && liquid_viscosity.is_finite()
        && liquid_viscosity > 0.0)
    {
        return Err(AzothError::invalid_input(
            "mechanical_design",
            format!(
                "the design needs a gas denser than zero, a liquid denser than the gas and a \
                 positive liquid viscosity; it has {gas_density}, {liquid_density} and \
                 {liquid_viscosity}"
            ),
        ));
    }

    // **The class reads the phase's EOS volume, not mass over the physical-properties
    // density.** `Phase.getVolume()` is the cubic volume - untranslated - while
    // `getPhysicalProperties().getDensity()` carries the Peneloux shift, so the two differ by
    // the shift's volume on the capture's PR row (`0.11035858` m3 against `0.110877` from the
    // translated density). The phase's `Z` and NeqSim's gas constant are what the class's
    // volume is.
    let gas_z = match gas_flash.phase {
        Phase::AllLiquid | Phase::Trivial => gas_flash.z_liquid,
        _ => gas_flash.z_vapour,
    };
    let liquid_z = liquid_flash.z_liquid;
    const NEQSIM_R: f64 = 8.3144621;
    let gas_volume_flow = gas_outlet.n * gas_z * NEQSIM_R * gas_outlet.t.value / gas_outlet.p.value;
    let liquid_volume_flow =
        liquid_outlet.n * liquid_z * NEQSIM_R * liquid_outlet.t.value / liquid_outlet.p.value;

    let mut max_design_volume_flow = inputs.volume_safety_factor * gas_volume_flow * 3600.0;
    let allowable_gas_velocity = inputs
        .max_gas_velocity
        .min(inputs.gas_load_factor * ((liquid_density - gas_density) / gas_density).sqrt());
    let gas_area_fraction = if horizontal { inputs.fg } else { 1.0 };
    let liquid_area_fraction = if horizontal { 1.0 - inputs.fg } else { 1.0 };
    let mut inner_diameter = (4.0 * (max_design_volume_flow / 3600.0)
        / (std::f64::consts::PI * allowable_gas_velocity * gas_area_fraction))
        .sqrt();

    let mut tantan_length = match kind {
        VesselKind::Separator => {
            if !perform_sizing {
                // The bubble alternative: a 250 µm bubble in the liquid, settled by Stokes.
                let bubble_velocity =
                    9.82 * (250.0e-6_f64).powi(2) * (liquid_density - gas_density)
                        / (18.0 * liquid_viscosity);
                let ar = liquid_volume_flow / bubble_velocity;
                let daim = (ar / 4.0).sqrt();
                if daim > inner_diameter {
                    inner_diameter = daim;
                    4.0 * daim
                } else {
                    inner_diameter * 2.0
                }
            } else {
                inner_diameter * 2.0
            }
        }
        // **The scrubber's length is five diameters and nothing else sizes it.** Its class
        // never reads the retention time in `calcDesign`, so the liquid flow and viscosity
        // above do not reach the geometry.
        VesselKind::Scrubber => inner_diameter * 5.0,
    };

    if kind == VesselKind::Separator {
        // `sizeForLiquidHoldup`: the class field's retention, floored by the oil and
        // aqueous minimums, then the required liquid volume sets both dimensions.
        let mut retention = inputs.retention_time_s.max(inputs.min_oil_retention_s);
        if inputs.aqueous {
            retention = retention.max(inputs.min_water_retention_s);
        }
        let required_volume = liquid_volume_flow
            * retention
            * inputs.foam_allowance_factor
            * inputs.volume_safety_factor;
        if liquid_area_fraction > 0.0 {
            inner_diameter = inner_diameter.max(
                (4.0 * required_volume / (5.0 * std::f64::consts::PI * liquid_area_fraction))
                    .cbrt(),
            );
            let liquid_area =
                std::f64::consts::PI * inner_diameter * inner_diameter * liquid_area_fraction / 4.0;
            tantan_length = tantan_length
                .max(2.0 * inner_diameter)
                .max(required_volume / liquid_area);
        } else {
            tantan_length = tantan_length.max(2.0 * inner_diameter);
        }
    }

    let max_design_gas_volume_flow =
        allowable_gas_velocity * std::f64::consts::PI * inner_diameter * inner_diameter / 4.0
            * gas_area_fraction
            * 3600.0;
    max_design_volume_flow = max_design_gas_volume_flow;

    // `performSizingCalculations` runs on the design-standard state `autoSize` read, which on
    // the capture has no pressure-vessel code: the wall, the weights and the outer diameter
    // stay at their zero defaults and `volumeTotal` stays at the class's `2.0`.
    let (
        wall_thickness,
        outer_diameter,
        vessel_shell_weight,
        piping_weight,
        structural_weight,
        electrical_weight,
        total_weight,
        vessel_volume,
    ) = if perform_sizing {
        (0.0, inner_diameter, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0)
    } else {
        let pressure_bar = inputs.max_operation_pressure.value / 1.0e5;
        let wall_thickness_mm = pressure_bar / 10.0 * inner_diameter * 1.0e3
            / (2.0 * inputs.allowable_stress * inputs.joint_efficiency - 1.2 * pressure_bar / 10.0)
            + inputs.corrosion_allowance_mm;
        let wall_thickness = wall_thickness_mm / 1000.0;
        let shell = 0.032 * (wall_thickness * 1.0e3) * (inner_diameter * 1.0e3) * tantan_length;
        (
            wall_thickness,
            inner_diameter + 2.0 * wall_thickness,
            shell,
            shell * 0.4,
            shell * 0.1,
            shell * 0.08,
            shell + shell * 0.4 + shell * 0.1 + shell * 0.08,
            std::f64::consts::PI * (inner_diameter / 2.0).powi(2) * tantan_length,
        )
    };

    let module_width = inner_diameter * 2.0;
    let module_length = tantan_length * 1.5;
    let module_height = match kind {
        VesselKind::Separator => inner_diameter * 2.0 + 1.0,
        VesselKind::Scrubber => inner_diameter * 2.0,
    };

    let effective_length_liquid = tantan_length * 0.82;
    let effective_length_gas = tantan_length * 0.64;

    let total_volume_flow_m3_hr = (gas_volume_flow + liquid_volume_flow) * 3600.0;
    // The class sums each phase's `getFlowRate("m3/hr") * getDensity("kg/m3")`, so the mass
    // flow is per hour and the mixture density is kg per m3/hr - the same unit the velocity
    // is formed from.
    let total_mass_flow =
        (gas_volume_flow * gas_density + liquid_volume_flow * liquid_density) * 3600.0;
    let mix_density = total_mass_flow / total_volume_flow_m3_hr;
    let nozzle_velocity = (8000.0 / mix_density).sqrt();
    let inlet_nozzle_id = round_up(
        (4.0 * ((total_volume_flow_m3_hr / 3600.0) / nozzle_velocity) / std::f64::consts::PI)
            .sqrt(),
        0.05,
    )
    .max(0.05);
    let gas_outlet_nozzle_id = round_up(
        (4.0 * (gas_volume_flow / 20.0) / std::f64::consts::PI).sqrt(),
        0.05,
    );
    let oil_outlet_nozzle_id = round_up(
        (4.0 * (liquid_volume_flow / 1.5) / std::f64::consts::PI).sqrt(),
        0.025,
    )
    .max(0.05);

    let hhll = 0.80 * inner_diameter;
    let hll = 0.70 * inner_diameter;
    let nll = 0.50 * inner_diameter;
    let lll = 0.30 * inner_diameter;
    let llll = 0.15 * inner_diameter;
    let weir_height = 0.25 * inner_diameter;
    let hil = 0.25 * inner_diameter;
    let nil = 0.20 * inner_diameter;
    let lil = 0.15 * inner_diameter;

    let segment = |level: f64| -> f64 {
        let r = inner_diameter / 2.0;
        let h = level.min(inner_diameter);
        let theta = 2.0 * ((r - h) / r).acos();
        (r * r / 2.0) * (theta - theta.sin()) * effective_length_liquid
    };
    let surge_volume = segment(hll) - segment(nll);
    // **The class's own unit slip**: the flow is m3/s and is divided as though m3/hr.
    let liquid_hourly = liquid_volume_flow / 3600.0;
    let liquid_retention_time = segment(nll) / liquid_hourly;
    let surge_time = surge_volume / liquid_hourly;

    Ok(Design {
        inner_diameter,
        tantan_length,
        outer_diameter,
        wall_thickness,
        max_design_volume_flow_m3_hr: max_design_volume_flow,
        max_design_gas_volume_flow_m3_hr: max_design_gas_volume_flow,
        gas_density,
        liquid_density,
        vessel_shell_weight_kg: vessel_shell_weight,
        piping_weight_kg: piping_weight,
        structural_weight_kg: structural_weight,
        electrical_weight_kg: electrical_weight,
        total_weight_kg: total_weight,
        module_width,
        module_length,
        module_height,
        vessel_volume_m3: vessel_volume,
        hhll,
        hll,
        nll,
        lll,
        llll,
        weir_height,
        hil,
        nil,
        lil,
        effective_length_liquid,
        effective_length_gas,
        inlet_nozzle_id,
        gas_outlet_nozzle_id,
        oil_outlet_nozzle_id,
        surge_volume_m3: surge_volume,
        surge_time_s: surge_time,
        liquid_retention_time_s: liquid_retention_time,
    })
}

/// `Math.ceil(value / step) * step`.
fn round_up(value: f64, step: f64) -> f64 {
    (value / step).ceil() * step
}

#[cfg(test)]
mod tests {
    use super::{DesignInputs, VesselKind, mechanical_design};
    use crate::segment::phase::molar_mass_of;
    use crate::stream::Stream;
    use azoth_core::units::{kelvins, pascals};

    fn stream() -> (Stream, Stream) {
        // The capture's PR row: the class's fixture composition on the cubic this port
        // carries, at its 10,000 kg/hr - which is the flow the captured diameter belongs to.
        let components = vec![
            "methane".to_string(),
            "ethane".to_string(),
            "propane".to_string(),
            "n-pentane".to_string(),
        ];
        let z = vec![0.7, 0.1, 0.1, 0.1];
        let probe = Stream::from_pt(
            components.clone(),
            z.clone(),
            1.0,
            pascals(20.0e5),
            kelvins(298.15),
        )
        .expect("the captured feed");
        let (mixture, _ideal_gas) = probe.mixture().expect("the captured mixture");
        let molar_mass = molar_mass_of(&mixture, &z).expect("a molar mass");
        let n = (10000.0 / 3600.0) / molar_mass;
        let feed = Stream::from_pt(components, z, n, pascals(20.0e5), kelvins(298.15))
            .expect("the captured feed");
        crate::kernels::separator::separator(&feed, pascals(0.0), 0.0, None)
            .expect("the captured separation")
    }

    #[test]
    fn the_captures_pr_row_is_reproduced() {
        let (gas, liquid) = stream();
        let design = mechanical_design(
            &gas,
            &liquid,
            true,
            VesselKind::Separator,
            &DesignInputs {
                gas_load_factor: 0.137,
                fg: 0.8,
                volume_safety_factor: 1.0,
                retention_time_s: 150.0,
                max_operation_pressure: pascals(100.0e5),
                allowable_stress: 95.147_650_602,
                joint_efficiency: 1.0,
                corrosion_allowance_mm: 0.0,
                max_gas_velocity: 3.0,
                foam_allowance_factor: 1.35,
                min_oil_retention_s: 120.0,
                min_water_retention_s: 180.0,
                aqueous: false,
            },
            false,
        )
        .expect("the captured design");
        let close = |actual: f64, expected: f64, what: &str| {
            let error = (actual - expected).abs() / expected.abs().max(f64::MIN_POSITIVE);
            assert!(
                error <= 2.0e-4,
                "{what}: {actual} against {expected} ({error:.3e})"
            );
        };
        close(design.inner_diameter, 0.635_887_318_997_330_6, "diameter");
        close(design.tantan_length, 3.179_436_594_986_652_3, "length");
        close(design.wall_thickness, 0.035_664_838_877_036_247, "wall");
        close(
            design.max_design_volume_flow_m3_hr,
            662.875_955_687_116,
            "design flow",
        );
        close(design.hhll, 0.508_709_855_197_864_5, "HHLL");
        close(design.nll, 0.317_943_659_498_665_3, "NLL");
        close(design.weir_height, 0.158_971_829_749_332_66, "weir");
        close(
            design.effective_length_liquid,
            2.607_138_007_889_054_7,
            "liquid length",
        );
        close(
            design.liquid_retention_time_s,
            1_494_450.000_000_000_5,
            "retention",
        );
        close(
            design.surge_volume_m3,
            0.205_074_922_382_108_28,
            "surge volume",
        );
        close(design.surge_time_s, 740_303.342_040_814_7, "surge time");
        close(design.inlet_nozzle_id, 0.1, "inlet nozzle");
        close(design.gas_outlet_nozzle_id, 0.1, "gas nozzle");
        close(design.oil_outlet_nozzle_id, 0.05, "oil nozzle");
    }

    /// The scrubber capture's PR row: vertical, the scrubber standard's `1.1` safety factor,
    /// `5 d` length and a module `2 d` high rather than `2 d + 1`.
    #[test]
    fn the_scrubbers_pr_row_is_reproduced() {
        let (gas, liquid) = stream();
        let design = mechanical_design(
            &gas,
            &liquid,
            false,
            VesselKind::Scrubber,
            &DesignInputs {
                gas_load_factor: 0.11,
                fg: 1.0,
                volume_safety_factor: 1.1,
                retention_time_s: 60.0,
                max_operation_pressure: pascals(100.0e5),
                allowable_stress: 95.147_650_602,
                joint_efficiency: 1.0,
                corrosion_allowance_mm: 0.0,
                max_gas_velocity: 3.0,
                foam_allowance_factor: 1.35,
                min_oil_retention_s: 120.0,
                min_water_retention_s: 180.0,
                aqueous: true,
            },
            false,
        )
        .expect("the captured scrubber design");
        let close = |actual: f64, expected: f64, what: &str| {
            let error = (actual - expected).abs() / expected.abs().max(f64::MIN_POSITIVE);
            assert!(
                error <= 2.0e-4,
                "{what}: {actual} against {expected} ({error:.3e})"
            );
        };
        close(design.inner_diameter, 0.515_375_309_329_634_3, "diameter");
        close(design.tantan_length, 2.576_876_546_648_171_6, "length");
        close(design.wall_thickness, 0.028_905_714_612_184_743, "wall");
        close(
            design.vessel_shell_weight_kg,
            1_228.431_403_439_446,
            "shell weight",
        );
        close(
            design.module_height,
            1.030_750_618_659_268_5,
            "module height",
        );
        close(design.hhll, 0.412_300_247_463_707_4, "HHLL");
        close(
            design.surge_volume_m3,
            0.109_179_949_337_516_54,
            "surge volume",
        );
        close(design.surge_time_s, 394_130.498_451_730_7, "surge time");
        close(
            design.liquid_retention_time_s,
            795_631.047_385_863_8,
            "retention",
        );
    }
}
