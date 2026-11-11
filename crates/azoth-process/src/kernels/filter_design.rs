//! `FilterMechanicalDesign.calcDesign` - the element count, vessel, boundary and weights.
//!
//! **The design is called by a caller, never by `Filter.run`.** It reads the inlet's volumetric
//! flow and density, the filter type's default face velocity and element geometry, and the
//! loading state the run just computed for its utilization checks. It writes the required element
//! count back to the equipment and, for a granular type, replaces the hydraulic media area with
//! the vessel cross-section - the port publishes that area so a caller can run again.
//!
//! **The pressure boundary is the class's own screening arithmetic, mixing bar and MPa as it
//! writes them**: `P = (designPressure_bara - 1.01325) * 0.1` MPa against an allowable stress in
//! MPa, a 6 mm shell floor and the head floored at the shell. It is ASME VIII Div. 1's shape, not
//! a code calculation, and the results are preliminary.
//!
//! Not here: the installation and maintenance cost arithmetic (USD, which this library has no
//! currency dimension for), the material grade string (stored and reported, read by no formula),
//! and the filter-type default pressure-drop model (the model id states its selector explicitly).

use std::f64::consts::PI;

use azoth_core::{AzothError, Result};

use crate::segment::phase::{system_mass_density, system_volumetric_flow};
use crate::stream::Stream;

/// `FilterType`: the eight constructions and the defaults their design reads.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FilterType {
    /// A replaceable depth or pleated cartridge element.
    Cartridge,
    /// A replaceable bag element.
    Bag,
    /// A Y-pattern screen or strainer.
    YStrainer,
    /// A basket screen with a larger dirt-holding area.
    BasketStrainer,
    /// A coalescing element for dispersed droplets.
    Coalescer,
    /// A fixed granular-media bed.
    GranularMedia,
    /// A granular bed designed for periodic backwash.
    BackwashableMedia,
    /// An activated-carbon or other replaceable sorbent bed.
    ActivatedCarbon,
}

/// The type-dependent design inputs.
#[derive(Debug, Clone, Copy)]
pub struct TypeDefaults {
    /// The preliminary maximum element face velocity, m/s.
    pub max_face_velocity_mps: f64,
    /// The effective filtration area per element, m2.
    pub element_area_m2: f64,
    /// The element envelope diameter, m.
    pub element_diameter_m: f64,
    /// The element length, m.
    pub element_length_m: f64,
    /// Whether the type is a granular bed rather than a cartridge bundle.
    pub granular: bool,
}

impl FilterType {
    /// The class's string spelling, as the model id's enum states it.
    #[must_use]
    pub fn parse(spelling: &str) -> Option<Self> {
        Some(match spelling {
            "cartridge" => Self::Cartridge,
            "bag" => Self::Bag,
            "y_strainer" => Self::YStrainer,
            "basket_strainer" => Self::BasketStrainer,
            "coalescer" => Self::Coalescer,
            "granular_media" => Self::GranularMedia,
            "backwashable_media" => Self::BackwashableMedia,
            "activated_carbon" => Self::ActivatedCarbon,
            _ => return None,
        })
    }

    /// `FilterType`'s constructor arguments, as the enum writes them.
    #[must_use]
    pub fn defaults(self) -> TypeDefaults {
        let (face, area, diameter, length, granular) = match self {
            Self::Cartridge => (0.15, 1.0, 0.15, 1.0, false),
            Self::Bag => (0.10, 0.5, 0.18, 0.8, false),
            Self::YStrainer => (1.5, 0.25, 0.20, 0.5, false),
            Self::BasketStrainer => (1.0, 0.5, 0.30, 0.8, false),
            Self::Coalescer => (0.08, 1.5, 0.15, 1.0, false),
            Self::GranularMedia | Self::BackwashableMedia | Self::ActivatedCarbon => {
                (0.02, 1.0, 0.0, 1.0, true)
            }
        };
        TypeDefaults {
            max_face_velocity_mps: face,
            element_area_m2: area,
            element_diameter_m: diameter,
            element_length_m: length,
            granular,
        }
    }
}

/// The design's settings, with the class's defaults already resolved by the model.
#[derive(Debug, Clone, Copy)]
pub struct DesignRequest {
    /// The type whose defaults fill the unstated settings.
    pub filter_type: FilterType,
    /// The element count the equipment carries before sizing.
    pub number_of_elements: u32,
    /// The user's face velocity, m/s; absent is the type's.
    pub max_face_velocity_mps: Option<f64>,
    /// The user's element area, m2; absent is the type's.
    pub element_area_m2: Option<f64>,
    /// The user's nozzle velocity, m/s; absent is the density switch.
    pub max_nozzle_velocity_mps: Option<f64>,
    /// The minimum vessel L/D ratio.
    pub min_ld_ratio: f64,
    /// The maximum vessel L/D ratio.
    pub max_ld_ratio: f64,
    /// The design pressure's multiplier over the operating pressure.
    pub design_pressure_margin: f64,
    /// The design pressure's additive margin, Pa.
    pub minimum_design_pressure_margin_pa: f64,
    /// The design temperature's margin, K.
    pub design_temperature_margin_k: f64,
    /// The allowable stress, Pa (the class states it in MPa).
    pub allowable_stress_pa: f64,
    /// The welded-joint efficiency.
    pub joint_efficiency: f64,
    /// The corrosion allowance, m.
    pub corrosion_allowance_m: f64,
    /// The element collapse/burst rating, Pa; positive.
    pub element_collapse_pressure_pa: f64,
    /// Whether element fabrication integrity has been recorded.
    pub element_integrity_verified: bool,
    /// Whether a beta curve or a positive nominal efficiency is stated.
    pub performance_data: bool,
    /// The granular bed's depth, m.
    pub media_bed_depth_m: f64,
    /// The terminal differential pressure, Pa.
    pub terminal_pressure_drop_pa: f64,
    /// The unrestricted drop the run computed, Pa.
    pub unrestricted_drop_pa: f64,
    /// The terminal-drop utilization the run computed.
    pub terminal_pressure_utilization: f64,
}

/// One preliminary design warning, in the class's order.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DesignWarning {
    /// The terminal differential pressure is exceeded.
    TerminalExceeded,
    /// The calculated differential pressure exceeds the element rating.
    CollapseExceeded,
    /// The terminal setting is at or above the collapse rating.
    TerminalAboveCollapse,
    /// The calculated face velocity exceeds the design limit.
    FaceVelocity,
    /// The selected nozzle exceeds the velocity limit.
    NozzleVelocity,
    /// Element fabrication integrity is not recorded.
    IntegrityNotVerified,
    /// No beta-ratio curve or positive nominal efficiency is configured.
    NoPerformanceData,
}

impl DesignWarning {
    /// The class's own wording.
    #[must_use]
    pub fn message(self) -> &'static str {
        match self {
            Self::TerminalExceeded => {
                "Terminal differential pressure is exceeded; replace or backwash the element"
            }
            Self::CollapseExceeded => {
                "Calculated differential pressure exceeds the element collapse/burst rating"
            }
            Self::TerminalAboveCollapse => {
                "Terminal replacement differential pressure must be below the element collapse \
                 rating"
            }
            Self::FaceVelocity => {
                "Calculated element face velocity exceeds the configured design limit"
            }
            Self::NozzleVelocity => {
                "Selected preliminary nozzle size exceeds the configured velocity limit"
            }
            Self::IntegrityNotVerified => {
                "Filter element fabrication integrity has not been recorded as verified"
            }
            Self::NoPerformanceData => {
                "No particle-size efficiency or beta-ratio test data have been configured"
            }
        }
    }
}

/// What `calcDesign` computed.
#[derive(Debug, Clone)]
pub struct FilterDesign {
    /// The design pressure, Pa.
    pub design_pressure_pa: f64,
    /// The design temperature, K.
    pub design_temperature_k: f64,
    /// The required element count.
    pub required_elements: f64,
    /// The vessel's inside diameter, m.
    pub inner_diameter_m: f64,
    /// The vessel's outside diameter, m.
    pub outer_diameter_m: f64,
    /// The tangent-to-tangent length, m.
    pub vessel_length_m: f64,
    /// The shell thickness, m.
    pub shell_thickness_m: f64,
    /// The head thickness, m.
    pub head_thickness_m: f64,
    /// The effective element area, m2.
    pub element_area_m2: f64,
    /// The calculated face velocity, m/s.
    pub calculated_face_velocity_mps: f64,
    /// The required nozzle inside diameter, m.
    pub required_nozzle_diameter_m: f64,
    /// The selected nominal nozzle diameter, m.
    pub selected_nozzle_diameter_m: f64,
    /// The velocity through the selected nozzle, m/s.
    pub calculated_nozzle_velocity_mps: f64,
    /// The empty vessel weight, kg.
    pub empty_vessel_weight_kg: f64,
    /// The total equipped weight, kg.
    pub total_equipped_weight_kg: f64,
    /// The terminal-drop utilization the run computed.
    pub terminal_pressure_utilization: f64,
    /// The collapse-rating utilization.
    pub collapse_pressure_utilization: f64,
    /// Whether both differential-pressure checks pass.
    pub differential_pressure_design_acceptable: bool,
    /// The vessel cross-section a granular design writes to the filter, m2.
    pub design_media_area_m2: Option<f64>,
    /// The preliminary warnings, in the class's order.
    pub warnings: Vec<DesignWarning>,
}

/// The nominal nozzle inside diameters the class selects from, mm.
pub const NOZZLE_DIAMETERS_MM: [f64; 21] = [
    15.0, 20.0, 25.0, 40.0, 50.0, 65.0, 80.0, 100.0, 150.0, 200.0, 250.0, 300.0, 350.0, 400.0,
    450.0, 500.0, 600.0, 750.0, 900.0, 1050.0, 1200.0,
];

/// `FilterMechanicalDesign.calcDesign`.
///
/// # Errors
/// [`azoth_core::AzothError::InvalidInput`] where the shell or head denominator is not positive,
/// which is the class's `IllegalStateException`, and whatever the flow or density lookup refuses.
pub fn design(feed: &Stream, request: &DesignRequest) -> Result<FilterDesign> {
    let defaults = request.filter_type.defaults();
    let operating_pressure_bara = feed.p.value / 1.0e5;
    let design_pressure_bar = (operating_pressure_bara * request.design_pressure_margin)
        .max(operating_pressure_bara + request.minimum_design_pressure_margin_pa / 1.0e5);
    let design_temperature_k = feed.t.value + request.design_temperature_margin_k;

    let measured_density = system_mass_density(feed).unwrap_or(50.0);
    let density = if measured_density.is_finite() && measured_density > 0.0 {
        measured_density
    } else {
        50.0
    };
    let mut volumetric = system_volumetric_flow(feed)?.max(0.0);
    if volumetric <= 0.0 {
        let mass_flow_kg_h = feed.n * feed.molar_mass()?.value * 3600.0;
        volumetric = (mass_flow_kg_h / 3600.0 / density).max(0.0);
    }

    let design_face_velocity = request
        .max_face_velocity_mps
        .unwrap_or(defaults.max_face_velocity_mps);
    let design_element_area = request.element_area_m2.unwrap_or(defaults.element_area_m2);
    let required_area = volumetric / design_face_velocity.max(1.0e-12);

    let (required_elements, mut inner_diameter, mut face_velocity) = if defaults.granular {
        let diameter = (4.0 * required_area / PI).sqrt().max(0.5);
        (1.0, diameter, volumetric / (PI * diameter * diameter / 4.0))
    } else {
        let count = (required_area / design_element_area.max(1.0e-12))
            .ceil()
            .min(f64::from(i32::MAX));
        let count = count.max(1.0).max(f64::from(request.number_of_elements));
        let bundle = count.sqrt() * defaults.element_diameter_m * 1.8;
        let diameter = (bundle + 0.2).max(0.5);
        (count, diameter, volumetric / (count * design_element_area))
    };

    let active_length = if defaults.granular {
        request.media_bed_depth_m
    } else {
        defaults.element_length_m
    };
    let mut vessel_length = active_length + inner_diameter;
    let ld_ratio = vessel_length / inner_diameter;
    if ld_ratio < request.min_ld_ratio {
        vessel_length = inner_diameter * request.min_ld_ratio;
    } else if ld_ratio > request.max_ld_ratio {
        inner_diameter = vessel_length / request.max_ld_ratio;
    }
    let design_media_area_m2 = if defaults.granular {
        let area = PI * inner_diameter * inner_diameter / 4.0;
        face_velocity = volumetric / area;
        Some(area)
    } else {
        None
    };

    let design_nozzle_velocity = request
        .max_nozzle_velocity_mps
        .unwrap_or(if density < 200.0 { 20.0 } else { 3.0 });
    let required_nozzle_diameter = if volumetric > 0.0 {
        (4.0 * volumetric / (PI * design_nozzle_velocity)).sqrt()
    } else {
        0.0
    };
    let selected_nozzle_diameter_mm = select_nozzle(required_nozzle_diameter * 1000.0);
    let selected_nozzle_diameter = selected_nozzle_diameter_mm / 1000.0;
    let selected_area = PI * selected_nozzle_diameter * selected_nozzle_diameter / 4.0;
    let calculated_nozzle_velocity = if selected_area > 0.0 {
        volumetric / selected_area
    } else {
        0.0
    };

    let design_pressure_mpa = (design_pressure_bar - 1.01325).max(0.0) * 0.1;
    let allowable_stress_mpa = request.allowable_stress_pa / 1.0e6;
    let shell_denominator =
        allowable_stress_mpa * request.joint_efficiency - 0.6 * design_pressure_mpa;
    if shell_denominator <= 0.0 {
        return Err(AzothError::invalid_input(
            "mechanical_design",
            "the shell design pressure exceeds the range of the configured material data",
        ));
    }
    let shell_thickness = (design_pressure_mpa * inner_diameter / 2.0 / shell_denominator
        + request.corrosion_allowance_m)
        .max(0.006);
    let head_denominator =
        2.0 * allowable_stress_mpa * request.joint_efficiency - 0.2 * design_pressure_mpa;
    if head_denominator <= 0.0 {
        return Err(AzothError::invalid_input(
            "mechanical_design",
            "the head design pressure exceeds the range of the configured material data",
        ));
    }
    let head_thickness = (design_pressure_mpa * inner_diameter / head_denominator
        + request.corrosion_allowance_m)
        .max(shell_thickness);
    let outer_diameter = inner_diameter + 2.0 * shell_thickness;

    let collapse_pressure_utilization =
        request.unrestricted_drop_pa / request.element_collapse_pressure_pa;
    let acceptable =
        request.terminal_pressure_utilization < 1.0 && collapse_pressure_utilization < 1.0;

    let mut warnings = Vec::new();
    if request.terminal_pressure_utilization >= 1.0 {
        warnings.push(DesignWarning::TerminalExceeded);
    }
    if collapse_pressure_utilization >= 1.0 {
        warnings.push(DesignWarning::CollapseExceeded);
    }
    if request.terminal_pressure_drop_pa >= request.element_collapse_pressure_pa {
        warnings.push(DesignWarning::TerminalAboveCollapse);
    }
    if face_velocity > design_face_velocity * 1.000_001 {
        warnings.push(DesignWarning::FaceVelocity);
    }
    if calculated_nozzle_velocity > design_nozzle_velocity * 1.000_001 {
        warnings.push(DesignWarning::NozzleVelocity);
    }
    if !request.element_integrity_verified {
        warnings.push(DesignWarning::IntegrityNotVerified);
    }
    if !request.performance_data {
        warnings.push(DesignWarning::NoPerformanceData);
    }

    let steel_density = 7850.0;
    let shell_area = PI * inner_diameter * vessel_length;
    let head_area = 2.0 * PI * (inner_diameter / 2.0).powi(2) * 1.084;
    let empty_vessel_weight = (shell_area + head_area) * shell_thickness * steel_density;
    let accessories = empty_vessel_weight * 0.30;
    let element_weight = if defaults.granular {
        let media_density = if request.filter_type == FilterType::ActivatedCarbon {
            500.0
        } else {
            1000.0
        };
        PI * inner_diameter * inner_diameter / 4.0 * request.media_bed_depth_m * media_density
    } else {
        required_elements * (defaults.element_length_m * 10.0).max(2.0)
    };
    let total_equipped_weight = empty_vessel_weight + accessories + element_weight;

    Ok(FilterDesign {
        design_pressure_pa: design_pressure_bar * 1.0e5,
        design_temperature_k,
        required_elements,
        inner_diameter_m: inner_diameter,
        outer_diameter_m: outer_diameter,
        vessel_length_m: vessel_length,
        shell_thickness_m: shell_thickness,
        head_thickness_m: head_thickness,
        element_area_m2: design_element_area,
        calculated_face_velocity_mps: face_velocity,
        required_nozzle_diameter_m: required_nozzle_diameter,
        selected_nozzle_diameter_m: selected_nozzle_diameter,
        calculated_nozzle_velocity_mps: calculated_nozzle_velocity,
        empty_vessel_weight_kg: empty_vessel_weight,
        total_equipped_weight_kg: total_equipped_weight,
        terminal_pressure_utilization: request.terminal_pressure_utilization,
        collapse_pressure_utilization,
        differential_pressure_design_acceptable: acceptable,
        design_media_area_m2,
        warnings,
    })
}

/// The next nominal nozzle diameter at or above the required one, mm; the required one above the
/// table.
fn select_nozzle(required_mm: f64) -> f64 {
    if required_mm <= 0.0 {
        return 0.0;
    }
    for nominal in NOZZLE_DIAMETERS_MM {
        if nominal >= required_mm {
            return nominal;
        }
    }
    required_mm
}

#[cfg(test)]
mod tests {
    use super::*;
    use azoth_core::units::{kelvins, pascals};

    fn gas_feed() -> Stream {
        Stream::from_pt(
            vec!["methane".to_string()],
            vec![1.0],
            173.145_781_822_463_25,
            pascals(5_000_000.0),
            kelvins(303.15),
        )
        .expect("the feed")
    }

    fn low_gas_feed() -> Stream {
        Stream::from_pt(
            vec!["methane".to_string()],
            vec![1.0],
            8.657_289_091_123_163,
            pascals(5_000_000.0),
            kelvins(303.15),
        )
        .expect("the feed")
    }

    fn liquid_feed() -> Stream {
        Stream::from_pt(
            vec!["n-butane".to_string()],
            vec![1.0],
            23.895_684_821_652_168,
            pascals(500_000.0),
            kelvins(300.0),
        )
        .expect("the feed")
    }

    fn request(filter_type: FilterType) -> DesignRequest {
        DesignRequest {
            filter_type,
            number_of_elements: 1,
            max_face_velocity_mps: None,
            element_area_m2: None,
            max_nozzle_velocity_mps: None,
            min_ld_ratio: 2.0,
            max_ld_ratio: 5.0,
            design_pressure_margin: 1.10,
            minimum_design_pressure_margin_pa: 350_000.0,
            design_temperature_margin_k: 25.0,
            allowable_stress_pa: 138.0e6,
            joint_efficiency: 0.85,
            corrosion_allowance_m: 0.003,
            element_collapse_pressure_pa: 1.0e6,
            element_integrity_verified: false,
            performance_data: false,
            media_bed_depth_m: 1.0,
            terminal_pressure_drop_pa: 500_000.0,
            unrestricted_drop_pa: 0.0,
            terminal_pressure_utilization: 0.0,
        }
    }

    fn close(actual: f64, expected: f64, tolerance: f64, what: &str) {
        let error = (actual - expected).abs() / expected.abs().max(f64::MIN_POSITIVE);
        assert!(
            error <= tolerance,
            "{what}: {actual} against {expected}, relative error {error:.3e} > {tolerance:.3e}"
        );
    }

    #[test]
    fn the_coalescer_default_matches_the_capture() {
        let mut request = request(FilterType::Coalescer);
        request.element_collapse_pressure_pa = 500_000.0;
        request.terminal_pressure_drop_pa = 100_000.0;
        request.unrestricted_drop_pa = 15_000.0;
        request.terminal_pressure_utilization = 0.15;
        request.element_integrity_verified = true;
        request.performance_data = true;
        let design = design(&gas_feed(), &request).expect("the design");
        assert_eq!(design.required_elements, 1.0);
        close(design.inner_diameter_m, 0.5, 1.0e-12, "inner diameter");
        close(design.vessel_length_m, 1.5, 1.0e-12, "vessel length");
        close(
            design.shell_thickness_m,
            0.014_832_889_206_146_602,
            1.0e-9,
            "shell",
        );
        close(design.selected_nozzle_diameter_m, 0.08, 1.0e-12, "nozzle");
        close(
            design.calculated_nozzle_velocity_mps,
            15.901_520_894_119_196,
            5.0e-4,
            "nozzle velocity",
        );
        close(
            design.empty_vessel_weight_kg,
            323.917_079_245_233_25,
            1.0e-9,
            "empty weight",
        );
        close(
            design.total_equipped_weight_kg,
            431.092_203_018_803_2,
            1.0e-9,
            "total weight",
        );
        assert!(design.warnings.is_empty());
    }

    #[test]
    fn the_basket_strainer_warns_in_the_class_order() {
        let mut request = request(FilterType::BasketStrainer);
        request.element_collapse_pressure_pa = 150_000.0;
        request.terminal_pressure_drop_pa = 100_000.0;
        request.unrestricted_drop_pa = 200_000.0;
        request.terminal_pressure_utilization = 2.0;
        let design = design(&gas_feed(), &request).expect("the design");
        close(design.inner_diameter_m, 0.74, 1.0e-12, "inner diameter");
        close(design.vessel_length_m, 1.54, 1.0e-12, "vessel length");
        assert_eq!(
            design.warnings,
            vec![
                DesignWarning::TerminalExceeded,
                DesignWarning::CollapseExceeded,
                DesignWarning::IntegrityNotVerified,
                DesignWarning::NoPerformanceData,
            ]
        );
        assert!(!design.differential_pressure_design_acceptable);
    }

    #[test]
    fn the_activated_carbon_uses_the_vessel_cross_section() {
        let mut request = request(FilterType::ActivatedCarbon);
        request.unrestricted_drop_pa = 4_476.258_265_553_568;
        request.terminal_pressure_utilization = 0.008_952_516_531_107_135;
        request.element_integrity_verified = true;
        request.performance_data = true;
        let design = design(&gas_feed(), &request).expect("the design");
        assert_eq!(design.required_elements, 1.0);
        close(
            design.inner_diameter_m,
            2.255_767_427_311_189_6,
            5.0e-4,
            "diameter",
        );
        close(
            design.vessel_length_m,
            4.511_534_854_622_379,
            5.0e-4,
            "length",
        );
        close(
            design.design_media_area_m2.expect("a media area"),
            3.996_488_097_749_557_6,
            5.0e-4,
            "media area",
        );
        assert!(design.warnings.is_empty());
    }

    #[test]
    fn the_long_bed_takes_the_maximum_ld_ratio() {
        let mut request = request(FilterType::GranularMedia);
        request.media_bed_depth_m = 4.0;
        request.unrestricted_drop_pa = 59.894_456_143_419_2;
        request.terminal_pressure_utilization = 1.197_889_122_868_384e-4;
        request.element_integrity_verified = true;
        let design = design(&low_gas_feed(), &request).expect("the design");
        close(
            design.inner_diameter_m,
            0.900_880_986_177_952_8,
            5.0e-4,
            "diameter",
        );
        close(
            design.vessel_length_m,
            4.504_404_930_889_764,
            5.0e-4,
            "length",
        );
        close(
            design.design_media_area_m2.expect("a media area"),
            0.637_418_586_795_286,
            5.0e-4,
            "media area",
        );
    }

    #[test]
    fn the_overrides_are_the_captures() {
        let mut request = request(FilterType::Cartridge);
        request.max_face_velocity_mps = Some(0.05);
        request.element_area_m2 = Some(0.5);
        request.max_nozzle_velocity_mps = Some(10.0);
        request.min_ld_ratio = 1.5;
        request.max_ld_ratio = 3.0;
        request.design_pressure_margin = 1.2;
        request.minimum_design_pressure_margin_pa = 500_000.0;
        request.design_temperature_margin_k = 15.0;
        request.allowable_stress_pa = 100.0e6;
        request.joint_efficiency = 0.7;
        request.corrosion_allowance_m = 0.006;
        request.element_integrity_verified = true;
        let design = design(&gas_feed(), &request).expect("the design");
        close(
            design.design_pressure_pa,
            6_000_000.0,
            1.0e-12,
            "design pressure",
        );
        close(
            design.design_temperature_k,
            318.15,
            1.0e-12,
            "design temperature",
        );
        assert_eq!(design.required_elements, 4.0);
        close(design.inner_diameter_m, 0.74, 1.0e-12, "inner diameter");
        close(design.vessel_length_m, 1.74, 1.0e-12, "vessel length");
        close(
            design.shell_thickness_m,
            0.038_839_055_716_983_23,
            1.0e-9,
            "shell",
        );
        close(design.selected_nozzle_diameter_m, 0.15, 1.0e-12, "nozzle");
        close(
            design.calculated_nozzle_velocity_mps,
            4.523_099_276_549_46,
            5.0e-4,
            "nozzle velocity",
        );
    }

    #[test]
    fn the_liquid_density_switches_the_nozzle_velocity() {
        let mut request = request(FilterType::Bag);
        request.element_integrity_verified = true;
        let design = design(&liquid_feed(), &request).expect("the design");
        close(design.inner_diameter_m, 0.524, 1.0e-12, "inner diameter");
        close(design.vessel_length_m, 1.324, 1.0e-12, "vessel length");
        close(design.shell_thickness_m, 0.006, 1.0e-12, "shell floor");
        close(design.selected_nozzle_diameter_m, 0.04, 1.0e-12, "nozzle");
        close(
            design.calculated_nozzle_velocity_mps,
            1.960_149_184_530_458_8,
            5.0e-4,
            "nozzle velocity",
        );
    }

    #[test]
    fn a_shell_beyond_the_material_range_is_refused() {
        let mut request = request(FilterType::Cartridge);
        request.allowable_stress_pa = 1.0e6;
        request.design_pressure_margin = 100.0;
        assert!(design(&gas_feed(), &request).is_err());
    }
}
