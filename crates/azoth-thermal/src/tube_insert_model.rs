//! `thermal.tube_insert_model` - twisted tape, wire matrix and coiled wire enhancement.
//!
//! Transcribed from `process/mechanicaldesign/heatexchanger/TubeInsertModel.java` at NeqSim
//! `118d9aff`: the enhancement and friction ratios of the three inserts, the performance
//! evaluation criterion, and the enhancement applied to a plain-tube state.
//!
//! Spec: `specs/calcs/thermal/tube_insert_model.toml`.

use azoth_core::units::{HeatTransfer, Pressure, pascals, watts_per_square_meter_kelvin};
use azoth_core::{AzothError, Result, apply_checks};
use std::f64::consts::PI;

use crate::results::TubeInsertModelResult;
use crate::spec_gen;

/// The insert type.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum InsertType {
    /// A plain tube.
    None,
    /// A twisted tape.
    TwistedTape,
    /// A wire matrix, HiTRAN-like.
    WireMatrix,
    /// A coiled wire.
    CoiledWire,
}

impl std::str::FromStr for InsertType {
    type Err = AzothError;

    fn from_str(value: &str) -> Result<Self> {
        match value {
            "none" => Ok(Self::None),
            "twisted_tape" => Ok(Self::TwistedTape),
            "wire_matrix" => Ok(Self::WireMatrix),
            "coiled_wire" => Ok(Self::CoiledWire),
            other => Err(AzothError::invalid_input(
                "insert",
                format!("{other} is not one of none, twisted_tape, wire_matrix or coiled_wire"),
            )),
        }
    }
}

/// The class's default parameters.
pub const TWIST_RATIO: f64 = 5.0;
/// The class's default tape thickness.
pub const TAPE_THICKNESS: f64 = 0.001;
/// The class's default tape width fraction.
pub const TAPE_WIDTH_FRACTION: f64 = 0.95;
/// The class's default wire diameter.
pub const WIRE_DIAMETER: f64 = 0.001;
/// The class's default wire coil pitch.
pub const WIRE_COIL_PITCH: f64 = 0.010;
/// The class's default matrix density.
pub const MATRIX_DENSITY: f64 = 0.5;
/// The class's default roughness ratio.
pub const ROUGHNESS_RATIO: f64 = 0.02;
/// The class's default helix angle, in degrees.
pub const HELIX_ANGLE: f64 = 45.0;

/// The composed request, in SI magnitudes.
#[derive(Debug, Clone)]
pub struct TubeInsertRequest {
    pub insert: InsertType,
    pub reynolds_number: f64,
    pub prandtl_number: f64,
    pub twist_ratio: Option<f64>,
    pub tape_thickness: Option<f64>,
    pub tape_width_fraction: Option<f64>,
    pub wire_diameter: Option<f64>,
    pub wire_coil_pitch: Option<f64>,
    pub matrix_density: Option<f64>,
    pub roughness_ratio: Option<f64>,
    pub helix_angle: Option<f64>,
    pub plain_heat_transfer_coefficient: f64,
    pub plain_pressure_drop: f64,
}

/// The class's blockage factor: `t/(t/0.05)` cancels the thickness, and the factor is clamped.
#[must_use]
#[allow(clippy::manual_clamp)] // `clamp` panics on NaN where the class's max/min propagate it
pub fn tape_blockage(thickness: f64, width: f64) -> f64 {
    let factor = PI / (PI - 4.0 * thickness * width / (thickness / 0.05));
    factor.min(1.5).max(1.0)
}

/// The swirl number, `sqrt(1 + (pi/(2y))**2)`.
#[must_use]
pub fn tape_swirl(twist_ratio: f64) -> f64 {
    (1.0 + (PI / (2.0 * twist_ratio)).powi(2)).sqrt()
}

/// The twisted tape's heat-transfer enhancement ratio.
#[must_use]
pub fn twisted_tape_heat(reynolds: f64, twist_ratio: f64, thickness: f64, width: f64) -> f64 {
    let y = twist_ratio.max(2.0);
    let blockage = tape_blockage(thickness, width);
    if reynolds < 2300.0 {
        return (1.0 + 0.769 / y).powf(1.0 / 3.0) * blockage;
    }
    let swirl = tape_swirl(y);
    let enhancement = swirl.powf(0.8) * blockage.powf(0.8);
    (enhancement * (1.0 + 2.0 / y) / swirl).max(1.0)
}

/// The twisted tape's friction-penalty ratio.
#[must_use]
pub fn twisted_tape_friction(reynolds: f64, twist_ratio: f64, thickness: f64, width: f64) -> f64 {
    let y = twist_ratio.max(2.0);
    let blockage = tape_blockage(thickness, width);
    if reynolds < 2300.0 {
        return (1.0 + 2.75 / y) * blockage.powf(1.29);
    }
    let swirl = tape_swirl(y);
    ((1.0 + 2.75 / y) * swirl.powf(1.8) / swirl).max(1.0)
}

/// The wire matrix's heat-transfer enhancement ratio.
#[must_use]
pub fn wire_matrix_heat(reynolds: f64, density: f64) -> f64 {
    if reynolds < 500.0 {
        return 1.0 + 5.0 * density;
    }
    if reynolds < 2300.0 {
        return 1.0 + 4.0 * density * (reynolds / 2300.0).powf(0.3);
    }
    (1.0 + 2.5 * density / (reynolds / 10_000.0).powf(0.2)).max(1.0)
}

/// The wire matrix's friction-penalty ratio.
#[must_use]
pub fn wire_matrix_friction(reynolds: f64, density: f64) -> f64 {
    let base = 1.0 + 8.0 * density;
    if reynolds < 2300.0 {
        return base * 1.5;
    }
    (base / (reynolds / 10_000.0).powf(0.15)).max(1.0)
}

/// The coiled wire's heat-transfer enhancement ratio.
#[must_use]
pub fn coiled_wire_heat(
    reynolds: f64,
    prandtl: f64,
    roughness: f64,
    pitch: f64,
    helix_angle: f64,
) -> f64 {
    if reynolds < 2300.0 {
        return 1.0 + 0.5 * roughness * 20.0;
    }
    let p_over_d = pitch / (pitch / roughness);
    let a = 2.64
        * reynolds.powf(0.036)
        * roughness.powf(0.212)
        * p_over_d.max(0.1).powf(-0.21)
        * (helix_angle / 90.0).powf(0.29)
        * prandtl.powf(-0.024);
    (1.0 + a.powi(7)).powf(1.0 / 7.0)
}

/// The coiled wire's friction-penalty ratio.
#[must_use]
pub fn coiled_wire_friction(reynolds: f64, roughness: f64, helix_angle: f64) -> f64 {
    if reynolds < 2300.0 {
        return 1.0 + 10.0 * roughness;
    }
    (1.0 + 29.1 * roughness.powf(0.67) * reynolds.powf(-0.06) * (helix_angle / 90.0).powf(0.49))
        .max(1.0)
}

/// The composed insert calculation.
///
/// # Errors
/// [`AzothError::OutOfRange`] for the spec's own bounds, which include a non-positive tape
/// thickness or coil pitch - the class's blockage and pitch ratios divide by them.
pub fn rate(request: &TubeInsertRequest) -> Result<TubeInsertModelResult> {
    let spec = &spec_gen::TUBE_INSERT_MODEL_SPEC;
    let mut warnings = Vec::new();
    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "reynolds_number" => Some(request.reynolds_number),
            "prandtl_number" => Some(request.prandtl_number),
            "twist_ratio" => request.twist_ratio,
            "tape_thickness" => request.tape_thickness,
            "tape_width_fraction" => request.tape_width_fraction,
            "wire_diameter" => request.wire_diameter,
            "wire_coil_pitch" => request.wire_coil_pitch,
            "matrix_density" => request.matrix_density,
            "roughness_ratio" => request.roughness_ratio,
            "helix_angle" => request.helix_angle,
            "plain_heat_transfer_coefficient" => Some(request.plain_heat_transfer_coefficient),
            "plain_pressure_drop" => Some(request.plain_pressure_drop),
            _ => None,
        },
        &mut warnings,
    )?;

    let twist_ratio = request.twist_ratio.unwrap_or(TWIST_RATIO);
    let thickness = request.tape_thickness.unwrap_or(TAPE_THICKNESS);
    let width = request.tape_width_fraction.unwrap_or(TAPE_WIDTH_FRACTION);
    let density = request.matrix_density.unwrap_or(MATRIX_DENSITY);
    let roughness = request.roughness_ratio.unwrap_or(ROUGHNESS_RATIO);
    let pitch = request.wire_coil_pitch.unwrap_or(WIRE_COIL_PITCH);
    let helix = request.helix_angle.unwrap_or(HELIX_ANGLE);
    // The wire diameter is carried by the class and read by no arithmetic.
    let _wire_diameter = request.wire_diameter.unwrap_or(WIRE_DIAMETER);

    let heat_ratio = match request.insert {
        InsertType::None => 1.0,
        InsertType::TwistedTape => {
            twisted_tape_heat(request.reynolds_number, twist_ratio, thickness, width)
        }
        InsertType::WireMatrix => wire_matrix_heat(request.reynolds_number, density),
        InsertType::CoiledWire => coiled_wire_heat(
            request.reynolds_number,
            request.prandtl_number,
            roughness,
            pitch,
            helix,
        ),
    };
    let friction_ratio = match request.insert {
        InsertType::None => 1.0,
        InsertType::TwistedTape => {
            twisted_tape_friction(request.reynolds_number, twist_ratio, thickness, width)
        }
        InsertType::WireMatrix => wire_matrix_friction(request.reynolds_number, density),
        InsertType::CoiledWire => coiled_wire_friction(request.reynolds_number, roughness, helix),
    };
    Ok(TubeInsertModelResult {
        heat_transfer_enhancement_ratio: heat_ratio,
        pressure_drop_penalty_ratio: friction_ratio,
        performance_evaluation_criteria: heat_ratio / friction_ratio.powf(1.0 / 3.0),
        enhanced_heat_transfer_coefficient: watts_per_square_meter_kelvin(
            request.plain_heat_transfer_coefficient * heat_ratio,
        ),
        enhanced_pressure_drop: pascals(request.plain_pressure_drop * friction_ratio),
        warnings,
    })
}

/// The flat entry point the binding and the case runner call.
///
/// # Errors
/// As [`rate`].
#[allow(clippy::too_many_arguments)]
pub fn tube_insert_model(
    insert: &str,
    reynolds_number: f64,
    prandtl_number: f64,
    twist_ratio: Option<f64>,
    tape_thickness: Option<f64>,
    tape_width_fraction: Option<f64>,
    wire_diameter: Option<f64>,
    wire_coil_pitch: Option<f64>,
    matrix_density: Option<f64>,
    roughness_ratio: Option<f64>,
    helix_angle: Option<f64>,
    plain_heat_transfer_coefficient: HeatTransfer,
    plain_pressure_drop: Pressure,
) -> Result<TubeInsertModelResult> {
    rate(&TubeInsertRequest {
        insert: insert.parse()?,
        reynolds_number,
        prandtl_number,
        twist_ratio,
        tape_thickness,
        tape_width_fraction,
        wire_diameter,
        wire_coil_pitch,
        matrix_density,
        roughness_ratio,
        helix_angle,
        plain_heat_transfer_coefficient: plain_heat_transfer_coefficient.value,
        plain_pressure_drop: plain_pressure_drop.value,
    })
}
