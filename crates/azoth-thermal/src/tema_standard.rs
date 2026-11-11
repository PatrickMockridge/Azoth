//! `thermal.tema_standard` - TEMA configuration, class limits and layout rules.
//!
//! Transcribed from `process/mechanicaldesign/heatexchanger/TEMAStandard.java` at NeqSim
//! `f0c7436c69`: the configuration cost and head predicates, the construction-class limits,
//! the tube-size, pitch and baffle tables, and the six numeric rules the design calculator
//! reads.
//!
//! Spec: `specs/calcs/thermal/tema_standard.toml`.

use azoth_core::units::{Length, meters};
use azoth_core::{AzothError, Result, apply_checks};

use crate::results::TemaStandardResult;
use crate::spec_gen;

/// The stationary head letter.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FrontHead {
    /// Channel with removable cover.
    A,
    /// Bonnet.
    B,
    /// Channel integral with the tubesheet.
    C,
    /// Channel integral with the tubesheet, large bore.
    N,
    /// Special high-pressure closure.
    D,
}

/// The shell letter.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Shell {
    /// One-pass shell.
    E,
    /// Two-pass shell with a longitudinal baffle.
    F,
    /// Split flow.
    G,
    /// Double split flow.
    H,
    /// Divided flow.
    J,
    /// Kettle reboiler.
    K,
    /// Cross flow.
    X,
}

/// The rear head letter.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RearHead {
    /// Fixed tubesheet like a B head.
    L,
    /// Fixed tubesheet like an A head.
    M,
    /// Fixed tubesheet like an N head.
    N,
    /// Outside packed floating head.
    P,
    /// Floating head with a backing device.
    S,
    /// Pull-through floating head.
    T,
    /// U-tube bundle.
    U,
    /// Externally sealed floating tubesheet.
    W,
}

/// The construction class.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TemaClass {
    /// Severe service.
    R,
    /// Moderate service.
    C,
    /// Light service.
    B,
}

/// The standard tube size.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TubeSize {
    /// 3/8 inch outside diameter.
    Tube3_8Inch,
    /// 1/2 inch outside diameter.
    Tube1_2Inch,
    /// 5/8 inch outside diameter.
    Tube5_8Inch,
    /// 3/4 inch outside diameter.
    Tube3_4Inch,
    /// 1 inch outside diameter.
    Tube1Inch,
}

/// The tube layout pattern.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PitchPattern {
    /// Triangular 30 degrees.
    Triangular30,
    /// Rotated triangular 60 degrees.
    Triangular60,
    /// Square 90 degrees.
    Square90,
    /// Rotated square 45 degrees.
    Square45,
}

/// The baffle type.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BaffleType {
    /// Single segmental.
    SingleSegmental,
    /// Double segmental.
    DoubleSegmental,
    /// Triple segmental.
    TripleSegmental,
    /// No tubes in the window.
    NoTubesInWindow,
    /// Disc and doughnut.
    DiscAndDoughnut,
    /// Rod baffles.
    RodBaffles,
}

macro_rules! parse_enum {
    ($name:ident, $field:literal, { $($text:literal => $variant:path),+ $(,)? }) => {
        impl std::str::FromStr for $name {
            type Err = AzothError;

            fn from_str(value: &str) -> Result<Self> {
                match value {
                    $($text => Ok($variant),)+
                    other => Err(AzothError::invalid_input(
                        $field,
                        format!("{other} is not one of the declared values"),
                    )),
                }
            }
        }
    };
}

parse_enum!(FrontHead, "front_head", {
    "A" => FrontHead::A, "B" => FrontHead::B, "C" => FrontHead::C,
    "N" => FrontHead::N, "D" => FrontHead::D,
});
parse_enum!(Shell, "shell", {
    "E" => Shell::E, "F" => Shell::F, "G" => Shell::G, "H" => Shell::H,
    "J" => Shell::J, "K" => Shell::K, "X" => Shell::X,
});
parse_enum!(RearHead, "rear_head", {
    "L" => RearHead::L, "M" => RearHead::M, "N" => RearHead::N, "P" => RearHead::P,
    "S" => RearHead::S, "T" => RearHead::T, "U" => RearHead::U, "W" => RearHead::W,
});
parse_enum!(TemaClass, "tema_class", {
    "R" => TemaClass::R, "C" => TemaClass::C, "B" => TemaClass::B,
});
parse_enum!(TubeSize, "tube_size", {
    "tube_3_8_inch" => TubeSize::Tube3_8Inch,
    "tube_1_2_inch" => TubeSize::Tube1_2Inch,
    "tube_5_8_inch" => TubeSize::Tube5_8Inch,
    "tube_3_4_inch" => TubeSize::Tube3_4Inch,
    "tube_1_inch" => TubeSize::Tube1Inch,
});
parse_enum!(PitchPattern, "pitch_pattern", {
    "triangular_30" => PitchPattern::Triangular30,
    "triangular_60" => PitchPattern::Triangular60,
    "square_90" => PitchPattern::Square90,
    "square_45" => PitchPattern::Square45,
});
parse_enum!(BaffleType, "baffle_type", {
    "single_segmental" => BaffleType::SingleSegmental,
    "double_segmental" => BaffleType::DoubleSegmental,
    "triple_segmental" => BaffleType::TripleSegmental,
    "no_tubes_in_window" => BaffleType::NoTubesInWindow,
    "disc_and_doughnut" => BaffleType::DiscAndDoughnut,
    "rod_baffles" => BaffleType::RodBaffles,
});

/// The composed request, in SI magnitudes.
#[derive(Debug, Clone)]
pub struct TemaStandardRequest {
    pub front_head: FrontHead,
    pub shell: Shell,
    pub rear_head: RearHead,
    pub tema_class: TemaClass,
    pub tube_size: TubeSize,
    pub pitch_pattern: PitchPattern,
    pub baffle_type: BaffleType,
    pub shell_inside_diameter: f64,
    pub tube_pitch: f64,
    pub tube_passes: f64,
    pub tube_material: Option<String>,
    pub needs_cleanable: bool,
    pub large_temperature_difference: bool,
    pub high_pressure: bool,
    pub hazardous: bool,
}

/// The configuration cost factor, adding the class's steps in its order.
#[must_use]
pub fn configuration_cost_factor(front: FrontHead, shell: Shell, rear: RearHead) -> f64 {
    let mut factor = 1.0;
    match front {
        FrontHead::D => factor += 0.4,
        FrontHead::A => factor += 0.15,
        FrontHead::C | FrontHead::N => factor += 0.1,
        FrontHead::B => {}
    }
    match shell {
        Shell::K => factor += 0.3,
        Shell::F => factor += 0.2,
        _ => {}
    }
    match rear {
        RearHead::T => factor += 0.25,
        RearHead::S => factor += 0.2,
        RearHead::P | RearHead::W => factor += 0.15,
        _ => {}
    }
    factor
}

/// The shell letter's pressure drop factor.
#[must_use]
pub fn shell_pressure_drop_factor(shell: Shell) -> f64 {
    match shell {
        Shell::E => 1.0,
        Shell::F => 0.8,
        Shell::G => 0.6,
        Shell::H => 0.5,
        Shell::J => 0.65,
        Shell::K => 0.7,
        Shell::X => 0.3,
    }
}

/// Whether the rear head is floating.
#[must_use]
pub fn rear_is_floating(rear: RearHead) -> bool {
    matches!(rear, RearHead::P | RearHead::S | RearHead::T | RearHead::W)
}

/// The class's three limit values: cost factor, minimum wall, minimum corrosion allowance.
#[must_use]
pub fn class_limits(class: TemaClass) -> (f64, f64, f64) {
    match class {
        TemaClass::R => (1.0, 0.0127, 0.00476),
        TemaClass::C => (0.8, 0.0095, 0.00318),
        TemaClass::B => (0.6, 0.00635, 0.00211),
    }
}

/// The tube size's outside diameter.
#[must_use]
pub fn tube_outside_diameter(size: TubeSize) -> f64 {
    match size {
        TubeSize::Tube3_8Inch => 0.009525,
        TubeSize::Tube1_2Inch => 0.0127,
        TubeSize::Tube5_8Inch => 0.015875,
        TubeSize::Tube3_4Inch => 0.01905,
        TubeSize::Tube1Inch => 0.0254,
    }
}

/// The pattern's layout angle, minimum pitch ratio and heat-transfer factor.
#[must_use]
pub fn pitch_limits(pattern: PitchPattern) -> (f64, f64, f64) {
    match pattern {
        PitchPattern::Triangular30 => (30.0, 1.25, 1.1),
        PitchPattern::Triangular60 => (60.0, 1.25, 1.15),
        PitchPattern::Square90 => (90.0, 1.25, 1.0),
        PitchPattern::Square45 => (45.0, 1.25, 1.05),
    }
}

/// The baffle type's heat-transfer and pressure-drop factors.
#[must_use]
pub fn baffle_factors(baffle: BaffleType) -> (f64, f64) {
    match baffle {
        BaffleType::SingleSegmental => (0.75, 1.0),
        BaffleType::DoubleSegmental => (0.5, 0.6),
        BaffleType::TripleSegmental => (0.35, 0.4),
        BaffleType::NoTubesInWindow => (0.6, 0.5),
        BaffleType::DiscAndDoughnut => (0.5, 0.55),
        BaffleType::RodBaffles => (0.2, 0.3),
    }
}

/// The class's recommended designation from the four application flags.
#[must_use]
pub fn recommend_configuration(
    needs_cleanable: bool,
    large_temperature_difference: bool,
    high_pressure: bool,
    hazardous: bool,
) -> &'static str {
    if hazardous {
        return if large_temperature_difference {
            "AEW"
        } else {
            "BEM"
        };
    }
    if high_pressure {
        return if large_temperature_difference {
            "NEU"
        } else {
            "NEN"
        };
    }
    if large_temperature_difference {
        return if needs_cleanable { "AES" } else { "AEU" };
    }
    if needs_cleanable { "AES" } else { "BEM" }
}

/// The minimum tube pitch, `OD * ratio`.
#[must_use]
pub fn min_tube_pitch(diameter: f64, pattern: PitchPattern) -> f64 {
    diameter * pitch_limits(pattern).1
}

/// The truncated tube-count estimate.
#[must_use]
pub fn estimate_tube_count(
    shell_inside_diameter: f64,
    tube_diameter: f64,
    tube_pitch: f64,
    pattern: PitchPattern,
    tube_passes: i32,
) -> f64 {
    let effective_diameter = shell_inside_diameter - 2.0 * tube_diameter;
    let pass_lane_factor = 1.0 - 0.05 * f64::from(tube_passes - 1);
    let effective_area =
        std::f64::consts::PI * (effective_diameter / 2.0).powi(2) * pass_lane_factor;
    let area_per_tube = match pattern {
        PitchPattern::Triangular30 | PitchPattern::Triangular60 => {
            tube_pitch * tube_pitch * 3.0f64.sqrt() / 2.0
        }
        PitchPattern::Square90 | PitchPattern::Square45 => tube_pitch * tube_pitch,
    };
    (effective_area / area_per_tube).trunc()
}

/// The material factor for the unsupported span, matched on the upper-cased name.
#[must_use]
pub fn material_factor(material: Option<&str>) -> f64 {
    let Some(name) = material else {
        return 1.0;
    };
    let upper = name.to_uppercase();
    if upper.contains("STAINLESS") || upper.contains("SS") {
        0.95
    } else if upper.contains("COPPER") || upper.contains("CU") {
        0.75
    } else if upper.contains("TITANIUM") || upper.contains("TI") {
        0.9
    } else {
        1.0
    }
}

/// The class's maximum unsupported tube span.
#[must_use]
pub fn max_unsupported_span(tube_diameter: f64, material: Option<&str>) -> f64 {
    let base = if tube_diameter <= 0.0127 {
        0.660
    } else if tube_diameter <= 0.01905 {
        0.889
    } else if tube_diameter <= 0.0254 {
        1.067
    } else {
        1.321
    };
    base * material_factor(material)
}

/// The minimum baffle spacing.
#[must_use]
pub fn min_baffle_spacing(shell_inside_diameter: f64, class: TemaClass) -> f64 {
    let floor = if class == TemaClass::B {
        0.0381
    } else {
        0.0508
    };
    (shell_inside_diameter / 5.0).max(floor)
}

/// The composed TEMA calculation.
///
/// # Errors
/// [`AzothError::InvalidInput`] for a non-integral pass count, and [`AzothError::OutOfRange`]
/// for the spec's own bounds.
pub fn rate(request: &TemaStandardRequest) -> Result<TemaStandardResult> {
    let spec = &spec_gen::TEMA_STANDARD_SPEC;
    let mut warnings = Vec::new();
    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "shell_inside_diameter" => Some(request.shell_inside_diameter),
            "tube_pitch" => Some(request.tube_pitch),
            "tube_passes" => Some(request.tube_passes),
            _ => None,
        },
        &mut warnings,
    )?;
    if !request.tube_passes.is_finite()
        || (request.tube_passes - request.tube_passes.round()).abs() > 0.0
        || request.tube_passes < 1.0
    {
        return Err(AzothError::invalid_input(
            "tube_passes",
            format!("{} is not a whole pass count", request.tube_passes),
        ));
    }
    let passes = request.tube_passes.round() as i32;
    let diameter = tube_outside_diameter(request.tube_size);
    let (angle, ratio, pitch_heat) = pitch_limits(request.pitch_pattern);
    let (baffle_heat, baffle_drop) = baffle_factors(request.baffle_type);
    let (class_cost, min_wall, min_corrosion) = class_limits(request.tema_class);
    let removable = rear_is_floating(request.rear_head) || request.rear_head == RearHead::U;
    Ok(TemaStandardResult {
        configuration_cost_factor: configuration_cost_factor(
            request.front_head,
            request.shell,
            request.rear_head,
        ),
        shell_pressure_drop_factor: shell_pressure_drop_factor(request.shell),
        bundle_removable: removable,
        good_thermal_expansion: removable,
        class_cost_factor: class_cost,
        min_tube_wall_thickness: meters(min_wall),
        min_corrosion_allowance: meters(min_corrosion),
        tube_outer_diameter: meters(diameter),
        tube_outer_diameter_inch: meters(diameter),
        layout_angle: angle,
        min_pitch_ratio: ratio,
        pitch_heat_transfer_factor: pitch_heat,
        baffle_heat_transfer_factor: baffle_heat,
        baffle_pressure_drop_factor: baffle_drop,
        min_tube_pitch: meters(min_tube_pitch(diameter, request.pitch_pattern)),
        estimated_tube_count: estimate_tube_count(
            request.shell_inside_diameter,
            diameter,
            request.tube_pitch,
            request.pitch_pattern,
            passes,
        ),
        max_unsupported_span: meters(max_unsupported_span(
            diameter,
            request.tube_material.as_deref(),
        )),
        min_baffle_spacing: meters(min_baffle_spacing(
            request.shell_inside_diameter,
            request.tema_class,
        )),
        max_baffle_spacing: meters(request.shell_inside_diameter),
        recommended_designation: recommend_configuration(
            request.needs_cleanable,
            request.large_temperature_difference,
            request.high_pressure,
            request.hazardous,
        )
        .to_string(),
        warnings,
    })
}

/// The flat entry point the binding and the case runner call.
///
/// # Errors
/// As [`rate`].
#[allow(clippy::too_many_arguments)]
pub fn tema_standard(
    front_head: &str,
    shell: &str,
    rear_head: &str,
    tema_class: &str,
    tube_size: &str,
    pitch_pattern: &str,
    baffle_type: &str,
    shell_inside_diameter: Length,
    tube_pitch: Length,
    tube_passes: f64,
    needs_cleanable: bool,
    large_temperature_difference: bool,
    high_pressure: bool,
    hazardous: bool,
    tube_material: Option<&str>,
) -> Result<TemaStandardResult> {
    rate(&TemaStandardRequest {
        front_head: front_head.parse()?,
        shell: shell.parse()?,
        rear_head: rear_head.parse()?,
        tema_class: tema_class.parse()?,
        tube_size: tube_size.parse()?,
        pitch_pattern: pitch_pattern.parse()?,
        baffle_type: baffle_type.parse()?,
        shell_inside_diameter: shell_inside_diameter.value,
        tube_pitch: tube_pitch.value,
        tube_passes,
        tube_material: tube_material.map(str::to_owned),
        needs_cleanable,
        large_temperature_difference,
        high_pressure,
        hazardous,
    })
}
