//! `thermal.shell_and_tube_vibration` - the tube-vibration screen.
//!
//! Transcribed from `process/mechanicaldesign/heatexchanger/VibrationAnalysis.java` at
//! NeqSim `118d9aff`. The screen fixes the natural frequency's end condition at
//! `pinned`, as `performScreening` does; the `fixed` and `clamped-pinned` branches are the
//! helper's and are measured as evidence rather than carried as inputs.
//!
//! Spec: `specs/calcs/thermal/shell_and_tube_vibration.toml`.

use azoth_core::units::{Length, MassDensity, Pressure, Velocity, meters_per_second, per_second};
use azoth_core::{AzothError, Result, apply_checks};

use crate::results::ShellAndTubeVibrationResult;
use crate::spec_gen;

/// Strouhal number the class uses when the pitch ratio is outside its correction band.
const STROUHAL_DEFAULT: f64 = 0.22;
/// Connors' constant for a triangular layout.
const CONNORS_TRIANGULAR: f64 = 3.3;
/// Connors' constant for a square layout.
const CONNORS_SQUARE: f64 = 3.4;
/// TEMA's safety factor on the critical velocity.
const FLUID_ELASTIC_SAFETY_FACTOR: f64 = 0.5;
/// The pinned-pinned eigenvalue constant.
const CN_PINNED: f64 = 9.87;

/// The screen's inputs, in SI magnitudes.
#[derive(Debug, Clone)]
pub struct VibrationRequest {
    pub tube_outside_diameter: f64,
    pub tube_inside_diameter: f64,
    pub unsupported_span: f64,
    pub tube_pitch: f64,
    pub tube_material_e: f64,
    pub tube_density: f64,
    pub crossflow_velocity: f64,
    pub shell_fluid_density: f64,
    pub tube_fluid_density: f64,
    pub shell_inside_diameter: f64,
    pub sonic_velocity: f64,
    pub damping_ratio: f64,
    pub triangular_pitch: bool,
}

/// The first natural frequency of a pinned tube span.
#[must_use]
pub fn natural_frequency(
    tube_od: f64,
    tube_id: f64,
    unsupported_span: f64,
    material_e: f64,
    tube_density: f64,
    tube_fluid_density: f64,
    shell_fluid_density: f64,
) -> f64 {
    if tube_od <= 0.0 || tube_id >= tube_od || unsupported_span <= 0.0 || material_e <= 0.0 {
        return 0.0;
    }
    let second_moment = std::f64::consts::PI / 64.0 * (tube_od.powi(4) - tube_id.powi(4));
    let tube_area = std::f64::consts::PI / 4.0 * (tube_od * tube_od - tube_id * tube_id);
    let core_area = std::f64::consts::PI / 4.0 * tube_id * tube_id;
    let metal_mass = tube_density * tube_area;
    let internal_mass = tube_fluid_density * core_area;
    let hydrodynamic_mass = shell_fluid_density * std::f64::consts::PI / 4.0 * tube_od * tube_od;
    let effective_mass = metal_mass + internal_mass + hydrodynamic_mass;
    if effective_mass <= 0.0 {
        return 0.0;
    }
    CN_PINNED / (2.0 * std::f64::consts::PI)
        * (material_e * second_moment / (effective_mass * unsupported_span.powi(4))).sqrt()
}

/// The Strouhal number, Fitzhugh's pitch correction inside `1 < r < 3`.
#[must_use]
pub fn strouhal_number(pitch_ratio: f64) -> f64 {
    if pitch_ratio > 1.0 && pitch_ratio < 3.0 {
        STROUHAL_DEFAULT * (1.0 + 0.2 / (pitch_ratio - 0.8))
    } else {
        STROUHAL_DEFAULT
    }
}

/// The vortex-shedding frequency.
#[must_use]
pub fn vortex_shedding_frequency(crossflow_velocity: f64, tube_od: f64, tube_pitch: f64) -> f64 {
    if tube_od <= 0.0 || crossflow_velocity <= 0.0 {
        return 0.0;
    }
    strouhal_number(tube_pitch / tube_od) * crossflow_velocity / tube_od
}

/// Connors' critical velocity; `f64::MAX` is the class's "not computable" sentinel.
#[must_use]
pub fn critical_velocity_connors(
    natural_frequency: f64,
    tube_od: f64,
    damping_ratio: f64,
    effective_mass_per_length: f64,
    shell_fluid_density: f64,
    triangular: bool,
) -> f64 {
    if natural_frequency <= 0.0 || tube_od <= 0.0 || shell_fluid_density <= 0.0 {
        // The class's Double.MAX_VALUE sentinel; the spec's bounds make this branch
        // unreachable (a positive OD and shell density, and a natural frequency the valid
        // geometry always gives).
        // numerics-ok: the class's Double.MAX_VALUE sentinel, unreachable on valid ranges.
        return f64::MAX;
    }
    let k = if triangular {
        CONNORS_TRIANGULAR
    } else {
        CONNORS_SQUARE
    };
    let mass_ratio = effective_mass_per_length / (shell_fluid_density * tube_od * tube_od);
    k * natural_frequency
        * tube_od
        * (2.0 * std::f64::consts::PI * damping_ratio * mass_ratio).sqrt()
}

/// The acoustic resonance frequency of the shell cavity.
#[must_use]
pub fn acoustic_frequency(shell_id: f64, sonic_velocity: f64, mode: u32) -> f64 {
    if shell_id <= 0.0 || sonic_velocity <= 0.0 || mode < 1 {
        return 0.0;
    }
    f64::from(mode) * sonic_velocity / (2.0 * shell_id)
}

/// The tube's effective mass per unit length, which Connors' criterion reads.
fn effective_mass(
    tube_od: f64,
    tube_id: f64,
    tube_density: f64,
    tube_fluid_density: f64,
    shell_fluid_density: f64,
) -> f64 {
    let tube_area = std::f64::consts::PI / 4.0 * (tube_od * tube_od - tube_id * tube_id);
    let core_area = std::f64::consts::PI / 4.0 * tube_id * tube_id;
    tube_density * tube_area
        + tube_fluid_density * core_area
        + shell_fluid_density * std::f64::consts::PI / 4.0 * tube_od * tube_od
}

/// `performScreening`.
///
/// # Errors
/// [`AzothError::InvalidInput`] for an inside diameter not smaller than the outside one, a
/// negative damping ratio or velocity, and [`AzothError::OutOfRange`] for the spec's own
/// bounds. A state whose critical velocity is not finite is refused rather than passed.
pub fn rate(request: &VibrationRequest) -> Result<ShellAndTubeVibrationResult> {
    let spec = &spec_gen::SHELL_AND_TUBE_VIBRATION_SPEC;
    let mut warnings = Vec::new();
    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "tube_outside_diameter" => Some(request.tube_outside_diameter),
            "tube_inside_diameter" => Some(request.tube_inside_diameter),
            "unsupported_span" => Some(request.unsupported_span),
            "tube_pitch" => Some(request.tube_pitch),
            "tube_material_e" => Some(request.tube_material_e),
            "tube_density" => Some(request.tube_density),
            "crossflow_velocity" => Some(request.crossflow_velocity),
            "shell_fluid_density" => Some(request.shell_fluid_density),
            "tube_fluid_density" => Some(request.tube_fluid_density),
            "shell_inside_diameter" => Some(request.shell_inside_diameter),
            "sonic_velocity" => Some(request.sonic_velocity),
            "damping_ratio" => Some(request.damping_ratio),
            _ => None,
        },
        &mut warnings,
    )?;
    validate(request)?;

    let mut result = screen(request);
    result.warnings = warnings;
    Ok(result)
}

/// `performScreening` without the input bounds and the guards, for the design calculator
/// that calls it as one stage of a longer pipeline.
#[must_use]
pub fn screen(request: &VibrationRequest) -> ShellAndTubeVibrationResult {
    let fn_value = natural_frequency(
        request.tube_outside_diameter,
        request.tube_inside_diameter,
        request.unsupported_span,
        request.tube_material_e,
        request.tube_density,
        request.tube_fluid_density,
        request.shell_fluid_density,
    );
    let fvs = vortex_shedding_frequency(
        request.crossflow_velocity,
        request.tube_outside_diameter,
        request.tube_pitch,
    );
    let vortex_ratio = if fn_value > 0.0 { fvs / fn_value } else { 0.0 };
    let vortex_critical = vortex_ratio > 0.8 && vortex_ratio < 1.2;

    let mass = effective_mass(
        request.tube_outside_diameter,
        request.tube_inside_diameter,
        request.tube_density,
        request.tube_fluid_density,
        request.shell_fluid_density,
    );
    let vcrit = critical_velocity_connors(
        fn_value,
        request.tube_outside_diameter,
        request.damping_ratio,
        mass,
        request.shell_fluid_density,
        request.triangular_pitch,
    );
    // numerics-ok: the same sentinel; an unreachable branch on the spec's valid ranges.
    let velocity_ratio = if vcrit > 0.0 && vcrit < f64::MAX {
        request.crossflow_velocity / (vcrit * FLUID_ELASTIC_SAFETY_FACTOR)
    } else {
        0.0
    };
    let fluid_elastic_critical = velocity_ratio > 1.0;

    let fac = acoustic_frequency(request.shell_inside_diameter, request.sonic_velocity, 1);
    let acoustic_ratio = if fac > 0.0 { fvs / fac } else { 0.0 };
    let acoustic_critical = acoustic_ratio > 0.8 && acoustic_ratio < 1.2;

    let passed = !vortex_critical && !fluid_elastic_critical && !acoustic_critical;
    ShellAndTubeVibrationResult {
        natural_frequency: per_second(fn_value),
        vortex_shedding_frequency: per_second(fvs),
        vortex_shedding_ratio: vortex_ratio,
        vortex_shedding_critical: vortex_critical,
        critical_velocity: meters_per_second(vcrit),
        velocity_ratio,
        fluid_elastic_critical,
        acoustic_frequency: per_second(fac),
        acoustic_ratio,
        acoustic_critical,
        passed,
        warnings: Vec::new(),
    }
}

/// The layout and damping guards the spec's bounds cannot express.
fn validate(request: &VibrationRequest) -> Result<()> {
    if request.tube_inside_diameter >= request.tube_outside_diameter {
        return Err(AzothError::invalid_input(
            "tube_inside_diameter",
            "the screen requires an inside diameter smaller than the outside diameter",
        ));
    }
    if request.damping_ratio < 0.0 {
        return Err(AzothError::invalid_input(
            "damping_ratio",
            "Connors' criterion takes the square root of the damping, so it cannot be negative",
        ));
    }
    if request.crossflow_velocity < 0.0 {
        return Err(AzothError::invalid_input(
            "crossflow_velocity",
            "a velocity cannot be negative",
        ));
    }
    Ok(())
}

/// The flat entry point the binding and the case runner call.
///
/// # Errors
/// As [`rate`].
#[allow(clippy::too_many_arguments)]
pub fn shell_and_tube_vibration(
    tube_outside_diameter: Length,
    tube_inside_diameter: Length,
    unsupported_span: Length,
    tube_pitch: Length,
    tube_material_e: Pressure,
    tube_density: MassDensity,
    crossflow_velocity: Velocity,
    shell_fluid_density: MassDensity,
    tube_fluid_density: MassDensity,
    shell_inside_diameter: Length,
    sonic_velocity: Velocity,
    damping_ratio: f64,
    triangular_pitch: bool,
) -> Result<ShellAndTubeVibrationResult> {
    rate(&VibrationRequest {
        tube_outside_diameter: tube_outside_diameter.value,
        tube_inside_diameter: tube_inside_diameter.value,
        unsupported_span: unsupported_span.value,
        tube_pitch: tube_pitch.value,
        tube_material_e: tube_material_e.value,
        tube_density: tube_density.value,
        crossflow_velocity: crossflow_velocity.value,
        shell_fluid_density: shell_fluid_density.value,
        tube_fluid_density: tube_fluid_density.value,
        shell_inside_diameter: shell_inside_diameter.value,
        sonic_velocity: sonic_velocity.value,
        damping_ratio,
        triangular_pitch,
    })
}
