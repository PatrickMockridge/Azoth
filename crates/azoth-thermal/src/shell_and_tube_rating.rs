//! `thermal.shell_and_tube_rating` - the shell-and-tube thermal-hydraulic rating.
//!
//! Transcribed from `process/mechanicaldesign/heatexchanger/ThermalDesignCalculator.java`
//! at NeqSim `f0c7436c69`. The port implements `calculateStrict()`'s refusal, not
//! `calculate()`'s silent zeros: a geometry the class would leave at zero is an error
//! here, because a zero coefficient is not an answer.
//!
//! ```text
//! h_tube  = Nu * k_tube / d_f
//! h_shell = h_ideal * Jc * Jl * Jb * Js * Jr          (Bell-Delaware)
//! 1/U     = 1/h_shell + R_f,shell + d_o ln(d_o/d_i)/(2 k_w)
//!           + (d_o/d_i) R_f,tube + d_o/(d_f h_tube)
//! ```
//!
//! Spec: `specs/calcs/thermal/shell_and_tube_rating.toml`.

use azoth_core::units::{
    Area, DynamicViscosity, Length, MassDensity, MassRate, SpecificHeatCapacity,
    ThermalConductivity, meters, meters_per_second, pascals, square_meters,
    watts_per_square_meter_kelvin,
};
use azoth_core::{AzothError, Result, apply_checks};

use crate::bell_delaware;
use crate::results::ShellAndTubeRatingResult;
use crate::spec_gen;

/// The shell-side analysis method.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ShellSideMethod {
    /// Kern's equivalent-diameter correlation.
    Kern,
    /// The Bell-Delaware correction-factor method.
    BellDelaware,
}

impl std::str::FromStr for ShellSideMethod {
    type Err = AzothError;

    fn from_str(value: &str) -> Result<Self> {
        match value {
            "kern" => Ok(Self::Kern),
            "bell_delaware" => Ok(Self::BellDelaware),
            other => Err(AzothError::invalid_input(
                "shell_side_method",
                format!("{other} is not one of kern or bell_delaware"),
            )),
        }
    }
}

/// The geometry and both fluids, in SI magnitudes.
///
/// Every field is a required input except the optional group at the end, which carries
/// the class's own default when absent.
#[derive(Debug, Clone)]
pub struct RatingRequest {
    pub tube_outside_diameter: f64,
    pub tube_inside_diameter: f64,
    pub tube_length: f64,
    pub tube_count: f64,
    pub tube_passes: f64,
    pub tube_pitch: f64,
    pub triangular_pitch: bool,
    pub shell_inside_diameter: f64,
    pub baffle_spacing: f64,
    pub baffle_count: f64,
    pub baffle_cut: f64,
    pub tube_wall_conductivity: f64,
    pub tube_density: f64,
    pub tube_viscosity: f64,
    pub tube_cp: f64,
    pub tube_conductivity: f64,
    pub tube_mass_flow: f64,
    pub shell_density: f64,
    pub shell_viscosity: f64,
    pub shell_cp: f64,
    pub shell_conductivity: f64,
    pub shell_mass_flow: f64,
    pub fouling_tube: f64,
    pub fouling_shell: f64,
    pub shell_side_method: ShellSideMethod,
    pub tube_fouling_layer_thickness: Option<f64>,
    pub tube_fouling_layer_conductivity: Option<f64>,
    pub shell_viscosity_wall: Option<f64>,
    pub has_sealing: Option<bool>,
    pub sealing_pairs: Option<f64>,
    pub tube_to_baffle_clearance: Option<f64>,
    pub shell_to_baffle_clearance: Option<f64>,
    pub bypass_area: Option<f64>,
}

/// The class's constructor defaults, for a caller that states only what differs.
impl Default for RatingRequest {
    fn default() -> Self {
        Self {
            tube_outside_diameter: 0.01905,
            tube_inside_diameter: 0.01483,
            tube_length: 6.096,
            tube_count: 100.0,
            tube_passes: 2.0,
            tube_pitch: 0.02381,
            triangular_pitch: true,
            shell_inside_diameter: 0.5,
            baffle_spacing: 0.2,
            baffle_count: 10.0,
            baffle_cut: 0.25,
            tube_wall_conductivity: 52.0,
            tube_density: 800.0,
            tube_viscosity: 5.0e-4,
            tube_cp: 2500.0,
            tube_conductivity: 0.13,
            tube_mass_flow: 10.0,
            shell_density: 1.5,
            shell_viscosity: 1.5e-5,
            shell_cp: 1100.0,
            shell_conductivity: 0.025,
            shell_mass_flow: 5.0,
            fouling_tube: 1.76e-4,
            fouling_shell: 1.76e-4,
            shell_side_method: ShellSideMethod::Kern,
            tube_fouling_layer_thickness: None,
            tube_fouling_layer_conductivity: None,
            shell_viscosity_wall: None,
            has_sealing: None,
            sealing_pairs: None,
            tube_to_baffle_clearance: None,
            shell_to_baffle_clearance: None,
            bypass_area: None,
        }
    }
}

/// The class's `calcDarcyFriction`: `64/Re` below 2300, the smooth Swamee-Jain form above.
#[must_use]
pub fn darcy_friction(re: f64) -> f64 {
    if re <= 0.0 {
        return 0.0;
    }
    if re < 2300.0 {
        return 64.0 / re;
    }
    let log_term = (5.74 / re.powf(0.9)).log10();
    0.25 / (log_term * log_term)
}

/// `HeatTransferCoefficientCalculator.calculateGnielinskiNusselt`.
#[must_use]
pub fn gnielinski_nusselt(re: f64, pr: f64, friction_factor: f64) -> f64 {
    if re < 3000.0 || pr < 0.5 {
        return laminar_nusselt();
    }
    let f = friction_factor;
    let numerator = (f / 8.0) * (re - 1000.0) * pr;
    let denominator = 1.0 + 12.7 * (f / 8.0).sqrt() * (pr.powf(2.0 / 3.0) - 1.0);
    (numerator / denominator).max(3.66)
}

/// `HeatTransferCalculator.calculateLaminarNusselt(true)`, constant wall temperature.
#[must_use]
pub fn laminar_nusselt() -> f64 {
    3.66
}

/// A count that must be a positive integer, as the class's `int` setters are.
fn integral_count(value: f64, field: &str) -> Result<usize> {
    if !value.is_finite() || value <= 0.0 || (value - value.round()).abs() > 0.0 {
        return Err(AzothError::invalid_input(
            field,
            format!("{value} is not a positive whole count"),
        ));
    }
    Ok(value.round() as usize)
}

fn optional_integral(value: Option<f64>, field: &str, default: usize) -> Result<usize> {
    match value {
        Some(value) if value < 0.0 || !value.is_finite() || (value - value.round()).abs() > 0.0 => {
            Err(AzothError::invalid_input(
                field,
                format!("{value} is not a non-negative whole count"),
            ))
        }
        Some(value) => Ok(value.round() as usize),
        None => Ok(default),
    }
}

/// The thermal-hydraulic rating.
///
/// # Errors
/// [`AzothError::InvalidInput`] for every state the class's `validateStrictInputs`
/// refuses, and [`AzothError::OutOfRange`] for the spec's own input bounds.
#[allow(clippy::too_many_lines)] // one transcription of the class's calculateStrict chain
pub fn rate(request: &RatingRequest) -> Result<ShellAndTubeRatingResult> {
    let spec = &spec_gen::SHELL_AND_TUBE_RATING_SPEC;
    let mut warnings = Vec::new();

    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "tube_outside_diameter" => Some(request.tube_outside_diameter),
            "tube_inside_diameter" => Some(request.tube_inside_diameter),
            "tube_length" => Some(request.tube_length),
            "tube_count" => Some(request.tube_count),
            "tube_passes" => Some(request.tube_passes),
            "tube_pitch" => Some(request.tube_pitch),
            "shell_inside_diameter" => Some(request.shell_inside_diameter),
            "baffle_spacing" => Some(request.baffle_spacing),
            "baffle_count" => Some(request.baffle_count),
            "baffle_cut" => Some(request.baffle_cut),
            "tube_wall_conductivity" => Some(request.tube_wall_conductivity),
            "tube_density" => Some(request.tube_density),
            "tube_viscosity" => Some(request.tube_viscosity),
            "tube_cp" => Some(request.tube_cp),
            "tube_conductivity" => Some(request.tube_conductivity),
            "tube_mass_flow" => Some(request.tube_mass_flow),
            "shell_density" => Some(request.shell_density),
            "shell_viscosity" => Some(request.shell_viscosity),
            "shell_cp" => Some(request.shell_cp),
            "shell_conductivity" => Some(request.shell_conductivity),
            "shell_mass_flow" => Some(request.shell_mass_flow),
            "fouling_tube" => Some(request.fouling_tube),
            "fouling_shell" => Some(request.fouling_shell),
            "tube_fouling_layer_thickness" => request.tube_fouling_layer_thickness,
            "tube_fouling_layer_conductivity" => request.tube_fouling_layer_conductivity,
            "sealing_pairs" => request.sealing_pairs,
            "tube_to_baffle_clearance" => request.tube_to_baffle_clearance,
            "shell_to_baffle_clearance" => request.shell_to_baffle_clearance,
            "bypass_area" => request.bypass_area,
            _ => None,
        },
        &mut warnings,
    )?;

    validate_strict(request)?;

    let mut result = evaluate(request)?;
    result.warnings = warnings;
    Ok(result)
}

/// The class's permissive arithmetic, without the input bounds and the strict guards.
///
/// `calculateZones` and the design calculator both call the class's `calculate()`, so their
/// kernels compose this rather than transcribing the rating a second time. `rate` is this
/// function plus the spec's input bounds, the strict validation and the warnings.
///
/// # Errors
/// [`AzothError::InvalidInput`] for a count that is not a positive whole number.
#[allow(clippy::too_many_lines)] // one transcription of the class's calculate chain
pub fn evaluate(request: &RatingRequest) -> Result<ShellAndTubeRatingResult> {
    let tube_count = integral_count(request.tube_count, "tube_count")?;
    let tube_passes = integral_count(request.tube_passes, "tube_passes")?;
    let baffle_count = integral_count(request.baffle_count, "baffle_count")?;
    let sealing_pairs = optional_integral(request.sealing_pairs, "sealing_pairs", 0)?;
    let has_sealing = request.has_sealing.unwrap_or(false);
    let deposit_thickness = request.tube_fouling_layer_thickness.unwrap_or(0.0);
    let deposit_conductivity = request.tube_fouling_layer_conductivity.unwrap_or(1.0);
    let tube_to_baffle = request.tube_to_baffle_clearance.unwrap_or(0.0004);
    let shell_to_baffle = request.shell_to_baffle_clearance.unwrap_or(0.003);
    let mut bypass = request.bypass_area.unwrap_or(0.0);

    let tube_od = request.tube_outside_diameter;
    let tube_id = request.tube_inside_diameter;
    let tube_length = request.tube_length;
    let tube_pitch = request.tube_pitch;
    let triangular = request.triangular_pitch;
    let shell_id = request.shell_inside_diameter;
    let baffle_spacing = request.baffle_spacing;
    let baffle_cut = request.baffle_cut;
    let tube_density = request.tube_density;
    let tube_viscosity = request.tube_viscosity;
    let tube_cp = request.tube_cp;
    let tube_conductivity = request.tube_conductivity;
    let tube_mass_flow = request.tube_mass_flow;
    let shell_density = request.shell_density;
    let shell_viscosity = request.shell_viscosity;
    let shell_cp = request.shell_cp;
    let shell_conductivity = request.shell_conductivity;
    let shell_mass_flow = request.shell_mass_flow;

    let effective_id = tube_id - 2.0 * deposit_thickness;
    let outside_area = std::f64::consts::PI * tube_od * tube_length * tube_count as f64;

    // Tube side: Gnielinski above 3000, the class's interpolation in the transition band,
    // and its laminar constant below 2300.
    let area_per_tube = std::f64::consts::PI / 4.0 * effective_id * effective_id;
    let tubes_per_pass = tube_count / tube_passes;
    let tube_flow_area = tubes_per_pass as f64 * area_per_tube;
    let tube_velocity = tube_mass_flow / (tube_density * tube_flow_area);
    let tube_re = tube_density * tube_velocity * effective_id / tube_viscosity;
    let tube_pr = tube_cp * tube_viscosity / tube_conductivity;
    let tube_nu = if tube_re > 3000.0 {
        gnielinski_nusselt(tube_re, tube_pr, darcy_friction(tube_re))
    } else if tube_re > 2300.0 {
        let laminar = laminar_nusselt();
        let turbulent = gnielinski_nusselt(3000.0, tube_pr, darcy_friction(3000.0));
        let fraction = (tube_re - 2300.0) / 700.0;
        laminar + fraction * (turbulent - laminar)
    } else {
        laminar_nusselt()
    };
    let tube_htc = tube_nu * tube_conductivity / effective_id;

    // Shell side. The wall viscosity is the bulk one when absent or negative, as the
    // class's `shellViscosityWall > 0` test says.
    let wall_viscosity = match request.shell_viscosity_wall {
        Some(value) if value > 0.0 => value,
        _ => shell_viscosity,
    };
    let shell_area = bell_delaware::crossflow_area(shell_id, baffle_spacing, tube_od, tube_pitch);
    let shell_mass_flux = shell_mass_flow / shell_area;
    let shell_velocity = shell_mass_flux / shell_density;
    let (shell_htc, shell_re, jc, jl, jb, js, jr) = match request.shell_side_method {
        ShellSideMethod::Kern => {
            let equivalent_diameter =
                bell_delaware::shell_equivalent_diameter(tube_od, tube_pitch, triangular);
            let re = shell_mass_flux * equivalent_diameter / shell_viscosity;
            let htc = bell_delaware::kern_shell_side_htc(
                shell_mass_flux,
                equivalent_diameter,
                shell_viscosity,
                shell_cp,
                shell_conductivity,
                wall_viscosity,
            );
            (htc, re, 1.0, 1.0, 1.0, 1.0, 1.0)
        }
        ShellSideMethod::BellDelaware => {
            let re = shell_mass_flux * tube_od / shell_viscosity;
            let ideal = bell_delaware::ideal_crossflow_htc(
                shell_mass_flux,
                tube_od,
                shell_viscosity,
                shell_cp,
                shell_conductivity,
                wall_viscosity,
                triangular,
            );
            let jc = bell_delaware::jc(baffle_cut);
            let jl = bell_delaware::jl(
                tube_to_baffle,
                shell_to_baffle,
                shell_area,
                tube_count,
                tube_od,
                baffle_spacing,
            );
            if bypass <= 0.0 {
                bypass = bell_delaware::bypass_area(shell_id, shell_id * 0.9, baffle_spacing);
            }
            let rows =
                bell_delaware::tube_rows_crossflow(shell_id, baffle_cut, tube_pitch, triangular);
            let jb = bell_delaware::jb(bypass, shell_area, has_sealing, sealing_pairs, rows);
            let js =
                bell_delaware::js(baffle_spacing, baffle_spacing, baffle_spacing, baffle_count);
            let jr = bell_delaware::jr(re, rows);
            let htc = bell_delaware::corrected_htc(ideal, jc, jl, jb, js, jr);
            (htc, re, jc, jl, jb, js, jr)
        }
    };

    // Overall coefficient on the outside area.
    let do_by_di = tube_od / tube_id;
    let r_shell = 1.0 / shell_htc;
    let r_foul_shell = request.fouling_shell;
    let r_wall = tube_od * do_by_di.ln() / (2.0 * request.tube_wall_conductivity);
    let r_foul_tube = do_by_di * request.fouling_tube;
    let r_tube = tube_od / effective_id / tube_htc;
    let mut total_resistance = r_shell + r_foul_shell + r_wall + r_foul_tube + r_tube;
    let deposit_resistance = if deposit_thickness > 0.0 {
        tube_od * (tube_id / effective_id).ln() / (2.0 * deposit_conductivity)
    } else {
        0.0
    };
    total_resistance += deposit_resistance;
    let overall_u = 1.0 / total_resistance;

    // Tube-side pressure drop: friction plus 2.5 velocity heads per pass.
    let friction = darcy_friction(tube_re);
    let tube_dp_friction =
        friction * tube_length / effective_id * tube_density * tube_velocity * tube_velocity / 2.0;
    let tube_dp_return = 2.5 * tube_density * tube_velocity * tube_velocity / 2.0;
    let tube_dp = tube_passes as f64 * (tube_dp_friction + tube_dp_return);

    // Shell-side pressure drop.
    let shell_dp = match request.shell_side_method {
        ShellSideMethod::Kern => {
            let equivalent_diameter =
                bell_delaware::shell_equivalent_diameter(tube_od, tube_pitch, triangular);
            bell_delaware::kern_shell_side_pressure_drop(
                shell_mass_flux,
                equivalent_diameter,
                shell_id,
                baffle_count,
                shell_density,
                shell_viscosity,
                wall_viscosity,
            )
        }
        ShellSideMethod::BellDelaware => {
            let rows =
                bell_delaware::tube_rows_crossflow(shell_id, baffle_cut, tube_pitch, triangular);
            let ideal_dp = bell_delaware::ideal_crossflow_pressure_drop(
                rows,
                shell_mass_flux,
                shell_density,
                shell_re,
                triangular,
            );
            let rl = bell_delaware::rl(
                tube_to_baffle,
                shell_to_baffle,
                shell_area,
                tube_count,
                tube_od,
                baffle_spacing,
            );
            if bypass <= 0.0 {
                bypass = bell_delaware::bypass_area(shell_id, shell_id * 0.9, baffle_spacing);
            }
            let rb = bell_delaware::rb(bypass, shell_area, has_sealing, sealing_pairs, rows);
            let window_area = std::f64::consts::PI / 4.0 * shell_id * shell_id * baffle_cut * 0.5;
            let rows_window = (rows / 3).max(1);
            let window_dp = bell_delaware::window_pressure_drop(
                window_area,
                shell_mass_flow,
                shell_density,
                rows_window,
            );
            let crossflow_dp = ideal_dp * (baffle_count as f64 - 1.0) * rl * rb;
            let windows_dp = baffle_count as f64 * window_dp * rl;
            let end_zone_dp = 2.0 * ideal_dp * (1.0 + rows as f64 / 5.0) * rb;
            crossflow_dp + windows_dp + end_zone_dp
        }
    };

    Ok(ShellAndTubeRatingResult {
        tube_side_heat_transfer_coefficient: watts_per_square_meter_kelvin(tube_htc),
        shell_side_heat_transfer_coefficient: watts_per_square_meter_kelvin(shell_htc),
        overall_heat_transfer_coefficient: watts_per_square_meter_kelvin(overall_u),
        tube_side_pressure_drop: pascals(tube_dp),
        shell_side_pressure_drop: pascals(shell_dp),
        tube_side_velocity: meters_per_second(tube_velocity),
        shell_side_velocity: meters_per_second(shell_velocity),
        tube_side_reynolds: tube_re,
        shell_side_reynolds: shell_re,
        effective_tube_inside_diameter: meters(effective_id),
        outside_heat_transfer_area: square_meters(outside_area),
        tube_fouling_layer_resistance_outside: deposit_resistance,
        jc,
        jl,
        jb,
        js,
        jr,
        warnings: Vec::new(),
    })
}

/// The flat entry point the binding and the spec's case runner call.
///
/// One parameter per declared input, in the spec's own order; the optional group is last.
///
/// # Errors
/// As [`rate`].
#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
pub fn shell_and_tube_rating(
    tube_outside_diameter: Length,
    tube_inside_diameter: Length,
    tube_length: Length,
    tube_count: f64,
    tube_passes: f64,
    tube_pitch: Length,
    triangular_pitch: bool,
    shell_inside_diameter: Length,
    baffle_spacing: Length,
    baffle_count: f64,
    baffle_cut: f64,
    tube_wall_conductivity: ThermalConductivity,
    tube_density: MassDensity,
    tube_viscosity: DynamicViscosity,
    tube_cp: SpecificHeatCapacity,
    tube_conductivity: ThermalConductivity,
    tube_mass_flow: MassRate,
    shell_density: MassDensity,
    shell_viscosity: DynamicViscosity,
    shell_cp: SpecificHeatCapacity,
    shell_conductivity: ThermalConductivity,
    shell_mass_flow: MassRate,
    fouling_tube: f64,
    fouling_shell: f64,
    shell_side_method: ShellSideMethod,
    tube_fouling_layer_thickness: Option<Length>,
    tube_fouling_layer_conductivity: Option<ThermalConductivity>,
    shell_viscosity_wall: Option<DynamicViscosity>,
    has_sealing: Option<bool>,
    sealing_pairs: Option<f64>,
    tube_to_baffle_clearance: Option<Length>,
    shell_to_baffle_clearance: Option<Length>,
    bypass_area: Option<Area>,
) -> Result<ShellAndTubeRatingResult> {
    rate(&RatingRequest {
        tube_outside_diameter: tube_outside_diameter.value,
        tube_inside_diameter: tube_inside_diameter.value,
        tube_length: tube_length.value,
        tube_count,
        tube_passes,
        tube_pitch: tube_pitch.value,
        triangular_pitch,
        shell_inside_diameter: shell_inside_diameter.value,
        baffle_spacing: baffle_spacing.value,
        baffle_count,
        baffle_cut,
        tube_wall_conductivity: tube_wall_conductivity.value,
        tube_density: tube_density.value,
        tube_viscosity: tube_viscosity.value,
        tube_cp: tube_cp.value,
        tube_conductivity: tube_conductivity.value,
        tube_mass_flow: tube_mass_flow.value,
        shell_density: shell_density.value,
        shell_viscosity: shell_viscosity.value,
        shell_cp: shell_cp.value,
        shell_conductivity: shell_conductivity.value,
        shell_mass_flow: shell_mass_flow.value,
        fouling_tube,
        fouling_shell,
        shell_side_method,
        tube_fouling_layer_thickness: tube_fouling_layer_thickness.map(|value| value.value),
        tube_fouling_layer_conductivity: tube_fouling_layer_conductivity.map(|value| value.value),
        shell_viscosity_wall: shell_viscosity_wall.map(|value| value.value),
        has_sealing,
        sealing_pairs,
        tube_to_baffle_clearance: tube_to_baffle_clearance.map(|value| value.value),
        shell_to_baffle_clearance: shell_to_baffle_clearance.map(|value| value.value),
        bypass_area: bypass_area.map(|value| value.value),
    })
}

/// The class's `validateStrictInputs`, in its own order.
pub(crate) fn validate_strict(request: &RatingRequest) -> Result<()> {
    let positives = [
        request.tube_outside_diameter,
        request.tube_inside_diameter,
        request.tube_length,
        request.tube_pitch,
        request.shell_inside_diameter,
        request.baffle_spacing,
        request.tube_wall_conductivity,
        request.tube_density,
        request.tube_viscosity,
        request.tube_cp,
        request.tube_conductivity,
        request.tube_mass_flow,
        request.shell_density,
        request.shell_viscosity,
        request.shell_cp,
        request.shell_conductivity,
        request.shell_mass_flow,
    ];
    if positives
        .iter()
        .any(|value| !value.is_finite() || *value <= 0.0)
    {
        return Err(AzothError::invalid_input(
            "geometry",
            "strict thermal rating requires positive finite geometry and fluid inputs",
        ));
    }
    if request.tube_outside_diameter <= request.tube_inside_diameter
        || request.tube_pitch <= request.tube_outside_diameter
        || request.shell_inside_diameter <= request.tube_outside_diameter
        || request.baffle_spacing > request.tube_length
    {
        return Err(AzothError::invalid_input(
            "geometry",
            "strict thermal rating requires OD > ID, pitch > OD, shell ID > OD, and baffle \
             spacing no greater than tube length",
        ));
    }
    let tube_count = integral_count(request.tube_count, "tube_count")?;
    let tube_passes = integral_count(request.tube_passes, "tube_passes")?;
    let baffle_count = integral_count(request.baffle_count, "baffle_count")?;
    if tube_count % tube_passes != 0 {
        return Err(AzothError::invalid_input(
            "tube_count",
            "strict thermal rating requires the tube count to divide evenly by the passes",
        ));
    }
    if baffle_count < 1 {
        return Err(AzothError::invalid_input(
            "baffle_count",
            "strict thermal rating requires at least one baffle",
        ));
    }
    if !request.baffle_cut.is_finite() || request.baffle_cut <= 0.0 || request.baffle_cut >= 0.5 {
        return Err(AzothError::invalid_input(
            "baffle_cut",
            "strict thermal rating requires a baffle cut between zero and one half",
        ));
    }
    if !request.fouling_tube.is_finite()
        || request.fouling_tube < 0.0
        || !request.fouling_shell.is_finite()
        || request.fouling_shell < 0.0
    {
        return Err(AzothError::invalid_input(
            "fouling",
            "strict thermal rating requires finite nonnegative fouling resistances",
        ));
    }
    if let Some(value) = request.shell_viscosity_wall {
        if !value.is_finite() || value == 0.0 {
            return Err(AzothError::invalid_input(
                "shell_viscosity_wall",
                "strict thermal rating requires a finite nonzero wall viscosity; negative \
                 selects the bulk viscosity",
            ));
        }
    }
    let thickness = request.tube_fouling_layer_thickness.unwrap_or(0.0);
    let conductivity = request.tube_fouling_layer_conductivity.unwrap_or(1.0);
    if !thickness.is_finite() || thickness < 0.0 {
        return Err(AzothError::invalid_input(
            "tube_fouling_layer_thickness",
            "the deposit thickness must be finite and nonnegative",
        ));
    }
    if !conductivity.is_finite() || conductivity <= 0.0 {
        return Err(AzothError::invalid_input(
            "tube_fouling_layer_conductivity",
            "the deposit conductivity must be finite and positive",
        ));
    }
    if request.tube_inside_diameter - 2.0 * thickness <= 0.0 {
        return Err(AzothError::invalid_input(
            "tube_fouling_layer_thickness",
            "the deposit must leave a positive tube flow diameter",
        ));
    }

    // The derived geometry the class also checks.
    let effective_id = request.tube_inside_diameter - 2.0 * thickness;
    let tube_flow_area = std::f64::consts::PI / 4.0
        * effective_id
        * effective_id
        * (tube_count as f64 / tube_passes as f64);
    let shell_area = bell_delaware::crossflow_area(
        request.shell_inside_diameter,
        request.baffle_spacing,
        request.tube_outside_diameter,
        request.tube_pitch,
    );
    let equivalent_diameter = bell_delaware::shell_equivalent_diameter(
        request.tube_outside_diameter,
        request.tube_pitch,
        request.triangular_pitch,
    );
    let outside_area = std::f64::consts::PI
        * request.tube_outside_diameter
        * request.tube_length
        * tube_count as f64;
    if [
        tube_flow_area,
        shell_area,
        equivalent_diameter,
        outside_area,
    ]
    .iter()
    .any(|value| !value.is_finite() || *value <= 0.0)
    {
        return Err(AzothError::invalid_input(
            "geometry",
            "strict thermal rating requires finite positive flow and heat-transfer geometry",
        ));
    }

    if request.shell_side_method == ShellSideMethod::BellDelaware {
        let tube_to_baffle = request.tube_to_baffle_clearance.unwrap_or(0.0004);
        let shell_to_baffle = request.shell_to_baffle_clearance.unwrap_or(0.003);
        let bypass = request.bypass_area.unwrap_or(0.0);
        let sealing_pairs = optional_integral(request.sealing_pairs, "sealing_pairs", 0)?;
        let has_sealing = request.has_sealing.unwrap_or(false);
        if !tube_to_baffle.is_finite()
            || tube_to_baffle < 0.0
            || tube_to_baffle >= request.tube_pitch - request.tube_outside_diameter
            || !shell_to_baffle.is_finite()
            || shell_to_baffle < 0.0
            || shell_to_baffle >= request.shell_inside_diameter
            || !bypass.is_finite()
            || bypass < 0.0
            || bypass > request.shell_inside_diameter * request.baffle_spacing
            || (has_sealing && sealing_pairs == 0)
        {
            return Err(AzothError::invalid_input(
                "bell_delaware",
                "strict Bell-Delaware rating requires valid clearances, bypass area and \
                 sealing-strip counts",
            ));
        }
    }
    Ok(())
}
