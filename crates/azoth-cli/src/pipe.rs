//! The `pipe` subcommand: pressure drop through a pipe with fittings.
//!
//! # Why this composes three calcs rather than being one
//!
//! Friction in the straight pipe and loss through fittings are computed by
//! different methods - Darcy-Weisbach over the pipe's own length, and the
//! equivalent-length method for the fittings - and adding them is a modelling
//! decision. Presenting the two contributions separately, as this report does,
//! means a user can see when one of them dominates and question it.
//!
//! # Assumptions this makes explicitly, because it is given less information
//! # than the calcs need
//!
//! The CLI takes a diameter and a flow, so it has to derive a velocity, and it
//! has no way to know the pipe's roughness. Both are surfaced in the output
//! rather than buried: the roughness in particular is a default the user should
//! override for anything that matters.

use azoth_core::units::{self, DynamicViscosity, Length, MassDensity, Velocity};
use azoth_core::{AzothError, Result, Warning, WarningCode};
use azoth_hydraulics as hyd;

/// Default absolute roughness, in metres: 0.046 mm, the usual figure quoted for
/// commercial steel.
///
/// A default rather than a value the library claims to know. Real pipe
/// roughness varies by material and by age - drawn tubing is an order of
/// magnitude smoother, cast iron several times rougher - so this is stated in
/// the output every time it is used.
pub const DEFAULT_ROUGHNESS_M: f64 = 4.6e-5;

/// Which explicit friction factor correlation to use.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FrictionMethod {
    /// The implicit Colebrook-White equation, solved iteratively.
    Colebrook,
    /// The explicit Swamee-Jain approximation to it.
    SwameeJain,
}

impl FrictionMethod {
    /// Parse the user's spelling.
    ///
    /// # Errors
    /// An unknown name is an error rather than a fallback, because silently
    /// choosing a correlation the user did not ask for changes the answer.
    pub fn parse(name: &str) -> Result<Self> {
        match name.trim().to_lowercase().as_str() {
            "colebrook" => Ok(Self::Colebrook),
            "swamee-jain" | "swamee_jain" | "swameejain" => Ok(Self::SwameeJain),
            other => Err(AzothError::invalid_input(
                "friction-method",
                format!("unknown method `{other}`; expected colebrook or swamee-jain"),
            )),
        }
    }

    /// Name used in the report.
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::Colebrook => "colebrook",
            Self::SwameeJain => "swamee-jain",
        }
    }
}

/// How the user expressed the flow rate.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FlowUnit {
    /// Cubic metres per hour.
    CubicMetresPerHour,
    /// Litres per second.
    LitresPerSecond,
    /// Kilograms per second.
    KilogramsPerSecond,
}

impl FlowUnit {
    /// Parse the user's spelling.
    ///
    /// # Errors
    /// An unknown name is an error rather than an assumed default: reading
    /// m3/h as L/s is a factor of 3600.
    pub fn parse(name: &str) -> Result<Self> {
        match name.trim().to_lowercase().as_str() {
            "m3/h" | "m3h" | "m³/h" => Ok(Self::CubicMetresPerHour),
            "l/s" | "ls" => Ok(Self::LitresPerSecond),
            "kg/s" | "kgs" => Ok(Self::KilogramsPerSecond),
            other => Err(AzothError::invalid_input(
                "flow-unit",
                format!("unknown flow unit `{other}`; expected m3/h, L/s or kg/s"),
            )),
        }
    }

    /// Name used in the report.
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::CubicMetresPerHour => "m3/h",
            Self::LitresPerSecond => "L/s",
            Self::KilogramsPerSecond => "kg/s",
        }
    }

    /// Convert a flow in this unit to cubic metres per second.
    fn to_cubic_metres_per_second(self, flow: f64, density: f64) -> f64 {
        match self {
            Self::CubicMetresPerHour => flow / 3600.0,
            Self::LitresPerSecond => flow / 1000.0,
            Self::KilogramsPerSecond => flow / density,
        }
    }
}

/// Everything the report needs, computed.
#[derive(Debug, Clone)]
pub struct PipeResult {
    /// Fluid name.
    pub fluid: String,
    /// Temperature the properties were taken at, in degrees Celsius.
    pub temperature_c: f64,
    /// Fluid density.
    pub density: MassDensity,
    /// Fluid dynamic viscosity.
    pub viscosity: DynamicViscosity,
    /// Pipe internal diameter.
    pub diameter: Length,
    /// Pipe length.
    pub length: Length,
    /// Absolute roughness used.
    pub roughness_m: f64,
    /// Relative roughness, `epsilon / D`.
    pub relative_roughness: f64,
    /// Volumetric flow rate, in m^3/s.
    pub flow_m3_s: f64,
    /// Bulk mean velocity.
    pub velocity: Velocity,
    /// Reynolds number.
    pub re: f64,
    /// Flow regime description.
    pub regime: azoth_core::FlowRegime,
    /// Friction factor.
    pub friction_factor: f64,
    /// Which correlation produced it.
    pub method: FrictionMethod,
    /// Iterations the solver took, when the method is implicit.
    pub iterations: Option<u32>,
    /// Fitting ids used.
    pub fittings: Vec<String>,
    /// Total resistance coefficient for those fittings.
    pub k_total: f64,
    /// Straight-pipe pressure drop.
    pub dp_straight: f64,
    /// Fitting pressure drop.
    pub dp_fittings: f64,
    /// Total pressure drop.
    pub dp_total: f64,
    /// Every caveat the contributing calcs raised, deduplicated.
    pub warnings: Vec<Warning>,
}

/// Compute the pressure drop for a pipe with fittings.
///
/// # Errors
/// Propagates whatever the underlying calcs reject: a non-positive diameter or
/// density, an unknown fitting id, an out-of-table temperature, and so on.
#[allow(clippy::too_many_arguments)] // one per physical input; grouping them would obscure the call sites
pub fn compute(
    fluid: &str,
    temperature_c: f64,
    flow: f64,
    flow_unit: FlowUnit,
    diameter_m: f64,
    length_m: f64,
    roughness_m: f64,
    fitting_ids: &[String],
    method: FrictionMethod,
) -> Result<PipeResult> {
    let table = hyd::provider_for(fluid)?;
    let density = table.density_at_celsius(temperature_c)?;
    let viscosity = table.dynamic_viscosity_at_celsius(temperature_c)?;

    let diameter = units::meters(diameter_m);
    let length = units::meters(length_m);

    // Velocity from the flow and the bore. A = pi D^2 / 4.
    let area = std::f64::consts::PI * diameter_m * diameter_m / 4.0;
    if area <= 0.0 {
        return Err(AzothError::out_of_range(
            "diameter",
            diameter_m,
            "must be positive: a zero bore has no flow area",
        ));
    }
    let flow_m3_s = flow_unit.to_cubic_metres_per_second(flow, density.value);
    let velocity = units::meters_per_second(flow_m3_s / area);

    let mut warnings = Vec::new();

    // Reynolds number. This calc raises on a non-positive density, diameter or
    // viscosity, which is what we want: those make the velocity meaningless too.
    let reynolds = hyd::reynolds_number(density, velocity, diameter, viscosity)?;
    warnings.extend(reynolds.warnings.iter().cloned());

    let relative_roughness = roughness_m / diameter_m;
    if relative_roughness > 0.05 {
        warnings.push(Warning::for_field(
            WarningCode::OutOfValidRange,
            "roughness",
            format!(
                "relative roughness {relative_roughness:.4} exceeds 0.05, above which \
                 the Moody/Colebrook framework is not meaningful. Check that the \
                 roughness ({roughness_m} m) and the diameter ({diameter_m} m) are both \
                 as intended."
            ),
        ));
    }

    // Friction factor.
    let (friction_factor, iterations) = match method {
        FrictionMethod::Colebrook => {
            let solved = hyd::friction_factor_colebrook(reynolds.re, relative_roughness)?;
            warnings.extend(solved.warnings.iter().cloned());
            (solved.f, Some(solved.iterations))
        }
        FrictionMethod::SwameeJain => {
            let solved = hyd::friction_factor_swamee_jain(reynolds.re, relative_roughness)?;
            warnings.extend(solved.warnings.iter().cloned());
            (solved.f, None)
        }
    };

    // Fittings. The basis is the *actual* friction factor at this Reynolds
    // number rather than Crane's f_T for fully turbulent flow - a common and
    // slightly more accurate variant, and the difference is stated here because
    // the two give different answers.
    let (k_total, dp_fittings) = if fitting_ids.is_empty() {
        (0.0, 0.0)
    } else {
        let refs: Vec<&str> = fitting_ids.iter().map(String::as_str).collect();
        let k = hyd::crane_k_factors(&refs, friction_factor)?;
        warnings.extend(k.warnings.iter().cloned());
        let velocity_head = density.value * velocity.value * velocity.value / 2.0;
        (k.k_total, k.k_total * velocity_head)
    };

    // Straight pipe. Viscosity is supplied even though darcy_weisbach does not
    // need it for the arithmetic: without it the calc cannot check the flow
    // regime, and would correctly report RANGE_CHECK_SKIPPED - which would be
    // misleading here, because this tool has just computed the Reynolds number
    // itself and knows the regime perfectly well.
    let straight = hyd::darcy_weisbach(
        friction_factor,
        length,
        diameter,
        density,
        velocity,
        Some(viscosity),
    )?;
    warnings.extend(straight.warnings.iter().cloned());
    let dp_straight = straight.dp.value;

    // Deduplicate. Each contributing calc checks its own inputs, so the same
    // condition can legitimately be reported more than once, and a warning list
    // that repeats itself teaches people to skim it.
    let mut seen = std::collections::HashSet::new();
    warnings.retain(|w| seen.insert((w.code, w.field.clone(), w.message.clone())));

    Ok(PipeResult {
        fluid: table.name.clone(),
        temperature_c,
        density,
        viscosity,
        diameter,
        length,
        roughness_m,
        relative_roughness,
        flow_m3_s,
        velocity,
        re: reynolds.re,
        regime: reynolds.regime,
        friction_factor,
        method,
        iterations,
        fittings: fitting_ids.to_vec(),
        k_total,
        dp_straight,
        dp_fittings,
        dp_total: dp_straight + dp_fittings,
        warnings,
    })
}
