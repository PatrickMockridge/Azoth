//! `thermal.fouling_model` - fixed, Kern-Seaton and Ebert-Panchal fouling.
//!
//! Transcribed from `process/mechanicaldesign/heatexchanger/FoulingModel.java` at NeqSim
//! `118d9aff`: the fixed TEMA resistance, the Kern-Seaton asymptote, and the Ebert-Panchal
//! threshold model with its analytical resistance, rate, threshold temperature and
//! threshold velocity.
//!
//! Spec: `specs/calcs/thermal/fouling_model.toml`.

use azoth_core::units::{
    DynamicViscosity, Length, MassDensity, MolarEnergy, ThermodynamicTemperature, Time, Velocity,
    hours, kelvins, meters_per_second, pascals,
};
use azoth_core::{AzothError, Result, apply_checks};

use crate::results::FoulingModelResult;
use crate::shell_and_tube_rating::darcy_friction;
use crate::spec_gen;

/// The class's universal gas constant, J/(mol*K).
const R_GAS: f64 = 8.314;
/// Seconds in an hour, the class's rate conversion.
const SECONDS_PER_HOUR: f64 = 3600.0;

/// The fouling model type.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FoulingModelType {
    /// A fixed TEMA resistance.
    Fixed,
    /// The Ebert-Panchal threshold model.
    EbertPanchal,
    /// The Kern-Seaton asymptote.
    KernSeaton,
}

impl std::str::FromStr for FoulingModelType {
    type Err = AzothError;

    fn from_str(value: &str) -> Result<Self> {
        match value {
            "fixed" => Ok(Self::Fixed),
            "ebert_panchal" => Ok(Self::EbertPanchal),
            "kern_seaton" => Ok(Self::KernSeaton),
            other => Err(AzothError::invalid_input(
                "model",
                format!("{other} is not one of fixed, ebert_panchal or kern_seaton"),
            )),
        }
    }
}

/// The composed request, in SI magnitudes.
#[derive(Debug, Clone)]
pub struct FoulingRequest {
    pub model: FoulingModelType,
    pub velocity: f64,
    pub fluid_density: f64,
    pub fluid_viscosity: f64,
    pub wall_temperature: f64,
    pub tube_inner_diameter: f64,
    pub alpha: f64,
    pub beta: f64,
    pub activation_energy: f64,
    pub gamma: f64,
    pub fixed_fouling_resistance: f64,
    pub rf_max: f64,
    pub time_constant: f64,
    pub time_hours: f64,
    pub target_fouling_resistance: f64,
    pub threshold_rate: f64,
    pub euler_steps: Option<f64>,
    pub euler_step_hours: Option<f64>,
}

/// The two derived condition values: Reynolds number and wall shear stress.
#[must_use]
pub fn conditions(velocity: f64, density: f64, viscosity: f64, diameter: f64) -> (f64, f64) {
    if viscosity <= 0.0 {
        return (0.0, 0.0);
    }
    let reynolds = density * velocity * diameter / viscosity;
    let friction = darcy_friction(reynolds);
    let shear = friction / 8.0 * density * velocity * velocity;
    (reynolds, shear)
}

/// The Ebert-Panchal deposition term, `alpha Re**beta exp(-E/(R Tw))`.
#[must_use]
pub fn deposition(
    alpha: f64,
    beta: f64,
    activation_energy: f64,
    reynolds: f64,
    wall_k: f64,
) -> f64 {
    if reynolds <= 0.0 || wall_k <= 0.0 {
        return 0.0;
    }
    alpha * reynolds.powf(beta) * (-activation_energy / (R_GAS * wall_k)).exp()
}

/// The instantaneous Ebert-Panchal rate at a resistance, per hour.
#[must_use]
#[allow(clippy::too_many_arguments)]
pub fn ebert_panchal_fouling_rate(
    alpha: f64,
    beta: f64,
    activation_energy: f64,
    gamma: f64,
    wall_shear_stress: f64,
    reynolds: f64,
    wall_k: f64,
    current_rf: f64,
) -> f64 {
    if wall_k <= 0.0 || reynolds <= 0.0 {
        return 0.0;
    }
    let deposition = deposition(alpha, beta, activation_energy, reynolds, wall_k);
    let removal = gamma * wall_shear_stress * current_rf;
    (deposition - removal) * SECONDS_PER_HOUR
}

/// The Ebert-Panchal resistance at a time, the analytical approach to the asymptote.
#[must_use]
#[allow(clippy::too_many_arguments)]
pub fn ebert_panchal_resistance(
    alpha: f64,
    beta: f64,
    activation_energy: f64,
    gamma: f64,
    wall_shear_stress: f64,
    reynolds: f64,
    wall_k: f64,
    time_hours: f64,
) -> f64 {
    if time_hours <= 0.0 {
        return 0.0;
    }
    let deposition = deposition(alpha, beta, activation_energy, reynolds, wall_k);
    let removal_coefficient = gamma * wall_shear_stress;
    if removal_coefficient <= 0.0 {
        return deposition * time_hours * SECONDS_PER_HOUR;
    }
    let asymptotic = deposition / removal_coefficient;
    let tau_hours = 1.0 / removal_coefficient / SECONDS_PER_HOUR;
    asymptotic * (1.0 - (-time_hours / tau_hours).exp())
}

/// The Kern-Seaton resistance at a time.
#[must_use]
pub fn kern_seaton_resistance(rf_max: f64, time_constant: f64, time_hours: f64) -> f64 {
    if time_hours <= 0.0 || rf_max <= 0.0 || time_constant <= 0.0 {
        return 0.0;
    }
    rf_max * (1.0 - (-time_hours / time_constant).exp())
}

/// The threshold wall temperature at a zero-resistance rate.
#[must_use]
pub fn threshold_temperature(
    alpha: f64,
    beta: f64,
    activation_energy: f64,
    reynolds: f64,
    threshold_rate: f64,
) -> f64 {
    if alpha <= 0.0 || threshold_rate <= 0.0 || reynolds <= 0.0 {
        return 0.0;
    }
    let ln_term = (threshold_rate / (alpha * reynolds.powf(beta))).ln();
    if ln_term >= 0.0 {
        // numerics-ok: the class's "any deposition exceeds the threshold" sentinel is
        // Double.MAX_VALUE, carried by both kernels and pinned in the capture.
        return f64::MAX;
    }
    -activation_energy / (R_GAS * ln_term)
}

/// The velocity at which the asymptote crosses the fixed resistance.
#[must_use]
#[allow(clippy::too_many_arguments)]
pub fn threshold_velocity(
    alpha: f64,
    beta: f64,
    activation_energy: f64,
    gamma: f64,
    density: f64,
    viscosity: f64,
    wall_k: f64,
    diameter: f64,
    fixed_resistance: f64,
) -> f64 {
    if wall_k <= 0.0 || density <= 0.0 || viscosity <= 0.0 || diameter <= 0.0 {
        return 0.0;
    }
    let mut low = 0.1;
    let mut high = 10.0;
    for _ in 0..50 {
        let mid = (low + high) / 2.0;
        let reynolds = density * mid * diameter / viscosity;
        let friction = darcy_friction(reynolds);
        let shear = friction / 8.0 * density * mid * mid;
        let deposition = deposition(alpha, beta, activation_energy, reynolds, wall_k);
        let removal = gamma * shear;
        let asymptotic = if removal > 0.0 {
            deposition / removal
        } else {
            // numerics-ok: the class's "no removal" sentinel is Double.MAX_VALUE; it is
            // above the fixed resistance, so the bisection moves its lower bound.
            f64::MAX
        };
        if asymptotic > fixed_resistance {
            low = mid;
        } else {
            high = mid;
        }
    }
    (low + high) / 2.0
}

/// The composed fouling calculation.
///
/// # Errors
/// [`AzothError::InvalidInput`] for a non-integral or negative Euler step count, and
/// [`AzothError::OutOfRange`] for the spec's own bounds.
pub fn rate(request: &FoulingRequest) -> Result<FoulingModelResult> {
    let spec = &spec_gen::FOULING_MODEL_SPEC;
    let mut warnings = Vec::new();
    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "velocity" => Some(request.velocity),
            "fluid_density" => Some(request.fluid_density),
            "fluid_viscosity" => Some(request.fluid_viscosity),
            "wall_temperature" => Some(request.wall_temperature),
            "tube_inner_diameter" => Some(request.tube_inner_diameter),
            "alpha" => Some(request.alpha),
            "activation_energy" => Some(request.activation_energy),
            "gamma" => Some(request.gamma),
            "fixed_fouling_resistance" => Some(request.fixed_fouling_resistance),
            "rf_max" => Some(request.rf_max),
            "time_constant" => Some(request.time_constant),
            "time_hours" => Some(request.time_hours),
            "target_fouling_resistance" => Some(request.target_fouling_resistance),
            "threshold_rate" => Some(request.threshold_rate),
            _ => None,
        },
        &mut warnings,
    )?;
    let steps = match request.euler_steps {
        Some(value) => {
            if !value.is_finite() || (value - value.round()).abs() > 0.0 || value < 0.0 {
                return Err(AzothError::invalid_input(
                    "euler_steps",
                    format!("{value} is not a whole step count"),
                ));
            }
            value.round() as i32
        }
        None => 0,
    };
    // The request carries SI seconds; the class's arithmetic is in hours.
    let time_constant = request.time_constant / SECONDS_PER_HOUR;
    let time_hours = request.time_hours / SECONDS_PER_HOUR;
    let step_hours = request.euler_step_hours.unwrap_or(0.0) / SECONDS_PER_HOUR;

    let (reynolds, shear) = conditions(
        request.velocity,
        request.fluid_density,
        request.fluid_viscosity,
        request.tube_inner_diameter,
    );
    let asymptotic = match request.model {
        FoulingModelType::Fixed => request.fixed_fouling_resistance,
        FoulingModelType::EbertPanchal => {
            let removal = request.gamma * shear;
            if removal > 0.0 {
                deposition(
                    request.alpha,
                    request.beta,
                    request.activation_energy,
                    reynolds,
                    request.wall_temperature,
                ) / removal
            } else {
                // numerics-ok: the class's "no removal" sentinel is Double.MAX_VALUE.
                f64::MAX
            }
        }
        FoulingModelType::KernSeaton => request.rf_max,
    };
    let resistance_at_time = match request.model {
        FoulingModelType::Fixed => request.fixed_fouling_resistance,
        FoulingModelType::EbertPanchal => ebert_panchal_resistance(
            request.alpha,
            request.beta,
            request.activation_energy,
            request.gamma,
            shear,
            reynolds,
            request.wall_temperature,
            time_hours,
        ),
        FoulingModelType::KernSeaton => {
            kern_seaton_resistance(request.rf_max, time_constant, time_hours)
        }
    };
    let mut current = 0.0;
    let mut operating = 0.0;
    for _ in 0..steps {
        if request.model == FoulingModelType::EbertPanchal {
            let fouling_rate = ebert_panchal_fouling_rate(
                request.alpha,
                request.beta,
                request.activation_energy,
                request.gamma,
                shear,
                reynolds,
                request.wall_temperature,
                current,
            );
            current = (current + fouling_rate * step_hours).max(0.0);
        }
        operating += step_hours;
    }
    let current_resistance = match request.model {
        FoulingModelType::Fixed => request.fixed_fouling_resistance,
        FoulingModelType::EbertPanchal => current,
        FoulingModelType::KernSeaton => {
            kern_seaton_resistance(request.rf_max, time_constant, operating)
        }
    };
    let time_to_target = if request.target_fouling_resistance >= asymptotic {
        -1.0
    } else {
        match request.model {
            FoulingModelType::Fixed => 0.0,
            FoulingModelType::EbertPanchal => {
                let removal = request.gamma * shear;
                let tau_hours = if removal > 0.0 {
                    1.0 / (removal * SECONDS_PER_HOUR)
                } else {
                    // numerics-ok: the class's "no removal" time constant is Double.MAX_VALUE;
                    // the product with the zero logarithm answers the class's -0.0.
                    f64::MAX
                };
                -tau_hours * (1.0 - request.target_fouling_resistance / asymptotic).ln()
            }
            FoulingModelType::KernSeaton => {
                -time_constant * (1.0 - request.target_fouling_resistance / request.rf_max).ln()
            }
        }
    };
    Ok(FoulingModelResult {
        fouling_resistance_at_time: resistance_at_time,
        current_fouling_resistance: current_resistance,
        asymptotic_fouling_resistance: asymptotic,
        threshold_wall_temperature: kelvins(threshold_temperature(
            request.alpha,
            request.beta,
            request.activation_energy,
            reynolds,
            request.threshold_rate,
        )),
        threshold_velocity: meters_per_second(threshold_velocity(
            request.alpha,
            request.beta,
            request.activation_energy,
            request.gamma,
            request.fluid_density,
            request.fluid_viscosity,
            request.wall_temperature,
            request.tube_inner_diameter,
            request.fixed_fouling_resistance,
        )),
        time_to_target: hours(time_to_target),
        operating_time: hours(operating),
        fouling_rate: ebert_panchal_fouling_rate(
            request.alpha,
            request.beta,
            request.activation_energy,
            request.gamma,
            shear,
            reynolds,
            request.wall_temperature,
            current_resistance,
        ),
        reynolds_number: reynolds,
        wall_shear_stress: pascals(shear),
        warnings,
    })
}

/// The flat entry point the binding and the case runner call.
///
/// # Errors
/// As [`rate`].
#[allow(clippy::too_many_arguments)]
pub fn fouling_model(
    model: &str,
    velocity: Velocity,
    fluid_density: MassDensity,
    fluid_viscosity: DynamicViscosity,
    wall_temperature: ThermodynamicTemperature,
    tube_inner_diameter: Length,
    alpha: f64,
    beta: f64,
    activation_energy: MolarEnergy,
    gamma: f64,
    fixed_fouling_resistance: f64,
    rf_max: f64,
    time_constant: Time,
    time_hours: Time,
    target_fouling_resistance: f64,
    threshold_rate: f64,
    euler_steps: Option<f64>,
    euler_step_hours: Option<Time>,
) -> Result<FoulingModelResult> {
    rate(&FoulingRequest {
        model: model.parse()?,
        velocity: velocity.value,
        fluid_density: fluid_density.value,
        fluid_viscosity: fluid_viscosity.value,
        wall_temperature: wall_temperature.value,
        tube_inner_diameter: tube_inner_diameter.value,
        alpha,
        beta,
        activation_energy: activation_energy.value,
        gamma,
        fixed_fouling_resistance,
        rf_max,
        time_constant: time_constant.value,
        time_hours: time_hours.value,
        target_fouling_resistance,
        threshold_rate,
        euler_steps,
        euler_step_hours: euler_step_hours.map(|value| value.value),
    })
}
