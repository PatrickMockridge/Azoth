//! `thermal.lmtd_correction_factor` - the Bowman-Mueller-Nagle correction.
//!
//! Transcribed from `process/mechanicaldesign/heatexchanger/LMTDcorrectionFactor.java` at
//! NeqSim `118d9aff`. The four temperatures enter as differences; `R`, `P` and `F_t` are
//! the class's own intermediates.
//!
//! Spec: `specs/calcs/thermal/lmtd_correction_factor.toml`.

use azoth_core::units::ThermodynamicTemperature;
use azoth_core::{AzothError, Result, apply_checks};

use crate::results::LmtdCorrectionFactorResult;
use crate::spec_gen;

/// The class's minimum acceptable correction factor.
pub const MIN_ACCEPTABLE_FT: f64 = 0.75;

/// The correction-factor request, in K.
#[derive(Debug, Clone)]
pub struct LmtdRequest {
    pub hot_inlet_temperature: f64,
    pub hot_outlet_temperature: f64,
    pub cold_inlet_temperature: f64,
    pub cold_outlet_temperature: f64,
    pub shell_passes: Option<f64>,
}

/// `R = (Th_in - Th_out)/(Tc_out - Tc_in)`, with a zero cold change answering 1.
#[must_use]
pub fn calc_r(hot_in: f64, hot_out: f64, cold_in: f64, cold_out: f64) -> f64 {
    let cold_change = cold_out - cold_in;
    if cold_change.abs() < 1e-10 {
        return 1.0;
    }
    (hot_in - hot_out) / cold_change
}

/// `P = (Tc_out - Tc_in)/(Th_in - Tc_in)`, with a zero inlet difference answering 0.
#[must_use]
pub fn calc_p(hot_in: f64, _hot_out: f64, cold_in: f64, cold_out: f64) -> f64 {
    let max_change = hot_in - cold_in;
    if max_change.abs() < 1e-10 {
        return 0.0;
    }
    (cold_out - cold_in) / max_change
}

/// The single-shell (1-2N) correction factor.
#[must_use]
pub fn calc_ft_single_shell(r: f64, p: f64) -> f64 {
    if p <= 0.0 || p >= 1.0 {
        return 1.0;
    }
    let w = (r * r + 1.0).sqrt();
    if (r - 1.0).abs() < 1e-6 {
        let sqrt_two = std::f64::consts::SQRT_2;
        let numerator = p * sqrt_two;
        let arg1 = 2.0 - p * (2.0 - sqrt_two);
        let arg2 = 2.0 - p * (2.0 + sqrt_two);
        if arg1 <= 0.0 || arg2 <= 0.0 || arg1 / arg2 <= 0.0 {
            return 0.0;
        }
        let denominator = (1.0 - p) * (arg1 / arg2).ln();
        if denominator.abs() < 1e-15 {
            return 1.0;
        }
        return numerator / denominator;
    }
    let rp = r * p;
    if rp >= 1.0 {
        return 0.0;
    }
    let log_arg1 = (1.0 - p) / (1.0 - rp);
    if log_arg1 <= 0.0 {
        return 0.0;
    }
    let arg2_num = 2.0 - p * (r + 1.0 - w);
    let arg2_den = 2.0 - p * (r + 1.0 + w);
    if arg2_num <= 0.0 || arg2_den <= 0.0 || arg2_num / arg2_den <= 0.0 {
        return 0.0;
    }
    let numerator = w * log_arg1.ln();
    let denominator = (r - 1.0) * (arg2_num / arg2_den).ln();
    if denominator.abs() < 1e-15 {
        return 1.0;
    }
    (numerator / denominator).clamp(0.0, 1.0)
}

/// The overall `P` converted to a per-shell `P1` for `n` shells in series.
#[must_use]
pub fn convert_p_to_per_shell(p: f64, r: f64, n: i32) -> f64 {
    if n <= 1 {
        return p;
    }
    if (r - 1.0).abs() < 1e-6 {
        return p / (f64::from(n) - f64::from(n - 1) * p);
    }
    let x = (1.0 - r * p) / (1.0 - p);
    if x <= 0.0 {
        return f64::NAN;
    }
    let x_to_one_over_n = x.powf(1.0 / f64::from(n));
    let denominator = x_to_one_over_n - r;
    if denominator.abs() < 1e-15 {
        return f64::NAN;
    }
    (x_to_one_over_n - 1.0) / denominator
}

/// The correction factor for `shell_passes` shell passes.
#[must_use]
pub fn calc_ft_from_rp(r: f64, p: f64, shell_passes: i32) -> f64 {
    if p <= 0.0 || p >= 1.0 {
        return 1.0;
    }
    if r <= 0.0 {
        return 1.0;
    }
    let passes = shell_passes.max(1);
    let mut per_shell = p;
    if passes > 1 {
        per_shell = convert_p_to_per_shell(p, r, passes);
        if per_shell <= 0.0 || per_shell >= 1.0 || per_shell.is_nan() {
            return 1.0;
        }
    }
    calc_ft_single_shell(r, per_shell)
}

/// The first pass count from one to six reaching 0.75, or -1.
#[must_use]
pub fn required_shell_passes(hot_in: f64, hot_out: f64, cold_in: f64, cold_out: f64) -> f64 {
    let r = calc_r(hot_in, hot_out, cold_in, cold_out);
    let p = calc_p(hot_in, hot_out, cold_in, cold_out);
    for passes in 1..=6 {
        if calc_ft_from_rp(r, p, passes) >= MIN_ACCEPTABLE_FT {
            return f64::from(passes);
        }
    }
    -1.0
}

/// The correction-factor calculation.
///
/// # Errors
/// [`AzothError::InvalidInput`] for a non-integral shell-pass count, and
/// [`AzothError::OutOfRange`] for the spec's temperature bounds.
pub fn rate(request: &LmtdRequest) -> Result<LmtdCorrectionFactorResult> {
    let spec = &spec_gen::LMTD_CORRECTION_FACTOR_SPEC;
    let mut warnings = Vec::new();
    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "hot_inlet_temperature" => Some(request.hot_inlet_temperature),
            "hot_outlet_temperature" => Some(request.hot_outlet_temperature),
            "cold_inlet_temperature" => Some(request.cold_inlet_temperature),
            "cold_outlet_temperature" => Some(request.cold_outlet_temperature),
            _ => None,
        },
        &mut warnings,
    )?;
    let passes = match request.shell_passes {
        Some(value) => {
            if !value.is_finite() || (value - value.round()).abs() > 0.0 {
                return Err(AzothError::invalid_input(
                    "shell_passes",
                    format!("{value} is not a whole shell-pass count"),
                ));
            }
            (value.round() as i32).max(1)
        }
        None => 1,
    };
    let r = calc_r(
        request.hot_inlet_temperature,
        request.hot_outlet_temperature,
        request.cold_inlet_temperature,
        request.cold_outlet_temperature,
    );
    let p = calc_p(
        request.hot_inlet_temperature,
        request.hot_outlet_temperature,
        request.cold_inlet_temperature,
        request.cold_outlet_temperature,
    );
    Ok(LmtdCorrectionFactorResult {
        r,
        p,
        ft: calc_ft_from_rp(r, p, passes),
        ft_1_shell_pass: calc_ft_from_rp(r, p, 1),
        ft_2_shell_passes: calc_ft_from_rp(r, p, 2),
        required_shell_passes: required_shell_passes(
            request.hot_inlet_temperature,
            request.hot_outlet_temperature,
            request.cold_inlet_temperature,
            request.cold_outlet_temperature,
        ),
        warnings,
    })
}

/// The flat entry point the binding and the case runner call.
///
/// # Errors
/// As [`rate`].
pub fn lmtd_correction_factor(
    hot_inlet_temperature: ThermodynamicTemperature,
    hot_outlet_temperature: ThermodynamicTemperature,
    cold_inlet_temperature: ThermodynamicTemperature,
    cold_outlet_temperature: ThermodynamicTemperature,
    shell_passes: Option<f64>,
) -> Result<LmtdCorrectionFactorResult> {
    rate(&LmtdRequest {
        hot_inlet_temperature: hot_inlet_temperature.value,
        hot_outlet_temperature: hot_outlet_temperature.value,
        cold_inlet_temperature: cold_inlet_temperature.value,
        cold_outlet_temperature: cold_outlet_temperature.value,
        shell_passes,
    })
}
