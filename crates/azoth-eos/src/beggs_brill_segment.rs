//! `eos.beggs_brill_segment` - NeqSim's `PipeBeggsAndBrills` segment correlation.
//!
//! One pipe segment's Beggs-and-Brill flow regime, liquid holdup with its inclination
//! correction, no-slip density and viscosity, Reynolds number, friction factor and pressure
//! drop - without the flash or the march, which are the model layer's.

use azoth_core::units::{
    DynamicViscosity, Length, MassDensity, Pressure, ThermodynamicTemperature, VolumeRate,
    meters_per_second, pascals,
};
use azoth_core::{AzothError, Result, apply_checks};

use crate::results::BeggsBrillSegmentResult;
use crate::spec_gen;

/// The class's field-unit constants: feet per metre, cubic feet per cubic metre, pounds per
/// cubic foot per kilogram per cubic metre, and bar per pound per (foot second squared).
const FEET_PER_METER: f64 = 3.2808399;
const FT3_PER_M3: f64 = FEET_PER_METER * FEET_PER_METER * FEET_PER_METER;
/// NeqSim's exact kg/m3 per lb/ft3, used to read a phase density the way the class does.
const KG_M3_PER_LB_FT3: f64 = 16.018_463_37;
/// The class's rounded factor when it writes a density back to the profile.
const PROFILE_LB_FT3_TO_KG_M3: f64 = 16.018_46;
const BAR_PER_IMPERIAL_DROP: f64 = 1.487_27e-5;
const CP_PER_PA_S: f64 = 1000.0;
const FROUDE_G: f64 = 32.174;
const HYDROSTATIC_G: f64 = 32.2;
const POUNDS_PER_SQUARE_INCH: f64 = 6_894.757_293_168;
const MINIMUM_SURFACE_TENSION_DYNE_CM: f64 = 1.0;

/// One pipe segment's Beggs-and-Brill state.
///
/// # Errors
/// * [`azoth_core::AzothError::OutOfRange`] if an input is outside the spec's hard range.
/// * [`azoth_core::AzothError::InvalidInput`] if neither phase flows.
///
/// # Example
/// ```
/// use azoth_core::units::{
///     cubic_meters_per_second, kelvins,
///     kilograms_per_cubic_meter, meters, pascal_seconds, pascals,
/// };
/// use azoth_eos::beggs_brill_segment::beggs_brill_segment;
///
/// let r = beggs_brill_segment(
///     cubic_meters_per_second(0.01796694168385081),
///     cubic_meters_per_second(1.42993961618865e-4),
///     kilograms_per_cubic_meter(69.58154569795788),
///     kilograms_per_cubic_meter(970.1200910245294),
///     pascal_seconds(1.5089488997328844e-5),
///     pascal_seconds(4.6535425555414334e-4),
///     kelvins(333.15),
///     pascals(9999960.703327656),
///     meters(0.1),
///     meters(1.0),
///     meters(0.0),
///     meters(1.0e-5),
/// )?;
/// assert_eq!(r.regime, "segregated");
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
#[allow(clippy::too_many_arguments)] // The signature is the spec's declared inputs.
pub fn beggs_brill_segment(
    gas_flow: VolumeRate,
    liquid_flow: VolumeRate,
    gas_density: MassDensity,
    liquid_density: MassDensity,
    gas_viscosity: DynamicViscosity,
    liquid_viscosity: DynamicViscosity,
    temperature: ThermodynamicTemperature,
    pressure: Pressure,
    diameter: Length,
    length: Length,
    elevation: Length,
    roughness: Length,
) -> Result<BeggsBrillSegmentResult> {
    let spec = &spec_gen::BEGGS_BRILL_SEGMENT_SPEC;
    let mut warnings = Vec::new();

    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "gas_flow" => Some(gas_flow.value),
            "liquid_flow" => Some(liquid_flow.value),
            "gas_density" => Some(gas_density.value),
            "liquid_density" => Some(liquid_density.value),
            "gas_viscosity" => Some(gas_viscosity.value),
            "liquid_viscosity" => Some(liquid_viscosity.value),
            "temperature" => Some(temperature.value),
            "pressure" => Some(pressure.value),
            "diameter" => Some(diameter.value),
            "length" => Some(length.value),
            "roughness" => Some(roughness.value),
            _ => None,
        },
        &mut warnings,
    )?;

    if gas_flow.value <= 0.0 && liquid_flow.value <= 0.0 {
        return Err(AzothError::invalid_input(
            "gas_flow",
            "at least one phase must flow; the superficial velocities are undefined otherwise",
        ));
    }

    // The class's own field units, which every correlation coefficient below is written in.
    let d_ft = diameter.value * FEET_PER_METER;
    let l_ft = length.value * FEET_PER_METER;
    let e_ft = elevation.value * FEET_PER_METER;
    let area = std::f64::consts::PI / 4.0 * d_ft * d_ft;
    let v_sl = liquid_flow.value * FT3_PER_M3 / area;
    let v_sg = gas_flow.value * FT3_PER_M3 / area;
    let v_mix = v_sl + v_sg;

    let rho_g = gas_density.value / KG_M3_PER_LB_FT3;
    let rho_l = liquid_density.value / KG_M3_PER_LB_FT3;
    let mu_g = gas_viscosity.value * CP_PER_PA_S;
    let mu_l = liquid_viscosity.value * CP_PER_PA_S;

    let (regime, _lambda, el, rho_mix, rho_ns, mu_ns, s_multiplier) = if liquid_flow.value <= 0.0 {
        ("single_phase", 0.0, 0.0, rho_g, rho_g, mu_g, 0.0)
    } else if gas_flow.value <= 0.0 {
        ("single_phase", 1.0, 1.0, rho_l, rho_l, mu_l, 0.0)
    } else {
        let lambda = v_sl / v_mix;
        let froude = v_mix * v_mix / (FROUDE_G * d_ft);
        let l1 = 316.0 * lambda.powf(0.302);
        let l2 = 0.000_925_2 * lambda.powf(-2.468_4);
        let l3 = 0.1 * lambda.powf(-1.451_6);
        let l4 = 0.5 * lambda.powf(-6.738);
        let regime = if (lambda < 0.01 && froude < l1) || (lambda >= 0.01 && froude < l2) {
            "segregated"
        } else if lambda >= 0.01 && froude >= l2 && froude <= l3 {
            "transition"
        } else if ((0.01..0.4).contains(&lambda) && froude > l3 && froude <= l1)
            || (lambda >= 0.4 && froude > l3 && froude <= l4)
        {
            "intermittent"
        } else if (lambda < 0.4 && froude >= l1) || (lambda >= 0.4 && froude > l4) {
            "distributed"
        } else {
            return Err(AzothError::invalid_input(
                "gas_flow",
                "the Beggs-and-Brill regime map has no regime for this state",
            ));
        };
        let a = (l3 - froude) / (l3 - l2);

        let nvl = liquid_velocity_number(rho_l, v_sl, temperature.value, pressure.value);
        let theta = (elevation.value / length.value).asin();
        let sin18 = (1.8 * theta).sin();

        let holdup = if regime == "transition" {
            a * inclination_corrected_holdup("segregated", lambda, froude, nvl, sin18, e_ft)
                + (1.0 - a)
                    * inclination_corrected_holdup("intermittent", lambda, froude, nvl, sin18, e_ft)
        } else {
            inclination_corrected_holdup(regime, lambda, froude, nvl, sin18, e_ft)
        };
        let el = lambda.max(holdup.min(1.0));
        let rho_mix = rho_l * el + rho_g * (1.0 - el);
        let rho_ns = rho_l * lambda + rho_g * (1.0 - lambda);
        let mu_ns = mu_l * lambda + mu_g * (1.0 - lambda);

        let y = lambda / (el * el);
        let s = if y > 1.0 && y < 1.2 {
            (2.2 * y - 1.2).ln()
        } else {
            let ln_y = y.ln();
            ln_y / (-0.0523 + 3.182 * ln_y - 0.8725 * ln_y * ln_y + 0.01853 * ln_y.powi(4))
        };
        (regime, lambda, el, rho_mix, rho_ns, mu_ns, s)
    };

    let reynolds = rho_ns * v_mix * d_ft * (16.0 / (3.28 * 3.28)) / (0.001 * mu_ns);
    let relative_roughness = roughness.value / diameter.value;
    let friction_factor = if reynolds.abs() < 1.0e-10 {
        0.0
    } else if reynolds.abs() < 2300.0 {
        64.0 / reynolds
    } else if reynolds.abs() < 4000.0 {
        let laminar = 64.0 / 2300.0;
        let turbulent =
            (1.0 / (-1.8 * ((relative_roughness / 3.7).powf(1.11) + 6.9 / 4000.0).log10())).powi(2);
        laminar + (turbulent - laminar) * (reynolds - 2300.0) / 1700.0
    } else {
        (1.0 / (-1.8 * ((relative_roughness / 3.7).powf(1.11) + 6.9 / reynolds).log10())).powi(2)
    };
    let two_phase = friction_factor * s_multiplier.exp();

    let hydrostatic_raw = rho_mix * HYDROSTATIC_G * e_ft;
    let friction_raw = two_phase * v_mix * v_mix * rho_ns * l_ft / (2.0 * d_ft);
    let to_pa = BAR_PER_IMPERIAL_DROP * 1.0e5;

    Ok(BeggsBrillSegmentResult {
        regime: regime.to_owned(),
        liquid_holdup: el,
        superficial_gas_velocity: meters_per_second(v_sg / FEET_PER_METER),
        superficial_liquid_velocity: meters_per_second(v_sl / FEET_PER_METER),
        superficial_mixture_velocity: meters_per_second(v_mix / FEET_PER_METER),
        no_slip_density: azoth_core::units::kilograms_per_cubic_meter(
            rho_ns * PROFILE_LB_FT3_TO_KG_M3,
        ),
        no_slip_viscosity: azoth_core::units::pascal_seconds(mu_ns / CP_PER_PA_S),
        reynolds,
        friction_factor,
        two_phase_friction_factor: two_phase,
        hydrostatic_drop: pascals(hydrostatic_raw * to_pa),
        friction_drop: pascals(friction_raw * to_pa),
        total_drop: pascals((hydrostatic_raw + friction_raw) * to_pa),
        warnings,
    })
}

fn liquid_velocity_number(
    liquid_density_lb_ft3: f64,
    liquid_velocity_ft_s: f64,
    temperature_k: f64,
    pressure_pa: f64,
) -> f64 {
    let specific_gravity = liquid_density_lb_ft3 / (1000.0 * 0.062_427_960_6);
    let api = 141.5 / specific_gravity - 131.5;
    let sigma68 = 39.0 - 0.2571 * api;
    let sigma100 = 37.5 - 0.2571 * api;
    let temperature_f = (temperature_k - 273.15) * 9.0 / 5.0 + 32.0;
    let mut sigma = if temperature_f > 100.0 {
        sigma100
    } else if temperature_f < 68.0 {
        sigma68
    } else {
        sigma68 + (temperature_f - 68.0) * (sigma100 - sigma68) / (100.0 - 68.0)
    };
    let pressure_psi = pressure_pa / POUNDS_PER_SQUARE_INCH;
    sigma *= 1.0 - 0.024 * pressure_psi.powf(0.45);
    sigma = if sigma.is_nan() {
        MINIMUM_SURFACE_TENSION_DYNE_CM
    } else {
        sigma.max(MINIMUM_SURFACE_TENSION_DYNE_CM)
    };
    1.938 * liquid_velocity_ft_s * (liquid_density_lb_ft3 / sigma).powf(0.25)
}

fn inclination_corrected_holdup(
    regime: &str,
    lambda: f64,
    froude: f64,
    nvl: f64,
    sin18: f64,
    elevation_ft: f64,
) -> f64 {
    let horizontal = match regime {
        "segregated" => 0.98 * lambda.powf(0.4846) / froude.powf(0.0868),
        "intermittent" => 0.845 * lambda.powf(0.5351) / froude.powf(0.0173),
        _ => 1.065 * lambda.powf(0.5824) / froude.powf(0.0609),
    };
    let log_arg = if elevation_ft > 0.0 {
        match regime {
            "segregated" => 0.011 * nvl.powf(3.539) / (lambda.powf(3.768) * froude.powf(1.614)),
            "intermittent" => 2.96 * lambda.powf(0.305) * froude.powf(0.0978) / nvl.powf(0.4473),
            // The published distributed-uphill coefficient is zero.
            _ => f64::NAN,
        }
    } else {
        4.70 * nvl.powf(0.1244) / (lambda.powf(0.3692) * froude.powf(0.5056))
    };
    let mut betta = if log_arg > 0.0 {
        (1.0 - lambda) * log_arg.ln()
    } else {
        0.0
    };
    if betta < 0.0 {
        betta = 0.0;
    }
    horizontal * (1.0 + betta * (sin18 - sin18.powi(3) / 3.0))
}
