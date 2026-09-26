//! `hydraulics.packing_hydraulics` - a packed bed's flooding, load, pressure drop and film
//! coefficients.
//!
//! Spec: `specs/calcs/hydraulics/packing_hydraulics.toml`, which records the four correlations,
//! the empirical bands and the clamps.
//!
//! **Three of the four are Onda's and the fourth is Leva's.** The wetted area, `kGa` and `kLa`
//! are Onda (1968); the bed's pressure drop is Leva's, with the liquid-loading correction; the
//! flooding velocity is the Eckert GPDC fit; and the HETP is a two-resistance combination of the
//! film coefficients. The calculator is `PackingHydraulicsCalculator`, which the rate-based
//! column's transport snapshot reads - `getWettedArea`, `getKGa`, `getKLa`,
//! `getPressureDropPerMeter` and `getPercentFlood` are the five numbers it takes.

use azoth_core::units::{
    DynamicViscosity, Length, MassDensity, MassRate, SurfaceTension, meters, pascals,
};
use azoth_core::{AzothError, Result, apply_checks};

use crate::packing::{PackingSpecification, packing_or_default};
use crate::results::PackingHydraulicsResult;
use crate::spec_gen;

/// The gravitational acceleration the class writes inline, in m/s**2.
const G: f64 = 9.81;

/// The reference liquid viscosity the flooding fit's correction uses, in Pa*s.
const WATER_VISCOSITY: f64 = 0.001;

/// The Eckert fit's three coefficients: `log10(Y) = c0 + c1 x + c2 x**2` with `x = log10(FLV)`.
const ECKERT_C0: f64 = -1.668;
const ECKERT_C1: f64 = -1.085;
const ECKERT_C2: f64 = -0.297;

/// The two bounds the flow parameter is clamped to before the fit is evaluated.
const MIN_FLOW_PARAMETER: f64 = 0.005;
const MAX_FLOW_PARAMETER: f64 = 5.0;

/// The Onda constants: the gas-film coefficient's `5.23`, the liquid's `0.0051`, and the
/// wetted-area expression's `-1.45`.
const ONDA_GAS_CONSTANT: f64 = 5.23;
const ONDA_LIQUID_CONSTANT: f64 = 0.0051;
const ONDA_WETTING_CONSTANT: f64 = -1.45;

/// The wetted-area ratio's own bounds, which the class clamps to.
const MIN_WETTED_RATIO: f64 = 0.2;

/// The minimum wetting rate per category, in m**3/(m**2 s).
const MIN_WETTING_RATE_STRUCTURED: f64 = 2.5e-5;
const MIN_WETTING_RATE_RANDOM: f64 = 5.0e-5;

/// The design window the verdict is stated over, in per cent of flood.
const DESIGN_MIN_FLOOD: f64 = 40.0;
const DESIGN_MAX_FLOOD: f64 = 80.0;

/// The floor the two-stripping-factor HETP is clamped up to, in metres.
const MIN_HETP_M: f64 = 0.1;

/// Everything the calculator reads, as one state.
#[derive(Debug, Clone, Copy)]
pub struct PackingState {
    /// The column's internal diameter.
    pub column_diameter: Length,
    /// The packed height, which the theoretical-stage count divides by the HETP.
    pub packed_height: f64,
    /// The vapour's mass flow.
    pub vapor_mass_flow: MassRate,
    /// The liquid's mass flow.
    pub liquid_mass_flow: MassRate,
    /// The vapour's mass density.
    pub vapor_density: MassDensity,
    /// The liquid's mass density.
    pub liquid_density: MassDensity,
    /// The vapour's dynamic viscosity.
    pub vapor_viscosity: DynamicViscosity,
    /// The liquid's dynamic viscosity.
    pub liquid_viscosity: DynamicViscosity,
    /// The interfacial surface tension, in N/m.
    pub surface_tension: SurfaceTension,
    /// The vapour's diffusivity, in m**2/s.
    pub vapor_diffusivity: f64,
    /// The liquid's diffusivity, in m**2/s.
    pub liquid_diffusivity: f64,
    /// The packing's relative hydraulic capacity factor, which the flooding velocity is scaled
    /// by. One is the class's default.
    pub hydraulic_capacity_factor: f64,
}

/// A packed bed's hydraulics, from the packing's geometry and the state it runs at.
///
/// **The packing resolves by name**, with the class's own fallback: a name the library does not
/// carry is `Pall-Ring-50`, and the compiled table's row for it **replaces the built-in** - the
/// file is loaded after the built-ins and both normalize to `pallring50`, so the plastic row's
/// geometry is the one that answers.
///
/// # Errors
/// * [`azoth_core::AzothError::OutOfRange`] if the diameter, either mass flow, either density or
///   either viscosity is not positive, or if the surface tension or a diffusivity is negative.
///
/// # Example
/// ```
/// use azoth_core::units::{
///     kilograms_per_cubic_meter, kilograms_per_second, meters, newtons_per_meter,
///     pascal_seconds, pascals,
/// };
/// use azoth_hydraulics::packing_hydraulics::{PackingState, packing_hydraulics};
///
/// let out = packing_hydraulics("Pall-Ring-50", PackingState {
///     column_diameter: meters(1.0),
///     packed_height: 5.0,
///     vapor_mass_flow: kilograms_per_second(0.35),
///     liquid_mass_flow: kilograms_per_second(3.5),
///     vapor_density: kilograms_per_cubic_meter(45.0),
///     liquid_density: kilograms_per_cubic_meter(990.0),
///     vapor_viscosity: pascal_seconds(1.8e-5),
///     liquid_viscosity: pascal_seconds(6.5e-4),
///     surface_tension: newtons_per_meter(0.072),
///     vapor_diffusivity: 3.3e-7,
///     liquid_diffusivity: 1.9e-9,
///     hydraulic_capacity_factor: 1.0,
/// })?;
/// assert!((out.wetted_area - 47.63208581918636).abs() < 1e-12);
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
pub fn packing_hydraulics(name: &str, state: PackingState) -> Result<PackingHydraulicsResult> {
    let spec = &spec_gen::PACKING_HYDRAULICS_SPEC;
    let mut warnings = Vec::new();

    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "column_diameter" => Some(state.column_diameter.value),
            "packed_height" => Some(state.packed_height),
            "vapor_mass_flow" => Some(state.vapor_mass_flow.value),
            "liquid_mass_flow" => Some(state.liquid_mass_flow.value),
            "vapor_density" => Some(state.vapor_density.value),
            "liquid_density" => Some(state.liquid_density.value),
            "vapor_viscosity" => Some(state.vapor_viscosity.value),
            "liquid_viscosity" => Some(state.liquid_viscosity.value),
            "surface_tension" => Some(state.surface_tension.value),
            "vapor_diffusivity" => Some(state.vapor_diffusivity),
            "liquid_diffusivity" => Some(state.liquid_diffusivity),
            "hydraulic_capacity_factor" => Some(state.hydraulic_capacity_factor),
            _ => None,
        },
        &mut warnings,
    )?;

    let packing = packing_or_default(name);
    let area =
        std::f64::consts::PI / 4.0 * state.column_diameter.value * state.column_diameter.value;
    if area <= 0.0 {
        return Err(AzothError::invalid_input(
            "column_diameter",
            "a column of no cross-section has no velocity to be in",
        ));
    }

    let flooding_velocity = flooding_velocity(&state, packing);
    let vapor_velocity = (state.vapor_mass_flow.value / state.vapor_density.value) / area;
    let f_factor = vapor_velocity * state.vapor_density.value.sqrt();
    let liquid_velocity = (state.liquid_mass_flow.value / state.liquid_density.value) / area;
    let percent_flood = if flooding_velocity > 0.0 {
        (vapor_velocity / flooding_velocity) * 100.0
    } else {
        0.0
    };

    let pressure_drop_per_meter = pressure_drop(&state, packing, vapor_velocity, area);
    let (wetted_area, k_ga, k_la) = mass_transfer(&state, packing, area, liquid_velocity);

    // ---- The HETP's three rungs, in the class's own order.
    let wetting_rate_minimum = if packing.category.eq_ignore_ascii_case("structured") {
        MIN_WETTING_RATE_STRUCTURED
    } else {
        MIN_WETTING_RATE_RANDOM
    };
    let wetting_rate = (state.liquid_mass_flow.value / state.liquid_density.value)
        / (area * packing.specific_surface_area);
    let wetting_ok = wetting_rate >= wetting_rate_minimum;

    let estimated = estimate_hetp(packing, state.column_diameter.value);
    let (htu_g, htu_l, htu_og, hetp) = if k_ga <= 0.0
        || k_la <= 0.0
        || state.vapor_mass_flow.value <= 0.0
        || state.liquid_mass_flow.value <= 0.0
    {
        (0.0, 0.0, 0.0, estimated)
    } else {
        let g_v = state.vapor_mass_flow.value / area;
        let g_l = state.liquid_mass_flow.value / area;
        let htu_g = g_v / (k_ga * state.vapor_density.value);
        let htu_l = g_l / k_la;
        // The stripping factor is one, which is the class's own constant, and the HETP is then
        // the limit of its formula.
        let htu_og = htu_g + htu_l;
        // The clamp is two-sided: up to the estimate's own two-fold band, and never below 0.1 m.
        let mut hetp = htu_og.max(MIN_HETP_M);
        if hetp > estimated * 2.0 {
            hetp = estimated * 2.0;
        }
        (htu_g, htu_l, htu_og, hetp)
    };
    let theoretical_stages = if hetp > 0.0 {
        state.packed_height / hetp
    } else {
        0.0
    };

    let design_ok = wetting_ok && (DESIGN_MIN_FLOOD..=DESIGN_MAX_FLOOD).contains(&percent_flood);

    apply_checks(
        spec.derived_checks(),
        |quantity| match quantity {
            "percent_flood" => Some(percent_flood),
            "wetted_area" => Some(wetted_area),
            _ => None,
        },
        &mut warnings,
    )?;

    Ok(PackingHydraulicsResult {
        packing_name: packing.name.clone(),
        packing_category: packing.category.clone(),
        specific_surface_area: packing.specific_surface_area,
        void_fraction: packing.void_fraction,
        packing_factor: packing.packing_factor,
        flooding_velocity,
        vapor_velocity,
        liquid_velocity,
        f_factor,
        percent_flood,
        pressure_drop_per_meter: pascals(pressure_drop_per_meter),
        total_pressure_drop: pascals(pressure_drop_per_meter * state.packed_height),
        wetted_area,
        k_ga,
        k_la,
        htu_g,
        htu_l,
        htu_og,
        hetp: meters(hetp),
        theoretical_stages,
        wetting_rate,
        minimum_wetting_rate: wetting_rate_minimum,
        wetting_ok,
        design_ok,
        warnings,
    })
}

/// `calculateFloodingVelocity`: the Eckert GPDC fit, solved for the velocity.
fn flooding_velocity(state: &PackingState, packing: &PackingSpecification) -> f64 {
    let mut flow_parameter = if state.vapor_mass_flow.value > 0.0 {
        (state.liquid_mass_flow.value / state.vapor_mass_flow.value)
            * (state.vapor_density.value / state.liquid_density.value).sqrt()
    } else {
        0.0
    };
    flow_parameter = flow_parameter.clamp(MIN_FLOW_PARAMETER, MAX_FLOW_PARAMETER);

    let x = flow_parameter.log10();
    let y_flood = 10f64.powf(ECKERT_C0 + ECKERT_C1 * x + ECKERT_C2 * x * x);

    let rho_factor =
        state.vapor_density.value / (G * (state.liquid_density.value - state.vapor_density.value));
    let mu_factor = (state.liquid_viscosity.value / WATER_VISCOSITY).powf(0.1);
    let velocity_squared = y_flood / (packing.packing_factor * rho_factor * mu_factor);
    velocity_squared.max(0.0).sqrt() * state.hydraulic_capacity_factor
}

/// `calculatePressureDrop`: Leva's dry form with its liquid-loading correction.
fn pressure_drop(
    state: &PackingState,
    packing: &PackingSpecification,
    vapor_velocity: f64,
    area: f64,
) -> f64 {
    let liquid_loading = state.liquid_mass_flow.value / area;
    let dry = packing.packing_factor * state.vapor_density.value * vapor_velocity * vapor_velocity
        / packing.void_fraction.powi(3)
        * 0.01;
    let leva_c = if packing.category.eq_ignore_ascii_case("structured") {
        0.010
    } else {
        0.015
    };
    (dry * 10f64.powf(leva_c * liquid_loading)).max(0.0)
}

/// `calculateMassTransfer`: Onda's wetted area and the two film coefficients.
fn mass_transfer(
    state: &PackingState,
    packing: &PackingSpecification,
    area: f64,
    liquid_velocity: f64,
) -> (f64, f64, f64) {
    let a = packing.specific_surface_area;
    let vapor_superficial = state.vapor_mass_flow.value / area;
    let liquid_superficial = state.liquid_mass_flow.value / area;
    if vapor_superficial <= 0.0 || liquid_superficial <= 0.0 {
        return (0.0, 0.0, 0.0);
    }

    // ---- The wetted area.
    let re_l = liquid_velocity * state.liquid_density.value / (state.liquid_viscosity.value * a);
    let fr_l = liquid_velocity * liquid_velocity * a / G;
    let we_l = liquid_velocity * liquid_velocity * state.liquid_density.value
        / (state.surface_tension.value * a);
    let sigma_ratio = packing.critical_surface_tension / state.surface_tension.value;
    let exponent = ONDA_WETTING_CONSTANT
        * sigma_ratio.powf(0.75)
        * re_l.powf(0.1)
        * fr_l.powf(-0.05)
        * we_l.powf(0.2);
    let wetted_area = (1.0 - exponent.exp()).clamp(MIN_WETTED_RATIO, 1.0) * a;

    let diameter = effective_packing_diameter(packing);

    // ---- The gas film coefficient.
    let re_v = vapor_superficial / (a * state.vapor_viscosity.value);
    let sc_v = state.vapor_viscosity.value / (state.vapor_density.value * state.vapor_diffusivity);
    let k_ga = ONDA_GAS_CONSTANT
        * (a * state.vapor_diffusivity / (diameter * diameter))
        * re_v.powf(0.7)
        * sc_v.powf(1.0 / 3.0);

    // ---- The liquid film coefficient.
    let re_l2 = liquid_superficial / (wetted_area * state.liquid_viscosity.value);
    let sc_l =
        state.liquid_viscosity.value / (state.liquid_density.value * state.liquid_diffusivity);
    let scale = (state.liquid_density.value / (state.liquid_viscosity.value * G)).powf(1.0 / 3.0);
    let k_la =
        ONDA_LIQUID_CONSTANT * re_l2.powf(2.0 / 3.0) * sc_l.powf(-0.5) * (a * diameter).powf(0.4)
            / scale
            * wetted_area;

    (wetted_area, k_ga, k_la)
}

/// `estimateHETP` and `getEffectivePackingDiameter`: the two empirical rules of thumb.
fn estimate_hetp(packing: &PackingSpecification, column_diameter: f64) -> f64 {
    if packing.category.eq_ignore_ascii_case("structured") {
        (100.0 / packing.specific_surface_area + 0.1).max(0.15)
    } else {
        let mut estimate = 0.12 + 0.012 * packing.nominal_size_mm;
        if column_diameter > 1.0 {
            estimate *= 1.0 + 0.1 * (column_diameter - 1.0);
        }
        estimate.max(0.15)
    }
}

/// `getEffectivePackingDiameter`: the nominal size where there is one, and `4 eps / a` otherwise.
fn effective_packing_diameter(packing: &PackingSpecification) -> f64 {
    if packing.nominal_size_mm > 0.0 {
        return packing.nominal_size_mm / 1000.0;
    }
    if packing.specific_surface_area > 0.0 && packing.void_fraction > 0.0 {
        return (4.0 * packing.void_fraction / packing.specific_surface_area).max(1.0e-4);
    }
    0.025
}
