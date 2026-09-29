//! `hydraulics.tray_hydraulics` - a tray's flooding, weeping, entrainment, pressure drop and
//! efficiency.
//!
//! Spec: `specs/calcs/hydraulics/tray_hydraulics.toml`, which records the five correlations, the
//! tray-type branches and the clamps.
//!
//! **This is `TrayHydraulicsCalculator`, and none of the three columns that would read it carries
//! its output.** `distillation_column`, `absorption_column` and `stripping_column` all declare
//! `max_allowable_gas_load_factor` and refuse it, and `getFsFactor()` is the quantity that would
//! answer it. The oracle is `validation/neqsim/captures/tray_hydraulics_probe.tsv`, eight states
//! driven through the class directly.

use azoth_core::units::{DynamicViscosity, Length, MassDensity, MassRate, SurfaceTension};
use azoth_core::{Result, apply_checks};

use crate::results::TrayHydraulicsResult;
use crate::spec_gen;

/// The gravitational acceleration the class writes inline, in m/s**2.
const G: f64 = 9.81;

/// The surface tension the flooding velocity's correction is referenced to, in N/m.
const REFERENCE_SURFACE_TENSION: f64 = 0.020;

/// The Fair capacity factor's correction exponents, both the class's own.
const FLV_CORRECTION_SCALE: f64 = 1.463;
const FLV_CORRECTION_POWER: f64 = 0.842;
/// The flow parameter the tabulated capacity factors are stated at, which normalises them.
const FLV_REFERENCE: f64 = 0.1;
/// The band the flow parameter is clamped to before the correction is taken.
const MIN_FLV: f64 = 0.01;
const MAX_FLV: f64 = 2.0;

/// The flow parameter's bands and Fair's `(a, b)` pairs for the entrainment fit.
const ENTRAINMENT_LOW_FLV: f64 = 0.02;
const ENTRAINMENT_HIGH_FLV: f64 = 0.2;
const ENTRAINMENT_LOW: (f64, f64) = (0.085, 4.2);
const ENTRAINMENT_MID: (f64, f64) = (0.05, 3.8);
const ENTRAINMENT_HIGH: (f64, f64) = (0.03, 3.5);
/// The entrainment above which the verdict turns.
const ENTRAINMENT_LIMIT: f64 = 0.1;

/// The orifice discharge coefficients, one per tray type.
const ORIFICE_SIEVE: f64 = 0.73;
const ORIFICE_VALVE: f64 = 0.80;
const ORIFICE_OTHER: f64 = 0.65;

/// The flooding velocity's tray-type factors.
const FLOOD_FACTOR_VALVE: f64 = 1.05;
const FLOOD_FACTOR_BUBBLE_CAP: f64 = 0.85;

/// The hole area's fallback fractions, for the types that do not state one.
const HOLE_FRACTION_VALVE: f64 = 0.14;
const HOLE_FRACTION_OTHER: f64 = 0.12;

/// The weir length's auto-derived fraction of the diameter.
const WEIR_LENGTH_FRACTION: f64 = 0.73;

/// The weeping minimum's two constants, from Sinnott's chart fit.
const WEEPING_CHART_BASE: f64 = 7.0;
const WEEPING_CHART_SLOPE: f64 = 0.27;
const WEEPING_HOLE_SLOPE: f64 = 0.90;
/// The reference hole diameter Sinnott's term is written against, in mm.
const WEEPING_REFERENCE_HOLE_MM: f64 = 25.4;
/// The numerator the class substitutes when Sinnott's term goes negative.
const WEEPING_NUMERATOR_FLOOR: f64 = 0.5;

/// O'Connell's efficiency constants and its two clamps.
const OCONNELL_BASE: f64 = 51.0;
const OCONNELL_SLOPE: f64 = 32.5;
const OCONNELL_ALPHA_MU_MIN: f64 = 0.1;
const OCONNELL_ALPHA_MU_MAX: f64 = 10.0;
const OCONNELL_EFFICIENCY_MIN: f64 = 10.0;
const OCONNELL_EFFICIENCY_MAX: f64 = 100.0;

/// The design window the flood is judged inside, in per cent.
const DESIGN_FLOOD_MIN: f64 = 50.0;
const DESIGN_FLOOD_MAX: f64 = 85.0;
/// The downcomer backup's allowance, as a fraction of `spacing + weir height`.
const BACKUP_LIMIT_FRACTION: f64 = 0.5;
/// The apron gap the downcomer loss is taken over, as a fraction of the weir length.
const APRON_GAP_FRACTION: f64 = 0.025;
/// The apron's own coefficient and its floor, both the class's.
const APRON_LOSS_COEFFICIENT: f64 = 166.0;
const APRON_AREA_FLOOR_M2: f64 = 0.001;

/// Everything `TrayHydraulicsCalculator.calculate` reads.
#[derive(Debug, Clone)]
pub struct TrayHydraulicsState {
    /// `sieve`, `valve` or `bubble-cap`, which four branches read.
    pub tray_type: String,
    /// The column's internal diameter.
    pub column_diameter: Length,
    /// The tray spacing.
    pub tray_spacing: Length,
    /// The weir height.
    pub weir_height: Length,
    /// The weir length. **At or below zero derives `0.73 D`**, which is the class's own rule.
    pub weir_length: Length,
    /// The downcomer area as a fraction of the total.
    pub downcommer_area_fraction: f64,
    /// The hole diameter. **The spec declares it in millimetres and the class works in them**,
    /// so the kernel converts back at the two places that read it.
    pub hole_diameter: Length,
    /// The hole area as a fraction of the active area, read for a `sieve` tray alone.
    pub hole_area_fraction: f64,
    /// The design flooding fraction. **Stored and read by nothing.**
    pub design_flood_fraction: f64,
    /// The vapour's mass flow.
    pub vapor_mass_flow: MassRate,
    /// The liquid's mass flow.
    pub liquid_mass_flow: MassRate,
    /// The vapour's mass density.
    pub vapor_density: MassDensity,
    /// The liquid's mass density.
    pub liquid_density: MassDensity,
    /// The liquid's dynamic viscosity.
    pub liquid_viscosity: DynamicViscosity,
    /// The interfacial surface tension.
    pub surface_tension: SurfaceTension,
    /// The relative volatility of the key components.
    pub relative_volatility: f64,
}

/// `TrayHydraulicsCalculator.calculate`: the ten steps, in the class's order.
///
/// # Errors
/// [`azoth_core::AzothError::OutOfRange`] for any declared input outside the range its own
/// division requires - see the spec's `[[valid_range]]` rows.
pub fn tray_hydraulics(state: TrayHydraulicsState) -> Result<TrayHydraulicsResult> {
    let spec = &spec_gen::TRAY_HYDRAULICS_SPEC;
    let mut warnings = Vec::new();

    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "column_diameter" => Some(state.column_diameter.value),
            "tray_spacing" => Some(state.tray_spacing.value),
            "weir_height" => Some(state.weir_height.value),
            "hole_diameter" => Some(state.hole_diameter.value * 1000.0),
            "hole_area_fraction" => Some(state.hole_area_fraction),
            "downcommer_area_fraction" => Some(state.downcommer_area_fraction),
            "vapor_mass_flow" => Some(state.vapor_mass_flow.value),
            "liquid_mass_flow" => Some(state.liquid_mass_flow.value),
            "vapor_density" => Some(state.vapor_density.value),
            "liquid_density" => Some(state.liquid_density.value),
            "liquid_viscosity" => Some(state.liquid_viscosity.value),
            "surface_tension" => Some(state.surface_tension.value),
            "relative_volatility" => Some(state.relative_volatility),
            _ => None,
        },
        &mut warnings,
    )?;

    let diameter = state.column_diameter.value;
    let spacing = state.tray_spacing.value;
    let weir_height_m = state.weir_height.value;
    let vapor_mass_flow = state.vapor_mass_flow.value;
    let liquid_mass_flow = state.liquid_mass_flow.value;
    let vapor_density = state.vapor_density.value;
    let liquid_density = state.liquid_density.value;
    let liquid_viscosity = state.liquid_viscosity.value;
    let surface_tension = state.surface_tension.value;
    let hole_diameter_mm = state.hole_diameter.value * 1000.0;
    let sieve = state.tray_type.eq_ignore_ascii_case("sieve");
    let valve = state.tray_type.eq_ignore_ascii_case("valve");
    let bubble_cap = state.tray_type.eq_ignore_ascii_case("bubble-cap");

    // ---- `calculateAreas`. `activeAreaFraction` is derived here and not an input.
    let total_area = std::f64::consts::FRAC_PI_4 * diameter * diameter;
    let downcommer_area = total_area * state.downcommer_area_fraction;
    let active_area_fraction = 1.0 - 2.0 * state.downcommer_area_fraction;
    let active_area = total_area * active_area_fraction;
    let hole_area = active_area
        * if sieve {
            state.hole_area_fraction
        } else if valve {
            HOLE_FRACTION_VALVE
        } else {
            HOLE_FRACTION_OTHER
        };
    // **The derived weir length is what every crest term reads**, so a stated one is not merely
    // reported back - it replaces the derived value everywhere below.
    let weir_length = if state.weir_length.value <= 0.0 {
        WEIR_LENGTH_FRACTION * diameter
    } else {
        state.weir_length.value
    };

    // ---- `calculateFloodingVelocity`.
    let mut flv = 0.0;
    if vapor_mass_flow > 0.0 && vapor_density > 0.0 && liquid_density > 0.0 {
        flv = (liquid_mass_flow / vapor_mass_flow) * (vapor_density / liquid_density).sqrt();
    }
    let k_base = capacity_factor(spacing, flv);
    let sigma_corrected = k_base * (surface_tension / REFERENCE_SURFACE_TENSION).powf(0.2);
    let tray_factor = if valve {
        FLOOD_FACTOR_VALVE
    } else if bubble_cap {
        FLOOD_FACTOR_BUBBLE_CAP
    } else {
        1.0
    };
    let k_final = sigma_corrected * tray_factor;
    let delta_rho = (liquid_density - vapor_density).max(0.1);
    let flooding_velocity = k_final * (delta_rho / vapor_density).sqrt();

    // ---- `calculateActualVelocity`. The net area is the total less **one** downcomer.
    let net_area = total_area - downcommer_area;
    let (actual_vapor_velocity, fs_factor, percent_flood) =
        if net_area <= 0.0 || vapor_density <= 0.0 {
            (0.0, 0.0, 0.0)
        } else {
            let vapor_volume_flow = vapor_mass_flow / vapor_density;
            let actual = vapor_volume_flow / net_area;
            let fs = actual * vapor_density.sqrt();
            let flood = if flooding_velocity > 0.0 {
                (actual / flooding_velocity) * 100.0
            } else {
                0.0
            };
            (actual, fs, flood)
        };

    // ---- `calculateWeepingCheck`.
    let liquid_vol_flow = liquid_mass_flow / liquid_density;
    let crest = |weir: f64| {
        if weir > 0.0 {
            750.0 * (liquid_vol_flow / weir).powf(2.0 / 3.0)
        } else {
            0.0
        }
    };
    let actual_hole_velocity = if hole_area > 0.0 && vapor_density > 0.0 {
        (vapor_mass_flow / vapor_density) / hole_area
    } else {
        0.0
    };
    // The class returns before the check for a bubble-cap tray, leaving `turndownRatio` at its
    // field initialiser - which is what the captured row's `0.0` is.
    let (minimum_vapor_velocity, weeping_ok, turndown_written) = if bubble_cap {
        (0.0, true, None)
    } else {
        let how_mm = crest(weir_length);
        let hw_plus_how = weir_height_m * 1000.0 + how_mm;
        let kw = if hw_plus_how < 5.0 {
            WEEPING_CHART_BASE
        } else {
            WEEPING_CHART_BASE + WEEPING_CHART_SLOPE * hw_plus_how.sqrt()
        };
        let mut numerator =
            kw - WEEPING_HOLE_SLOPE * (WEEPING_REFERENCE_HOLE_MM - hole_diameter_mm);
        if numerator < 0.0 {
            numerator = WEEPING_NUMERATOR_FLOOR;
        }
        let u_min_hole = numerator / vapor_density.max(0.01).sqrt();
        let minimum = if hole_area > 0.0 {
            u_min_hole * hole_area / (total_area - downcommer_area)
        } else {
            0.0
        };
        let ratio = if u_min_hole > 0.0 {
            actual_hole_velocity / u_min_hole
        } else {
            1.0
        };
        (minimum, actual_hole_velocity >= u_min_hole, Some(ratio))
    };

    // ---- `calculateEntrainment`.
    let entrainment = if percent_flood / 100.0 <= 0.0 {
        0.0
    } else {
        let mut flv_entrainment = 0.0;
        if vapor_mass_flow > 0.0 {
            flv_entrainment =
                (liquid_mass_flow / vapor_mass_flow) * (vapor_density / liquid_density).sqrt();
        }
        let (a, b) = if flv_entrainment < ENTRAINMENT_LOW_FLV {
            ENTRAINMENT_LOW
        } else if flv_entrainment < ENTRAINMENT_HIGH_FLV {
            ENTRAINMENT_MID
        } else {
            ENTRAINMENT_HIGH
        };
        (a * (percent_flood / 100.0).powf(b)).min(1.0)
    };
    let entrainment_ok = entrainment < ENTRAINMENT_LIMIT;

    // ---- `calculatePressureDrop`.
    let orifice_coefficient = if sieve {
        ORIFICE_SIEVE
    } else if valve {
        ORIFICE_VALVE
    } else {
        ORIFICE_OTHER
    };
    let dry_mm = if orifice_coefficient > 0.0 && liquid_density > 0.0 {
        51.0 * (actual_hole_velocity / orifice_coefficient).powi(2)
            * (vapor_density / liquid_density)
    } else {
        0.0
    };
    let dry = dry_mm * liquid_density * G / 1000.0;
    let h_liquid_mm = weir_height_m * 1000.0 + crest(weir_length);
    let liquid_head = h_liquid_mm * liquid_density * G / 1000.0;
    let residual_mm = if hole_diameter_mm > 0.0 && liquid_density > 0.0 {
        6.0 * surface_tension * 1000.0 / (liquid_density * G * hole_diameter_mm / 1000.0)
    } else {
        0.0
    };
    let residual = residual_mm * liquid_density * G / 1000.0;
    let total = dry + liquid_head + residual;

    // ---- `calculateDowncommerBackup`.
    let ht_liquid = total / (liquid_density * G);
    let apron_area = (weir_length * APRON_GAP_FRACTION).max(APRON_AREA_FLOOR_M2);
    let apron_flow = liquid_mass_flow / liquid_density.max(1.0);
    let apron = APRON_LOSS_COEFFICIENT * (apron_flow / apron_area).powi(2) / 1000.0;
    let liquid_in_downcomer = (weir_height_m * 1000.0 + crest(weir_length)) / 1000.0;
    let downcommer_backup = ht_liquid + liquid_in_downcomer + apron;
    let max_allowed = spacing + weir_height_m;
    let downcommer_backup_fraction = if max_allowed > 0.0 {
        downcommer_backup / max_allowed
    } else {
        0.0
    };
    let downcommer_backup_ok = downcommer_backup < BACKUP_LIMIT_FRACTION * max_allowed;

    // ---- `calculateTrayEfficiency`. O'Connell takes the viscosity in centipoise.
    let alpha_mu = (state.relative_volatility * liquid_viscosity * 1000.0)
        .clamp(OCONNELL_ALPHA_MU_MIN, OCONNELL_ALPHA_MU_MAX);
    let eo_percent = (OCONNELL_BASE - OCONNELL_SLOPE * alpha_mu.log10())
        .clamp(OCONNELL_EFFICIENCY_MIN, OCONNELL_EFFICIENCY_MAX);
    let tray_efficiency = eo_percent / 100.0;

    // ---- `calculateTurndownRatio`. **The same number as the weeping check's ratio** wherever
    // both are defined, so it only bites where `minimum_vapor_velocity` is zero - which is a
    // bubble-cap tray, and there the field keeps its initialiser.
    let turndown_ratio = match turndown_written {
        Some(_) if actual_vapor_velocity > 0.0 && minimum_vapor_velocity > 0.0 => {
            actual_vapor_velocity / minimum_vapor_velocity
        }
        Some(ratio) => ratio,
        None => 0.0,
    };

    // ---- `assessDesign`. Five conditions.
    let design_ok = weeping_ok
        && entrainment_ok
        && downcommer_backup_ok
        && (DESIGN_FLOOD_MIN..=DESIGN_FLOOD_MAX).contains(&percent_flood);

    Ok(TrayHydraulicsResult {
        flooding_velocity,
        actual_vapor_velocity,
        percent_flood,
        minimum_vapor_velocity,
        fs_factor,
        weeping_ok,
        entrainment,
        entrainment_ok,
        downcommer_backup: azoth_core::units::meters(downcommer_backup),
        downcommer_backup_fraction,
        downcommer_backup_ok,
        total_tray_pressure_drop: azoth_core::units::pascals(total),
        total_tray_pressure_drop_mbar: total / 100.0,
        dry_tray_pressure_drop: azoth_core::units::pascals(dry),
        liquid_head_pressure_drop: azoth_core::units::pascals(liquid_head),
        residual_head_pressure_drop: azoth_core::units::pascals(residual),
        tray_efficiency,
        turndown_ratio,
        calculated_weir_length: azoth_core::units::meters(weir_length),
        active_area: azoth_core::units::square_meters(active_area),
        total_area: azoth_core::units::square_meters(total_area),
        hole_area: azoth_core::units::square_meters(hole_area),
        downcommer_area: azoth_core::units::square_meters(downcommer_area),
        design_ok,
        warnings,
    })
}

/// `getCapacityFactor`: the tabulated Fair factors interpolated over spacing, then corrected.
///
/// The tabulation is the class's own - `0.15` m to `0.91` m in six rows, `0.025` to `0.105` -
/// with its linear extrapolation past the last row.
fn capacity_factor(spacing: f64, flv: f64) -> f64 {
    let table = [
        (0.15, 0.025),
        (0.23, 0.040),
        (0.30, 0.050),
        (0.46, 0.070),
        (0.61, 0.085),
        (0.91, 0.105),
    ];
    let k_at_flv_01 = if spacing <= table[0].0 {
        table[0].1
    } else if spacing > table[table.len() - 1].0 {
        table[table.len() - 1].1 + (spacing - table[table.len() - 1].0) * 0.02
    } else {
        let mut value = table[0].1;
        for pair in table.windows(2) {
            let (low, high) = (pair[0], pair[1]);
            if spacing <= high.0 {
                value = low.1 + (spacing - low.0) / (high.0 - low.0) * (high.1 - low.1);
                break;
            }
        }
        value
    };
    // **The clamp is the class's, and it bites**: a flow parameter below `0.01` or above `2` is
    // read at the bound, so the correction is bounded too.
    let flv_clamped = flv.clamp(MIN_FLV, MAX_FLV);
    let correction = (-FLV_CORRECTION_SCALE * flv_clamped.powf(FLV_CORRECTION_POWER)).exp();
    let normalizer = (-FLV_CORRECTION_SCALE * FLV_REFERENCE.powf(FLV_CORRECTION_POWER)).exp();
    k_at_flv_01 * correction / normalizer
}
