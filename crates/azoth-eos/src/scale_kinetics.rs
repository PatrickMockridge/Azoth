//! `eos.scale_kinetics` - NeqSim's `ScaleKinetics`.
//!
//! A classical-nucleation induction time, a parabolic surface-growth rate and a
//! mass-transport rate over the supersaturation a saturation index states, with the
//! rate-controlling regime reported.

use azoth_core::units::{MolarConcentration, meters_per_second, seconds};
use azoth_core::{Result, apply_checks};

use crate::results::ScaleKineticsResult;
use crate::spec_gen;

/// The class's year and calcite molar volume, in seconds and cubic metres per mole.
const SECONDS_PER_YEAR: f64 = 3.1536e7;
const CALCITE_MOLAR_VOLUME: f64 = 37.0e-6;

/// The screening kinetics of a supersaturated solution.
///
/// # Errors
/// * [`azoth_core::AzothError::OutOfRange`] if an input is outside the spec's hard range.
///
/// # Example
/// ```
/// use azoth_core::units::moles_per_cubic_meter;
/// use azoth_eos::scale_kinetics::scale_kinetics;
///
/// let r = scale_kinetics(
///     1.0,               // SI
///     3.6,               // 0.001 h
///     3.0,               // B
///     1.5854895991882293e-10, // 5 mm/yr
///     2.0,               // order
///     1.0e-5,            // m/s
///     moles_per_cubic_meter(0.0),
///     moles_per_cubic_meter(0.0),
/// )?;
/// assert_eq!(r.limiting_regime, "REACTION");
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
#[allow(clippy::too_many_arguments)] // The signature is the spec's declared inputs.
pub fn scale_kinetics(
    saturation_index: f64,
    nucleation_constant: f64,
    nucleation_shape_factor: f64,
    surface_rate_constant: f64,
    surface_reaction_order: f64,
    mass_transfer_coefficient: f64,
    bulk_concentration: MolarConcentration,
    equilibrium_concentration: MolarConcentration,
) -> Result<ScaleKineticsResult> {
    let spec = &spec_gen::SCALE_KINETICS_SPEC;
    let mut warnings = Vec::new();

    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "nucleation_constant" => Some(nucleation_constant),
            "surface_rate_constant" => Some(surface_rate_constant),
            "mass_transfer_coefficient" => Some(mass_transfer_coefficient),
            "bulk_concentration" => Some(bulk_concentration.value),
            "equilibrium_concentration" => Some(equilibrium_concentration.value),
            _ => None,
        },
        &mut warnings,
    )?;

    // The setters' clamps: a negative duration, rate or concentration is zero.
    let nucleation_constant = nucleation_constant.max(0.0);
    let nucleation_shape_factor = nucleation_shape_factor.max(0.0);
    let surface_rate_constant = surface_rate_constant.max(0.0);
    let mass_transfer_coefficient = mass_transfer_coefficient.max(0.0);
    let bulk = bulk_concentration.value.max(0.0);
    let equilibrium = equilibrium_concentration.value.max(0.0);

    let supersaturation_ratio = 10.0f64.powf(saturation_index);

    let (
        induction_time,
        surface_limited_growth,
        transport_limited_growth,
        effective_growth_rate,
        limiting_regime,
    ) = if supersaturation_ratio <= 1.0 {
        (f64::INFINITY, 0.0, 0.0, 0.0, "NONE".to_owned())
    } else {
        let ln_s = supersaturation_ratio.ln();
        let induction_hours =
            nucleation_constant / 3600.0 * (nucleation_shape_factor / (ln_s * ln_s)).exp();

        let rate_mm_yr = surface_rate_constant * 1000.0 * SECONDS_PER_YEAR;
        let surface_mm_yr = rate_mm_yr * (supersaturation_ratio - 1.0).powf(surface_reaction_order);

        // The class's flux is `k_m * (c_bulk - c_eq)` in mol/(m2 s), with the driving
        // concentration in mol/m3, and its wall recession is that times the molar volume.
        let driving_mol_m3 = (bulk - equilibrium).max(0.0);
        let transport_m_s = mass_transfer_coefficient * driving_mol_m3 * CALCITE_MOLAR_VOLUME;
        let transport_mm_yr = transport_m_s * 1000.0 * SECONDS_PER_YEAR;

        let (effective_mm_yr, regime) =
            if transport_mm_yr <= 0.0 || surface_mm_yr <= transport_mm_yr {
                (surface_mm_yr, "REACTION")
            } else {
                (transport_mm_yr, "TRANSPORT")
            };

        (
            induction_hours * 3600.0,
            surface_mm_yr * 1e-3 / SECONDS_PER_YEAR,
            transport_mm_yr * 1e-3 / SECONDS_PER_YEAR,
            effective_mm_yr * 1e-3 / SECONDS_PER_YEAR,
            regime.to_owned(),
        )
    };

    Ok(ScaleKineticsResult {
        supersaturation_ratio,
        induction_time: seconds(induction_time),
        surface_limited_growth: meters_per_second(surface_limited_growth),
        transport_limited_growth: meters_per_second(transport_limited_growth),
        effective_growth_rate: meters_per_second(effective_growth_rate),
        limiting_regime,
        warnings,
    })
}
