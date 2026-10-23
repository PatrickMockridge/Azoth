//! `eos.stream_scale_analyzer` - NeqSim's `StreamScaleAnalyzer`.
//!
//! The coupled multi-mineral equilibrium at one litre of water, with each mineral's
//! precipitated mass converted to a rate by the aqueous throughput, and the dominant
//! mineral reported.

use azoth_core::Result;
use azoth_core::apply_checks;
use azoth_core::units::{
    MassDensity, Pressure, ThermodynamicTemperature, VolumeRate, kilograms_per_second, liters,
};

use crate::model_gen;
use crate::multi_mineral_scale_equilibrium::multi_mineral_scale_equilibrium;
use crate::results::{StreamScaleAnalyzerResult, StreamScaleMineralRow};

/// The scale mass rates of an aqueous stream.
///
/// # Errors
/// * [`azoth_core::AzothError::OutOfRange`] if an input is outside a declared range.
/// * [`azoth_core::AzothError::SolverNotConverged`] is not returned: a capped solve is a
///   warning beside its last iterate, as the equilibrium reports it.
///
/// # Example
/// ```
/// use azoth_core::units::{
///     cubic_meters_per_second, kelvins, kilograms_per_cubic_meter, pascals,
/// };
/// use azoth_eos::stream_scale_analyzer::stream_scale_analyzer;
///
/// let r = stream_scale_analyzer(
///     kilograms_per_cubic_meter(1.0),
///     kilograms_per_cubic_meter(1.0),
///     kilograms_per_cubic_meter(3.0),
///     kilograms_per_cubic_meter(0.0),
///     kilograms_per_cubic_meter(0.0),
///     kilograms_per_cubic_meter(0.0),
///     kilograms_per_cubic_meter(0.2),
///     kilograms_per_cubic_meter(1.5),
///     kilograms_per_cubic_meter(60.0),
///     kelvins(343.15),
///     pascals(1.0e7),
///     pascals(1.0e5),
///     0.0,
///     false,
///     false,
///     "davies",
///     cubic_meters_per_second(100000.0 / 86400.0 / 1000.0),
///     None,
///     None,
///     None,
/// )?;
/// assert_eq!(r.dominant_scale, "BaSO4");
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
#[allow(clippy::too_many_arguments)] // The signature is the spec's declared inputs.
pub fn stream_scale_analyzer(
    calcium: MassDensity,
    barium: MassDensity,
    strontium: MassDensity,
    iron: MassDensity,
    magnesium: MassDensity,
    sodium: MassDensity,
    bicarbonate: MassDensity,
    sulphate: MassDensity,
    total_dissolved_solids: MassDensity,
    temperature: ThermodynamicTemperature,
    pressure: Pressure,
    co2_partial_pressure: Pressure,
    ph: f64,
    auto_ph: bool,
    second_order_pressure_correction: bool,
    activity_model: &str,
    water_flow: VolumeRate,
    bdot_ion_size: Option<f64>,
    bdot_b: Option<f64>,
    pitzer_ionic_strength_molal: Option<f64>,
) -> Result<StreamScaleAnalyzerResult> {
    let spec = &model_gen::STREAM_SCALE_ANALYZER_SPEC;
    let mut warnings = Vec::new();

    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "temperature" => Some(temperature.value),
            "pressure" => Some(pressure.value),
            "co2_partial_pressure" => Some(co2_partial_pressure.value),
            "water_flow" => Some(water_flow.value),
            _ => None,
        },
        &mut warnings,
    )?;

    // The analyser forces the equilibrium's water volume to one litre; everything else is
    // the equilibrium's own state, and its warnings name the ions and conditions.
    let equilibrium = multi_mineral_scale_equilibrium(
        calcium,
        barium,
        strontium,
        iron,
        magnesium,
        sodium,
        bicarbonate,
        sulphate,
        total_dissolved_solids,
        temperature,
        pressure,
        co2_partial_pressure,
        ph,
        auto_ph,
        second_order_pressure_correction,
        activity_model,
        bdot_ion_size,
        bdot_b,
        pitzer_ionic_strength_molal,
        Some(liters(1.0)),
        None,
        None,
        None,
    )?;
    // A skipped bound is the inner model's own optional-input check; the outer model neither
    // declares those bounds nor exposes two of the inputs, so it cannot report them.
    warnings.extend(
        equilibrium
            .warnings
            .iter()
            .filter(|warning| warning.code != azoth_core::WarningCode::RangeCheckSkipped)
            .cloned(),
    );

    let minerals: Vec<StreamScaleMineralRow> = equilibrium
        .minerals
        .iter()
        .map(|mineral| StreamScaleMineralRow {
            mineral: mineral.name.clone(),
            initial_si: mineral.initial_saturation_index,
            final_si: mineral.final_saturation_index,
            precipitated_mass: mineral.precipitated_mass,
            scale_rate: kilograms_per_second(mineral.precipitated_mass.value * water_flow.value),
        })
        .collect();

    let mut dominant = "none".to_owned();
    let mut best = 0.0;
    for row in &minerals {
        if row.precipitated_mass.value > best {
            best = row.precipitated_mass.value;
            dominant.clone_from(&row.mineral);
        }
    }

    Ok(StreamScaleAnalyzerResult {
        minerals,
        total_scale_mass: equilibrium.total_mass,
        total_scale_rate: kilograms_per_second(equilibrium.total_mass.value * water_flow.value),
        dominant_scale: dominant,
        has_scaling: equilibrium.total_mass.value > 0.0,
        warnings,
    })
}
