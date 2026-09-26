//! `process.shortcut_distillation_column` - the FUG column as a registered id.
//!
//! Spec: `specs/models/process/shortcut_distillation_column.toml`. The arithmetic is
//! [`crate::kernels::shortcut_distillation_column`]; what is here is the boundary a case and
//! a cross-impl test address.
//!
//! **The first `procedure` in this namespace.** Every other `process.*` model is `direct`,
//! so this is the first one the generator emits an `[algorithm]` static for - the Underwood
//! bisection, with the class's own tolerance and iteration cap. Nothing on the process side
//! reads [`azoth_core::ModelSpec::algorithm`]; the block is the contract the two
//! implementations are held to, and the case is what holds them.

use azoth_core::units::{MolarEnergy, Power, Pressure, ThermodynamicTemperature, joules_per_mole};
use azoth_core::{CalcResult, Result, Warning, apply_checks};
use serde::Serialize;

use crate::executor::json::{scalar, warnings as wire_warnings};
use crate::kernels::shortcut_distillation_column as kernel;
use crate::model_gen;
use crate::stream::Stream;

/// Result of `process.shortcut_distillation_column`.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct ShortcutDistillationColumnResult {
    /// Distillate molar flow, mol/s.
    pub distillate_n: f64,
    /// Distillate composition.
    pub distillate_z: Vec<f64>,
    /// Distillate pressure.
    #[serde(serialize_with = "scalar")]
    pub distillate_p: Pressure,
    /// Distillate temperature.
    #[serde(serialize_with = "scalar")]
    pub distillate_t: ThermodynamicTemperature,
    /// Distillate molar enthalpy.
    #[serde(serialize_with = "scalar")]
    pub distillate_h: MolarEnergy,
    /// Bottoms molar flow, mol/s.
    pub bottoms_n: f64,
    /// Bottoms composition.
    pub bottoms_z: Vec<f64>,
    /// Bottoms pressure.
    #[serde(serialize_with = "scalar")]
    pub bottoms_p: Pressure,
    /// Bottoms temperature.
    #[serde(serialize_with = "scalar")]
    pub bottoms_t: ThermodynamicTemperature,
    /// Bottoms molar enthalpy.
    #[serde(serialize_with = "scalar")]
    pub bottoms_h: MolarEnergy,
    /// Fenske's minimum stages.
    pub minimum_stages: f64,
    /// Underwood's minimum reflux ratio.
    pub minimum_reflux_ratio: f64,
    /// Molokanov's actual stage count.
    pub actual_stages: f64,
    /// The reflux ratio the correlation was evaluated at.
    pub actual_reflux_ratio: f64,
    /// The feed stage counted from the top.
    pub feed_tray_number: i64,
    /// The condenser duty the class reports.
    #[serde(serialize_with = "scalar")]
    pub condenser_duty: Power,
    /// The reboiler duty the class reports.
    #[serde(serialize_with = "scalar")]
    pub reboiler_duty: Power,
    /// `alpha_LK/HK`.
    pub relative_volatility: f64,
    /// Caveats.
    #[serde(serialize_with = "wire_warnings")]
    pub warnings: Vec<Warning>,
}

impl ShortcutDistillationColumnResult {
    /// The result of one kernel call.
    ///
    /// **The warnings are the caller's, because the two callers have different ones.** A case
    /// runs the spec's checks through [`apply_checks`]; a flowsheet's equivalent is the
    /// *checker's*, which reports through the envelope's diagnostics rather than through a result.
    #[must_use]
    pub fn of(outcome: &kernel::ShortcutColumn, warnings: Vec<Warning>) -> Self {
        Self {
            distillate_n: outcome.distillate.n,
            distillate_z: outcome.distillate.z.clone(),
            distillate_p: outcome.distillate.p,
            distillate_t: outcome.distillate.t,
            distillate_h: joules_per_mole(outcome.distillate.h.value),
            bottoms_n: outcome.bottoms.n,
            bottoms_z: outcome.bottoms.z.clone(),
            bottoms_p: outcome.bottoms.p,
            bottoms_t: outcome.bottoms.t,
            bottoms_h: joules_per_mole(outcome.bottoms.h.value),
            minimum_stages: outcome.minimum_stages,
            minimum_reflux_ratio: outcome.minimum_reflux_ratio,
            actual_stages: outcome.actual_stages,
            actual_reflux_ratio: outcome.actual_reflux_ratio,
            feed_tray_number: outcome.feed_tray_number,
            condenser_duty: outcome.condenser_duty,
            reboiler_duty: outcome.reboiler_duty,
            relative_volatility: outcome.relative_volatility,
            warnings,
        }
    }
}

impl CalcResult for ShortcutDistillationColumnResult {
    const CALC_ID: &'static str = "process.shortcut_distillation_column";
    const FIELDS: &'static [&'static str] = &[
        "distillate_n",
        "distillate_z",
        "distillate_p",
        "distillate_t",
        "distillate_h",
        "bottoms_n",
        "bottoms_z",
        "bottoms_p",
        "bottoms_t",
        "bottoms_h",
        "minimum_stages",
        "minimum_reflux_ratio",
        "actual_stages",
        "actual_reflux_ratio",
        "feed_tray_number",
        "condenser_duty",
        "reboiler_duty",
        "relative_volatility",
        "warnings",
    ];

    fn warnings(&self) -> &[Warning] {
        &self.warnings
    }
}

/// Solve a shortcut distillation column.
///
/// # Errors
/// [`azoth_core::AzothError::InvalidInput`] for a key that is not a component, a relative
/// volatility at or below one, a reflux multiplier at or below one, a split that empties a
/// product, and whatever the feed's flash refuses.
#[allow(clippy::too_many_arguments)] // one parameter per declared input, and there are eleven
pub fn shortcut_distillation_column(
    components: &[String],
    feed_n: f64,
    feed_z: &[f64],
    feed_p: Pressure,
    feed_t: ThermodynamicTemperature,
    light_key: &str,
    heavy_key: &str,
    light_key_recovery_distillate: f64,
    heavy_key_recovery_bottoms: f64,
    reflux_ratio_multiplier: f64,
    condenser_pressure: Option<Pressure>,
    reboiler_pressure: Option<Pressure>,
) -> Result<ShortcutDistillationColumnResult> {
    let spec = &model_gen::SHORTCUT_DISTILLATION_COLUMN_SPEC;
    let mut warnings = Vec::new();
    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "reflux_ratio_multiplier" => Some(reflux_ratio_multiplier),
            "light_key_recovery_distillate" => Some(light_key_recovery_distillate),
            "heavy_key_recovery_bottoms" => Some(heavy_key_recovery_bottoms),
            "feed_t" => Some(feed_t.value),
            _ => None,
        },
        &mut warnings,
    )?;

    let feed = Stream::from_pt(components.to_vec(), feed_z.to_vec(), feed_n, feed_p, feed_t)?;
    let out = kernel(
        &feed,
        light_key,
        heavy_key,
        light_key_recovery_distillate,
        heavy_key_recovery_bottoms,
        reflux_ratio_multiplier,
        condenser_pressure,
        reboiler_pressure,
    )?;

    Ok(ShortcutDistillationColumnResult::of(&out, warnings))
}
