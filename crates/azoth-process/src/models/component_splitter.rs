//! `process.component_splitter` - the component splitter's kernel as a registered id.
//!
//! Spec: `specs/models/process/component_splitter.toml`. The arithmetic is
//! [`crate::kernels::component_splitter`]; what is here is the boundary a case and a
//! cross-impl test address.
//!
//! **Two named outlets rather than a `many` one**, because the class fixes the count at two
//! and the palette entry's `products` is corrected to match in the commit that carries this
//! kernel.

use azoth_core::units::{MolarEnergy, Pressure, ThermodynamicTemperature, joules_per_mole};
use azoth_core::{CalcResult, Result, Warning, apply_checks};

use crate::kernels::component_splitter as kernel;
use crate::model_gen;
use crate::stream::Stream;

/// Result of `process.component_splitter`.
#[derive(Debug, Clone, PartialEq)]
pub struct ComponentSplitterResult {
    /// Overhead molar flow, mol/s.
    pub overhead_n: f64,
    /// Overhead composition.
    pub overhead_z: Vec<f64>,
    /// Overhead pressure.
    pub overhead_p: Pressure,
    /// Overhead temperature.
    pub overhead_t: ThermodynamicTemperature,
    /// Overhead molar enthalpy, the state's at its own composition.
    pub overhead_h: MolarEnergy,
    /// Bottoms molar flow, mol/s.
    pub bottoms_n: f64,
    /// Bottoms composition.
    pub bottoms_z: Vec<f64>,
    /// Bottoms pressure.
    pub bottoms_p: Pressure,
    /// Bottoms temperature.
    pub bottoms_t: ThermodynamicTemperature,
    /// Bottoms molar enthalpy, the state's at its own composition.
    pub bottoms_h: MolarEnergy,
    /// Caveats.
    pub warnings: Vec<Warning>,
}

impl CalcResult for ComponentSplitterResult {
    const CALC_ID: &'static str = "process.component_splitter";
    const FIELDS: &'static [&'static str] = &[
        "overhead_n",
        "overhead_z",
        "overhead_p",
        "overhead_t",
        "overhead_h",
        "bottoms_n",
        "bottoms_z",
        "bottoms_p",
        "bottoms_t",
        "bottoms_h",
        "warnings",
    ];

    fn warnings(&self) -> &[Warning] {
        &self.warnings
    }
}

/// Divide a stream between two outlets component by component.
///
/// # Errors
/// [`azoth_core::AzothError::InvalidInput`] if the factors are not one per component, or if
/// any is outside `[0, 1]`, and whatever the databank or a flash refuses.
pub fn component_splitter(
    components: &[String],
    feed_n: f64,
    feed_z: &[f64],
    feed_p: Pressure,
    feed_t: ThermodynamicTemperature,
    split_factors: &[f64],
) -> Result<ComponentSplitterResult> {
    let spec = &model_gen::COMPONENT_SPLITTER_SPEC;
    let mut warnings = Vec::new();
    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "feed_t" => Some(feed_t.value),
            "split_factors" => split_factors.first().copied(),
            _ => None,
        },
        &mut warnings,
    )?;

    let feed = Stream::from_pt(components.to_vec(), feed_z.to_vec(), feed_n, feed_p, feed_t)?;
    let (overhead, bottoms) = kernel(&feed, split_factors)?;

    Ok(ComponentSplitterResult {
        overhead_n: overhead.n,
        overhead_z: overhead.z,
        overhead_p: overhead.p,
        overhead_t: overhead.t,
        overhead_h: joules_per_mole(overhead.h.value),
        bottoms_n: bottoms.n,
        bottoms_z: bottoms.z,
        bottoms_p: bottoms.p,
        bottoms_t: bottoms.t,
        bottoms_h: joules_per_mole(bottoms.h.value),
        warnings,
    })
}
