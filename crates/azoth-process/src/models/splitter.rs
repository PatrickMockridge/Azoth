//! `process.splitter` - the splitter's kernel as a registered id.
//!
//! Spec: `specs/models/process/splitter.toml`. The arithmetic is
//! [`crate::kernels::splitter`]; what is here is the boundary a case and a cross-impl
//! test address.
//!
//! **This is the first id with a `many` outlet port**, so its result is the shape the
//! module doc describes: one field per record field, each carrying one entry per outlet,
//! and `z` a matrix with one row per outlet.

use azoth_core::units::{MolarEnergy, Pressure, ThermodynamicTemperature, joules_per_mole};
use azoth_core::{CalcResult, Result, Warning, apply_checks};

use crate::kernels::splitter as kernel;
use crate::model_gen;
use crate::stream::Stream;

/// Result of `process.splitter`.
#[derive(Debug, Clone, PartialEq)]
pub struct SplitterResult {
    /// Molar flow of each outlet, mol/s.
    pub products_n: Vec<f64>,
    /// Composition of each outlet, one row per outlet.
    pub products_z: Vec<Vec<f64>>,
    /// Pressure of each outlet, which a splitter does not change.
    pub products_p: Vec<Pressure>,
    /// Temperature of each outlet, which a splitter does not change.
    pub products_t: Vec<ThermodynamicTemperature>,
    /// Molar enthalpy of each outlet, which is the feed's.
    pub products_h: Vec<MolarEnergy>,
    /// Caveats.
    pub warnings: Vec<Warning>,
}

impl CalcResult for SplitterResult {
    const CALC_ID: &'static str = "process.splitter";
    const FIELDS: &'static [&'static str] = &[
        "products_n",
        "products_z",
        "products_p",
        "products_t",
        "products_h",
        "warnings",
    ];

    fn warnings(&self) -> &[Warning] {
        &self.warnings
    }
}

/// Split a stream into several with the same state, in proportion to `split_factors`.
///
/// # Errors
/// [`azoth_core::AzothError::InvalidInput`] if there is no outlet or the factors do not
/// sum to a positive value, and whatever the databank refuses.
pub fn splitter(
    components: &[String],
    feed_n: f64,
    feed_z: &[f64],
    feed_p: Pressure,
    feed_t: ThermodynamicTemperature,
    split_factors: &[f64],
) -> Result<SplitterResult> {
    let spec = &model_gen::SPLITTER_SPEC;
    let mut warnings = Vec::new();
    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "feed_t" => Some(feed_t.value),
            _ => None,
        },
        &mut warnings,
    )?;

    let feed = Stream::from_pt(components.to_vec(), feed_z.to_vec(), feed_n, feed_p, feed_t)?;
    let out = kernel(&feed, split_factors)?;

    Ok(SplitterResult {
        products_n: out.iter().map(|s| s.n).collect(),
        products_z: out.iter().map(|s| s.z.clone()).collect(),
        products_p: out.iter().map(|s| s.p).collect(),
        products_t: out.iter().map(|s| s.t).collect(),
        products_h: out.iter().map(|s| joules_per_mole(s.h.value)).collect(),
        warnings,
    })
}
