//! `process.manifold` - the manifold's kernel as a registered id.
//!
//! Spec: `specs/models/process/manifold.toml`. The arithmetic is
//! [`crate::kernels::manifold`], which composes [`crate::kernels::mixer`] and
//! [`crate::kernels::splitter`] rather than restating either.
//!
//! **`many` at both ends**, which no other id in this tier is: the feeds cross as the same
//! vectors a `many` inlet writes and the outlets as the ones a `many` outlet does.

use azoth_core::units::{MolarEnergy, Pressure, ThermodynamicTemperature, joules_per_mole};
use azoth_core::{AzothError, CalcResult, Result, Warning, apply_checks};

use crate::kernels::manifold::manifold as kernel;
use crate::model_gen;
use crate::stream::Stream;

/// Result of `process.manifold`.
#[derive(Debug, Clone, PartialEq)]
pub struct ManifoldResult {
    /// Molar flow out, mol/s, one entry per outlet: the mixture's times its fraction.
    pub products_n: Vec<f64>,
    /// Outlet compositions, one row per outlet.
    pub products_z: Vec<Vec<f64>>,
    /// Outlet pressures, Pa.
    pub products_p: Vec<Pressure>,
    /// Outlet temperatures, K.
    pub products_t: Vec<ThermodynamicTemperature>,
    /// Outlet molar enthalpies, J/mol.
    pub products_h: Vec<MolarEnergy>,
    /// Caveats.
    pub warnings: Vec<Warning>,
}

impl CalcResult for ManifoldResult {
    const CALC_ID: &'static str = "process.manifold";
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

/// Join a manifold's feeds and divide the mixture between its outlets.
///
/// # Errors
/// [`azoth_core::AzothError::InvalidInput`] if there is no feed, if the feeds' shapes
/// disagree, or if the factors are not positive, and whatever the databank or a flash
/// refuses.
#[allow(clippy::too_many_arguments)] // one parameter per declared input, and there are six
pub fn manifold(
    components: &[String],
    feed_n: &[f64],
    feed_z: &[Vec<f64>],
    feed_p: &[Pressure],
    feed_t: &[ThermodynamicTemperature],
    split_factors: &[f64],
) -> Result<ManifoldResult> {
    let spec = &model_gen::MANIFOLD_SPEC;
    let mut warnings = Vec::new();
    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "feed_t" => feed_t.first().map(|t| t.value),
            "split_factors" => split_factors.first().copied(),
            _ => None,
        },
        &mut warnings,
    )?;

    let shape = |what: &str, got: usize| -> Result<()> {
        if got == feed_n.len() {
            Ok(())
        } else {
            Err(AzothError::invalid_input(
                what,
                format!(
                    "a manifold's feeds are one port, so every field has one entry per feed: \
                     `feed_n` has {} and `{what}` has {got}",
                    feed_n.len()
                ),
            ))
        }
    };
    shape("feed_z", feed_z.len())?;
    shape("feed_p", feed_p.len())?;
    shape("feed_t", feed_t.len())?;

    let feeds: Vec<Stream> = (0..feed_n.len())
        .map(|i| {
            Stream::from_pt(
                components.to_vec(),
                feed_z[i].clone(),
                feed_n[i],
                feed_p[i],
                feed_t[i],
            )
        })
        .collect::<Result<_>>()?;

    let products = kernel(&feeds, split_factors)?;

    Ok(ManifoldResult {
        products_n: products.iter().map(|out| out.n).collect(),
        products_z: products.iter().map(|out| out.z.clone()).collect(),
        products_p: products.iter().map(|out| out.p).collect(),
        products_t: products.iter().map(|out| out.t).collect(),
        products_h: products
            .iter()
            .map(|out| joules_per_mole(out.h.value))
            .collect(),
        warnings,
    })
}
