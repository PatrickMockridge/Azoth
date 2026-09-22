//! `process.mixer` - the mixer's kernel as a registered id.
//!
//! Spec: `specs/models/process/mixer.toml`. The arithmetic is
//! [`crate::kernels::mixer`]; what is here is the boundary a case and a cross-impl test
//! address.
//!
//! **The first id with a `many` inlet port.** Its feeds cross as the same vectors a
//! `many` outlet writes, one entry per inlet — and each inlet is four of the record's five
//! fields, because the kernel builds it with [`Stream::from_pt`] and the enthalpy is what
//! that computes.

use azoth_core::units::{MolarEnergy, Pressure, ThermodynamicTemperature, joules_per_mole};
use azoth_core::{AzothError, CalcResult, Result, Warning, apply_checks};

use crate::kernels::mixer as kernel;
use crate::model_gen;
use crate::stream::Stream;

/// Result of `process.mixer`.
#[derive(Debug, Clone, PartialEq)]
pub struct MixerResult {
    /// Molar flow out, mol/s: the feeds' sum.
    pub product_n: f64,
    /// Outlet composition, the flow-weighted average of the feeds'.
    pub product_z: Vec<f64>,
    /// Outlet pressure: `outlet_pressure` when given, else the lowest feed pressure.
    pub product_p: Pressure,
    /// Outlet temperature, solved from the joined enthalpy at the outlet pressure.
    pub product_t: ThermodynamicTemperature,
    /// Outlet molar enthalpy, the flow-weighted average of the feeds'.
    pub product_h: MolarEnergy,
    /// Caveats.
    pub warnings: Vec<Warning>,
}

impl CalcResult for MixerResult {
    const CALC_ID: &'static str = "process.mixer";
    const FIELDS: &'static [&'static str] = &[
        "product_n",
        "product_z",
        "product_p",
        "product_t",
        "product_h",
        "warnings",
    ];

    fn warnings(&self) -> &[Warning] {
        &self.warnings
    }
}

/// Join several streams into one, conserving molar flow and enthalpy.
///
/// # Errors
/// [`azoth_core::AzothError::InvalidInput`] if there is no feed, if the feeds' shapes
/// disagree, or if `outlet_pressure` is not positive, and whatever the databank refuses.
#[allow(clippy::too_many_arguments)] // one parameter per declared input, and there are six
pub fn mixer(
    components: &[String],
    feed_n: &[f64],
    feed_z: &[Vec<f64>],
    feed_p: &[Pressure],
    feed_t: &[ThermodynamicTemperature],
    outlet_pressure: Option<Pressure>,
) -> Result<MixerResult> {
    let spec = &model_gen::MIXER_SPEC;
    let mut warnings = Vec::new();
    // One call over every declared input, for the reason `reactions.chemical_equilibrium`
    // gives. `outlet_pressure` is optional, so `None` here is a skipped check and not a
    // pass - which is the whole reason the bound is declared on it rather than nowhere.
    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "outlet_pressure" => outlet_pressure.map(|p| p.value),
            _ => None,
        },
        &mut warnings,
    )?;

    // The feeds are one fluid, so the component list is declared once and every feed
    // carries a composition over it. A short feed is a caller error rather than something
    // to zip past, because a mixer that silently dropped an inlet would answer a smaller
    // question with a plausible number.
    let shape = |what: &str, got: usize| -> Result<()> {
        if got == feed_n.len() {
            Ok(())
        } else {
            Err(AzothError::invalid_input(
                what,
                format!(
                    "a mixer's feeds are one port, so every field has one entry per feed: \
                     `feed_n` has {} and `{what}` has {got}",
                    feed_n.len()
                ),
            ))
        }
    };
    shape("feed_z", feed_z.len())?;
    shape("feed_p", feed_p.len())?;
    shape("feed_t", feed_t.len())?;

    let inlets: Vec<Stream> = (0..feed_n.len())
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

    let product = kernel(&inlets, outlet_pressure)?;

    Ok(MixerResult {
        product_n: product.n,
        product_z: product.z,
        product_p: product.p,
        product_t: product.t,
        product_h: joules_per_mole(product.h.value),
        warnings,
    })
}
