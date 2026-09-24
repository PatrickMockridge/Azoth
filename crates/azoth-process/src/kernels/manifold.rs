//! `unit_ops.manifold` - a mixer and a splitter composed, with a low-flow rule between.

use azoth_core::{AzothError, Result};

use crate::kernels::{mixer, splitter};
use crate::stream::Stream;

/// The low-flow threshold `Mixer.mixStream` compares against, in kg/hr.
///
/// `ProcessEquipmentBaseClass.DEFAULT_MINIMUM_FLOW`, and **the palette declares no way to
/// change it**, so this is the whole of the rule here. It is inert for every flow a case can
/// state: `1e-20` kg/hr is `1e-24` mol/s for a light hydrocarbon.
const DEFAULT_MINIMUM_FLOW_KG_PER_HOUR: f64 = 1e-20;

/// Join a manifold's feeds, then divide the mixture between its outlets.
///
/// **`Manifold.run` is four statements**, and the third and fourth are the composition:
/// `propagateMinimumFlow()` pushes the threshold onto both children, `localmixer.run()`
/// joins the feeds, `refreshLocalSplitter()` attaches the splitter's inlet to the mixture,
/// and `localsplitter.run()` divides it. Nothing else in the class is the manifold's own, so
/// this composes [`mixer`] and [`splitter`] rather than restating either.
///
/// **The one rule that is new is the low-flow filter**, and it belongs to the mixer: an
/// inlet whose kg/hr is at or below the threshold is not mixed. Measured, a feed stated at
/// zero flow leaves the mixture as the other feed alone - the captured third row, where the
/// product's flow and composition are the second feed's exactly.
///
/// # Errors
/// [`azoth_core::AzothError::InvalidInput`] if there is no feed, if the feeds' shapes
/// disagree, or if the factors are not positive, and whatever the mixer or the splitter
/// refuses.
pub fn manifold(feeds: &[Stream], split_factors: &[f64]) -> Result<Vec<Stream>> {
    let mut active: Vec<Stream> = Vec::with_capacity(feeds.len());
    for feed in feeds {
        // `mass_flow` in kg/s against a threshold NeqSim states in kg/hr, so the conversion
        // is written out rather than folded into the constant.
        if feed.mass_flow()? * 3600.0 > DEFAULT_MINIMUM_FLOW_KG_PER_HOUR {
            active.push(feed.clone());
        }
    }
    if active.is_empty() {
        return Err(AzothError::invalid_input(
            "feeds",
            "a manifold with no feed above the low-flow threshold has nothing to divide",
        ));
    }

    let product = mixer(&active, None)?;
    splitter(&product, split_factors)
}
