//! `unit_ops.component_splitter` - a per-component split, with each outlet flashed.

use azoth_core::units::pascals;
use azoth_core::{AzothError, Result};

use crate::stream::Stream;

/// Divide a stream between two outlets **component by component**, flashing each.
///
/// **The factor is per component, not per outlet.** `ComponentSplitter.run` loops
/// `for i in 0..2` - an overhead and a bottoms - and reads `splitFactor[k]` as the fraction
/// of component *k* routed to the first, the second taking the remainder. So the two outlets
/// carry different *compositions* rather than the same state at different flows, which is
/// what separates this from `unit_ops.splitter`.
///
/// **Each outlet is flashed**, at the feed's own temperature and pressure:
/// `run` sets the component moles on an empty fluid, calls `init(0)` and runs a `TPflash`.
/// A component split is a material balance, and the flash is what makes each outlet a
/// *state* - so this is the one id in the splitter family whose outlets can land on
/// different phases.
///
/// **The class's `mass` basis is a no-op and this ports one basis.** The mass branch
/// computes `feedMoles * M * fraction / M`, which is `feedMoles * fraction` - the molar
/// split - so a caller who asks for the mass basis gets the same answer.
///
/// # Errors
/// [`azoth_core::AzothError::InvalidInput`] if the factors are not one per component, or if
/// any is outside `[0, 1]`, and whatever the flash refuses.
pub fn component_splitter(feed: &Stream, split_factors: &[f64]) -> Result<(Stream, Stream)> {
    if split_factors.len() != feed.z.len() {
        return Err(AzothError::invalid_input(
            "split_factors",
            format!(
                "a component split routes each component, so the factors are one per \
                 component: {} factors for {} components",
                split_factors.len(),
                feed.z.len()
            ),
        ));
    }
    for (index, factor) in split_factors.iter().enumerate() {
        if !(0.0..=1.0).contains(factor) {
            return Err(AzothError::invalid_input(
                "split_factors",
                format!(
                    "a component's split fraction is a fraction, and component {index}'s is \
                     {factor}"
                ),
            ));
        }
    }

    // The component molar flows the split divides, and the two sets of shares.
    let overhead: Vec<f64> = feed
        .z
        .iter()
        .zip(split_factors)
        .map(|(zi, factor)| feed.n * zi * factor)
        .collect();
    let bottoms: Vec<f64> = feed
        .z
        .iter()
        .zip(split_factors)
        .map(|(zi, factor)| feed.n * zi * (1.0 - factor))
        .collect();

    let outlet = |amounts: &[f64]| -> Result<Stream> {
        let n: f64 = amounts.iter().sum();
        if n <= 0.0 {
            // `run` sets the system without flashing when the total is zero, and a stream of
            // nothing is the feed's record at zero flow - which `from_pt` still resolves.
            return Stream::from_pt(
                feed.components.clone(),
                feed.z.clone(),
                0.0,
                pascals(feed.p.value),
                feed.t,
            );
        }
        let z: Vec<f64> = amounts.iter().map(|moles| moles / n).collect();
        Stream::from_pt(feed.components.clone(), z, n, pascals(feed.p.value), feed.t)
    };

    Ok((outlet(&overhead)?, outlet(&bottoms)?))
}
