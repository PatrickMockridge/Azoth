//! `unit_ops.splitter` - one stream, split into several with the same state.

use azoth_core::AzothError;

use crate::stream::Stream;

/// Split a feed into several outlets with the same composition, pressure,
/// temperature and molar enthalpy, scaled by the given fractions.
///
/// The fractions are normalised, so `[0.3, 0.7]` and `[3.0, 7.0]` mean the same
/// split. Nothing else changes - a splitter does not flash, mix or lose pressure.
pub fn splitter(feed: &Stream, fractions: &[f64]) -> azoth_core::Result<Vec<Stream>> {
    if fractions.is_empty() {
        return Err(AzothError::invalid_input(
            "fractions",
            "a splitter needs at least one outlet",
        ));
    }
    let total: f64 = fractions.iter().sum();
    if total <= 0.0 {
        return Err(AzothError::invalid_input(
            "fractions",
            "the split fractions must sum to a positive value",
        ));
    }

    Ok(fractions
        .iter()
        .map(|&f| Stream {
            components: feed.components.clone(),
            z: feed.z.clone(),
            n: feed.n * f / total,
            p: feed.p,
            t: feed.t,
            h: feed.h,
        })
        .collect())
}
