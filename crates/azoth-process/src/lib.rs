//! Unit operations: the process layer.
//!
//! A **unit operation** takes one or more streams and returns one or more streams,
//! changing their state by a physical rule: a separator splits a feed into a gas and a
//! liquid at the same temperature and pressure, a valve drops the pressure at constant
//! enthalpy, a mixer blends several feeds into one. This crate is where those live, and
//! it is the fourth level of the composition in `docs/src/spec.md` S4 - above a
//! calculation and a model, below a flowsheet.
//!
//! # A unit operation is a model, and nothing about it is new
//!
//! Each one has a spec under `specs/models/process/`, a spec-declared procedure, two
//! implementations and a worked example. What differs from `azoth-eos` is only what the
//! procedure is *over*: a mixture state that carries a molar flow and a composition.
//!
//! That is the whole of the port's Pareto argument, and it is measurable rather than
//! asserted. NeqSim's `Separator.run` is 109 lines of a 4,511-line file and its physics
//! is one `TPflash` followed by a phase split; its `ThrottlingValve.run` is 91 lines of
//! 1,894 and its physics is an isenthalpic flash. The rest of each file is performance
//! charts, entrainment models, geometry and mechanical design, none of which is here.
//!
//! # Why this crate depends on another
//!
//! [`azoth-eos`] is a sibling domain, and the rule stated in `docs/src/spec.md` S3 is
//! that a domain depends only on [`azoth_core`]. This crate is the exception that rule
//! has to make room for: a unit operation *is* a flash call plus arithmetic, so a
//! process layer that could not call the flashes would not be a process layer. S3
//! records the tier - cross-domain composition lives above the domains - and this is
//! the first member of it.
//!
//! # Molar flow is a bare `f64`
//!
//! Every unit operation takes its flow in mol/s as a plain number rather than a `uom`
//! quantity, which is the one place in this workspace where a dimensional quantity is
//! not carried by a type. The reason is that `uom` has no molar-flow quantity to use.
//! `mol/s` is in the units vocabulary, its conversion path is the identity because it
//! is already SI base, and `crates/azoth-core/src/units.rs` records the same reasoning
//! where a reader will meet it.

use azoth_core::{AzothError, ModelSpec, Result};

/// The algorithm a model's spec fixes, which its `kind` says it has.
///
/// Re-exported from [`azoth_eos`] rather than copied. The check is about the *spec*
/// rather than about any one model, and a second copy is a second place for it to
/// drift - which is the same argument that put it in one place there.
pub use azoth_eos::algorithm_of;

pub mod compressor;
pub mod expander;
pub mod heater;
pub mod isentropic;
pub mod mixer;
pub mod model_gen;
pub mod pump;
pub mod results;
pub mod separator;
pub mod splitter;
pub mod throttling_valve;

pub use compressor::compressor;
pub use expander::expander;
pub use heater::heater;
pub use mixer::mixer;
pub use pump::pump;
pub use results::{
    CompressorResult, ExpanderResult, HeaterResult, MixerResult, PumpResult, SeparatorResult,
    SplitterResult, ThrottlingValveResult,
};
pub use separator::separator;
pub use splitter::splitter;
pub use throttling_valve::throttling_valve;

/// Reported when a flash reports a phase that contradicts its own vapour fraction.
///
/// A model that says `two_phase` and reports no `beta` has contradicted itself, and
/// there is no defensible way to route the streams: `beta` absent means "the feed is
/// single phase, and this model does not know which", so a unit operation that split on
/// it would be inventing a number. Returned rather than panicked, per this workspace's
/// no-panic rule.
pub(crate) fn beta_of(model: &ModelSpec, beta: Option<f64>) -> Result<f64> {
    beta.ok_or_else(|| AzothError::InvalidInput {
        field: "beta".to_string(),
        reason: format!(
            "`{}` reports a two-phase split with no vapour fraction, so there is nothing \
             to split the feed by",
            model.id
        ),
    })
}
