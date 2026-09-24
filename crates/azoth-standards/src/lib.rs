//! Standard and regulatory calculations.
//!
//! The namespace NeqSim files under `standards/`, whose first member is ISO 6976: the
//! calorific values, density, relative density and Wobbe index of a natural gas from its
//! composition. It is here rather than in `azoth-eos` because it is not a property the
//! equation of state knows - it is a *standard's* tabulated constants, summed - and
//! because `unit_ops.flare`'s duty reads it through `Stream.LCV()`.
//!
//! * [`iso6976`] - calorific values and density from a composition
//! * [`iso6976_constants`] - the standard's per-component table, compiled from NeqSim

pub mod iso6976;
pub mod iso6976_constants;
pub mod model_gen;
pub mod results;

pub use iso6976::iso6976;
pub use results::Iso6976Result;
