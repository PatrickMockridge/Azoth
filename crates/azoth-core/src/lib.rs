//! Core types for `azoth`, an open and validated library of chemical
//! engineering calculations.
//!
//! This crate deliberately contains no engineering calculations. It holds the
//! vocabulary they share - units, warnings, errors, results, range checks - so
//! that every calculation in the workspace reports failures, caveats and units
//! the same way.
//!
//! [`solver`] is the one module that is *machinery* rather than vocabulary, and
//! it is here for the same reason: a named solution scheme is shared the way a
//! unit is, because two implementations of one equation must run the same
//! scheme to agree numerically. Like [`range`], it is not a calculation - it has
//! no equation and no spec - and it grows only when the spec schema's
//! `solver.kind` grows.
//!
//! The two ideas worth understanding before reading anything else:
//!
//! 1. **Warnings are not errors.** A value outside the range in which a
//!    correlation was validated is still a value, and refusing to return it
//!    would be less useful than returning it with a warning. What the library
//!    must never do is return it *silently*. See [`warning`].
//!
//! 2. **A check that could not run is not a check that passed.** When an
//!    optional input is missing, the range check that depends on it emits
//!    [`WarningCode::RangeCheckSkipped`] rather than quietly succeeding. See
//!    [`range`].
//!
//! Every calculation also carries a [`CalcResult::CALC_ID`] matching its spec
//! file and a [`CalcResult::FIELDS`] list that is asserted against the Python
//! result dataclass, so the two implementations cannot drift apart unnoticed.

pub mod error;
pub mod range;
pub mod result;
pub mod solver;
pub mod spec;
pub mod units;
pub mod warning;

pub use error::{AzothError, Result, Warned};
pub use range::{Band, RangeCheck, Severity, apply_checks};
pub use result::{CalcResult, FlowRegime};
pub use solver::{
    Convergence, CubicRootsOutcome, SolverKind, SolverOutcome, cubic_roots, fixed_point,
    require_converged, require_cubic_converged,
};
pub use spec::{CalcSpec, SolverSpec, SpecCheck, TestCase};
pub use warning::{Warning, WarningCode};

/// Version of the `azoth` library.
pub const VERSION: &str = env!("CARGO_PKG_VERSION");
