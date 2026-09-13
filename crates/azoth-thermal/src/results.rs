//! Result types for the thermal calculations.
//!
//! The same contract as the hydraulics crate's results: one struct per
//! calculation, field names identical to the Python result dataclass and listed in
//! [`CalcResult::FIELDS`], with a test asserting the three agree. See
//! `crates/azoth-core/src/result.rs` for why the duplication is deliberate.

use azoth_core::{CalcResult, Warning};
use uom::si::f64::Power;

/// Result of `thermal.conduction_plane_wall`.
#[derive(Debug, Clone, PartialEq)]
pub struct ConductionPlaneWallResult {
    /// Heat flow rate through the wall. Signed: it follows the sign of the
    /// temperature difference, and there is no warning for a negative value
    /// because a negative value is a direction rather than a problem.
    pub q: Power,
    /// Caveats.
    pub warnings: Vec<Warning>,
}

impl CalcResult for ConductionPlaneWallResult {
    const CALC_ID: &'static str = "thermal.conduction_plane_wall";
    const FIELDS: &'static [&'static str] = &["q", "warnings"];

    fn warnings(&self) -> &[Warning] {
        &self.warnings
    }
}
