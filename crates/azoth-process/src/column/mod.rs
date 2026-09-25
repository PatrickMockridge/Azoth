//! The distillation column's internals: the stage and the two ends.
//!
//! **Not a palette surface, and this paragraph is why.** `unit_ops.distillation_column` is
//! one entry, and its parameters are the column's; these are the pieces that entry's kernel
//! is built from - a tray, a condenser and a reboiler. They have no entry, no model id and no
//! palette declaration, because a stage is not a unit operation a flowsheet connects: the
//! column's two products are. The owner's decision was that the stages come first,
//! internally, and that separately-declared tray/condenser/reboiler palette entries are a
//! later expert feature.
//!
//! So a reader looking for `unit_ops.tray` finds nothing, which is the correct answer, and
//! this is where it is stated rather than left as an absence.
//!
//! **Each of the three is oracled on its own**, against
//! `validation/neqsim/captures/process_column_tray.tsv` and its siblings: every one of them
//! is a `SimulationInterface` and runs standalone, so the probe drives one stage directly and
//! prints both of its outlets. That is what lets D2 pin a stage's arithmetic before a column
//! exists to contain it.

pub mod block_tridiagonal;
pub mod condenser;
pub mod naphtali_sandholm;
pub mod reboiler;
pub mod tray;

pub use condenser::{CondenserMode, CondenserOutcome, condenser};
pub use reboiler::{ReboilerMode, ReboilerOutcome, reboiler};
pub use tray::{SideDraws, TrayOutcome, tray};

/// The vapour and liquid shares of a flash.
///
/// **The phase decides, and `beta` does not.** NeqSim's outlet getters look for a phase
/// *type*: a subcooled stage has no gas phase and reports a vapour of zero flow while its
/// `getBeta()` still answers `1.0`, so reading the split off `beta` would give that stage all
/// of its flow as vapour - the opposite of the measurement.
///
/// A trivial solution is refused: every K-value straddled one, so nothing proves which single
/// phase is there, and NeqSim reads a phase type this model does not derive.
pub(crate) fn phase_fractions(
    phase: azoth_eos::Phase,
    beta: Option<f64>,
) -> azoth_core::Result<(f64, f64)> {
    use azoth_core::AzothError;
    use azoth_eos::Phase;
    match phase {
        Phase::TwoPhase => beta.map(|b| (b, 1.0 - b)).ok_or_else(|| {
            AzothError::invalid_input("flash", "a two-phase answer with no vapour fraction")
        }),
        Phase::AllVapour => Ok((1.0, 0.0)),
        Phase::AllLiquid => Ok((0.0, 1.0)),
        Phase::Trivial => Err(AzothError::invalid_input(
            "flash",
            "the flash is a trivial solution - every K-value straddled one, so nothing here \
             proves which single phase is present",
        )),
    }
}
