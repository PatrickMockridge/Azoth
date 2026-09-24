//! `unit_ops.stripping_column` - the tray stripper, which is the absorber renamed.
//!
//! `StrippingColumn extends AbsorptionColumn` and its ninety lines add **nothing** to the
//! arithmetic: `addStrippingGasStream` is `addGasInStream`, `addRichLiquidStream` is
//! `addSolventInStream`, and the two product getters are the base's under other names. The
//! class's own javadoc states the reason - absorption and stripping are one set of
//! counter-current equilibrium-stage equations, and the thermodynamic driving force sets the
//! direction of transfer - so the stripping gas enters stage 0 and the rich liquid the top
//! stage, exactly where an absorber's gas and solvent enter.
//!
//! So this module is a rename and not a second kernel: the equations, the tray and the solver
//! are [`super::absorption_column`]'s, and a reader looking for a second implementation finds
//! the same one under two names, which is what the class is.

pub use super::absorption_column::{
    AbsorberOutcome as StripperOutcome, AbsorberSetup as StripperSetup,
    absorption_column as stripping_column,
};
