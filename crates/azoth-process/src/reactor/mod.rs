//! The reactors' internals.
//!
//! **Not a palette surface.** `unit_ops.plug_flow_reactor` is one entry and its parameters are
//! the reactor's; [`stepper`] is the axial march `PlugFlowReactor.run` inlines, which has no
//! entry, no model id and no declaration of its own. It is here rather than in
//! `azoth-reactions` because it is not a reaction - it is the loop that walks a reactor's
//! length, and it is the first integrator anywhere in this workspace (`azoth-reactions`'
//! `linalg` is LU, rank and null-space work).
//!
//! **The class hands it an ordinary differential equation, not a paper.** `PlugFlowReactor`
//! carries an `IntegrationMethod { EULER, RK4 }` enum, defaults to RK4, and writes both
//! schemes out inline. So this module is the *class's* arithmetic and not a better one: the
//! callers of [`stepper::march`] get the same two schemes, the same fixed step
//! `length / numberOfSteps`, and the same absence of step-size control.

pub mod catalyst_bed;
pub mod gibbs_database;
pub mod gibbs_solver;
pub mod kinetic_reaction;
pub mod stepper;

pub use catalyst_bed::CatalystBed;
pub use gibbs_database::{GibbsDatabase, GibbsSpecies};
pub use gibbs_solver::{GibbsSettings, GibbsState, solve};
pub use kinetic_reaction::{KineticReaction, RateBasis, RateType};
pub use stepper::{Scheme, march};
