//! Equations of state.
//!
//! The third namespace, and the first that is not a correlation over pipe or wall
//! geometry. It exists because a cubic equation of state decomposes the same way
//! every other calculation here does - into constitutive coefficients, each with a
//! published source and a worked example - and the decomposition is worth having
//! explicitly, since a wrong coefficient is invisible downstream.
//!
//! * [`pr_kappa`] - the Peng-Robinson alpha-function coefficient
//! * [`pr_alpha_ab`] - the alpha function and the reduced attraction parameters
//! * [`pr_z_factor`] - the compressibility factor, the cubic's real roots
//!
//! # Why the coefficients come first
//!
//! The obvious shape for this namespace is one calc per equation of state, taking
//! a component and returning a Z factor. That shape hides the constants: a
//! transposed digit in the `omega**2` coefficient produces an equation that still
//! runs, still converges and is slightly wrong everywhere, and no check on the Z
//! factor can see it, because the Z factor is computed *from* the coefficient.
//!
//! Splitting the constitutive coefficients out gives each one its own spec, its own
//! worked example and its own cross-language test. It also makes a modification
//! cheap: Peng-Robinson-Stryjek-Vera changes the temperature dependence of
//! `kappa` and nothing else, so under this decomposition it is one new calc rather
//! than a fork of the whole chain.
//!
//! # Dimensionless by construction
//!
//! Everything here works in reduced variables - `Tr`, `Pr`, and the dimensionless
//! coefficients that follow from them. `A = 0.45724 * alpha * Pr / Tr**2` needs no
//! gas constant, no pressure unit and no temperature unit, and the vapour-liquid
//! mixing rule is linear in `A` and `B` exactly as it is in `a` and `b`. That is
//! why this crate has no `uom` dependency: the dimensional conversion belongs at
//! the boundary, in one place, where it can be tested once.

pub mod pr_alpha_ab;
pub mod pr_kappa;
pub mod pr_z_factor;
pub mod results;
pub mod spec_gen;

pub use pr_alpha_ab::{OMEGA_A, OMEGA_B, pr_alpha_ab};
pub use pr_kappa::pr_kappa;
pub use pr_z_factor::pr_z_factor;
pub use results::{PrAlphaAbResult, PrKappaResult, PrZFactorResult, RootStructure};
