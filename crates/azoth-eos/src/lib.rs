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
//! * [`pr_departure`] - fugacity coefficient and departure functions
//! * [`pr_molar_volume`] - molar volume, the one dimensional quantity here
//! * [`prsv_kappa`] - the Stryjek-Vera coefficient, for the same alpha function
//! * [`vdw1f_mix_binary`] - van der Waals one-fluid mixing for a binary
//! * [`rachford_rice_binary`] - the vapour fraction that solves Rachford-Rice
//! * [`pr_mass_density`] - mass density, and the end of the path to something useful
//!
//! # Models
//!
//! The models below are *models* rather than calculations: their specs fix procedures
//! rather than equations, they live under `specs/models/`, and `tools/gen_models.py`
//! generates them into [`model_gen`] rather than into `spec_gen`. Each composes the
//! kernels above and adds a search over them.
//!
//! * [`pure_saturation`] - the pressure at which a pure component's two roots have
//!   equal fugacity, by bisection
//! * [`pt_flash`] - the two-phase split of a mixture at a fixed state, by successive
//!   substitution
//! * [`stability_test`] - whether a feed is stable as a single phase, by the
//!   tangent-plane criterion
//! * [`bubble_pressure`] - the pressure at which a liquid first gives off vapour
//! * [`dew_pressure`] - the pressure at which a vapour first condenses
//!
//! The last four share the mixture fugacity coefficient in [`mixture`], which no
//! registered calculation covers, and the two phase-boundary models share their
//! iteration in [`phase_boundary`].
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

use azoth_core::{AzothError, ModelAlgorithm, ModelSpec, Result};

/// The algorithm a model's spec fixes, which its `kind` says it has.
///
/// The schema requires an `algorithm` block for a `procedure` and forbids one for a
/// `direct` model, so a missing one where a procedure is expected means the generated
/// table and the schema disagree - a generator defect rather than a caller error.
/// Returned rather than panicked, per this crate's no-panic rule.
///
/// Here rather than in each of the three models that need it, because three copies of
/// a check is three places for it to drift, and the check is about the *spec* rather
/// than about any one model.
pub(crate) fn algorithm_of(spec: &ModelSpec) -> Result<&'static ModelAlgorithm> {
    spec.algorithm.ok_or_else(|| AzothError::InvalidInput {
        field: "algorithm".to_string(),
        reason: format!(
            "model `{}` is declared `{}` but its spec carries no algorithm, so there is \
             no procedure to run",
            spec.id, spec.kind
        ),
    })
}

pub mod bubble_pressure;
pub mod critical_point;
pub mod dew_pressure;
pub mod ideal_gas_cp;
pub mod mixture;
pub mod model_gen;
pub mod molar_enthalpy_entropy;
pub mod ph_flash;
pub mod phase_boundary;
pub mod pr_alpha_ab;
pub mod pr_departure;
pub mod pr_kappa;
pub mod pr_mass_density;
pub mod pr_molar_volume;
pub mod pr_z_factor;
pub mod prsv_kappa;
pub mod pt_flash;
pub mod pure_saturation;
pub mod rachford_rice_binary;
pub mod results;
pub mod spec_gen;
pub mod stability_test;
pub mod vdw1f_mix_binary;

pub use bubble_pressure::bubble_pressure;
pub use critical_point::critical_point;
pub use dew_pressure::dew_pressure;
pub use ideal_gas_cp::{REFERENCE_TEMPERATURE, ideal_gas_cp};
pub use mixture::{Component, Mixture, PhaseState, ReducedParameters, RootSide};
pub use molar_enthalpy_entropy::{IdealGasModel, molar_enthalpy_entropy};
pub use ph_flash::ph_flash;
pub use phase_boundary::{Incipient, PhaseBoundary, phase_boundary_pressure};
pub use pr_alpha_ab::{OMEGA_A, OMEGA_B, pr_alpha_ab};
pub use pr_departure::pr_departure;
pub use pr_kappa::pr_kappa;
pub use pr_mass_density::pr_mass_density;
pub use pr_molar_volume::{MOLAR_GAS_CONSTANT, pr_molar_volume};
pub use pr_z_factor::pr_z_factor;
pub use prsv_kappa::prsv_kappa;
pub use pt_flash::pt_flash;
pub use pure_saturation::pure_saturation;
pub use rachford_rice_binary::rachford_rice_binary;
pub use results::{
    BubblePressureResult, CriticalPointResult, DewPressureResult, IdealGasCpResult,
    MolarEnthalpyEntropyResult, Phase, PrAlphaAbResult, PrDepartureResult, PrKappaResult,
    PrMassDensityResult, PrMolarVolumeResult, PrZFactorResult, PrsvKappaResult, PtFlashResult,
    PureSaturationResult, RachfordRiceBinaryResult, RootStructure, StabilityTestResult,
    StabilityVerdict, Vdw1fMixBinaryResult,
};
pub use stability_test::stability_test;
pub use vdw1f_mix_binary::vdw1f_mix_binary;
