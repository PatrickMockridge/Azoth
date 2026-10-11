//! Petroleum-fraction characterisation.
//!
//! The front door for a fluid that has no databank row. A crude is an *assay* — a list of
//! cuts, each a molar mass and a density — and everything a cubic needs about a cut follows
//! from those two numbers by correlation. NeqSim's `thermo/characterization/` is that
//! subsystem, and this crate ports it.
//!
//! It is a namespace of its own rather than a family under `eos`, mirroring NeqSim's own
//! package split: what is here is not an equation of state and is not reached from one. A cut
//! is characterised first and *then* handed to a mixture, which is the caller's step.
//!
//! * [`tbp_cut_properties`] — a cut's critical properties, by any of NeqSim's ten models

pub mod assay_mass_fractions;
pub mod characterise_plus_fraction;
pub mod characterize_to_reference;
pub mod lumping;
pub mod model_gen;
pub mod oil_assay_derivation;
pub mod pedersen_plus_split;
pub mod results;
pub mod set_heavy_tbp_to_plus;
pub mod spec_gen;
pub mod tbp_closure;
pub mod tbp_cut_boundaries;
pub mod tbp_cut_properties;
pub mod tbp_density;
pub mod tbp_grouping;
pub mod unported;
pub mod unported_gen;
pub mod whitson_gamma_split;

pub use assay_mass_fractions::{AssayBasis, AssayCut, assay_mass_fractions};
pub use characterise_plus_fraction::{PlusModel, characterise_plus_fraction};
pub use characterize_to_reference::{
    ReferenceRow, SourceRow, characterize_to_reference,
};
pub use lumping::lumping;
pub use oil_assay_derivation::{AssayCutDerivation, oil_assay_derivation};
pub use pedersen_plus_split::pedersen_plus_split;
pub use results::{
    AssayMassFractionsResult, CharacterisePlusFractionResult, CharacterizeToReferenceResult,
    LumpingResult, OilAssayDerivationResult, PedersenPlusSplitResult, SetHeavyTbpToPlusResult,
    TbpClosureResult, TbpCutBoundariesResult, TbpCutPropertiesResult, TbpDensityResult,
    TbpGroupingResult,
    WhitsonGammaSplitResult,
};
pub use set_heavy_tbp_to_plus::set_heavy_tbp_to_plus;
pub use tbp_closure::{TbpClosureKind, tbp_closure};
pub use tbp_cut_boundaries::tbp_cut_boundaries;
pub use tbp_cut_properties::{TbpModel, tbp_cut_properties};
pub use tbp_density::tbp_density;
pub use tbp_grouping::tbp_grouping;
pub use whitson_gamma_split::{WhitsonDensityModel, gamma, p0_p1, whitson_gamma_split};
