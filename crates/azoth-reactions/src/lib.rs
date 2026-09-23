//! Chemical reaction equilibrium.
//!
//! NeqSim's `chemicalreactions/` builds an element matrix from the substances a
//! fluid carries, reduces the reaction set to a linearly independent basis,
//! solves for each component's reference potential from `sum(nu_i mu_i) = -RT
//! ln K`, and minimises the Gibbs energy subject to the element balances. This
//! crate is that chain, one id per step, with the reference potentials supplied
//! by the caller rather than read from a live system.
//!
//! * [`databank`] - the element table, the stoichiometry, and the reaction rows
//! * [`diis`] - the Pulay accelerator the multiphase RAND solve drives
//! * [`equilibrium_constant`] - `ln K`, its derivative and the heat of reaction
//! * [`reference_potentials`] - the independent basis, and the potentials from it
//! * [`chemical_equilibrium`] - the Smith-Missen reactive solve
//! * [`reactive_phase`] - which phase a fluid solves its reactions in
//! * [`reactive_phase_equilibrium`] - that solve as an operation on one phase
//! * [`linalg`] - the exact rank and the LU solve those two rest on
//! * [`reactive_tp_flash`] - the reactive flash stack as a model, over `azoth-eos`

pub mod chemical_equilibrium;
pub mod databank;
pub mod diis;
pub mod equilibrium_constant;
pub mod formula_matrix;
pub mod linalg;
pub mod lp_seed;
pub mod model_gen;
pub mod rand_solver;
pub mod reactive_flash;
pub mod reactive_phase;
pub mod reactive_phase_equilibrium;
pub mod reactive_stability;
pub mod reactive_tp_flash;
pub mod reference_potentials;
pub mod spec_gen;

pub use chemical_equilibrium::{ChemicalEquilibriumResult, chemical_equilibrium};
pub use equilibrium_constant::{EquilibriumConstantResult, equilibrium_constant};
pub use reactive_phase::{is_reactive_phase, reactive_phase_index};
pub use reactive_phase_equilibrium::{ReactivePhaseEquilibriumResult, reactive_phase_equilibrium};
pub use reference_potentials::{ReferencePotentialsResult, reference_potentials};
