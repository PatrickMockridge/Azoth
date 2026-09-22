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
//! * [`equilibrium_constant`] - `ln K`, its derivative and the heat of reaction
//! * [`reference_potentials`] - the independent basis, and the potentials from it
//! * [`chemical_equilibrium`] - the Smith-Missen reactive solve
//! * [`linalg`] - the exact rank and the LU solve those two rest on

pub mod chemical_equilibrium;
pub mod databank;
pub mod equilibrium_constant;
pub mod linalg;
pub mod model_gen;
pub mod reference_potentials;
pub mod spec_gen;

pub use chemical_equilibrium::{ChemicalEquilibriumResult, chemical_equilibrium};
pub use equilibrium_constant::{EquilibriumConstantResult, equilibrium_constant};
pub use reference_potentials::{ReferencePotentialsResult, reference_potentials};
