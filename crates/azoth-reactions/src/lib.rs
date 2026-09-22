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

pub mod databank;
