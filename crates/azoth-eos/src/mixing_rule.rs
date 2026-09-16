//! The mixing rule that combines per-component attraction and repulsion parameters
//! into a mixture's.
//!
//! NeqSim's `EosMixingRuleHandler` is a strategy object with one class per rule. This
//! mirrors that seam: each rule is a variant here, and the temperature dependence of
//! the binary interaction parameters is resolved once per state by
//! [`MixingRule::effective_kij`] rather than per iteration.

/// How the binary interaction parameters enter the mixture.
#[derive(Debug, Clone, PartialEq)]
pub enum MixingRule {
    /// Constant `kij` - the classic van der Waals one-fluid rule, and NeqSim's
    /// `ClassicSRK`.
    Classic {
        /// The symmetric `N x N` interaction matrix, flattened row-major, zero diagonal.
        kij: Vec<f64>,
    },
    /// Temperature-dependent `kij`, NeqSim's `ClassicSRKT`.
    ///
    /// `kij(T) = kij + kij_t / T` when `inverse_temperature`, else
    /// `kij(T) = kij + kij_t * (T / 273.15 - 1)`. The reference temperature of 273.15 K
    /// is NeqSim's, copied verbatim.
    ClassicT {
        /// The constant part, as in [`MixingRule::Classic`].
        kij: Vec<f64>,
        /// The temperature coefficient, one entry per interaction.
        kij_t: Vec<f64>,
        /// `true` for the `kij_t / T` form; `false` for the linear-in-`T` form.
        inverse_temperature: bool,
    },
}

impl MixingRule {
    /// The constant part of the interaction parameter between `i` and `j`.
    ///
    /// The base `kij`, before any temperature correction. Exposed for callers such as
    /// the critical-point construction that evaluate at a fixed temperature and need
    /// the matrix the mixture was built from.
    #[must_use]
    pub fn base_kij(&self, i: usize, j: usize, n: usize) -> f64 {
        let kij = match self {
            MixingRule::Classic { kij } | MixingRule::ClassicT { kij, .. } => kij,
        };
        kij[i * n + j]
    }

    /// The interaction matrix at an absolute temperature, flattened row-major.
    ///
    /// For [`MixingRule::Classic`] this is the matrix itself; for
    /// [`MixingRule::ClassicT`] the temperature correction is applied. Resolved once
    /// per `reduced_parameters` call, so the mixing rule's temperature dependence does
    /// not re-evaluate on every flash iteration.
    #[must_use]
    pub fn effective_kij(&self, t_kelvin: f64) -> Vec<f64> {
        match self {
            MixingRule::Classic { kij } => kij.clone(),
            MixingRule::ClassicT {
                kij,
                kij_t,
                inverse_temperature,
            } => kij
                .iter()
                .zip(kij_t)
                .map(|(k0, kt)| {
                    if *inverse_temperature {
                        k0 + kt / t_kelvin
                    } else {
                        k0 + kt * (t_kelvin / 273.15 - 1.0)
                    }
                })
                .collect(),
        }
    }
}
