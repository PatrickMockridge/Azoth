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
    /// Per-pair temperature-dependent `kij`, NeqSim's `ClassicSRKT2`.
    ///
    /// Each interaction picks its own form - `kij(T) = kij + kij_t * T` where its flag
    /// is false, `kij(T) = kij + kij_t / T` where it is true. The linear form is plain
    /// `T`, not [`MixingRule::ClassicT`]'s `T / 273.15 - 1`.
    ClassicT2 {
        /// The constant part, as in [`MixingRule::Classic`].
        kij: Vec<f64>,
        /// The temperature coefficient, one entry per interaction.
        kij_t: Vec<f64>,
        /// One flag per interaction: `true` for the `kij_t / T` form, `false` for
        /// `kij_t * T`.
        inverse_temperature: Vec<bool>,
    },
    /// The salinity-dependent Soreide-Whitson aqueous rule, NeqSim's
    /// `WhitsonSoreideMixingRule`.
    ///
    /// For the water-rich phase the water-gas interactions are replaced by the salinity
    /// correlations; every other pair keeps the base `kij`. The rule is *asymmetric* -
    /// the correlation applies only where the second component is water, exactly as
    /// NeqSim writes it - and phase-dependent, so it is resolved by
    /// [`MixingRule::phase_kij`] rather than [`MixingRule::effective_kij`].
    SoreideWhitson {
        /// The base interaction matrix, as in [`MixingRule::Classic`].
        kij: Vec<f64>,
        /// One role per component, in component order.
        roles: Vec<SoreideWhitsonRole>,
        /// Equivalent NaCl molality, in mol/kg water.
        salinity: f64,
    },
    /// The Huron-Vidal rule, NeqSim's `SRKHuronVidal2`.
    ///
    /// The attraction parameter mixes the pure components' `a_i/b_i` with the excess
    /// Gibbs energy of the co-volume-weighted NRTL in [`crate::hv_ge`], so this rule
    /// resolves its `a_mix` and fugacity in the mixture layer rather than through a
    /// `kij` matrix. `kij` is the interaction matrix the NRTL's non-`hv_pairs` read.
    HuronVidal {
        /// The base interaction matrix, as in [`MixingRule::Classic`].
        kij: Vec<f64>,
        /// The fitted NRTL energy `Dij`, in Kelvin.
        hv_gij: Vec<f64>,
        /// The temperature coefficient `DijT`.
        hv_gij_t: Vec<f64>,
        /// The fitted non-randomness `alpha`.
        hv_alpha: Vec<f64>,
        /// One flag per interaction: `true` for the fitted NRTL pair, `false` for the
        /// cubic's own excess energy.
        hv_pairs: Vec<bool>,
    },
    /// The Wong-Sandler rule, NeqSim's `WongSandlerMixingRule`.
    ///
    /// Like [`MixingRule::HuronVidal`] but with a GE-dependent `b_mix` instead of the
    /// classic co-volume sum. NeqSim reads the cached, DijT-free activity coefficients
    /// for this rule, so it resolves separately from [`MixingRule::HuronVidal`].
    WongSandler {
        /// The interaction matrix the rule's `b_mix` and the NRTL's classic pairs read.
        kij: Vec<f64>,
        /// The fitted NRTL energy `Dij`, in Kelvin.
        hv_gij: Vec<f64>,
        /// The fitted non-randomness `alpha`.
        hv_alpha: Vec<f64>,
        /// One flag per interaction: `true` for the fitted NRTL pair, `false` for the
        /// cubic's own excess energy.
        hv_pairs: Vec<bool>,
    },
}

/// A component's role in the Soreide-Whitson aqueous correlation.
///
/// Resolved from the component name at construction - azoth's `Component` carries no
/// name, so the caller that named it also assigns the role.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SoreideWhitsonRole {
    /// Water, `H2O`.
    Water,
    /// Nitrogen, `N2`.
    Nitrogen,
    /// Carbon dioxide, `CO2`.
    CarbonDioxide,
    /// A hydrocarbon, or any other gas NeqSim's `else` branch routes here.
    Hydrocarbon,
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
            MixingRule::Classic { kij }
            | MixingRule::ClassicT { kij, .. }
            | MixingRule::ClassicT2 { kij, .. }
            | MixingRule::SoreideWhitson { kij, .. }
            | MixingRule::HuronVidal { kij, .. }
            | MixingRule::WongSandler { kij, .. } => kij,
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
            MixingRule::ClassicT2 {
                kij,
                kij_t,
                inverse_temperature,
            } => kij
                .iter()
                .zip(kij_t)
                .zip(inverse_temperature)
                .map(|((k0, kt), inverse)| {
                    if *inverse {
                        k0 + kt / t_kelvin
                    } else {
                        k0 + kt * t_kelvin
                    }
                })
                .collect(),
            MixingRule::SoreideWhitson { kij, .. }
            | MixingRule::HuronVidal { kij, .. }
            | MixingRule::WongSandler { kij, .. } => kij.clone(),
        }
    }

    /// The interaction matrix at a phase's composition, flattened row-major.
    ///
    /// The phase-independent rules resolve their matrix once per state in
    /// [`MixingRule::effective_kij`] and this simply returns it. The Soreide-Whitson
    /// rule is phase-dependent - its salinity correlation applies only to the
    /// water-rich phase - so it is computed here, where the composition is known.
    ///
    /// `base` is the state's resolved matrix, `reduced_temperatures` the per-component
    /// `T/Tc`, and `acentric_factors` the per-component `omega`.
    #[must_use]
    pub fn phase_kij(
        &self,
        base: &[f64],
        reduced_temperatures: &[f64],
        acentric_factors: &[f64],
        x: &[f64],
    ) -> Vec<f64> {
        let MixingRule::SoreideWhitson {
            roles, salinity, ..
        } = self
        else {
            return base.to_vec();
        };
        let n = x.len();
        // NeqSim's `water.x > 0.8` gate: the aqueous correlation applies only to the
        // water-rich phase; every other phase keeps the base matrix.
        let is_aqueous = roles
            .iter()
            .zip(x)
            .any(|(&role, &xi)| role == SoreideWhitsonRole::Water && xi > 0.8);
        if !is_aqueous {
            return base.to_vec();
        }
        let s = *salinity;
        (0..n)
            .flat_map(|i| {
                (0..n).map(move |j| {
                    if roles[j] != SoreideWhitsonRole::Water {
                        return base[i * n + j];
                    }
                    let tr = reduced_temperatures[i];
                    match roles[i] {
                        SoreideWhitsonRole::Nitrogen => {
                            0.997
                                * (-1.70235 * (1.0 + 0.025587 * s.powf(0.75))
                                    + 0.44338 * (1.0 + 0.08126 * s.powf(0.75)) * tr)
                        }
                        SoreideWhitsonRole::CarbonDioxide => {
                            // NeqSim's ladder has a 0.8 branch above 3.5 mol/kg, but the
                            // 0.9 branch above 2.0 fires first, so it is unreachable.
                            let multip_k = if s > 2.0 { 0.9 } else { 1.0 };
                            multip_k
                                * 0.989
                                * (-0.31092 * (1.0 + 0.15587 * s.powf(0.75))
                                    + 0.2358 * (1.0 + 0.17837 * s.powf(0.98)) * tr
                                    - 21.2566 * (-6.7222_f64.powf(tr) - s).exp())
                        }
                        SoreideWhitsonRole::Water => 0.0,
                        SoreideWhitsonRole::Hydrocarbon => {
                            let c0 = 0.017407;
                            let c1 = 0.033516;
                            let c2 = 0.011478;
                            let a0 = 1.112 - 1.7369 * acentric_factors[i].powf(-0.1);
                            let a1 = 1.1001 + 0.83 * acentric_factors[i];
                            let a2 = -0.15742 - 1.0988 * acentric_factors[i];
                            0.777
                                * ((1.0 + c0 * s) * a0
                                    + (1.0 + c1 * s) * a1 * tr
                                    + (1.0 + c2 * s) * a2 * tr * tr)
                        }
                    }
                })
            })
            .collect()
    }
}
