//! The mixing rule that combines per-component attraction and repulsion parameters
//! into a mixture's.
//!
//! NeqSim's `EosMixingRuleHandler` is a strategy object with one class per rule. This
//! mirrors that seam: each rule is a variant here, and the temperature dependence of
//! the binary interaction parameters is resolved once per state by
//! [`MixingRule::effective_kij`] rather than per iteration.

use crate::databank::UnifacUmrpruParameters;

/// NeqSim's `hwfc` for the `UNIFAC_UMRPRU` GE model: `EosMixingRuleHandler.init` sets
/// `-1/0.53` for that model and the cubic's own constant for every other.
///
/// **A property of the pairing, not of the cubic.** The same handler computes the
/// Huron-Vidal rule's coefficient from the cubic's `delta` parameters, so a UMR rule
/// written with the cubic's `hv_constant()` would be out by a factor that no pure fluid
/// and no single-component comparison can show.
pub const UMR_HWFC: f64 = -1.0 / 0.53;

/// The UMR attraction coefficients `qPure_i + hwfc ln gamma_i`.
///
/// `EosMixingRuleHandler.init`'s one line, and the whole of what the rule adds to the
/// cubic: the mixture's `alpha_mix` is `sum_i x_i` of this. `qpure[i]` is NeqSim's
/// `qPure[i] = a_i^T/(b_i R T)`.
#[must_use]
pub fn umr_ader(qpure: &[f64], ln_gamma: &[f64]) -> Vec<f64> {
    qpure
        .iter()
        .zip(ln_gamma)
        .map(|(&q, &lg)| q + UMR_HWFC * lg)
        .collect()
}

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
    /// classic co-volume sum, and with the temperature coefficient read from the pairs'
    /// own `WSGIJT`/`WSGJIT` rather than from the Huron-Vidal ones. NeqSim builds a
    /// `PhaseGENRTLmodifiedHV` from `NRTLDijT` for this rule, so the two differ in which
    /// column feeds `DijT` as well as in the mixing; `Mixture::ws_ader` carries the
    /// measurement that settles which.
    WongSandler {
        /// The interaction matrix the rule's `b_mix` and the NRTL's classic pairs read.
        kij: Vec<f64>,
        /// The fitted NRTL energy `Dij`, in Kelvin.
        hv_gij: Vec<f64>,
        /// The fitted temperature coefficient `DijT`, in Kelvin per Kelvin.
        hv_gij_t: Vec<f64>,
        /// The fitted non-randomness `alpha`.
        hv_alpha: Vec<f64>,
        /// One flag per interaction: `true` for the fitted NRTL pair, `false` for the
        /// cubic's own excess energy.
        hv_pairs: Vec<bool>,
    },
    /// The UMR universal mixing rule, NeqSim's `SRKHuronVidal2` driven by
    /// `PhaseGEUnifacUMRPRU`.
    ///
    /// The same shape as [`MixingRule::HuronVidal`] - an attraction coefficient per
    /// component that the mixture averages into `A = n B R T alpha_mix`, with the
    /// co-volume left at the linear sum - but the excess Gibbs energy is UNIFAC's rather
    /// than NRTL's, and its coefficient is NeqSim's `hwfc = -1/0.53` rather than the
    /// cubic's own. `EosMixingRuleHandler.init` writes
    /// `alpha_mix = sum_i x_i (a_i^T/(b_i R T) + hwfc ln gamma_i)` for exactly the
    /// `UNIFAC_UMRPRU` GE model, and `calcA` is `n B R T alpha_mix`.
    ///
    /// **The UNIFAC parameters travel with the rule rather than being looked up here.**
    /// The rule is built once per mixture and the interaction matrix is re-evaluated per
    /// state; resolving the tables inside would put a databank read on the flash's inner
    /// loop, and the rule would no longer be a value a caller can construct and inspect.
    Umr {
        /// The interaction matrix the cubic's classic pairs read. NeqSim's UMR-CPA sets
        /// the UNIFAC groups and reads no `kij` column, so this is zero for it; it is
        /// here because the rule is still an equation-of-state mixing rule and a caller
        /// may pair it with one.
        kij: Vec<f64>,
        /// The UNIFAC-UMR-PRU group basis and interaction tables.
        unifac: UnifacUmrpruParameters,
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
            | MixingRule::WongSandler { kij, .. }
            | MixingRule::Umr { kij, .. } => kij,
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
            | MixingRule::WongSandler { kij, .. }
            | MixingRule::Umr { kij, .. } => kij.clone(),
        }
    }

    /// `T * d kij / dT`, one entry per interaction, at a temperature.
    ///
    /// Needed wherever the fugacity's temperature derivative is, because `A_ij` carries
    /// the interaction parameter and a rule whose `kij` moves with temperature moves it.
    /// The rules that do not are zero here rather than absent: a caller differentiating
    /// `A_ij` adds this term unconditionally, and a `match` that omitted it for three of
    /// the five rules is the shape of mistake this accessor exists to prevent.
    #[must_use]
    pub fn t_d_effective_kij(&self, t_kelvin: f64) -> Vec<f64> {
        match self {
            MixingRule::ClassicT {
                kij,
                kij_t,
                inverse_temperature,
            } => kij
                .iter()
                .zip(kij_t)
                .map(|(_k0, kt)| {
                    if *inverse_temperature {
                        -kt / t_kelvin
                    } else {
                        kt * t_kelvin / 273.15
                    }
                })
                .collect(),
            MixingRule::ClassicT2 {
                kij_t,
                inverse_temperature,
                ..
            } => kij_t
                .iter()
                .zip(inverse_temperature)
                .map(|(kt, inverse)| {
                    if *inverse {
                        -kt / t_kelvin
                    } else {
                        kt * t_kelvin
                    }
                })
                .collect(),
            _ => vec![0.0; self.interaction_count()],
        }
    }

    /// How many interactions the rule's matrix has, diagonal included.
    #[must_use]
    pub fn interaction_count(&self) -> usize {
        match self {
            MixingRule::Classic { kij }
            | MixingRule::ClassicT { kij, .. }
            | MixingRule::ClassicT2 { kij, .. }
            | MixingRule::SoreideWhitson { kij, .. }
            | MixingRule::HuronVidal { kij, .. }
            | MixingRule::WongSandler { kij, .. }
            | MixingRule::Umr { kij, .. } => kij.len(),
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
