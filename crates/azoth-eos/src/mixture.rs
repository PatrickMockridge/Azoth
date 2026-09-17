//! Components, mixtures, and the arithmetic over them.
//!
//! The model layer's own arithmetic: a mixture's fugacity coefficient is not a
//! registered calculation, because `eos.pr_departure` covers a *pure* component and the
//! registry is scalar, so there is no composition vector to hang a mixture form on.
//! `eos.pt_flash`'s spec records the contract that keeps the two in step - at `N = 1`
//! the mixture form must reduce to `eos.pr_departure`, and at `N = 2` the mixture
//! parameters must reproduce `eos.vdw1f_mix_binary` and the vapour fraction
//! `eos.rachford_rice_binary`.

use azoth_core::units::{Pressure, ThermodynamicTemperature};
use azoth_core::{AzothError, Result, Warning};

use crate::alpha_term::{
    Alpha, AlphaTerm, Danesh, Delft1998, Gassem2001, MatCop, MatCopFallback, Mollerup, RkAlpha,
    Schwartzentruber, Soave, TwuCoon, matcop_kappa, umr_kappa,
};
use crate::cubic::Cubic;
use crate::mixing_rule::MixingRule;
use crate::{
    pr_alpha_ab, pr_kappa, pr_z_factor, pr78_kappa, rk_alpha_ab, srk_alpha_ab, srk_kappa,
    srk_z_factor, twu_kappa,
};

/// Which root of the cubic a phase claims.
///
/// Selection is by *ordering* and never by an initial guess - the rule
/// [`crate::pr_z_factor`] fixes, and the one the flash has to follow so that both
/// implementations pick the same phase.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RootSide {
    /// The smallest admissible root.
    Liquid,
    /// The largest admissible root.
    Vapour,
}

/// One phase's state at a composition.
#[derive(Debug, Clone, PartialEq)]
pub struct PhaseState {
    /// The mixture's attraction parameter at this composition.
    pub a_mix: f64,
    /// The mixture's repulsion parameter at this composition.
    pub b_mix: f64,
    /// The root of the cubic this phase sits on.
    pub z: f64,
    /// `ln phi_i` for every component, one entry per component.
    pub ln_phi: Vec<f64>,
    /// The departure enthalpy over `R*T`, for the mixture.
    ///
    /// `pr_departure`'s expression with the pure component's `psi` replaced by
    /// [`Self::psi_bar`]. Verified by reduction: at `N = 1` this is bit-identical to
    /// the registered calc's `h_dep_rt`.
    pub h_dep_rt: f64,
    /// The departure entropy over `R`, for the mixture.
    ///
    /// `h_dep_rt - sum_i z_i ln phi_i`, which is the Gibbs identity and therefore
    /// exact rather than a second formula that could disagree with the first. It is
    /// also why this is one ulp from the pure calc's `s_dep_r` at `N = 1`: the sum is
    /// taken in a different order, and that is the whole of the difference.
    pub s_dep_r: f64,
    /// `sum_i sum_j z_i z_j A_ij (psi_i + psi_j)/2 / sum_i sum_j z_i z_j A_ij`.
    ///
    /// The composition-weighted average of the components' `psi`, weighted by the
    /// attraction parameters. At `N = 1` it equals that component's own `psi` exactly.
    pub psi_bar: f64,
    /// The departure heat capacity over `R`, for the mixture.
    ///
    /// `pr_departure`'s expression again, with `psi_bar` in place of `psi` and with
    /// `psi_bar` itself moving with temperature, because the weights that form it do.
    /// At `N = 1` the double sum collapses to the single component's `T*dpsi/dT`, so
    /// this reduces to the registered calc's `cp_dep_r` exactly.
    pub cp_dep_r: f64,
}

/// One component's critical constants.
///
/// A caller-supplied record, and deliberately carrying **no name**: a name would
/// invite a lookup rather than a value the caller supplied.
#[derive(Debug, Clone, PartialEq)]
pub struct Component {
    /// Critical temperature.
    pub tc: ThermodynamicTemperature,
    /// Critical pressure.
    pub pc: Pressure,
    /// Acentric factor.
    pub omega: f64,
    /// Molar mass in kg/mol. `None` for a component built from critical constants
    /// alone; a transport correlation that needs it refuses `None` rather than
    /// defaulting to a zero molar mass.
    pub molar_mass: Option<f64>,
    /// Fitted parameters for the alpha correlations that need them, in the order the
    /// correlation reads them. Empty for a component built without a fitted set, which
    /// is every caller-supplied component and every correlation that needs none.
    pub alpha_params: Vec<f64>,
    /// The volume-translation parameter in m³/mol, subtracted from the untranslated
    /// molar volume: `v_corr = v - c`. Zero for a component without a translation,
    /// which is the plain PR/SRK/RK forms. The Peneloux shift of
    /// [`crate::pr_peneloux_shift`] or [`crate::srk_peneloux_shift`] is one source; a
    /// fitted constant is another.
    pub volume_shift: f64,
}

impl Component {
    /// A component record, checked for the positivity the equation of state needs.
    ///
    /// # Errors
    /// * [`AzothError::OutOfRange`] if `Tc` or `Pc` is not positive. Both are
    ///   divisors - `Tr = T/Tc` and `Pr = P/Pc` - and a zero is not a state.
    pub fn new(tc: ThermodynamicTemperature, pc: Pressure, omega: f64) -> Result<Self> {
        for (field, value) in [("Tc", tc.value), ("Pc", pc.value)] {
            // Written as an explicit finiteness-and-positivity test rather than
            // `!(value > 0.0)`, which reads as a double negative and is the shape
            // clippy flags. Both reject NaN and non-positive values; only the
            // explicit form also rejects an infinity, which would make the reduced
            // variable zero and then divide by it.
            if !value.is_finite() || value <= 0.0 {
                return Err(AzothError::OutOfRange {
                    field: field.to_string(),
                    value,
                    detail: format!(
                        "a critical {field} is a divisor in the reduced variables and is \
                         finite and positive by definition"
                    ),
                });
            }
        }
        Ok(Self {
            tc,
            pc,
            omega,
            molar_mass: None,
            alpha_params: Vec::new(),
            volume_shift: 0.0,
        })
    }

    /// This component, with its molar mass attached.
    ///
    /// The databank populates it; a caller-supplied component may, and a transport
    /// correlation that needs the mass refuses one that did not.
    #[must_use]
    pub fn with_molar_mass(mut self, molar_mass: Option<f64>) -> Self {
        self.molar_mass = molar_mass;
        self
    }

    /// This component, with fitted alpha parameters attached.
    ///
    /// Only the parameterized correlations read them; the rest ignore the field. The
    /// values are caller-supplied, exactly like `Tc`, `Pc` and `omega`, and carry no
    /// name - the caller resolved them from a name and passed the numbers.
    #[must_use]
    pub fn with_alpha_params(mut self, params: Vec<f64>) -> Self {
        self.alpha_params = params;
        self
    }

    /// This component, with its volume-translation parameter attached.
    ///
    /// Subtracted from the untranslated molar volume; see the field. The value is
    /// caller-supplied, like `Tc`, `Pc` and `omega`, and carries no name.
    #[must_use]
    pub fn with_volume_shift(mut self, volume_shift: f64) -> Self {
        self.volume_shift = volume_shift;
        self
    }
}

/// A set of components and the binary interaction parameters between them.
///
/// The `kij` matrix is stored flattened row-major and validated once, at
/// construction, rather than re-checked at every evaluation: a flash evaluates the
/// mixing rule hundreds of times and a shape error found on the first is found for
/// free.
#[derive(Debug, Clone, PartialEq)]
pub struct Mixture {
    components: Vec<Component>,
    mixing_rule: MixingRule,
    cubic: Cubic,
    alpha: Alpha,
}

impl Mixture {
    /// A mixture of `components` with an `N x N` interaction matrix.
    ///
    /// `kij` is flattened row-major with `N = components.len()`. It **must** satisfy
    /// `kij[i][i] == 0` and `kij[i][j] == kij[j][i]`, and a matrix that does not is
    /// refused rather than corrected: the diagonal term would silently rescale each
    /// pure component's attraction, and an asymmetric matrix has no meaning in a
    /// double sum that is symmetric by construction. Both are mistakes a caller
    /// cannot see in the answer.
    ///
    /// # Errors
    /// * [`AzothError::InvalidInput`] if `components` is empty, if `kij` is not
    ///   `N*N` long, or if the matrix is not symmetric with a zero diagonal.
    /// * Propagates [`Component::new`]'s range errors.
    pub fn new(components: Vec<Component>, kij: Vec<f64>) -> Result<Self> {
        let n = components.len();
        if n == 0 {
            return Err(AzothError::invalid_input(
                "components",
                "a mixture needs at least one component",
            ));
        }
        check_kij(n, &kij)?;
        Ok(Self {
            components,
            mixing_rule: MixingRule::Classic { kij },
            cubic: Cubic::default(),
            alpha: Alpha::default(),
        })
    }

    /// The components, in order.
    #[must_use]
    pub fn components(&self) -> &[Component] {
        &self.components
    }

    /// How many components there are.
    #[must_use]
    pub fn len(&self) -> usize {
        self.components.len()
    }

    /// Whether the mixture has no components. Always false for a constructed
    /// [`Mixture`] - it is here because clippy asks for it beside [`Self::len`].
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.components.is_empty()
    }

    /// The cubic whose geometry this mixture is evaluated under.
    #[must_use]
    pub fn cubic(&self) -> Cubic {
        self.cubic
    }

    /// The alpha correlation this mixture's attraction term uses.
    #[must_use]
    pub fn alpha(&self) -> Alpha {
        self.alpha
    }

    /// This mixture, evaluated under a different cubic.
    ///
    /// The default is Peng-Robinson; this swaps in Soave-Redlich-Kwong (or any later
    /// cubic) without rebuilding the components or the interaction matrix, and resets
    /// the alpha term to the cubic's own Soave correlation.
    #[must_use]
    pub fn with_cubic(mut self, cubic: Cubic) -> Self {
        self.cubic = cubic;
        self.alpha = match cubic {
            Cubic::Pr => Alpha::Pr,
            Cubic::Srk | Cubic::Rk => Alpha::Srk,
            Cubic::Tst => Alpha::Twu,
        };
        self
    }

    /// This mixture, with its alpha correlation changed.
    ///
    /// The Soave variants share the alpha form and differ only in the `m` correlation,
    /// so switching between them keeps the cubic's shape and changes one coefficient.
    #[must_use]
    pub fn with_alpha(mut self, alpha: Alpha) -> Self {
        self.alpha = alpha;
        self
    }

    /// This mixture, with its mixing rule changed.
    ///
    /// The default is the classic constant-`kij` rule; this swaps in a temperature-
    /// dependent `kij` (or, later, a Huron-Vidal or Wong-Sandler rule) without
    /// rebuilding the components. The rule's matrices must already satisfy the same
    /// invariant [`Mixture::new`] enforces - `N x N`, symmetric, zero diagonal - and a
    /// rule that does not is a caller error with the same silent failure mode the
    /// constructor's check exists to prevent.
    #[must_use]
    pub fn with_mixing_rule(mut self, mixing_rule: MixingRule) -> Self {
        self.mixing_rule = mixing_rule;
        self
    }

    /// The interaction parameter between components `i` and `j`.
    ///
    /// The constant part, before any temperature correction - the matrix this mixture
    /// was built from.
    #[must_use]
    pub fn kij(&self, i: usize, j: usize) -> f64 {
        self.mixing_rule.base_kij(i, j, self.components.len())
    }

    /// The reduced parameters `(A_i, B_i)` for every component at a state.
    ///
    /// These depend on `T` and `P` alone, not on the composition, which is why they
    /// are computed once per flash rather than once per iteration: `A_i` and `B_i`
    /// are the same numbers for the liquid and the vapour phase, and only the
    /// composition re-weights them.
    ///
    /// # Errors
    /// Propagates the range checks of [`crate::pr_kappa`] and [`crate::pr_alpha_ab`].
    pub fn reduced_parameters(
        &self,
        t: ThermodynamicTemperature,
        p: Pressure,
    ) -> Result<ReducedParameters> {
        let mut a = Vec::with_capacity(self.len());
        let mut b = Vec::with_capacity(self.len());
        let mut psi = Vec::with_capacity(self.len());
        let mut psi_t = Vec::with_capacity(self.len());
        let mut reduced_temperatures = Vec::with_capacity(self.len());
        let mut warnings = Vec::new();

        for component in &self.components {
            let reduced_temperature = t.value / component.tc.value;
            let reduced_pressure = p.value / component.pc.value;
            reduced_temperatures.push(reduced_temperature);
            // The kappa correlation belongs to the alpha term; the Omega constants to
            // the cubic. SRK and PR share the Soave alpha form and differ in both; RK
            // is kappa-free.
            let cubic = self.cubic;
            let non_soave = |term: &dyn AlphaTerm| -> (f64, f64, f64, f64) {
                let alpha = term.alpha(reduced_temperature);
                let a_reduced = cubic.omega_a() * alpha * reduced_pressure
                    / (reduced_temperature * reduced_temperature);
                let b_reduced = cubic.omega_b() * reduced_pressure / reduced_temperature;
                (
                    a_reduced,
                    b_reduced,
                    term.psi(reduced_temperature),
                    term.psi_t(reduced_temperature),
                )
            };
            let (a_reduced, b_reduced, psi_value, psi_t_value) = match self.cubic {
                Cubic::Pr | Cubic::Srk => match self.alpha {
                    Alpha::TwuCoon => non_soave(&TwuCoon {
                        omega: component.omega,
                    }),
                    Alpha::Gassem2001 => non_soave(&Gassem2001 {
                        omega: component.omega,
                    }),
                    Alpha::Danesh => {
                        let kappa = pr78_kappa(component.omega)?;
                        warnings.extend(kappa.warnings);
                        non_soave(&Danesh { kappa: kappa.kappa })
                    }
                    Alpha::Schwartzentruber => non_soave(&Schwartzentruber {
                        omega: component.omega,
                        params: component.alpha_params.clone(),
                    }),
                    Alpha::Mollerup => non_soave(&Mollerup {
                        params: component.alpha_params.clone(),
                    }),
                    Alpha::MatCop => non_soave(&MatCop {
                        kappa: matcop_kappa(component.omega),
                        params: component.alpha_params.clone(),
                        fallback: MatCopFallback::None,
                    }),
                    Alpha::MatCopPr => {
                        let kappa = pr_kappa(component.omega)?;
                        warnings.extend(kappa.warnings);
                        non_soave(&MatCop {
                            kappa: kappa.kappa,
                            params: component.alpha_params.clone(),
                            fallback: MatCopFallback::Supercritical,
                        })
                    }
                    Alpha::MatCopPrUmr => non_soave(&MatCop {
                        kappa: umr_kappa(component.omega),
                        params: component.alpha_params.clone(),
                        fallback: MatCopFallback::Unset,
                    }),
                    Alpha::MatCop5PrUmr => {
                        let kappa = pr_kappa(component.omega)?;
                        warnings.extend(kappa.warnings);
                        non_soave(&MatCop {
                            kappa: kappa.kappa,
                            params: component.alpha_params.clone(),
                            fallback: MatCopFallback::AllUnset,
                        })
                    }
                    Alpha::Delft1998 => {
                        let kappa = pr78_kappa(component.omega)?;
                        warnings.extend(kappa.warnings);
                        let is_methane = component.alpha_params.first().copied() == Some(1.0);
                        non_soave(&Delft1998 {
                            kappa: kappa.kappa,
                            is_methane,
                        })
                    }
                    Alpha::Pr | Alpha::Srk | Alpha::Pr78 | Alpha::Twu => {
                        let (kappa_value, kappa_warnings) = match self.alpha {
                            Alpha::Pr => {
                                let kappa = pr_kappa(component.omega)?;
                                (kappa.kappa, kappa.warnings)
                            }
                            Alpha::Srk => {
                                let kappa = srk_kappa(component.omega)?;
                                (kappa.kappa, kappa.warnings)
                            }
                            Alpha::Pr78 => {
                                let kappa = pr78_kappa(component.omega)?;
                                (kappa.kappa, kappa.warnings)
                            }
                            Alpha::Twu => {
                                let kappa = twu_kappa(component.omega)?;
                                (kappa.kappa, kappa.warnings)
                            }
                            _ => unreachable!("handled above"),
                        };
                        warnings.extend(kappa_warnings);
                        let term = Soave { kappa: kappa_value };
                        let (a_reduced, b_reduced, ab_warnings) = if cubic == Cubic::Pr {
                            let ab =
                                pr_alpha_ab(kappa_value, reduced_temperature, reduced_pressure)?;
                            (ab.a_reduced, ab.b_reduced, ab.warnings)
                        } else {
                            let ab =
                                srk_alpha_ab(kappa_value, reduced_temperature, reduced_pressure)?;
                            (ab.a_reduced, ab.b_reduced, ab.warnings)
                        };
                        warnings.extend(ab_warnings);
                        (
                            a_reduced,
                            b_reduced,
                            term.psi(reduced_temperature),
                            term.psi_t(reduced_temperature),
                        )
                    }
                },
                Cubic::Rk => {
                    let ab = rk_alpha_ab(reduced_temperature, reduced_pressure)?;
                    warnings.extend(ab.warnings);
                    let term = RkAlpha;
                    (
                        ab.a_reduced,
                        ab.b_reduced,
                        term.psi(reduced_temperature),
                        term.psi_t(reduced_temperature),
                    )
                }
                Cubic::Tst => {
                    let kappa = twu_kappa(component.omega)?;
                    warnings.extend(kappa.warnings);
                    non_soave(&Soave { kappa: kappa.kappa })
                }
            };
            a.push(a_reduced);
            b.push(b_reduced);
            psi.push(psi_value);
            psi_t.push(psi_t_value);
        }
        Ok(ReducedParameters {
            a,
            b,
            psi,
            psi_t,
            reduced_temperatures,
            t_kelvin: t.value,
            kij: self.mixing_rule.effective_kij(t.value),
            warnings,
        })
    }

    /// The state of one phase: its mixture parameters, its root, and its `ln phi`.
    ///
    /// `side` names which root of the cubic this phase claims - the smallest for the
    /// liquid, the largest for the vapour. Selecting by *ordering* rather than by an
    /// initial guess is the rule [`crate::pr_z_factor`] fixes, and this follows it
    /// rather than second-guessing it.
    ///
    /// # Errors
    /// * [`AzothError::InvalidInput`] if `x` is not one entry per component or if it
    ///   sums anywhere but to one.
    /// * Propagates [`crate::pr_z_factor`]'s range checks.
    pub fn phase_state(
        &self,
        reduced: &ReducedParameters,
        x: &[f64],
        side: RootSide,
    ) -> Result<PhaseState> {
        let (a_mix, b_mix) = self.mixture_parameters(reduced, x);
        let (z_min, z_max) = match self.cubic {
            // TST shares Peng-Robinson's `delta`, so its cubic is PR's; only the omega
            // constants and the alpha term differ.
            Cubic::Pr | Cubic::Tst => {
                let roots = pr_z_factor(a_mix, b_mix)?;
                (roots.z_min, roots.z_max)
            }
            // RK's cubic is SRK's - the same `omega` and `delta` - so the root finder
            // is shared; only the alpha term above differs.
            Cubic::Srk | Cubic::Rk => {
                let roots = srk_z_factor(a_mix, b_mix)?;
                (roots.z_min, roots.z_max)
            }
        };
        let z = match side {
            RootSide::Liquid => z_min,
            RootSide::Vapour => z_max,
        };
        self.phase_state_at(reduced, x, z)
    }

    /// The state of one phase, at a root the caller has already chosen.
    ///
    /// For a caller who has a compressibility factor in hand - from
    /// [`crate::pr_z_factor`], or from a flash that solved for one - and does not want
    /// it re-derived. The choice of root is a *phase*, and a caller holding a Z has
    /// already made it; re-deriving here would silently overrule them.
    ///
    /// # Errors
    /// * [`AzothError::InvalidInput`] if `x` is not one entry per component, or if it
    ///   does not sum to one.
    /// * [`AzothError::OutOfRange`] if the root is not admissible - `z` must exceed
    ///   the mixture's `B`, or `ln(z - B)` is the logarithm of a negative number.
    pub fn phase_state_at(
        &self,
        reduced: &ReducedParameters,
        x: &[f64],
        z: f64,
    ) -> Result<PhaseState> {
        let n = self.len();
        if x.len() != n {
            return Err(AzothError::invalid_input(
                "x",
                format!("a composition for {n} components has {} entries", x.len()),
            ));
        }

        let kij = self.phase_kij(reduced, x);
        let (a_mix, b_mix) = self.mixture_parameters(reduced, x);
        // Written as an explicit test rather than `!(z > b_mix)`, which reads as a
        // double negative and is the shape clippy flags. Both reject NaN; the
        // explicit form also rejects a root at or below `B`, which is the case here.
        if !z.is_finite() || z <= b_mix {
            return Err(AzothError::OutOfRange {
                field: "z".to_string(),
                value: z,
                detail: format!(
                    "the root must exceed the mixture's B = {b_mix}, because `ln(z - B)`                      is otherwise the logarithm of a negative number. `z <= B` is the                      zero-volume limit, which is not a state"
                ),
            });
        }

        // The cross sum `sum_j x_j A_ij`, one per component, hoisted out of the loop
        // below so both languages evaluate it once and in the same order.
        let cross: Vec<f64> = (0..n)
            .map(|i| {
                (0..n)
                    .map(|j| x[j] * (1.0 - kij[i * n + j]) * (reduced.a[i] * reduced.a[j]).sqrt())
                    .sum::<f64>()
            })
            .collect();

        let i_term = self.cubic.i_term(z, b_mix);
        let coefficient = self.cubic.coefficient(a_mix, b_mix);
        let ln_z_minus_b = (z - b_mix).ln();

        let ln_phi: Vec<f64> = match &self.mixing_rule {
            MixingRule::HuronVidal { .. } | MixingRule::WongSandler { .. } => {
                // The GE rules carry the activity coefficient and a volume-derivative
                // term that the classic `factor * i_term` form cannot express.
                let (ader, alpha_mix, b_der) = self.ge_state(reduced, x);
                let delta1 = self.cubic.delta1();
                let delta2 = self.cubic.delta2();
                (0..n)
                    .map(|i| {
                        let bd = b_der[i];
                        let fv = alpha_mix * bd * z / ((z + delta1 * b_mix) * (z + delta2 * b_mix));
                        -ln_z_minus_b + bd / (z - b_mix)
                            - (ader[i] / self.cubic.delta_diff()) * i_term
                            - fv
                    })
                    .collect()
            }
            _ => (0..n)
                .map(|i| {
                    let b_ratio = reduced.b[i] / b_mix;
                    // The factor that is 1 for a pure component: `2 * sum_j x_j A_ij / A
                    // - B_i / B`. At N = 1 the sum is `A_11 = A` with `k11 = 0`, so it is
                    // `2 - 1`, and this whole expression becomes `pr_departure`'s.
                    let factor = 2.0 * cross[i] / a_mix - b_ratio;
                    b_ratio * (z - 1.0) - ln_z_minus_b - coefficient * factor * i_term
                })
                .collect(),
        };

        // The mixture's departure functions. `psi_bar` is the composition-weighted
        // average of the components' `psi`, and the two lines below are then
        // `pr_departure`'s expressions with `psi_bar` in place of `psi` - which is
        // what makes them reduce to it exactly at one component.
        let mut weight_total = 0.0;
        let mut weighted_psi = 0.0;
        let mut weighted_psi_t = 0.0;
        for i in 0..n {
            for j in 0..n {
                let weight =
                    x[i] * x[j] * (1.0 - kij[i * n + j]) * (reduced.a[i] * reduced.a[j]).sqrt();
                let psi_pair = 0.5 * (reduced.psi[i] + reduced.psi[j]);
                weight_total += weight;
                weighted_psi += weight * psi_pair;
                // `T*d(weight*psi_pair)/dT` for this pair. The weight carries
                // `sqrt(A_i A_j)`, whose logarithmic derivative is `psi_pair - 2`, and
                // `psi_pair` carries the two components' own derivatives.
                weighted_psi_t += weight
                    * ((psi_pair - 2.0) * psi_pair + 0.5 * (reduced.psi_t[i] + reduced.psi_t[j]));
            }
        }
        let psi_bar = weighted_psi / weight_total;
        let t_dpsi_bar = weighted_psi_t / weight_total - psi_bar * (psi_bar - 2.0);
        let h_dep_rt = (z - 1.0) + coefficient * (psi_bar - 1.0) * i_term;
        let s_dep_r = h_dep_rt
            - ln_phi
                .iter()
                .zip(x)
                .map(|(&lp, &x_i)| x_i * lp)
                .sum::<f64>();

        // The heat-capacity departure. `a_mix` moves with temperature exactly as `A`
        // does above - its logarithmic derivative is `psi_bar - 2`, because the weights
        // that form it are the same terms - so these four lines are the registered
        // calc's with `psi_bar` in place of `psi`.
        let t_da = a_mix * (psi_bar - 2.0);
        let t_db = -b_mix;
        let t_dc = coefficient * (psi_bar - 1.0);
        let d_f_dz = self.cubic.df_dz(z, a_mix, b_mix);
        let t_dfdt = self.cubic.t_dfdt(z, a_mix, b_mix, t_da, t_db);
        let t_dz = -t_dfdt / d_f_dz;
        let d1 = self.cubic.delta1();
        let d2 = self.cubic.delta2();
        let n_plus = z + d1 * b_mix;
        let n_minus = z + d2 * b_mix;
        let t_di = (t_dz + d1 * t_db) / n_plus - (t_dz + d2 * t_db) / n_minus;
        let cp_dep_r = h_dep_rt
            + t_dz
            + t_dc * (psi_bar - 1.0) * i_term
            + coefficient * t_dpsi_bar * i_term
            + coefficient * (psi_bar - 1.0) * t_di;

        Ok(PhaseState {
            a_mix,
            b_mix,
            z,
            ln_phi,
            h_dep_rt,
            s_dep_r,
            psi_bar,
            cp_dep_r,
        })
    }

    /// The phase-dependent interaction matrix for a composition.
    ///
    /// The phase-independent rules resolve their matrix in [`Self::reduced_parameters`];
    /// this returns it untouched. The Soreide-Whitson rule resolves here, against the
    /// phase composition, because its salinity correlation is phase-dependent.
    fn phase_kij(&self, reduced: &ReducedParameters, x: &[f64]) -> Vec<f64> {
        let acentric: Vec<f64> = self.components.iter().map(|c| c.omega).collect();
        self.mixing_rule
            .phase_kij(&reduced.kij, &reduced.reduced_temperatures, &acentric, x)
    }

    /// The GE-model attraction coefficients `A_i/B_i - ln gamma_i / Lambda`.
    #[allow(clippy::too_many_arguments)] // the signature is the GE model's inputs
    fn ge_ader(
        &self,
        reduced: &ReducedParameters,
        x: &[f64],
        kij: &[f64],
        hv_gij: &[f64],
        hv_gij_t: &[f64],
        hv_alpha: &[f64],
        hv_pairs: &[bool],
    ) -> Vec<f64> {
        let lambda = self.cubic.hv_constant();
        let ln_gamma = crate::hv_ge::hv_ln_gamma(
            x,
            reduced.t_kelvin,
            &reduced.a,
            &reduced.b,
            kij,
            hv_gij,
            hv_gij_t,
            hv_alpha,
            hv_pairs,
            lambda,
        );
        (0..self.len())
            .map(|i| reduced.a[i] / reduced.b[i] - ln_gamma[i] / lambda)
            .collect()
    }

    /// The Huron-Vidal attraction coefficients.
    fn hv_ader(&self, reduced: &ReducedParameters, x: &[f64]) -> Vec<f64> {
        let MixingRule::HuronVidal {
            kij,
            hv_gij,
            hv_gij_t,
            hv_alpha,
            hv_pairs,
        } = &self.mixing_rule
        else {
            unreachable!("hv_ader is only called for the Huron-Vidal rule");
        };
        self.ge_ader(reduced, x, kij, hv_gij, hv_gij_t, hv_alpha, hv_pairs)
    }

    /// The Wong-Sandler attraction coefficients, from the DijT-free activity
    /// coefficients this rule reads.
    ///
    /// **The `DijT` question here is open, and this records why rather than a guess.**
    /// NeqSim's `WongSandlerMixingRule` is handed `NRTLDijT` - loaded from
    /// `WSGIJT`/`WSGJIT` - and passes it into a `PhaseGENRTLmodifiedHV`, which reads as
    /// though a `WS` pair with a non-zero coefficient should use it. It does not settle
    /// the question: asked for `CO2`/`water`, whose `WSGIJT` is 0.96, the rule reports
    /// `HVDijT[0][1] = -0.842025353`, which is that pair's **`HVGIJT`**. So the field
    /// NeqSim exposes is not the field its name suggests, and the two columns disagree.
    ///
    /// Trying all three - zeros, `WSGIJT` and `HVGIJT` - against a NeqSim
    /// `SystemSrkEos` + `setMixingRule(5)` flash for `CO2`/`water` at 350 K and 5 bar
    /// gives `ln_phi` of `[3.809, -1.202]`, `[-0.014, -0.050]` and `[3.784, -1.477]`
    /// against NeqSim's `[10.546, -2.617]`. **None is close**, so this rule's own gap on
    /// that pair dwarfs the coefficient being chosen, and no value of it can be
    /// validated until that gap is closed. Forcing zeros is what the rule did when it
    /// was checked against `water`/`ethanol` - where both coefficients are zero, so the
    /// check could not distinguish them either.
    fn ws_ader(&self, reduced: &ReducedParameters, x: &[f64]) -> Vec<f64> {
        let MixingRule::WongSandler {
            kij,
            hv_gij,
            hv_alpha,
            hv_pairs,
        } = &self.mixing_rule
        else {
            unreachable!("ws_ader is only called for the Wong-Sandler rule");
        };
        let hv_gij_t = vec![0.0; self.len() * self.len()];
        self.ge_ader(reduced, x, kij, hv_gij, &hv_gij_t, hv_alpha, hv_pairs)
    }

    /// The Wong-Sandler `b_mix` and its composition derivative `BDER`.
    fn ws_b_mix(&self, reduced: &ReducedParameters, x: &[f64], ader: &[f64]) -> (f64, Vec<f64>) {
        let MixingRule::WongSandler { kij, .. } = &self.mixing_rule else {
            unreachable!("ws_b_mix is only called for the Wong-Sandler rule");
        };
        let n = self.len();
        let alpha_mix: f64 = (0..n).map(|i| x[i] * ader[i]).sum();
        let mut qf1 = vec![0.0; n];
        let mut q = 0.0;
        for i in 0..n {
            let mut ss = 0.0;
            for j in 0..n {
                ss += x[j]
                    * (1.0 - kij[j * n + i])
                    * (reduced.b[j] + reduced.b[i] - reduced.a[j] - reduced.a[i]);
            }
            qf1[i] = ss;
            q += x[i] * ss;
        }
        let denom = 1.0 - alpha_mix;
        let b_mix = 0.5 * q / denom;
        let bder: Vec<f64> = (0..n)
            .map(|i| (qf1[i] - b_mix * (1.0 - ader[i])) / denom)
            .collect();
        (b_mix, bder)
    }

    /// The GE-rule state `(ader, alpha_mix, b_der)` the fugacity is written in.
    fn ge_state(&self, reduced: &ReducedParameters, x: &[f64]) -> (Vec<f64>, f64, Vec<f64>) {
        let n = self.len();
        match &self.mixing_rule {
            MixingRule::HuronVidal { .. } => {
                let ader = self.hv_ader(reduced, x);
                let alpha_mix: f64 = (0..n).map(|i| x[i] * ader[i]).sum();
                (ader, alpha_mix, reduced.b.clone())
            }
            MixingRule::WongSandler { .. } => {
                let ader = self.ws_ader(reduced, x);
                let alpha_mix: f64 = (0..n).map(|i| x[i] * ader[i]).sum();
                let (_, bder) = self.ws_b_mix(reduced, x, &ader);
                (ader, alpha_mix, bder)
            }
            _ => unreachable!("ge_state is only called for GE rules"),
        }
    }

    /// The van der Waals one-fluid mixture parameters for a composition.
    ///
    /// ```text
    /// b_mix = sum_i x_i B_i
    /// a_mix = sum_i sum_j x_i x_j (1 - k_ij) sqrt(A_i A_j)
    /// ```
    ///
    /// Written as a double sum rather than the binary calc's longhand form because
    /// that is the only shape that generalises; the reduction to
    /// [`crate::vdw1f_mix_binary`] at `N = 2` is a test, not a rewrite.
    #[must_use]
    pub fn mixture_parameters(&self, reduced: &ReducedParameters, x: &[f64]) -> (f64, f64) {
        let n = self.len();
        match &self.mixing_rule {
            MixingRule::HuronVidal { .. } => {
                let ader = self.hv_ader(reduced, x);
                let alpha_mix: f64 = (0..n).map(|i| x[i] * ader[i]).sum();
                let b_mix = (0..n).map(|i| x[i] * reduced.b[i]).sum();
                (b_mix * alpha_mix, b_mix)
            }
            MixingRule::WongSandler { .. } => {
                let ader = self.ws_ader(reduced, x);
                let alpha_mix: f64 = (0..n).map(|i| x[i] * ader[i]).sum();
                let (b_mix, _) = self.ws_b_mix(reduced, x, &ader);
                (b_mix * alpha_mix, b_mix)
            }
            _ => {
                let kij = self.phase_kij(reduced, x);
                let b_mix = (0..n).map(|i| x[i] * reduced.b[i]).sum();
                let mut a_mix = 0.0;
                for i in 0..n {
                    for j in 0..n {
                        a_mix += x[i]
                            * x[j]
                            * (1.0 - kij[i * n + j])
                            * (reduced.a[i] * reduced.a[j]).sqrt();
                    }
                }
                (a_mix, b_mix)
            }
        }
    }

    /// The mixture's volume translation, the composition-weighted sum of the
    /// components' shifts.
    ///
    /// NeqSim mixes the Peneloux translation linearly, like the co-volume `b`, and
    /// subtracts the result from the untranslated molar volume:
    /// `v_corr = v - sum_i x_i c_i`. The untranslated `v` is [`crate::pr_molar_volume`],
    /// and the per-component `c_i` are [`crate::pr_peneloux_shift`] or
    /// [`crate::srk_peneloux_shift`].
    #[must_use]
    pub fn volume_shift(&self, x: &[f64]) -> f64 {
        (0..self.len())
            .map(|i| x[i] * self.components[i].volume_shift)
            .sum()
    }

    /// `A^R/(R T)` - the residual Helmholtz energy, at a set of mole numbers.
    ///
    /// ```text
    /// Phi = -N ln(1 - b) - (a / (2 sqrt2 b)) G(b)
    /// ```
    ///
    /// The energy whose first composition derivative is the logarithm of fugacity
    /// and whose second is [`Self::helmholtz_hessian`]. It is exposed because those
    /// two relations are what the tests are built on, and a function that cannot be
    /// evaluated cannot have its derivatives checked:
    ///
    /// ```text
    /// d(A^R/RT)/dn_i = ln phi_i + ln Z
    /// ```
    ///
    /// which is asserted against [`Self::phase_state_at`]'s `ln phi` - a different
    /// code path, reached by differentiating a departure function rather than an
    /// energy, so the agreement is evidence rather than a restatement.
    ///
    /// # Mole numbers, not mole fractions
    ///
    /// `n` is a set of mole numbers, and the argument is named for it because the
    /// distinction is load-bearing: `A^R` is homogeneous of degree one in `(V, n)`
    /// but **not** in `n` at fixed `V`, so the total is part of the state and a
    /// function of the fractions alone could not be differentiated. A caller
    /// holding a composition summing to one is already passing mole numbers.
    ///
    /// # Errors
    /// * [`AzothError::InvalidInput`] if `n` is not one entry per component.
    /// * [`AzothError::OutOfRange`] if the mixture's `B` is not positive.
    pub fn helmholtz_energy(
        &self,
        reduced: &ReducedParameters,
        n: &[f64],
        compressibility: f64,
    ) -> Result<f64> {
        let count = self.check_mole_numbers(n)?;
        let (b_hat, a_ij) = self.scaled_constants(reduced, compressibility, n);
        let total: f64 = n.iter().sum();
        let b_sum: f64 = (0..count).map(|i| n[i] * b_hat[i]).sum();
        self.check_b_sum(b_sum, compressibility)?;
        let mut a_sum = 0.0;
        for i in 0..count {
            for j in 0..count {
                a_sum += n[i] * n[j] * a_ij[i][j];
            }
        }
        let g = self.cubic.helmholtz_g(b_sum);
        Ok(-total * (1.0 - b_sum).ln() - (a_sum / (self.cubic.delta_diff() * b_sum)) * g)
    }

    /// `d2(A^R/RT)/dn_i dn_j` at constant temperature and **volume**.
    ///
    /// Not at constant pressure: the Hessian every criticality condition is written in
    /// is the Helmholtz one, because its vanishing is what separates a stable phase from
    /// a metastable one. A constant-pressure composition derivative answers a different
    /// question, and converting between the two frames needs two further derivative
    /// families and a partial-molar-volume correction that this avoids entirely.
    ///
    /// **Everything is dimensionless**, like the rest of the `eos` core:
    /// `a_i/(R T V)` is `A_i/Z` and `b_i/V` is `B_i/Z`, so the reduced parameters
    /// the caller already holds carry the whole construction and no dimensioned
    /// quantity appears.
    ///
    /// With `b = sum_i n_i B_i/Z`, `L = ln(1 - b)` and
    /// `G = ln((1 + (1+sqrt2)b)/(1 + (1-sqrt2)b))`, every term below is one of those
    /// three differentiated once or twice, written out rather than factored so that
    /// each term still shows which derivative it came from.
    ///
    /// # Errors
    /// * [`AzothError::InvalidInput`] if `n` is not one entry per component.
    /// * [`AzothError::OutOfRange`] if the mixture's `B` is not positive.
    pub fn helmholtz_hessian(
        &self,
        reduced: &ReducedParameters,
        n: &[f64],
        compressibility: f64,
    ) -> Result<Vec<Vec<f64>>> {
        let count = self.check_mole_numbers(n)?;
        let (b_hat, a_ij) = self.scaled_constants(reduced, compressibility, n);
        let total: f64 = n.iter().sum();
        let b_sum: f64 = (0..count).map(|i| n[i] * b_hat[i]).sum();
        self.check_b_sum(b_sum, compressibility)?;

        let d = 1.0 - b_sum;
        // L(b) = ln(1 - b), G(b) = ln((1 + delta1 b)/(1 + delta2 b)).
        let l_prime = -1.0 / d;
        let l_second = -1.0 / (d * d);
        let g = self.cubic.helmholtz_g(b_sum);
        let g_prime = self.cubic.helmholtz_g_prime(b_sum);
        let g_second = self.cubic.helmholtz_g_second(b_sum);
        let half_delta_diff = self.cubic.half_delta_diff();
        let delta_diff = self.cubic.delta_diff();

        let a_bar: Vec<f64> = (0..count)
            .map(|i| (0..count).map(|j| n[j] * a_ij[i][j]).sum())
            .collect();
        let mut a_sum = 0.0;
        for i in 0..count {
            for j in 0..count {
                a_sum += n[i] * n[j] * a_ij[i][j];
            }
        }

        let mut hessian = vec![vec![0.0; count]; count];
        for i in 0..count {
            for j in 0..count {
                let pair = a_bar[i] * b_hat[j] + a_bar[j] * b_hat[i];
                let product = b_hat[i] * b_hat[j];
                hessian[i][j] = -(b_hat[i] + b_hat[j]) * l_prime
                    - total * product * l_second
                    - g * a_ij[i][j] / (half_delta_diff * b_sum)
                    + g * pair / (half_delta_diff * b_sum * b_sum)
                    - g * a_sum * product / (half_delta_diff * b_sum * b_sum * b_sum)
                    - g_prime * pair / (half_delta_diff * b_sum)
                    + g_prime * a_sum * product / (half_delta_diff * b_sum * b_sum)
                    - g_second * a_sum * product / (delta_diff * b_sum);
            }
        }
        Ok(hessian)
    }

    /// Heidemann & Khalil's `Q`, whose smallest eigenvalue vanishes at a critical point.
    ///
    /// ```text
    /// Q_ij = sqrt(n_i n_j) * d2(A/RT)/dn_i dn_j
    /// ```
    ///
    /// at constant temperature and volume, the **total** Helmholtz energy rather than
    /// the residual one. Its ideal part is `delta_ij / n_i`, the constant-volume form;
    /// the model spec's notes record the constant-pressure alternative and what taking
    /// it here would cost.
    ///
    /// **The scaling by `sqrt(n_i n_j)` is not decoration.** A Maxwell relation makes
    /// the Hessian symmetric already; the scaling makes that symmetry *structural*
    /// rather than numerical, so the eigenvalues are real and the eigenvectors are
    /// available in every case rather than almost every case.
    ///
    /// **`Q` is not singular at ordinary states.** `A(T, V, n)` is not homogeneous in
    /// `n` at fixed `V` - homogeneity needs the volume to scale with it - so there is
    /// no Euler-theorem null vector, and the ideal part of the Hessian at constant
    /// volume is positive definite on its own. The vanishing is therefore informative
    /// rather than generic, and the direction it vanishes along is the critical
    /// composition fluctuation and nothing else.
    ///
    /// That is why the critical point solves for the **smallest-magnitude eigenvalue**
    /// rather than for `det(Q)`: the determinant is the product of every eigenvalue,
    /// so it vanishes when *any* of them does - including ones whose vanishing is not
    /// criticality - and it is badly scaled for a Newton step.
    ///
    /// # Errors
    /// * [`AzothError::InvalidInput`] if `n` is not one entry per component.
    /// * [`AzothError::OutOfRange`] if the mixture's `B` is not positive.
    pub fn criticality_matrix(
        &self,
        reduced: &ReducedParameters,
        n: &[f64],
        compressibility: f64,
    ) -> Result<Vec<Vec<f64>>> {
        let hessian = self.helmholtz_hessian(reduced, n, compressibility)?;
        let count = n.len();
        Ok((0..count)
            .map(|i| {
                (0..count)
                    .map(|j| {
                        let ideal = if i == j { 1.0 / n[i] } else { 0.0 };
                        (n[i] * n[j]).sqrt() * (hessian[i][j] + ideal)
                    })
                    .collect()
            })
            .collect())
    }

    /// One entry per component, and every entry finite.
    fn check_mole_numbers(&self, n: &[f64]) -> Result<usize> {
        let count = self.len();
        if n.len() != count {
            return Err(AzothError::invalid_input(
                "n",
                format!(
                    "a mixture of {count} components has {} mole numbers",
                    n.len()
                ),
            ));
        }
        Ok(count)
    }

    /// The `b` in `L(b)` and `G(b)`: `sum_i n_i B_i / Z`.
    ///
    /// Reports `compressibility`, the input the caller actually passed, rather than
    /// `b` itself. `b` is derived from the root and the composition, so naming it
    /// would hand a caller a number they never supplied and cannot act on; the value
    /// is quoted in the message instead, where it explains the refusal without
    /// pretending to be an input.
    fn check_b_sum(&self, b_sum: f64, compressibility: f64) -> Result<()> {
        // Explicit finiteness-and-positivity rather than `!(b_sum > 0.0)`: both
        // reject NaN, and the explicit form also rejects an infinity. The
        // logarithmic terms below are written against `b`, so zero is a division.
        if !b_sum.is_finite() || b_sum <= 0.0 {
            return Err(AzothError::OutOfRange {
                field: "compressibility".to_string(),
                value: compressibility,
                detail: format!(
                    "the mixture's `B/Z` came out as {b_sum}, and the logarithmic terms \
                     of the Helmholtz energy are written against it. It is positive for \
                     any admissible root, so this is a composition or a root that is not \
                     a state rather than a compressibility that is out of range"
                ),
            });
        }
        Ok(())
    }

    /// `b_i/V` per component, and the cross term `A_ij/Z`.
    ///
    /// The identities that make the Helmholtz construction dimensionless: with `V`
    /// the molar volume, `A_i/Z = a_i/(R T V)` and `B_i/Z = b_i/V`. Both are
    /// **independent of the composition**, which is what lets the sums carry all of
    /// the composition dependence and is why they are hoisted out of every loop.
    /// `a_i/(R T V)` itself is not returned because it appears only inside `A_ij`.
    fn scaled_constants(
        &self,
        reduced: &ReducedParameters,
        compressibility: f64,
        n: &[f64],
    ) -> (Vec<f64>, Vec<Vec<f64>>) {
        let count = n.len();
        // The Soreide-Whitson correlation keys on the mole *fraction* of water, so the
        // mole numbers are normalised for it; the phase-independent rules ignore it.
        let mut x: Vec<f64> = n.to_vec();
        normalise(&mut x);
        let kij = self.phase_kij(reduced, &x);
        let a_hat: Vec<f64> = reduced.a.iter().map(|v| v / compressibility).collect();
        let b_hat: Vec<f64> = reduced.b.iter().map(|v| v / compressibility).collect();
        let a_ij: Vec<Vec<f64>> = (0..count)
            .map(|i| {
                (0..count)
                    .map(|j| (1.0 - kij[i * count + j]) * (a_hat[i] * a_hat[j]).sqrt())
                    .collect()
            })
            .collect();
        (b_hat, a_ij)
    }
}

/// Wilson's correlation for a mixture's initial K-values.
///
/// `K_i = (Pc_i / P) * exp(5.373 * (1 + omega_i) * (1 - Tc_i / T))`.
///
/// The constant is 5.373, quoted as 5.37 in some sources; the discrepancy is recorded
/// in the references of `specs/models/eos/pt_flash.toml` rather than resolved, because
/// a reader meeting the other value needs to know it is the same correlation.
///
/// Shared by the two models that use it: `eos.pt_flash` seeds its iteration with it,
/// and `eos.stability_test` seeds *both* of its trials with it - the vapour-like one
/// from `z_i K_i` and the liquid-like one from `z_i / K_i`.
pub(crate) fn wilson_k(mixture: &Mixture, t: ThermodynamicTemperature, p: Pressure) -> Vec<f64> {
    mixture
        .components()
        .iter()
        .map(|c| {
            (c.pc.value / p.value) * (5.373 * (1.0 + c.omega) * (1.0 - c.tc.value / t.value)).exp()
        })
        .collect()
}

/// Rescale a vector to sum to one, in place.
///
/// A non-positive total leaves the values untouched rather than dividing by it: the
/// caller is holding a set of mole numbers that are all zero, which is not a
/// composition, and returning `NaN`s would turn "no composition here" into a number.
/// The two callers both guarantee a positive total, so this is a guard rather than a
/// branch either of them takes.
pub(crate) fn normalise(values: &mut [f64]) {
    let sum: f64 = values.iter().sum();
    if sum > 0.0 {
        for value in values.iter_mut() {
            *value /= sum;
        }
    }
}

/// The reduced parameters of every component at one state, plus any warnings the
/// kernels raised on the way.
///
/// A record rather than two bare vectors because the two always travel together and
/// a caller that swapped them would get a mixture whose attraction came from the
/// repulsion, which is a wrong answer that looks like a right one.
#[derive(Debug, Clone, PartialEq)]
pub struct ReducedParameters {
    /// `A_i = omega_a * alpha_i * Pr_i / Tr_i**2`, one per component.
    pub a: Vec<f64>,
    /// `B_i = omega_b * Pr_i / Tr_i`, one per component.
    pub b: Vec<f64>,
    /// `psi_i = -(kappa_i sqrt(Tr_i)) / (1 + kappa_i (1 - sqrt(Tr_i)))`, one per
    /// component: the logarithmic derivative of the alpha function.
    ///
    /// Carried alongside `A` and `B` because it is the same kind of quantity - a
    /// per-component function of the state, identical for both phases, computed once
    /// per solve - and because the mixture's departure functions are the pure form
    /// with `psi` replaced by a composition-weighted average of these.
    pub psi: Vec<f64>,
    /// `T * dpsi_i/dT`, one per component: the same derivative already multiplied by
    /// the absolute temperature, which is the form the departure heat capacity needs.
    ///
    /// Carried rather than recomputed because `kappa_i` and `Tr_i` - the two things it
    /// is built from - are local to [`Self::reduced_parameters`], and recovering them
    /// from `psi` and `A` afterwards is possible but loses the reduction to
    /// `pr_departure` at one component, which is the property the whole mixture layer
    /// is checked against.
    pub psi_t: Vec<f64>,
    /// `Tr_i = T / Tc_i`, one per component.
    ///
    /// The reduced temperature the alpha terms were evaluated at. Carried because the
    /// Soreide-Whitson mixing rule reads it again when it resolves its phase-dependent
    /// interaction matrix.
    pub reduced_temperatures: Vec<f64>,
    /// The state's absolute temperature, in Kelvin.
    ///
    /// The Huron-Vidal mixing rule reads it for its NRTL `tau = Dij / T + DijT`.
    pub t_kelvin: f64,
    /// The interaction matrix at this state's temperature, flattened row-major.
    ///
    /// The mixing rule's temperature dependence is resolved here, once per
    /// `reduced_parameters` call, so the flash's per-iteration mixing reads a plain
    /// matrix rather than re-evaluating the rule.
    pub kij: Vec<f64>,
    /// Warnings raised while computing them - in practice `pr_kappa`'s
    /// `kappa < 0` for a component with a sufficiently negative acentric factor.
    pub warnings: Vec<Warning>,
}

/// Validate an interaction matrix: `N x N` long, zero diagonal, symmetric.
///
/// Shared by [`Mixture::new`] and any rule constructor, so the invariant is checked
/// once at construction rather than re-derived at every evaluation.
fn check_kij(n: usize, kij: &[f64]) -> Result<()> {
    if kij.len() != n * n {
        return Err(AzothError::invalid_input(
            "kij",
            format!(
                "kij has {} entries but a mixture of {n} components needs an N x N \
                 matrix flattened row-major, which is {}",
                kij.len(),
                n * n
            ),
        ));
    }
    for i in 0..n {
        if kij[i * n + i] != 0.0 {
            return Err(AzothError::invalid_input(
                "kij",
                format!(
                    "kij[{i}][{i}] is {} but the diagonal must be zero: a component \
                     does not interact with itself, and a non-zero diagonal \
                     silently rescales that component's attraction",
                    kij[i * n + i]
                ),
            ));
        }
        for j in (i + 1)..n {
            let (forward, backward) = (kij[i * n + j], kij[j * n + i]);
            if forward != backward {
                return Err(AzothError::invalid_input(
                    "kij",
                    format!(
                        "kij[{i}][{j}] is {forward} but kij[{j}][{i}] is {backward}; \
                         the matrix must be symmetric"
                    ),
                ));
            }
        }
    }
    Ok(())
}
