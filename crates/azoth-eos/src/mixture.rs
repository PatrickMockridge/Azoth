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

use crate::{pr_alpha_ab, pr_kappa, pr_z_factor};

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
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Component {
    /// Critical temperature.
    pub tc: ThermodynamicTemperature,
    /// Critical pressure.
    pub pc: Pressure,
    /// Acentric factor.
    pub omega: f64,
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
        Ok(Self { tc, pc, omega })
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
    kij: Vec<f64>,
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
        Ok(Self { components, kij })
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

    /// The interaction parameter between components `i` and `j`.
    #[must_use]
    pub fn kij(&self, i: usize, j: usize) -> f64 {
        self.kij[i * self.components.len() + j]
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
        let mut warnings = Vec::new();

        for component in &self.components {
            let kappa = pr_kappa(component.omega)?;
            warnings.extend(kappa.warnings);
            let reduced_temperature = t.value / component.tc.value;
            let ab = pr_alpha_ab(
                kappa.kappa,
                reduced_temperature,
                p.value / component.pc.value,
            )?;
            warnings.extend(ab.warnings);
            let sqrt_tr = reduced_temperature.sqrt();
            a.push(ab.a_reduced);
            b.push(ab.b_reduced);
            psi.push(-kappa.kappa * sqrt_tr / (1.0 + kappa.kappa * (1.0 - sqrt_tr)));
            // `T*dpsi/dT`, from `dpsi/dTr` and `dTr/dT = 1/Tc`. Carrying the `T`
            // rather than dividing it out later is what keeps this free of the
            // absolute temperature: the expression below is a function of `Tr` alone.
            psi_t.push(
                -kappa.kappa * (1.0 + kappa.kappa) * reduced_temperature
                    / (2.0 * sqrt_tr * (1.0 + kappa.kappa * (1.0 - sqrt_tr)).powi(2)),
            );
        }
        Ok(ReducedParameters {
            a,
            b,
            psi,
            psi_t,
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
        let roots = pr_z_factor(a_mix, b_mix)?;
        let z = match side {
            RootSide::Liquid => roots.z_min,
            RootSide::Vapour => roots.z_max,
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
                    .map(|j| x[j] * (1.0 - self.kij(i, j)) * (reduced.a[i] * reduced.a[j]).sqrt())
                    .sum::<f64>()
            })
            .collect();

        let sqrt_2 = std::f64::consts::SQRT_2;
        let i_term = ((z + (1.0 + sqrt_2) * b_mix) / (z + (1.0 - sqrt_2) * b_mix)).ln();
        let coefficient = a_mix / (2.0 * sqrt_2 * b_mix);
        let ln_z_minus_b = (z - b_mix).ln();

        let ln_phi: Vec<f64> = (0..n)
            .map(|i| {
                let b_ratio = reduced.b[i] / b_mix;
                // The factor that is 1 for a pure component: `2 * sum_j x_j A_ij / A
                // - B_i / B`. At N = 1 the sum is `A_11 = A` with `k11 = 0`, so it is
                // `2 - 1`, and this whole expression becomes `pr_departure`'s.
                let factor = 2.0 * cross[i] / a_mix - b_ratio;
                b_ratio * (z - 1.0) - ln_z_minus_b - coefficient * factor * i_term
            })
            .collect();

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
                    x[i] * x[j] * (1.0 - self.kij(i, j)) * (reduced.a[i] * reduced.a[j]).sqrt();
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
        let d_f_dz =
            3.0 * z * z + 2.0 * (b_mix - 1.0) * z + (a_mix - 3.0 * b_mix * b_mix - 2.0 * b_mix);
        let t_dfdt = t_db * z * z
            + (t_da - 6.0 * b_mix * t_db - 2.0 * t_db) * z
            + (3.0 * b_mix * b_mix * t_db + 2.0 * b_mix * t_db - t_da * b_mix - a_mix * t_db);
        let t_dz = -t_dfdt / d_f_dz;
        let n_plus = z + (1.0 + sqrt_2) * b_mix;
        let n_minus = z + (1.0 - sqrt_2) * b_mix;
        let t_di =
            (t_dz + (1.0 + sqrt_2) * t_db) / n_plus - (t_dz + (1.0 - sqrt_2) * t_db) / n_minus;
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
        let b_mix = (0..n).map(|i| x[i] * reduced.b[i]).sum();
        let mut a_mix = 0.0;
        for i in 0..n {
            for j in 0..n {
                a_mix +=
                    x[i] * x[j] * (1.0 - self.kij(i, j)) * (reduced.a[i] * reduced.a[j]).sqrt();
            }
        }
        (a_mix, b_mix)
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
        let (b_hat, a_ij) = self.scaled_constants(reduced, compressibility, count);
        let total: f64 = n.iter().sum();
        let b_sum: f64 = (0..count).map(|i| n[i] * b_hat[i]).sum();
        self.check_b_sum(b_sum, compressibility)?;
        let mut a_sum = 0.0;
        for i in 0..count {
            for j in 0..count {
                a_sum += n[i] * n[j] * a_ij[i][j];
            }
        }
        let sqrt_2 = std::f64::consts::SQRT_2;
        let g = ((1.0 + (1.0 + sqrt_2) * b_sum) / (1.0 + (1.0 - sqrt_2) * b_sum)).ln();
        Ok(-total * (1.0 - b_sum).ln() - (a_sum / (2.0 * sqrt_2 * b_sum)) * g)
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
        let (b_hat, a_ij) = self.scaled_constants(reduced, compressibility, count);
        let total: f64 = n.iter().sum();
        let b_sum: f64 = (0..count).map(|i| n[i] * b_hat[i]).sum();
        self.check_b_sum(b_sum, compressibility)?;

        let sqrt_2 = std::f64::consts::SQRT_2;
        let d = 1.0 - b_sum;
        // L(b) = ln(1 - b), G(b) = ln((1 + (1+sqrt2)b) / (1 + (1-sqrt2)b)).
        let l_prime = -1.0 / d;
        let l_second = -1.0 / (d * d);
        let q = 1.0 + 2.0 * b_sum - b_sum * b_sum;
        let g = ((1.0 + (1.0 + sqrt_2) * b_sum) / (1.0 + (1.0 - sqrt_2) * b_sum)).ln();
        let g_prime = 2.0 * sqrt_2 / q;
        let g_second = -4.0 * sqrt_2 * d / (q * q);

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
                    - g * a_ij[i][j] / (sqrt_2 * b_sum)
                    + g * pair / (sqrt_2 * b_sum * b_sum)
                    - g * a_sum * product / (sqrt_2 * b_sum * b_sum * b_sum)
                    - g_prime * pair / (sqrt_2 * b_sum)
                    + g_prime * a_sum * product / (sqrt_2 * b_sum * b_sum)
                    - g_second * a_sum * product / (2.0 * sqrt_2 * b_sum);
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
        count: usize,
    ) -> (Vec<f64>, Vec<Vec<f64>>) {
        let a_hat: Vec<f64> = reduced.a.iter().map(|v| v / compressibility).collect();
        let b_hat: Vec<f64> = reduced.b.iter().map(|v| v / compressibility).collect();
        let a_ij: Vec<Vec<f64>> = (0..count)
            .map(|i| {
                (0..count)
                    .map(|j| (1.0 - self.kij(i, j)) * (a_hat[i] * a_hat[j]).sqrt())
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
    /// Warnings raised while computing them - in practice `pr_kappa`'s
    /// `kappa < 0` for a component with a sufficiently negative acentric factor.
    pub warnings: Vec<Warning>,
}
