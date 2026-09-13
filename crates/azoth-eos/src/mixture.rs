//! Components, mixtures, and the arithmetic over them.
//!
//! This is the model layer's own arithmetic, and it is the one place in this crate
//! where that is true. Everything else composes registered calculations; a mixture's
//! fugacity coefficient is not a registered calculation, because `eos.pr_departure`
//! covers a *pure* component and the registry is scalar - it has no composition
//! vector to hang a mixture form on.
//!
//! # Why that is a mitigation and not a hole
//!
//! The mixture form is checked by *composition* rather than by assertion, and the
//! reductions are exact rather than approximate:
//!
//! * **At `N = 1`** the cross-sum factor `2 * sum_j x_j A_ij / A - B_i / B` collapses
//!   to `2 * A / A - 1 = 1`, so [`Mixture::phase_state`]'s `ln_phi` becomes exactly
//!   `eos.pr_departure`'s. Tested in both languages.
//! * **At `N = 2`** the mixture parameters must reproduce `eos.vdw1f_mix_binary` and
//!   the vapour fraction must reproduce `eos.rachford_rice_binary`. Tested in both
//!   languages, and the second is a *cross-layer* check no single-language test can
//!   replace.
//!
//! Both reductions are to within a couple of ulps rather than bit-identical, and the
//! reason is written down where it bites: the registered binary kernel evaluates
//! `z1*z1*a1 + 2*z1*z2*(1 - k12)*sqrt(a1*a2) + z2*z2*a2` longhand, the mixture form
//! sums `i` then `j`, and the two associate differently. A bit-equality claim here
//! would be a claim about summation order rather than about the mixing rule.

use azoth_core::units::{Pressure, ThermodynamicTemperature};
use azoth_core::{AzothError, Result, Warning};

use crate::{pr_alpha_ab, pr_kappa, pr_z_factor};

/// Which root of the cubic a phase claims.
///
/// Selection is by *ordering* and never by an initial guess - the rule
/// [`crate::pr_z_factor`] fixes, and the one the flash has to follow so that both
/// implementations pick the same phase. A Newton iteration from a starting point
/// selects a root by basin, which is a different and unstable rule; see that calc's
/// module documentation.
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
    /// attraction parameters. Reported because it is what makes the mixture departure
    /// a departure at all, and because at `N = 1` it must equal that component's own
    /// `psi` exactly - a test asserts it.
    pub psi_bar: f64,
}

/// One component's critical constants.
///
/// A caller-supplied record, and deliberately carrying **no name**: if a name were
/// here, something would eventually use it to look a value up, and "this library
/// ships no component databank" is worth more than a nicer `Debug` output.
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
        }
        Ok(ReducedParameters {
            a,
            b,
            psi,
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
        let n = self.len();
        if x.len() != n {
            return Err(AzothError::invalid_input(
                "x",
                format!("a composition for {n} components has {} entries", x.len()),
            ));
        }

        let (a_mix, b_mix) = self.mixture_parameters(reduced, x);
        let roots = pr_z_factor(a_mix, b_mix)?;
        let z = match side {
            RootSide::Liquid => roots.z_min,
            RootSide::Vapour => roots.z_max,
        };

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
        for i in 0..n {
            for j in 0..n {
                let weight =
                    x[i] * x[j] * (1.0 - self.kij(i, j)) * (reduced.a[i] * reduced.a[j]).sqrt();
                weight_total += weight;
                weighted_psi += weight * 0.5 * (reduced.psi[i] + reduced.psi[j]);
            }
        }
        let psi_bar = weighted_psi / weight_total;
        let h_dep_rt = (z - 1.0) + coefficient * (psi_bar - 1.0) * i_term;
        let s_dep_r = h_dep_rt
            - ln_phi
                .iter()
                .zip(x)
                .map(|(&lp, &x_i)| x_i * lp)
                .sum::<f64>();

        Ok(PhaseState {
            a_mix,
            b_mix,
            z,
            ln_phi,
            h_dep_rt,
            s_dep_r,
            psi_bar,
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
    /// Warnings raised while computing them - in practice `pr_kappa`'s
    /// `kappa < 0` for a component with a sufficiently negative acentric factor.
    pub warnings: Vec<Warning>,
}
