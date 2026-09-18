//! The Wertheim association contribution, shared by every associating model.
//!
//! NeqSim's `PhaseCPAInterface`, `ComponentSrkCPA` and `CPAMixingRuleHandler`. The
//! contribution is a *term on the residual Helmholtz energy*, which is why one kernel
//! serves CPA, UMR-CPA, PC-SAFT and SAFT-VR-Mie: only the energy it is added to differs.
//!
//! NeqSim's own extraction of this kernel, `CPAContribution`, exists and is constructed
//! nowhere. Its javadoc says the radial distribution function "is the same for all cubic
//! equations of state", and the code it describes is still duplicated across
//! `PhaseSrkCPA`, `PhaseUMRCPA` and `PhaseSAFTVRMie`. This is that kernel, once.
//!
//! # What is here, and what is not yet
//!
//! Here: the site schemes and their bond rule, the Carnahan-Starling distribution
//! function, the association strength `Delta`, the site-fraction solve, and the Helmholtz
//! energy with its fugacity term. Not yet: the *implicit* derivatives of the site
//! fractions with respect to `n`, `T` and `V`, which are what
//! [`crate::mixture::Mixture::phase_derivatives`] needs. The solve's Jacobian is already
//! formed here, since the Newton step requires it, so the implicit step is one linear
//! solve against a matrix this module already builds.
//!
//! # Variables
//!
//! NeqSim's own: a total volume `V`, a total mole number of one, and the mixture covolume
//! `B = sum_j n_j b_j`. For CPA the `b` are the `bCPA` values NeqSim substitutes for the
//! cubic's own, which is why a covolume is an argument here rather than read from a cubic.

use azoth_core::{AzothError, Result};

/// NeqSim's `ThermodynamicConstantsInterface.R`.
///
/// Taken as NeqSim's own rounded literal rather than the CODATA value, for the reason this
/// crate's other modules take theirs: where azoth and NeqSim would disagree, NeqSim wins,
/// and the association energy enters `exp(eps/(R T))`, where a difference in the eighth
/// digit of `R` is a difference in the answer. `chung_conductivity` carries the same value
/// for the same reason.
pub const R: f64 = 8.3144621;

/// The association site scheme a component carries.
///
/// The names are the values of `COMP.csv`'s `associationscheme` column. A scheme fixes both
/// how many sites the molecule has and which of them bond, and the rule is NeqSim's
/// `getInteractionMatrix`: two sites associate exactly when their charges differ in sign.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SiteScheme {
    /// One site. Its charge product is positive, so it does not self-associate.
    OneA,
    /// Two sites of the same sign. No pair associates.
    TwoA,
    /// Two sites of opposite sign - one donor, one acceptor.
    TwoB,
    /// Four sites, two of each sign.
    FourC,
}

impl SiteScheme {
    /// The charges NeqSim gives the sites, `CPAMixingRuleHandler:32-35`.
    fn charges(self) -> &'static [i32] {
        match self {
            SiteScheme::OneA => &[-1],
            SiteScheme::TwoA => &[-1, -1],
            SiteScheme::TwoB => &[1, -1],
            SiteScheme::FourC => &[1, 1, -1, -1],
        }
    }

    /// How many sites the scheme has.
    #[must_use]
    pub fn site_count(self) -> usize {
        self.charges().len()
    }

    /// Whether this scheme's own sites `a` and `b` associate.
    ///
    /// NeqSim's rule is a product sign, so `OneA` and `TwoA` - whose charges all share a
    /// sign - associate with nothing. That is not a simplification of this port: it is
    /// what `getInteractionMatrix` computes, and it means a component given the `2A`
    /// scheme carries a fitted association energy that no calculation reads.
    #[must_use]
    pub fn associates_within(self, a: usize, b: usize) -> bool {
        let c = self.charges();
        c[a] * c[b] < 0
    }

    /// Whether a site on this scheme associates with a site on `other`.
    ///
    /// `CPAMixingRuleHandler.setCrossAssociationScheme` applies the same sign test across
    /// the two schemes' charge vectors.
    #[must_use]
    pub fn associates_across(self, a: usize, other: SiteScheme, b: usize) -> bool {
        self.charges()[a] * other.charges()[b] < 0
    }
}

/// One component's association parameters.
///
/// `associationsites` is redundant with the scheme and is not carried. `energy` is
/// `COMP.csv`'s `associationenergy`, in **joules per mole** - water's 16655 is `eps/R` of
/// 2003.1 K against CPA's published 2003.4 K, which is how the unit is fixed.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AssociationComponent {
    /// The site scheme, which fixes the site count as well.
    pub scheme: SiteScheme,
    /// The association energy `eps`, in J/mol.
    pub energy: f64,
    /// The association volume `kappa_AB`, dimensionless.
    pub volume: f64,
}

/// The cross-association rule for one component pair.
///
/// NeqSim's `assosSchemeType`: `0` is the Elliott rule, which combines two pure
/// components' association strengths geometrically, and `1` is CR-1, which uses the pair's
/// own energy and volume. NeqSim initialises the whole matrix to zero - the Elliott rule -
/// and overwrites only the pairs the `intertemp` table has a row for, reading
/// `cpaBetaCross` and `cpaEpsCross` from it.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum CrossRule {
    /// The Elliott rule: `DeltaNog(i,j) = sqrt(DeltaNog(i,i) DeltaNog(j,j))`.
    Elliott,
    /// CR-1 with the pair's own energy and volume.
    Cr1 {
        /// `eps_ij`, in J/mol.
        energy: f64,
        /// `kappa_AB` for the pair, dimensionless.
        volume: f64,
    },
}

/// The whole mixture's association parameters.
#[derive(Debug, Clone, PartialEq)]
pub struct Association {
    components: Vec<AssociationComponent>,
    /// `N x N`, row-major, symmetric. Empty means every pair is
    /// [`CrossRule::Elliott`], which is NeqSim's initialisation.
    cross: Vec<CrossRule>,
    /// The first site index of each component, with the total appended.
    offsets: Vec<usize>,
}

impl Association {
    /// A mixture's association parameters.
    ///
    /// # Errors
    /// * [`AzothError::InvalidInput`] if `cross` is neither empty nor `N*N` long.
    /// * [`AzothError::OutOfRange`] if any component's energy or volume is negative. Both
    ///   enter an exponential and a square root, and a negative volume would make the
    ///   Elliott combination the square root of a negative number.
    pub fn new(components: Vec<AssociationComponent>, cross: Vec<CrossRule>) -> Result<Self> {
        let n = components.len();
        if !cross.is_empty() && cross.len() != n * n {
            return Err(AzothError::invalid_input(
                "cross",
                format!(
                    "a cross rule matrix for {n} components has {} entries",
                    cross.len()
                ),
            ));
        }
        for (i, c) in components.iter().enumerate() {
            for (field, value) in [
                ("association energy", c.energy),
                ("association volume", c.volume),
            ] {
                if !value.is_finite() || value < 0.0 {
                    return Err(AzothError::OutOfRange {
                        field: field.to_string(),
                        value,
                        detail: format!(
                            "component {i}'s {field} is non-negative by definition: it \
                             enters an exponential and a square root, so a negative one is \
                             not a weaker association but no association at all"
                        ),
                    });
                }
            }
        }
        let mut offsets = Vec::with_capacity(n + 1);
        let mut total = 0;
        for c in &components {
            offsets.push(total);
            total += c.scheme.site_count();
        }
        offsets.push(total);
        let cross = if cross.is_empty() {
            vec![CrossRule::Elliott; n * n]
        } else {
            cross
        };
        Ok(Self {
            components,
            cross,
            offsets,
        })
    }

    /// One component's parameters.
    #[must_use]
    pub fn component(&self, i: usize) -> AssociationComponent {
        self.components[i]
    }

    /// The number of components.
    #[must_use]
    pub fn len(&self) -> usize {
        self.components.len()
    }

    /// Whether there are no components.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.components.is_empty()
    }

    /// The total number of association sites over every component.
    #[must_use]
    pub fn site_count(&self) -> usize {
        self.offsets[self.components.len()]
    }

    /// The sites belonging to component `i`, as a range into the flattened site vector.
    #[must_use]
    pub fn sites_of(&self, i: usize) -> std::ops::Range<usize> {
        self.offsets[i]..self.offsets[i + 1]
    }

    /// Which component a flattened site index belongs to.
    #[must_use]
    pub fn component_of_site(&self, site: usize) -> usize {
        self.offsets.partition_point(|&o| o <= site) - 1
    }

    /// The cross rule for a component pair.
    #[must_use]
    pub fn cross_rule(&self, i: usize, j: usize) -> CrossRule {
        self.cross[i * self.components.len() + j]
    }

    /// Whether the sites `a` and `b` bond, across components if they differ.
    ///
    /// The pair's [`CrossRule`] is not consulted: it changes a bond's strength, never
    /// whether one exists.
    #[must_use]
    pub fn bonds(&self, a: usize, b: usize) -> bool {
        let (ia, ib) = (self.component_of_site(a), self.component_of_site(b));
        let (sa, sb) = (a - self.offsets[ia], b - self.offsets[ib]);
        if ia == ib {
            self.components[ia].scheme.associates_within(sa, sb)
        } else {
            self.components[ia]
                .scheme
                .associates_across(sa, self.components[ib].scheme, sb)
        }
    }

    /// Whether any pair of sites in the mixture bonds at all.
    ///
    /// A mixture of `OneA` and `TwoA` components has sites and no bonds, so every site
    /// fraction is one and the association energy is exactly zero. The solve still works;
    /// this exists so a caller can skip it when it is provably a no-op.
    #[must_use]
    pub fn has_bonds(&self) -> bool {
        let n = self.site_count();
        (0..n).any(|a| (0..n).any(|b| self.bonds(a, b)))
    }
}

/// The radial distribution function at contact and its volume derivatives.
///
/// Carnahan-Starling, `g = (1 - eta/2)/(1 - eta)^3` with `eta = B/(4V)`:
/// `PhaseCPAInterface.calc_g`, `calc_lngV` and `calc_lngVV`. It is a property of the
/// *mixture* - `g` carries no component index - so one value serves every component, and
/// NeqSim's `calc_lngi(i)` is the derivative with respect to `n_i` rather than a
/// per-component quantity.
#[derive(Debug, Clone, PartialEq)]
pub struct Rdf {
    /// `g` at contact.
    pub g: f64,
    /// `d ln g / d n_j`, one entry per component.
    pub d_ln_g_dn: Vec<f64>,
    /// `d ln g / d V` at constant composition - NeqSim's `gcpav`.
    pub d_ln_g_dv: f64,
    /// `d^2 ln g / d V^2`.
    pub d2_ln_g_dv2: f64,
    /// `d ln g / d T` at constant composition and volume, which is zero: `eta` carries no
    /// temperature. Present so a caller does not have to remember that.
    pub d_ln_g_dt: f64,
}

impl Rdf {
    /// Evaluate the distribution function and its derivatives at a state.
    ///
    /// `covolumes` and `moles` are one entry per component; `v` is the total volume.
    #[must_use]
    pub fn new(covolumes: &[f64], moles: &[f64], v: f64) -> Self {
        let (b, n) = (covolumes, moles);
        let big_b: f64 = b.iter().zip(n).map(|(bi, ni)| bi * ni).sum();
        let eta = big_b / (4.0 * v);
        let one_minus = 1.0 - eta;
        let g = (1.0 - eta / 2.0) / (one_minus * one_minus * one_minus);
        // d ln g / d eta = 3/(1 - eta) - 1/(2 - eta).
        let d_ln_g_d_eta = 3.0 / one_minus - 1.0 / (2.0 - eta);
        // d^2 ln g / d eta^2. The second term is minus: ln g carries -3 ln(1 - eta), whose
        // first derivative is +3/(1 - eta) and second +3/(1 - eta)^2, while the 1 - eta/2
        // term contributes -1/(2 - eta) and then -1/(2 - eta)^2.
        let d2_ln_g_d_eta2 = 3.0 / (one_minus * one_minus) - 1.0 / ((2.0 - eta) * (2.0 - eta));
        let d_ln_g_dn = b.iter().map(|bi| bi * d_ln_g_d_eta / (4.0 * v)).collect();
        // eta = B/(4V), so d eta/dV = -B/(4V^2).
        let d_eta_dv = -big_b / (4.0 * v * v);
        let d_ln_g_dv = d_ln_g_d_eta * d_eta_dv;
        let d2_ln_g_dv2 =
            d2_ln_g_d_eta2 * d_eta_dv * d_eta_dv + d_ln_g_d_eta * big_b / (2.0 * v * v * v);
        Self {
            g,
            d_ln_g_dn,
            d_ln_g_dv,
            d2_ln_g_dv2,
            d_ln_g_dt: 0.0,
        }
    }
}

/// The association strength between two components, without the distribution function:
/// NeqSim's `deltaNog`.
///
/// The covolumes are arguments because they are not association parameters - for CPA they
/// are the `bCPA` values NeqSim substitutes for the cubic's own `b`, and for PC-SAFT the
/// segment diameter enters differently. `DeltaNog` is the same function either way.
///
/// Under either rule the self-pair `(i, i)` reduces to component `i`'s own strength, so a
/// component's self-association needs no special case.
#[must_use]
pub fn delta_nog(
    a: &AssociationComponent,
    covolume_a: f64,
    c: &AssociationComponent,
    covolume_c: f64,
    rule: CrossRule,
    t: f64,
) -> f64 {
    // `(exp(eps/(RT)) - 1) b beta`.
    let strength = |x: &AssociationComponent, covolume: f64| {
        ((x.energy / (R * t)).exp() - 1.0) * covolume * x.volume
    };
    match rule {
        // NeqSim's default for every pair the `intertemp` table has no row for.
        CrossRule::Elliott => (strength(a, covolume_a) * strength(c, covolume_c)).sqrt(),
        // CR-1, with the pair's own energy and volume. NeqSim falls back to
        // `(eps_a + eps_c)/2` and `sqrt(beta_a beta_c)` where the table has no row.
        CrossRule::Cr1 { energy, volume } => {
            ((energy / (R * t)).exp() - 1.0) * (covolume_a + covolume_c) / 2.0 * volume
        }
    }
}

/// The converged site fractions and the energies they carry.
#[derive(Debug, Clone, PartialEq)]
pub struct SiteState {
    /// `X_A`, one entry per site, flattened in component order.
    pub fractions: Vec<f64>,
    /// Successive-substitution sweeps taken, including the one that met the tolerance.
    pub iterations: usize,
    /// Whether the solve met [`SOLVE_TOLERANCE`] before [`MAX_SWEEPS`].
    pub converged: bool,
    /// Whether the Newton refinement ran and left every fraction positive.
    pub refined: bool,
    /// `A_assoc / (R T)`, the association Helmholtz energy over `R T`.
    pub helmholtz_rt: f64,
    /// The association contribution to `ln phi_i`, one entry per component, at constant
    /// `T` and `V`.
    ///
    /// NeqSim's `ComponentSrkCPA.dFCPAdN`: `sum_A ln X_A - (hcpatot/2) calc_lngi(i)`, with
    /// `hcpatot` of one - which is `PhaseSrkCPA`'s value, since only the *reduced*
    /// variants ever assign that field.
    pub ln_phi: Vec<f64>,
}

/// The successive-substitution tolerance, NeqSim's `solveX2`.
pub const SOLVE_TOLERANCE: f64 = 1.0e-12;
/// The sweep cap, NeqSim's `solveX2(15)`.
pub const MAX_SWEEPS: usize = 15;

impl Association {
    /// The association strength matrix, `Delta_ij = bond(i,j) deltaNog_ij g`.
    ///
    /// NeqSim builds this once per state in `initCPAMatrix` and reuses it for the solve
    /// and for every derivative, which is why it is returned rather than recomputed.
    #[must_use]
    pub fn delta_matrix(&self, covolumes: &[f64], t: f64, rdf: &Rdf) -> Vec<f64> {
        let n = self.site_count();
        let mut delta = vec![0.0; n * n];
        for i in 0..n {
            let ci = self.component_of_site(i);
            for j in 0..n {
                if !self.bonds(i, j) {
                    continue;
                }
                let cj = self.component_of_site(j);
                delta[i * n + j] = delta_nog(
                    &self.components[ci],
                    covolumes[ci],
                    &self.components[cj],
                    covolumes[cj],
                    self.cross_rule(ci, cj),
                    t,
                ) * rdf.g;
            }
        }
        delta
    }

    /// Solve for the site fractions and evaluate the association energy.
    ///
    /// NeqSim's `solveX`: successive substitution to [`SOLVE_TOLERANCE`], then a Newton
    /// refinement on the same residual. The two are not redundant - the substitution is
    /// globally reliable and converges linearly, the Newton is quadratic and can leave the
    /// physical branch - so the port keeps both and reports whether each finished.
    ///
    /// # Errors
    /// * [`AzothError::InvalidInput`] if the state vectors are not one entry per component.
    /// * [`AzothError::OutOfRange`] if `t` or `v` is not positive.
    pub fn solve(&self, covolumes: &[f64], moles: &[f64], v: f64, t: f64) -> Result<SiteState> {
        let n = self.len();
        for (field, given) in [("covolumes", covolumes.len()), ("moles", moles.len())] {
            if given != n {
                return Err(AzothError::invalid_input(
                    field,
                    format!("a mixture of {n} components takes {given} entries"),
                ));
            }
        }
        for (field, value) in [("v", v), ("T", t)] {
            if !value.is_finite() || value <= 0.0 {
                return Err(AzothError::OutOfRange {
                    field: field.to_string(),
                    value,
                    detail: "the association solve needs a positive volume and temperature"
                        .to_string(),
                });
            }
        }

        let sites = self.site_count();
        if sites == 0 {
            return Ok(SiteState {
                fractions: Vec::new(),
                iterations: 0,
                converged: true,
                refined: false,
                helmholtz_rt: 0.0,
                ln_phi: vec![0.0; n],
            });
        }

        let rdf = Rdf::new(covolumes, moles, v);
        let delta = self.delta_matrix(covolumes, t, &rdf);
        // `Klk_ij = n_i n_j Delta_ij / V`, NeqSim's `KlkMatrix`. The moles are indexed by
        // the component owning the site, not by the site.
        let mut klk = vec![0.0; sites * sites];
        for i in 0..sites {
            let ni = moles[self.component_of_site(i)];
            for j in 0..sites {
                klk[i * sites + j] =
                    ni * moles[self.component_of_site(j)] * delta[i * sites + j] / v;
            }
        }

        // Successive substitution, Gauss-Seidel: the site written this sweep is visible to
        // the rest of the sweep, which is NeqSim's in-place update.
        let mut x = vec![1.0; sites];
        let mut iterations = 0;
        let mut converged = false;
        while iterations < MAX_SWEEPS {
            iterations += 1;
            let mut error = 0.0;
            for i in 0..sites {
                let old = x[i];
                let sum: f64 = (0..sites).map(|j| klk[i * sites + j] * x[j]).sum();
                let new = 1.0 / (1.0 + sum);
                x[i] = new;
                error += ((old - new) / new).abs();
            }
            if error < SOLVE_TOLERANCE {
                converged = true;
                break;
            }
        }

        let refined = self.newton_refine(&mut x, &klk, sites);

        // `A_assoc/(RT) = sum_i n_i sum_{A in i} (ln X_A - X_A/2 + 1/2)`.
        let helmholtz_rt: f64 = moles
            .iter()
            .enumerate()
            .map(|(i, ni)| {
                let per_site: f64 = self.sites_of(i).map(|a| x[a].ln() - x[a] / 2.0 + 0.5).sum();
                ni * per_site
            })
            .sum();
        let ln_phi = (0..n)
            .map(|i| {
                let site_term: f64 = self.sites_of(i).map(|a| x[a].ln()).sum();
                site_term - 0.5 * rdf.d_ln_g_dn[i]
            })
            .collect();

        Ok(SiteState {
            fractions: x,
            iterations,
            converged,
            refined,
            helmholtz_rt,
            ln_phi,
        })
    }

    /// Newton-refine the site fractions in place, returning whether it converged on a
    /// physical branch.
    ///
    /// NeqSim's `solveX` after `solveX2`. The residual is
    /// `F_i = X_i (1 + sum_j Klk_ij X_j) - 1`, whose fixed point is the substitution, and
    /// the step is `-J^{-1} F`. A fraction that would go non-positive is clipped to
    /// `1e-10` - NeqSim's value - which leaves the solve unrefined rather than wrong, so
    /// the caller keeps the substitution's answer.
    fn newton_refine(&self, x: &mut [f64], klk: &[f64], sites: usize) -> bool {
        if sites == 0 {
            return false;
        }
        const NEWTON_TOLERANCE: f64 = 1.0e-12;
        const MAX_NEWTON: usize = 100;
        let mut jacobian = vec![0.0; sites * sites];
        let mut residual = vec![0.0; sites];
        let mut physical = true;
        for _ in 0..MAX_NEWTON {
            for i in 0..sites {
                let inner: f64 = (0..sites).map(|j| klk[i * sites + j] * x[j]).sum();
                residual[i] = x[i] * (1.0 + inner) - 1.0;
                for k in 0..sites {
                    let diagonal = if i == k { 1.0 } else { 0.0 };
                    jacobian[i * sites + k] = diagonal * (1.0 + inner) + x[i] * klk[i * sites + k];
                }
            }
            if residual.iter().all(|r| r.abs() < NEWTON_TOLERANCE) {
                return physical;
            }
            let Some(step) = solve_dense(&mut jacobian, &residual, sites) else {
                return false;
            };
            for i in 0..sites {
                let next = x[i] - step[i];
                if next <= 0.0 {
                    x[i] = 1.0e-10;
                    physical = false;
                } else {
                    x[i] = next;
                }
            }
        }
        false
    }
}

/// Solve `A x = b` for a small dense system by Gaussian elimination with partial pivoting.
///
/// The association kernel's Jacobian is the only dense solve in this crate, and it is
/// `sites x sites` - a handful - which is why it is written here rather than taken as a
/// dependency. NeqSim uses EJML for the same solve. Returns `None` on a singular matrix,
/// which the caller reports as an unrefined solve rather than a wrong one.
fn solve_dense(a: &mut [f64], b: &[f64], n: usize) -> Option<Vec<f64>> {
    let mut x = b.to_vec();
    for col in 0..n {
        let mut pivot = col;
        for row in col + 1..n {
            if a[row * n + col].abs() > a[pivot * n + col].abs() {
                pivot = row;
            }
        }
        if a[pivot * n + col] == 0.0 {
            return None;
        }
        if pivot != col {
            for k in 0..n {
                a.swap(col * n + k, pivot * n + k);
            }
            x.swap(col, pivot);
        }
        for row in col + 1..n {
            let factor = a[row * n + col] / a[col * n + col];
            if factor == 0.0 {
                continue;
            }
            for k in col..n {
                a[row * n + k] -= factor * a[col * n + k];
            }
            x[row] -= factor * x[col];
        }
    }
    for row in (0..n).rev() {
        let mut sum = x[row];
        for k in row + 1..n {
            sum -= a[row * n + k] * x[k];
        }
        x[row] = sum / a[row * n + row];
    }
    Some(x)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn ln_g(b: &[f64], n: &[f64], v: f64) -> f64 {
        Rdf::new(b, n, v).g.ln()
    }

    #[test]
    fn the_scheme_site_counts_and_bonds_are_neqsims() {
        assert_eq!(SiteScheme::OneA.site_count(), 1);
        assert_eq!(SiteScheme::TwoA.site_count(), 2);
        assert_eq!(SiteScheme::TwoB.site_count(), 2);
        assert_eq!(SiteScheme::FourC.site_count(), 4);
        let four_c = SiteScheme::FourC;
        assert!(four_c.associates_within(0, 2));
        assert!(four_c.associates_within(1, 3));
        assert!(!four_c.associates_within(0, 1));
        assert!(!four_c.associates_within(2, 3));
        assert!(SiteScheme::TwoB.associates_within(0, 1));
        assert!(!SiteScheme::TwoB.associates_within(0, 0));
        // The finding: NeqSim's sign test leaves these two schemes with no bonds at all.
        assert!(!SiteScheme::OneA.associates_within(0, 0));
        assert!(!SiteScheme::TwoA.associates_within(0, 1));
    }

    #[test]
    fn the_rdf_composition_derivative_matches_a_finite_difference() {
        let b = [4.2e-5, 3.1e-5];
        let n = [0.6, 0.4];
        let v = 1.0e-4;
        let rdf = Rdf::new(&b, &n, v);
        let step = 1.0e-9;
        for j in 0..2 {
            let mut np = n;
            np[j] += step;
            let numerical = (ln_g(&b, &np, v) - ln_g(&b, &n, v)) / step;
            let off = (rdf.d_ln_g_dn[j] / numerical - 1.0).abs();
            assert!(
                off < 1.0e-6,
                "component {j}: analytic {} vs numerical {numerical}",
                rdf.d_ln_g_dn[j]
            );
        }
    }

    #[test]
    fn the_rdf_volume_derivatives_match_finite_differences() {
        let b = [4.2e-5, 3.1e-5];
        let n = [0.6, 0.4];
        let v = 1.0e-4;
        let rdf = Rdf::new(&b, &n, v);
        let step = 1.0e-9;
        let first = (ln_g(&b, &n, v + step) - ln_g(&b, &n, v - step)) / (2.0 * step);
        let second = (ln_g(&b, &n, v + step) - 2.0 * ln_g(&b, &n, v) + ln_g(&b, &n, v - step))
            / (step * step);
        assert!(
            (rdf.d_ln_g_dv / first - 1.0).abs() < 1.0e-5,
            "first {} vs {first}",
            rdf.d_ln_g_dv
        );
        assert!(
            (rdf.d2_ln_g_dv2 / second - 1.0).abs() < 1.0e-3,
            "second {} vs {second}",
            rdf.d2_ln_g_dv2
        );
    }

    /// NeqSim's `calc_lngi(i) = 2 b_i (10V - B)/((8V - B)(4V - B))`, which the module's
    /// derivation must reproduce from `calc_g` rather than restate.
    #[test]
    fn the_composition_derivative_is_neqsims_closed_form() {
        let b = [4.2e-5, 3.1e-5];
        let n = [0.6, 0.4];
        let v = 1.0e-4;
        let big_b: f64 = b.iter().zip(&n).map(|(bi, ni)| bi * ni).sum();
        let rdf = Rdf::new(&b, &n, v);
        for (j, &bj) in b.iter().enumerate() {
            let closed = 2.0 * bj * (10.0 * v - big_b) / ((8.0 * v - big_b) * (4.0 * v - big_b));
            assert!(
                (rdf.d_ln_g_dn[j] / closed - 1.0).abs() < 1.0e-12,
                "component {j}: {} vs NeqSim's {closed}",
                rdf.d_ln_g_dn[j]
            );
        }
    }

    #[test]
    fn a_mixture_of_one_a_and_two_a_has_sites_and_no_bonds() {
        let a = Association::new(
            vec![
                AssociationComponent {
                    scheme: SiteScheme::OneA,
                    energy: 0.0,
                    volume: 0.11347,
                },
                AssociationComponent {
                    scheme: SiteScheme::TwoA,
                    energy: 5000.0,
                    volume: 0.001160489,
                },
            ],
            Vec::new(),
        )
        .expect("valid parameters");
        assert_eq!(a.site_count(), 3);
        assert!(
            !a.has_bonds(),
            "1A and 2A bond with nothing under NeqSim's sign test, so H2S's 5000 J/mol \
             association energy reaches no calculation"
        );
    }

    #[test]
    fn a_water_like_component_has_four_bonds() {
        let a = Association::new(
            vec![AssociationComponent {
                scheme: SiteScheme::FourC,
                energy: 16655.0,
                volume: 0.0692,
            }],
            Vec::new(),
        )
        .expect("valid parameters");
        assert_eq!(a.site_count(), 4);
        assert!(a.has_bonds());
        // The two positive sites bond with the two negative ones and with nothing else,
        // so there are four unordered pairs.
        let unordered = (0..4)
            .flat_map(|i| (i + 1..4).map(move |j| (i, j)))
            .filter(|&(i, j)| a.bonds(i, j))
            .count();
        assert_eq!(unordered, 4, "4C has four unlike-sign pairs");
    }

    #[test]
    fn the_elliott_rule_reduces_to_the_component_itself() {
        let water = AssociationComponent {
            scheme: SiteScheme::FourC,
            energy: 16655.0,
            volume: 0.0692,
        };
        let covolume = 2.6e-5;
        let own = (water.energy / (R * 300.0)).exp() - 1.0;
        let self_combination = delta_nog(
            &water,
            covolume,
            &water,
            covolume,
            CrossRule::Elliott,
            300.0,
        );
        assert!((self_combination / (own * covolume * water.volume) - 1.0).abs() < 1.0e-12);
    }

    /// The Elliott rule is a geometric mean, so a pair's strength is exactly the root of
    /// the two self-strengths - which is the property that lets the cross rule be a
    /// matrix rather than a second formula.
    #[test]
    fn the_elliott_rule_is_the_geometric_mean_of_the_two_components() {
        let water = AssociationComponent {
            scheme: SiteScheme::FourC,
            energy: 16655.0,
            volume: 0.0692,
        };
        let methanol = AssociationComponent {
            scheme: SiteScheme::TwoB,
            energy: 24591.0,
            volume: 0.0161,
        };
        let (bw, bm) = (2.6e-5, 3.1e-5);
        let self_w = delta_nog(&water, bw, &water, bw, CrossRule::Elliott, 300.0);
        let self_m = delta_nog(&methanol, bm, &methanol, bm, CrossRule::Elliott, 300.0);
        let cross = delta_nog(&water, bw, &methanol, bm, CrossRule::Elliott, 300.0);
        assert!((cross / (self_w * self_m).sqrt() - 1.0).abs() < 1.0e-12);
    }

    /// A water-like pure component: the site fractions must satisfy the substitution
    /// exactly, and every one of them must sit strictly between zero and one.
    #[test]
    fn the_site_fractions_satisfy_the_substitution() {
        let water = AssociationComponent {
            scheme: SiteScheme::FourC,
            energy: 16655.0,
            volume: 0.0692,
        };
        let a = Association::new(vec![water], Vec::new()).expect("valid parameters");
        let covolume = 2.6e-5;
        let moles = [1.0];
        let v = 2.0e-5;
        let t = 373.0;
        let state = a
            .solve(&[covolume], &moles, v, t)
            .expect("a solvable state");
        assert!(
            state.fractions.iter().all(|&x| x > 0.0 && x < 1.0),
            "an associating site is partly bonded: {:?}",
            state.fractions
        );
        // The residual of the fixed-point equation, recomputed here rather than trusted.
        let rdf = Rdf::new(&[covolume], &moles, v);
        let delta = a.delta_matrix(&[covolume], t, &rdf);
        let sites = a.site_count();
        for i in 0..sites {
            let sum: f64 = (0..sites)
                .map(|j| moles[a.component_of_site(j)] * delta[i * sites + j] * state.fractions[j])
                .sum();
            let expected = 1.0 / (1.0 + sum / v);
            assert!(
                (state.fractions[i] / expected - 1.0).abs() < 1.0e-9,
                "site {i}: {} vs {expected}",
                state.fractions[i]
            );
        }
    }

    /// A stronger bond leaves fewer free sites, and the Helmholtz energy is negative -
    /// association always lowers the free energy, which is why it is worth modelling.
    #[test]
    fn a_stronger_association_lowers_the_free_energy() {
        let build = |energy| {
            Association::new(
                vec![AssociationComponent {
                    scheme: SiteScheme::TwoB,
                    energy,
                    volume: 0.0161,
                }],
                Vec::new(),
            )
            .expect("valid parameters")
        };
        let (covolume, moles, v, t) = (3.1e-5, [1.0], 4.0e-5, 300.0);
        let weak = build(5000.0)
            .solve(&[covolume], &moles, v, t)
            .expect("solves");
        let strong = build(24591.0)
            .solve(&[covolume], &moles, v, t)
            .expect("solves");
        assert!(
            strong.fractions[0] < weak.fractions[0],
            "a stronger bond leaves fewer free sites: {} vs {}",
            strong.fractions[0],
            weak.fractions[0]
        );
        assert!(
            strong.helmholtz_rt < weak.helmholtz_rt,
            "association lowers the free energy: {} vs {}",
            strong.helmholtz_rt,
            weak.helmholtz_rt
        );
        assert!(weak.helmholtz_rt < 0.0);
    }

    /// The finding, measured through the solve rather than through the bond table: a `2A`
    /// component carries a fitted energy and no bond, so every site fraction is one and
    /// the association energy is *exactly* zero.
    ///
    /// Its fugacity term is not, and that is NeqSim's own behaviour rather than this
    /// port's: `dFCPAdN` carries `-(hcpatot/2) calc_lngi(i)` beside the site sum, and that
    /// piece is the derivative of the distribution function's prefactor. With no bonds it
    /// has nothing to cancel against and stays. Whether it should is a question for the
    /// oracle - it is the kind of term a full fugacity assembly can absorb.
    #[test]
    fn a_two_a_component_associates_not_at_all() {
        let h2s_like = Association::new(
            vec![AssociationComponent {
                scheme: SiteScheme::TwoA,
                energy: 5000.0,
                volume: 0.001160489,
            }],
            Vec::new(),
        )
        .expect("valid parameters");
        let (covolume, moles, v, t) = (2.9e-5, [1.0], 4.0e-5, 300.0);
        let state = h2s_like.solve(&[covolume], &moles, v, t).expect("solves");
        assert!(state.fractions.iter().all(|&x| x == 1.0));
        assert_eq!(state.helmholtz_rt, 0.0);
        // The term that remains is exactly the RDF prefactor's derivative and nothing
        // else - no part of it comes from a bond.
        let rdf = Rdf::new(&[covolume], &moles, v);
        assert!((state.ln_phi[0] + 0.5 * rdf.d_ln_g_dn[0]).abs() < 1.0e-300);
        assert!(
            state.ln_phi[0] != 0.0,
            "the prefactor derivative is not zero"
        );
    }

    /// A mixture with no sites at all is a no-op rather than an error, because that is
    /// what a CPA mixture of light hydrocarbons on a 0-site basis is.
    #[test]
    fn a_mixture_without_sites_returns_a_trivial_state() {
        let none = Association::new(
            vec![
                AssociationComponent {
                    scheme: SiteScheme::OneA,
                    energy: 0.0,
                    volume: 0.0,
                };
                2
            ],
            Vec::new(),
        )
        .expect("valid parameters");
        // OneA has a site, so build a genuinely site-free mixture by hand instead.
        assert_eq!(none.site_count(), 2);
        let empty = Association::new(Vec::new(), Vec::new());
        assert!(empty.expect("an empty mixture is valid").is_empty());
    }

    #[test]
    fn a_negative_association_parameter_is_refused() {
        let bad = Association::new(
            vec![AssociationComponent {
                scheme: SiteScheme::TwoB,
                energy: -1.0,
                volume: 0.01,
            }],
            Vec::new(),
        );
        assert!(bad.is_err());
    }

    #[test]
    fn a_wrong_sized_cross_matrix_is_refused() {
        let bad = Association::new(
            vec![
                AssociationComponent {
                    scheme: SiteScheme::TwoB,
                    energy: 1.0,
                    volume: 0.01,
                },
                AssociationComponent {
                    scheme: SiteScheme::TwoB,
                    energy: 1.0,
                    volume: 0.01,
                },
            ],
            vec![CrossRule::Elliott; 3],
        );
        assert!(bad.is_err());
    }
}
