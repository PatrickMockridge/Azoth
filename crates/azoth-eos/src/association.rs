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
//! # What is here
//!
//! The site schemes and their bond rule, the Carnahan-Starling distribution function, the
//! association strength `Delta`, the site-fraction solve, the Helmholtz energy with its
//! fugacity term, and the *implicit* derivatives of the site fractions with respect to
//! `n`, `T` and `V` - which are what [`crate::mixture::Mixture::phase_derivatives`] needs,
//! and which are implicit because the site fractions solve a nonlinear system rather than
//! being a formula.
//!
//! The tests come in two kinds and both are needed. One kind differentiates the solve's own
//! output, which checks the kernel against itself. The other compares against NeqSim's
//! numbers, printed by `validation/neqsim/CpaProbe.java`, which is the only check that the
//! kernel is NeqSim's and not merely self-consistent - and it is the one that settled
//! `hCPA`. Neither catches what the other does: the scaling error in the fixed point was
//! invisible to the second kind until the first existed to fail, and the `hCPA` error was
//! invisible to the first kind entirely.
//!
//! [`crate::mixture::Mixture`] consumes it: `with_association` turns the contribution on,
//! `reduced_parameters` substitutes the fitted `a` and `b`, and `phase_state_at` adds the
//! fugacity and enthalpy terms.
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
    /// No sites: a component the table gives no scheme, and the stand-in a mixture
    /// carries for every non-associating member so the site offsets line up.
    ///
    /// Not a value `COMP.csv` holds - its marker for the same thing is `0`, which
    /// [`SiteScheme::from_databank_name`] answers `None` to. This is the kernel's own
    /// spelling of it, because a mixture needs one entry per component.
    NonAssociating,
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
            SiteScheme::NonAssociating => &[],
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

    /// The scheme a databank name denotes.
    ///
    /// `None` for `"0"`, which is the table's marker for a component with no scheme at all,
    /// and for any name this library does not carry. A site *count* cannot stand in for
    /// this: the table carries `1A` at zero sites and both `2A` and `2B` at two, and
    /// NeqSim's `setAssociationScheme` switches on the name.
    #[must_use]
    pub fn from_databank_name(name: &str) -> Option<Self> {
        match name.trim() {
            "1A" => Some(Self::OneA),
            "2A" => Some(Self::TwoA),
            "2B" => Some(Self::TwoB),
            "4C" => Some(Self::FourC),
            _ => None,
        }
    }

    /// Whether this scheme's sites bond with each other at all.
    ///
    /// False for [`SiteScheme::OneA`] and [`SiteScheme::TwoA`], whose charge vectors carry
    /// one sign, so NeqSim's product test is positive for every pair and the interaction
    /// matrix is all zeros. Such a component never self-associates - though it still
    /// *cross*-associates with an oppositely-charged partner, which is why this is a
    /// property of the scheme and not a verdict on the component.
    #[must_use]
    pub fn self_bonds(self) -> bool {
        let n = self.site_count();
        (0..n).any(|a| (0..n).any(|b| self.associates_within(a, b)))
    }
}

/// The association parameters a databank row carries, for one component.
///
/// The two fitted cubic sets are in **NeqSim's internal scale** — `a` in
/// `Pa m**6/mol**2 x 1e5` and `b` in `m**3/mol x 1e5`, from `Component.java:526-531` —
/// because the compiled table carries what the source table carries. The conversion belongs
/// with the model that reads them, and the manifest marks them `neqsim-internal` rather than
/// the `dimensionless` it used to claim. Water is the check: `b_srk` is 1.4515 internal,
/// i.e. `1.4515e-5 m**3/mol`, against the cubic's own `2.11e-5`.
#[derive(Debug, Clone, PartialEq)]
pub struct AssociationRecord {
    /// The site scheme, which the site count cannot reconstruct.
    pub scheme: SiteScheme,
    /// The site count the table states. Carried beside the scheme because the two disagree
    /// upstream - `1A` appears at zero sites - and a silent resolution would hide that.
    pub sites: u32,
    /// The association energy `eps`, in J/mol.
    pub energy: f64,
    /// `kappa_AB` for the SRK family.
    pub volume_srk: f64,
    /// The fitted attraction for SRK-CPA, in NeqSim's internal scale.
    pub a_srk: f64,
    /// The fitted covolume for SRK-CPA, in NeqSim's internal scale.
    pub b_srk: f64,
    /// The SRK alpha correlation's `m`.
    pub m_srk: f64,
    /// `kappa_AB` for the PR family.
    pub volume_pr: f64,
    /// The fitted attraction for PR-CPA, in NeqSim's internal scale.
    pub a_pr: f64,
    /// The fitted covolume for PR-CPA, in NeqSim's internal scale.
    pub b_pr: f64,
    /// The PR alpha correlation's `m`.
    pub m_pr: f64,
    /// `racketZCPA`, the Rackett compressibility NeqSim's CPA volume correction reads.
    pub racket_z: f64,
    /// `volcorrCPA_T`, the CPA volume-translation coefficient.
    ///
    /// NeqSim's `SystemSrkCPA` calls `useVolumeCorrection(true)`, and
    /// `ComponentSrk.getVolumeCorrection` is `0.40768 (0.29441 - Z_RA) R Tc/Pc`, so a
    /// component's root carries a translation this library does not yet apply. It is
    /// **not** the cause of the root gap recorded in `tests/databank.rs` - for water the
    /// shift is `-2.5148594e-5` internal, two orders of magnitude short - but it is real
    /// and it is carried here for when it is applied.
    pub volume_correction: f64,
}

/// The factor from NeqSim's internal `a` and `b` to SI, `Component.java:526-531`.
pub const NEQSIM_INTERNAL_TO_SI: f64 = 1.0e-5;

/// Which cubic family's fitted parameter set to read.
///
/// A record carries both because the two differ and neither is derived from the other.
/// Water's `kappa_AB` is 0.0692 for SRK against 0.046473789 for PR, and its fitted
/// covolume is 1.4515 against 1.456360879 - so choosing is a selection, not a
/// conversion, and reading the wrong family is a different fluid.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AssociationCubic {
    /// The SRK family - `aCPA_SRK`, `bCPA_SRK`, `mCPA_SRK`, `associationboundingvolume_SRK`.
    Srk,
    /// The PR family - the `_PR` columns.
    Pr,
}

/// A component that carries no association, for a mixture's shape.
///
/// A mixture's [`Association`] needs one entry per component so the site offsets line
/// up, and a non-associating member contributes none.
pub const NON_ASSOCIATING: AssociationComponent = AssociationComponent {
    scheme: SiteScheme::NonAssociating,
    energy: 0.0,
    volume: 0.0,
};

impl AssociationCubic {
    /// The family a cubic's geometry belongs to.
    ///
    /// By geometry and not by name: `Cubic::Tst` shares Peng-Robinson's `omega` and
    /// `delta`, and `Cubic::Rk` shares Soave's, so each reads the family it is shaped
    /// like. Only `Srk` and `Pr` are families NeqSim builds a CPA model on - the other
    /// two would be this library's own combination - which is why every cubic maps and
    /// none is refused.
    #[must_use]
    pub fn of(cubic: crate::cubic::Cubic) -> Self {
        match cubic {
            crate::cubic::Cubic::Srk | crate::cubic::Cubic::Rk => Self::Srk,
            crate::cubic::Cubic::Pr | crate::cubic::Cubic::Tst => Self::Pr,
        }
    }
}

impl AssociationRecord {
    /// Whether this record carries a usable fitted set at one cubic family.
    ///
    /// **NeqSim's own guard, and load-bearing.** `ComponentSrkCPA` substitutes the
    /// fitted values only `if (Math.abs(aCPA) > 1e-6)`, in its internal units - and the
    /// table has rows where that matters in both directions: CO2 names the `2A` scheme
    /// with every fitted value zero, and H2S and benzene carry an SRK set and a PR set
    /// that is entirely zero. Substituting unconditionally gives those a covolume of
    /// zero, which makes the reduced pressure `NaN` rather than a different answer.
    #[must_use]
    pub fn has_fitted_set(&self, cubic: AssociationCubic) -> bool {
        let internal = match cubic {
            AssociationCubic::Srk => self.a_srk,
            AssociationCubic::Pr => self.a_pr,
        };
        internal.abs() > 1.0e-6
    }
    /// The kernel's per-component record at one cubic family.
    ///
    /// This is the bridge from the databank's shape to the kernel's, and the only place
    /// the cubic-specific selection happens.
    #[must_use]
    pub fn at(&self, cubic: AssociationCubic) -> AssociationComponent {
        AssociationComponent {
            scheme: self.scheme,
            energy: self.energy,
            volume: match cubic {
                AssociationCubic::Srk => self.volume_srk,
                AssociationCubic::Pr => self.volume_pr,
            },
        }
    }

    /// The fitted attraction `a0` at one family, in SI - `Pa m**6/mol**2`.
    #[must_use]
    pub fn attraction(&self, cubic: AssociationCubic) -> f64 {
        let internal = match cubic {
            AssociationCubic::Srk => self.a_srk,
            AssociationCubic::Pr => self.a_pr,
        };
        internal * NEQSIM_INTERNAL_TO_SI
    }

    /// The fitted covolume at one family, in SI - `m**3/mol`.
    ///
    /// **Not the cubic's own.** NeqSim's `ComponentSrkCPA` substitutes this for the
    /// `b` a cubic derives from `Tc` and `Pc`, and the difference is large: water's is
    /// `1.4515e-5` against `0.08664 R Tc/Pc`'s `2.11e-5`.
    #[must_use]
    pub fn covolume(&self, cubic: AssociationCubic) -> f64 {
        let internal = match cubic {
            AssociationCubic::Srk => self.b_srk,
            AssociationCubic::Pr => self.b_pr,
        };
        internal * NEQSIM_INTERNAL_TO_SI
    }

    /// The fitted Soave alpha coefficient `m` at one family.
    ///
    /// NeqSim's `ComponentSrkCPA` calls `getAttractiveTerm().setm(mCPA)`, so the alpha
    /// function is Soave's form with a *fitted* coefficient rather than the one
    /// `0.480 + 1.574 omega - 0.176 omega**2` would give.
    #[must_use]
    pub fn alpha_m(&self, cubic: AssociationCubic) -> f64 {
        match cubic {
            AssociationCubic::Srk => self.m_srk,
            AssociationCubic::Pr => self.m_pr,
        }
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
    /// The packing fraction `eta = B/(4V)`.
    pub eta: f64,
    /// `d^2 ln g / d eta^2`, which a second derivative with respect to a *pair* of mole
    /// numbers needs: `d^2 ln g / dn_i dn_j = b_i b_j (d^2 ln g/d eta^2)/(16 V^2)`.
    pub d2_ln_g_d_eta2: f64,
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
            eta,
            d2_ln_g_d_eta2,
        }
    }

    /// `d ln g / d eta`, recomputed from the stored packing fraction.
    fn d_ln_g_d_eta(&self) -> f64 {
        3.0 / (1.0 - self.eta) - 1.0 / (2.0 - self.eta)
    }

    /// `d^2 ln g / dn_i dn_j`, from the packing fraction's second derivative.
    ///
    /// `eta` is linear in every mole number with coefficient `b_i/(4V)`, so the chain rule
    /// contributes `b_i b_j/(16 V^2)` and nothing else.
    #[must_use]
    pub fn d2_ln_g_dn_dn(&self, b_i: f64, b_j: f64, v: f64) -> f64 {
        b_i * b_j * self.d2_ln_g_d_eta2 / (16.0 * v * v)
    }

    /// `d^2 ln g / dn_i dV`, which the association fugacity's volume derivative needs.
    ///
    /// `d ln g/dn_i` is `b_i (d ln g/d eta)/(4V)`, and `eta` moves with `V` as `-B/(4V^2)`.
    #[must_use]
    pub fn d2_ln_g_dn_dv(&self, b_i: f64, sum_b: f64, v: f64) -> f64 {
        -b_i * (sum_b * self.d2_ln_g_d_eta2 / (16.0 * v * v * v)
            + self.d_ln_g_d_eta() / (4.0 * v * v))
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
    /// `hCPA`, the unbonded site count weighted by moles: `sum_i n_i sum_{A in i}(1 - X_A)`.
    ///
    /// NeqSim's `PhaseCPAInterface.calc_hCPA`, and the factor multiplying the distribution
    /// function's derivative in every association term. It is *not* one, which is what
    /// comparing against `validation/neqsim/CpaProbe.java` settled: the port first carried
    /// the class field's initialiser and was 4% out on `dFCPAdN`.
    pub unbonded_sites: f64,
    /// `A_assoc / (R T)`, the association Helmholtz energy over `R T`.
    pub helmholtz_rt: f64,
    /// `d (A_assoc/RT) / dV` at constant composition and temperature.
    pub d_helmholtz_dv: f64,
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
                unbonded_sites: 0.0,
                helmholtz_rt: 0.0,
                d_helmholtz_dv: 0.0,
                ln_phi: vec![0.0; n],
            });
        }

        let rdf = Rdf::new(covolumes, moles, v);
        let delta = self.delta_matrix(covolumes, t, &rdf);
        // `K_ij = n_j Delta_ij / V`, so that `sum_j K_ij X_j` is `S_i` and the fixed point
        // `X_i = 1/(1 + S_i)` is NeqSim's `solveX2`.
        //
        // NeqSim's *Newton* path builds `Klk_ij = n_i n_j Delta_ij / V` instead, but pairs
        // it with the `n_i`-scaled residual `n_i (1/X_i - 1) - sum_j Klk_ij X_j`, which
        // vanishes at the same point. Taking that matrix without its scaling gives each
        // site fraction a spurious factor of its component's mole number, which is what a
        // finite difference against the substitution catches and nothing else does.
        let mut klk = vec![0.0; sites * sites];
        for i in 0..sites {
            for j in 0..sites {
                klk[i * sites + j] = moles[self.component_of_site(j)] * delta[i * sites + j] / v;
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
        // `hCPA = sum_i n_i sum_{A in i} (1 - X_A)`: the number of unbonded sites, weighted
        // by the moles that carry them. NeqSim's `calc_hCPA`, assigned in
        // `PhaseSrkCPA.calcPressure` - not its initialiser, which is a `1.0` the class
        // overwrites before any association term is read.
        let unbonded_sites: f64 = moles
            .iter()
            .enumerate()
            .map(|(i, ni)| {
                let per_site: f64 = self.sites_of(i).map(|a| 1.0 - x[a]).sum();
                ni * per_site
            })
            .sum();
        let ln_phi = (0..n)
            .map(|i| {
                let site_term: f64 = self.sites_of(i).map(|a| x[a].ln()).sum();
                site_term - 0.5 * unbonded_sites * rdf.d_ln_g_dn[i]
            })
            .collect();
        // NeqSim's `dFCPAdV = (hCPA/(2V)) (1 - V d ln g/dV)`.
        let d_helmholtz_dv = unbonded_sites * (1.0 - v * rdf.d_ln_g_dv) / (2.0 * v);

        Ok(SiteState {
            fractions: x,
            iterations,
            converged,
            refined,
            unbonded_sites,
            helmholtz_rt,
            d_helmholtz_dv,
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
            // One right-hand side, so a fresh factorization each step - the Jacobian moves
            // with `X` and reusing it would solve the previous step's system.
            let mut step = vec![residual.clone()];
            if !solve_many(&jacobian, &mut step, sites) {
                return false;
            }
            let step = &step[0];
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

/// The logarithmic temperature derivative of a component's Boltzmann factor,
/// `(exp(e/(RT)) - 1)`.
///
/// Zero at zero energy rather than a division by zero: a component with no association
/// energy has a strength of zero, and the derivative of `log 0` is not the question being
/// asked. NeqSim guards the same ratio with a `1e-50` threshold.
fn boltzmann_d_ln_dt(energy: f64, t: f64) -> f64 {
    let ratio = energy / (R * t);
    let boltzmann = ratio.exp() - 1.0;
    if boltzmann.abs() < 1.0e-50 {
        return 0.0;
    }
    -ratio * ratio.exp() / (t * boltzmann)
}

/// The logarithmic temperature derivative of `DeltaNog` for a pair.
fn delta_nog_d_ln_dt(
    a: &AssociationComponent,
    c: &AssociationComponent,
    rule: CrossRule,
    t: f64,
) -> f64 {
    match rule {
        // `sqrt(K_a K_c)`, so the log-derivative is the mean of the two components'.
        CrossRule::Elliott => {
            0.5 * (boltzmann_d_ln_dt(a.energy, t) + boltzmann_d_ln_dt(c.energy, t))
        }
        CrossRule::Cr1 { energy, .. } => boltzmann_d_ln_dt(energy, t),
    }
}

/// The site fractions' implicit derivatives, and the association fugacity derivatives.
///
/// Every one of these is an *implicit* derivative: the site fractions solve a nonlinear
/// system, so their derivatives are one linear solve against that system's Jacobian, which
/// [`Association::solve`] has already formed for its Newton step. Their being implicit is
/// the whole reason this is a kernel rather than a formula.
#[derive(Debug, Clone, PartialEq)]
pub struct SiteDerivatives {
    /// `d X_A / d n_a`, indexed `[site][component]`, at constant `T` and `V`.
    pub d_fractions_dn: Vec<Vec<f64>>,
    /// `d X_A / d T` at constant composition and volume.
    pub d_fractions_dt: Vec<f64>,
    /// `d X_A / d V` at constant composition and temperature.
    pub d_fractions_dv: Vec<f64>,
    /// `d (ln phi_i^assoc) / d n_a`, at constant `T` and `V`. Indexed `[i][a]`.
    pub d_ln_phi_dn: Vec<Vec<f64>>,
    /// `d (ln phi_i^assoc) / d T` at constant composition and volume.
    pub d_ln_phi_dt: Vec<f64>,
    /// `d (ln phi_i^assoc) / d V` at constant composition and temperature.
    pub d_ln_phi_dv: Vec<f64>,
    /// `d (A_assoc/(R T)) / d T` at constant composition and volume.
    ///
    /// Not a fugacity derivative and not derivable from the others: the enthalpy
    /// departure needs `A - T dA/dT`, and the fugacity coefficient is `dA/dn` alone.
    pub d_helmholtz_dt: f64,
    /// `d^2 (A_assoc/(R T)) / dV^2` at constant composition and temperature.
    pub d2_helmholtz_dv2: f64,
    /// `d^2 X_A / dV^2`, one per site.
    pub d2_fractions_dv2: Vec<f64>,
}

impl Association {
    /// The implicit derivatives of the site fractions and of the association fugacity.
    ///
    /// `state` must be the converged [`SiteState`] at this `(covolumes, moles, v, t)`; the
    /// residual is
    ///
    /// ```text
    /// F_i = X_i (1 + S_i) - 1,   S_i = (1/V) sum_k m_k Delta_ik X_k
    /// ```
    ///
    /// whose Jacobian is `dF_i/dX_k = delta_ik (1 + S_i) + X_i m_k Delta_ik / V`. Each
    /// parameter's derivative is `-J^{-1} dF/d(parameter)` with `X` held, and the three
    /// right-hand sides differ only in which part of `S_i` moves:
    ///
    /// * `n_a` moves the mole number in the sum and `g` through `B`;
    /// * `T` moves `DeltaNog` alone, since `g` carries no temperature;
    /// * `V` moves the explicit `1/V` and `g` through `eta`.
    ///
    /// # Errors
    /// * [`AzothError::InvalidInput`] if the vectors are not one entry per component, or if
    ///   `state` is not for this mixture.
    /// * [`AzothError::OutOfRange`] if `v` or `t` is not positive.
    pub fn derivatives(
        &self,
        covolumes: &[f64],
        moles: &[f64],
        v: f64,
        t: f64,
        state: &SiteState,
    ) -> Result<SiteDerivatives> {
        let n = self.len();
        for (field, given) in [("covolumes", covolumes.len()), ("moles", moles.len())] {
            if given != n {
                return Err(AzothError::invalid_input(
                    field,
                    format!("a mixture of {n} components takes {given} entries"),
                ));
            }
        }
        let sites = self.site_count();
        if state.fractions.len() != sites {
            return Err(AzothError::invalid_input(
                "state",
                format!(
                    "the site state carries {} fractions for a mixture of {sites} sites",
                    state.fractions.len()
                ),
            ));
        }
        if sites == 0 {
            return Ok(SiteDerivatives {
                d_fractions_dn: Vec::new(),
                d_fractions_dt: Vec::new(),
                d_fractions_dv: Vec::new(),
                d_ln_phi_dn: vec![vec![0.0; n]; n],
                d_ln_phi_dt: vec![0.0; n],
                d_ln_phi_dv: vec![0.0; n],
                d_helmholtz_dt: 0.0,
                d2_helmholtz_dv2: 0.0,
                d2_fractions_dv2: Vec::new(),
            });
        }

        let rdf = Rdf::new(covolumes, moles, v);
        let delta = self.delta_matrix(covolumes, t, &rdf);
        let x = &state.fractions;
        let sum_b: f64 = covolumes.iter().zip(moles).map(|(b, m)| b * m).sum();

        // `S_i`, and the Jacobian `dF_i/dX_k`.
        let mut s = vec![0.0; sites];
        let mut jacobian = vec![0.0; sites * sites];
        for i in 0..sites {
            for k in 0..sites {
                let m_k = moles[self.component_of_site(k)];
                s[i] += m_k * delta[i * sites + k] * x[k];
            }
            s[i] /= v;
            for k in 0..sites {
                let m_k = moles[self.component_of_site(k)];
                let diagonal = if i == k { 1.0 } else { 0.0 };
                jacobian[i * sites + k] =
                    diagonal * (1.0 + s[i]) + x[i] * m_k * delta[i * sites + k] / v;
            }
        }

        // The three right-hand sides, negated because each is `J X' = -dF/d(parameter)`.
        let mut rhs = vec![vec![0.0; sites]; n + 2];
        for (a, row) in rhs.iter_mut().enumerate().take(n) {
            for i in 0..sites {
                let explicit: f64 = (0..sites)
                    .filter(|&k| self.component_of_site(k) == a)
                    .map(|k| delta[i * sites + k] * x[k])
                    .sum();
                row[i] = -x[i] * (explicit / v + rdf.d_ln_g_dn[a] * s[i]);
            }
        }
        for i in 0..sites {
            let ci = self.component_of_site(i);
            let mut total = 0.0;
            for k in 0..sites {
                let ck = self.component_of_site(k);
                total += moles[ck]
                    * delta[i * sites + k]
                    * delta_nog_d_ln_dt(
                        &self.components[ci],
                        &self.components[ck],
                        self.cross_rule(ci, ck),
                        t,
                    )
                    * x[k];
            }
            rhs[n][i] = -x[i] * total / v;
            rhs[n + 1][i] = -x[i] * s[i] * (rdf.d_ln_g_dv - 1.0 / v);
        }
        if !solve_many(&jacobian, &mut rhs, sites) {
            return Err(AzothError::invalid_input(
                "site fractions",
                "the site-fraction Jacobian is singular at this state, so the association \
                 derivatives are not defined here",
            ));
        }

        // `dX/dn_a` is column `a` of the solved right-hand sides.
        let d_fractions_dn = (0..sites)
            .map(|i| (0..n).map(|a| rhs[a][i]).collect())
            .collect();
        let d_fractions_dt = rhs[n].clone();
        let d_fractions_dv = rhs[n + 1].clone();

        // `ln phi_i^assoc = sum_{A in i} ln X_A - (h/2) d ln g/dn_i`, with
        // `h = sum_k n_k sum_{B in k}(1 - X_B)`. `h` moves with every parameter too, so
        // each derivative is the site sum plus a product rule over `h` and the
        // distribution function's second derivative.
        let h = state.unbonded_sites;
        let d_h_dn: Vec<f64> = (0..n)
            .map(|a| {
                let explicit: f64 = self.sites_of(a).map(|site| 1.0 - x[site]).sum();
                let through_x: f64 = (0..sites)
                    .map(|site| moles[self.component_of_site(site)] * rhs[a][site])
                    .sum();
                explicit - through_x
            })
            .collect();
        let d_h_dt: f64 = -(0..sites)
            .map(|site| moles[self.component_of_site(site)] * rhs[n][site])
            .sum::<f64>();
        let d_h_dv: f64 = -(0..sites)
            .map(|site| moles[self.component_of_site(site)] * rhs[n + 1][site])
            .sum::<f64>();
        // `d(A/(RT))/dT = sum_i n_i sum_{A in i} (1/X_A - 1/2) X_A^(T)`, which is the
        // enthalpy departure's term and not a fugacity one.
        let d_helmholtz_dt: f64 = (0..sites)
            .map(|site| {
                let m = moles[self.component_of_site(site)];
                m * (1.0 / x[site] - 0.5) * rhs[n][site]
            })
            .sum();

        // --- second order -----------------------------------------------------------
        //
        // The site fractions are implicit twice over. Differentiating
        // `F_A(X(V), V) = 0` twice gives
        //
        //   J X^{(VV)} = -(F_{VV} + J_V X^{(V)}) - H_{AVC} X^{(V)}_C
        //                 - H_{ABC} X^{(V)}_B X^{(V)}_C
        //
        // and the two `H` terms are the ones a *partial* `J_V` drops. `J_V` in the first
        // bracket is `dJ/dV` with `X` **held**, and the total derivative is
        // `dJ/dV = J_V + H_{ABC} X^{(V)}_C`, so the Hessian enters twice: once contracted
        // with `F`'s mixed `V,X` partial and once bilinear in `X^{(V)}`. `H_{ABC}` is
        // `(1/V)(d_AB m_C Delta_AC + d_AC m_B Delta_AB)`, because `S_A` is linear in `X`.
        //
        // Dropping them is not a small error - it is a factor of 1.8 to 2.0 - and it is
        // invisible to any check that compares the identity against itself, because a
        // consistently partial derivation satisfies a consistently partial identity.
        let mut jacobian_v = vec![0.0; sites * sites];
        for a in 0..sites {
            for b in 0..sites {
                let drho = rdf.d_ln_g_dv - 1.0 / v;
                let diagonal = if a == b { 1.0 } else { 0.0 };
                jacobian_v[a * sites + b] = diagonal * s[a] * drho
                    + x[a] * moles[self.component_of_site(b)] * delta[a * sites + b] * drho / v;
            }
        }
        let x_v = &rhs[n + 1];
        let mut rhs2 = vec![vec![0.0; sites]; n + 1];
        for a in 0..sites {
            let drho = rdf.d_ln_g_dv - 1.0 / v;
            let f_vv = x[a] * (s[a] * drho * drho + s[a] * (rdf.d2_ln_g_dv2 + 1.0 / (v * v)));
            // `J_V X^{(V)}`, plus the two Hessian contractions.
            let mut total = 0.0;
            for b in 0..sites {
                total += jacobian_v[a * sites + b] * x_v[b];
            }
            // `H_{AVC} X^{(V)}_C`, with `dF_{A,V}/dX_C = drho (X_A m_C Delta_AC / V + d_AC S_A)`.
            for c in 0..sites {
                let h_avc = drho
                    * (x[a] * moles[self.component_of_site(c)] * delta[a * sites + c] / v
                        + if a == c { s[a] } else { 0.0 });
                total += h_avc * x_v[c];
            }
            // `H_{ABC} X^{(V)}_B X^{(V)}_C`.
            for b in 0..sites {
                for c in 0..sites {
                    let m_b = moles[self.component_of_site(b)];
                    let m_c = moles[self.component_of_site(c)];
                    let h_abc = (if a == b {
                        m_c * delta[a * sites + c]
                    } else {
                        0.0
                    } + if a == c {
                        m_b * delta[a * sites + b]
                    } else {
                        0.0
                    }) / v;
                    total += h_abc * x_v[b] * x_v[c];
                }
            }
            rhs2[n][a] = -(f_vv + total);
        }
        if !solve_many(&jacobian, &mut rhs2, sites) {
            return Err(AzothError::invalid_input(
                "site fractions",
                "the site-fraction Jacobian is singular, so the association's second \
                 derivatives are not defined at this state",
            ));
        }

        let d2_helmholtz_dv2: f64 = (0..sites)
            .map(|a| {
                let m = moles[self.component_of_site(a)];
                let weight = 1.0 / x[a] - 0.5;
                m * (-x_v[a] * x_v[a] / (x[a] * x[a]) + weight * rhs2[n][a])
            })
            .sum();

        let mut d_ln_phi_dn = vec![vec![0.0; n]; n];
        let mut d_ln_phi_dt = vec![0.0; n];
        let mut d_ln_phi_dv = vec![0.0; n];
        for i in 0..n {
            for a in 0..n {
                let site_sum: f64 = self.sites_of(i).map(|site| rhs[a][site] / x[site]).sum();
                d_ln_phi_dn[i][a] = site_sum
                    - 0.5
                        * (d_h_dn[a] * rdf.d_ln_g_dn[i]
                            + h * rdf.d2_ln_g_dn_dn(covolumes[i], covolumes[a], v));
            }
            let site_sum_t: f64 = self.sites_of(i).map(|site| rhs[n][site] / x[site]).sum();
            d_ln_phi_dt[i] = site_sum_t - 0.5 * d_h_dt * rdf.d_ln_g_dn[i];
            let site_sum_v: f64 = self
                .sites_of(i)
                .map(|site| rhs[n + 1][site] / x[site])
                .sum();
            d_ln_phi_dv[i] = site_sum_v
                - 0.5 * (d_h_dv * rdf.d_ln_g_dn[i] + h * rdf.d2_ln_g_dn_dv(covolumes[i], sum_b, v));
        }

        Ok(SiteDerivatives {
            d_fractions_dn,
            d_fractions_dt,
            d_fractions_dv,
            d_ln_phi_dn,
            d_ln_phi_dt,
            d_ln_phi_dv,
            d_helmholtz_dt,
            d2_helmholtz_dv2,
            d2_fractions_dv2: rhs2[n].clone(),
        })
    }
}

/// Solve `A x = b` for several right-hand sides at once, in place.
///
/// The association kernel's Jacobian is the only dense solve in this crate, and it is
/// `sites x sites` - a handful - which is why it is written here rather than taken as a
/// dependency. NeqSim uses EJML for the same solve. The factorization is shared across the
/// right-hand sides because the implicit derivatives are one solve each against the *same*
/// matrix, which is also why `A` is not consumed: elimination runs on a copy.
///
/// Returns `false` on a singular matrix, leaving `rhs` untouched, so a caller reports an
/// undefined derivative rather than a wrong one.
fn solve_many(a: &[f64], rhs: &mut [Vec<f64>], n: usize) -> bool {
    let mut lu = a.to_vec();
    for col in 0..n {
        let mut pivot = col;
        for row in col + 1..n {
            if lu[row * n + col].abs() > lu[pivot * n + col].abs() {
                pivot = row;
            }
        }
        if lu[pivot * n + col] == 0.0 {
            return false;
        }
        if pivot != col {
            for k in 0..n {
                lu.swap(col * n + k, pivot * n + k);
            }
            for x in rhs.iter_mut() {
                x.swap(col, pivot);
            }
        }
        for row in col + 1..n {
            let factor = lu[row * n + col] / lu[col * n + col];
            if factor == 0.0 {
                continue;
            }
            for k in col..n {
                lu[row * n + k] -= factor * lu[col * n + k];
            }
            for x in rhs.iter_mut() {
                x[row] -= factor * x[col];
            }
        }
    }
    for x in rhs.iter_mut() {
        for row in (0..n).rev() {
            let mut sum = x[row];
            for k in row + 1..n {
                sum -= lu[row * n + k] * x[k];
            }
            x[row] = sum / lu[row * n + row];
        }
    }
    true
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
    /// component carries a fitted energy and no bond, so every site fraction is one, the
    /// association energy is exactly zero, and - once `h` is the unbonded site count
    /// rather than the class field's initialiser - so is the fugacity term, because there
    /// are no unbonded sites to weight the distribution function's derivative with.
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
        assert_eq!(state.unbonded_sites, 0.0);
        assert_eq!(state.ln_phi[0], 0.0);
        assert_eq!(state.d_helmholtz_dv, 0.0);
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

    /// A water/methanol mixture, which is the smallest case with a cross-association
    /// between two different schemes and therefore exercises the Elliott combination.
    fn water_and_methanol() -> (Association, [f64; 2], [f64; 2], f64, f64) {
        let a = Association::new(
            vec![
                AssociationComponent {
                    scheme: SiteScheme::FourC,
                    energy: 16655.0,
                    volume: 0.0692,
                },
                AssociationComponent {
                    scheme: SiteScheme::TwoB,
                    energy: 24591.0,
                    volume: 0.0161,
                },
            ],
            Vec::new(),
        )
        .expect("valid parameters");
        (a, [2.6e-5, 3.1e-5], [0.6, 0.4], 6.0e-5, 350.0)
    }

    /// Every composition derivative, against a finite difference of the solve itself.
    ///
    /// This is the check the module exists to pass: the site fractions are the solution of
    /// a nonlinear system, so their derivatives are implicit, and nothing but a numerical
    /// derivative of the solve can confirm the linear algebra that produces them.
    #[test]
    fn the_site_fraction_composition_derivatives_match_the_solve() {
        let (a, b, n, v, t) = water_and_methanol();
        let state = a.solve(&b, &n, v, t).expect("solves");
        let d = a.derivatives(&b, &n, v, t, &state).expect("differentiable");
        let step = 1.0e-9;
        for comp in 0..2 {
            let mut np = n;
            np[comp] += step;
            let shifted = a.solve(&b, &np, v, t).expect("solves");
            for site in 0..a.site_count() {
                let numerical = (shifted.fractions[site] - state.fractions[site]) / step;
                let analytic = d.d_fractions_dn[site][comp];
                assert!(
                    (analytic / numerical - 1.0).abs() < 1.0e-5,
                    "site {site} w.r.t. n_{comp}: analytic {analytic} vs numerical {numerical}"
                );
            }
        }
    }

    #[test]
    fn the_site_fraction_temperature_derivative_matches_the_solve() {
        let (a, b, n, v, t) = water_and_methanol();
        let state = a.solve(&b, &n, v, t).expect("solves");
        let d = a.derivatives(&b, &n, v, t, &state).expect("differentiable");
        let step = 1.0e-4;
        let hot = a.solve(&b, &n, v, t + step).expect("solves");
        let cold = a.solve(&b, &n, v, t - step).expect("solves");
        for site in 0..a.site_count() {
            let numerical = (hot.fractions[site] - cold.fractions[site]) / (2.0 * step);
            let analytic = d.d_fractions_dt[site];
            assert!(
                (analytic / numerical - 1.0).abs() < 1.0e-5,
                "site {site} w.r.t. T: analytic {analytic} vs numerical {numerical}"
            );
        }
    }

    #[test]
    fn the_site_fraction_volume_derivative_matches_the_solve() {
        let (a, b, n, v, t) = water_and_methanol();
        let state = a.solve(&b, &n, v, t).expect("solves");
        let d = a.derivatives(&b, &n, v, t, &state).expect("differentiable");
        let step = 1.0e-12;
        let big = a.solve(&b, &n, v + step, t).expect("solves");
        let small = a.solve(&b, &n, v - step, t).expect("solves");
        for site in 0..a.site_count() {
            let numerical = (big.fractions[site] - small.fractions[site]) / (2.0 * step);
            let analytic = d.d_fractions_dv[site];
            assert!(
                (analytic / numerical - 1.0).abs() < 1.0e-4,
                "site {site} w.r.t. V: analytic {analytic} vs numerical {numerical}"
            );
        }
    }

    /// The fugacity derivatives, against a finite difference of `ln phi`'s own definition
    /// recomputed at the shifted state - not against a stored value.
    #[test]
    fn the_fugacity_composition_derivatives_match_the_definition() {
        let (a, b, n, v, t) = water_and_methanol();
        let state = a.solve(&b, &n, v, t).expect("solves");
        let d = a.derivatives(&b, &n, v, t, &state).expect("differentiable");
        let step = 1.0e-9;
        for comp in 0..2 {
            let mut np = n;
            np[comp] += step;
            let shifted = a.solve(&b, &np, v, t).expect("solves");
            for i in 0..2 {
                let numerical = (shifted.ln_phi[i] - state.ln_phi[i]) / step;
                let analytic = d.d_ln_phi_dn[i][comp];
                assert!(
                    (analytic / numerical - 1.0).abs() < 1.0e-4,
                    "ln phi_{i} w.r.t. n_{comp}: analytic {analytic} vs numerical {numerical}"
                );
            }
        }
    }

    /// The two remaining fugacity derivatives, against the same finite difference.
    ///
    /// `ln phi`'s definition is recomputed at the shifted state rather than read from a
    /// stored value, so this compares the analytic derivative with the model, not with a
    /// second copy of the answer.
    #[test]
    fn the_fugacity_temperature_and_volume_derivatives_match_the_definition() {
        let (a, b, n, v, t) = water_and_methanol();
        let state = a.solve(&b, &n, v, t).expect("solves");
        let d = a.derivatives(&b, &n, v, t, &state).expect("differentiable");

        let dt = 1.0e-4;
        let hot = a.solve(&b, &n, v, t + dt).expect("solves");
        let cold = a.solve(&b, &n, v, t - dt).expect("solves");
        for i in 0..2 {
            let numerical = (hot.ln_phi[i] - cold.ln_phi[i]) / (2.0 * dt);
            assert!(
                (d.d_ln_phi_dt[i] / numerical - 1.0).abs() < 1.0e-5,
                "ln phi_{i} w.r.t. T: analytic {} vs numerical {numerical}",
                d.d_ln_phi_dt[i]
            );
        }

        let dv = 1.0e-12;
        let big = a.solve(&b, &n, v + dv, t).expect("solves");
        let small = a.solve(&b, &n, v - dv, t).expect("solves");
        for i in 0..2 {
            let numerical = (big.ln_phi[i] - small.ln_phi[i]) / (2.0 * dv);
            assert!(
                (d.d_ln_phi_dv[i] / numerical - 1.0).abs() < 1.0e-4,
                "ln phi_{i} w.r.t. V: analytic {} vs numerical {numerical}",
                d.d_ln_phi_dv[i]
            );
        }
    }

    /// The kernel against NeqSim's own numbers.
    ///
    /// `validation/neqsim/CpaProbe.java` drives `SystemSrkCPA` for water/methanol 0.6/0.4 at
    /// 300 K and 100 bar and prints every quantity below. The covolumes it prints are the
    /// databank's `bcpa_srk` column verbatim - water 1.4515, methanol 3.0978 - and the
    /// total volume is NeqSim's internal scaled one. That scaling is legitimate here
    /// because the kernel is invariant under a common rescaling of `V` and `b`: `DeltaNog`
    /// carries one factor of `b` and `S_i` divides by one factor of `V`.
    ///
    /// This is the check the module's other tests cannot make. They differentiate the
    /// kernel's own solve, so they prove it is self-consistent; this one proves it is
    /// *NeqSim's*. It is what caught `hCPA` being the unbonded site count rather than the
    /// initialiser of the class field.
    #[test]
    fn the_kernel_reproduces_neqsims_cpa_probe() {
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
        let a = Association::new(vec![water, methanol], Vec::new()).expect("valid parameters");
        let covolumes = [1.4515, 3.0978];
        let moles = [0.6, 0.4];
        let v = 2.620_326_748_288_29;
        let t = 300.0;

        let rdf = Rdf::new(&covolumes, &moles, v);
        // 1e-10 rather than 1e-15, and the gap is arithmetic order rather than formula:
        // NeqSim forms the packing fraction as `(B/n)/(4 (V/n))` where this forms `B/(4V)`,
        // and writes the distribution function as `(2 - eta)/(2 (1 - eta)^3)` where this
        // writes the identical `(1 - eta/2)/(1 - eta)^3`. Both reorderings are worth a few
        // units in the twelfth digit and nothing more, which is what the measured 1.2e-12
        // on `g` is.
        let close = |got: f64, want: f64, what: &str| {
            assert!(
                (got / want - 1.0).abs() < 1.0e-10,
                "{what}: azoth {got} vs NeqSim {want}"
            );
        };
        close(rdf.g, 1.765_205_648_571_49, "g at contact");
        close(rdf.d_ln_g_dn[0], 0.443_178_855_724_964, "d ln g/dn_water");
        close(
            rdf.d_ln_g_dn[1],
            0.945_834_970_213_430,
            "d ln g/dn_methanol",
        );
        close(rdf.d_ln_g_dv, -0.245_862_964_206_216, "d ln g/dV");

        let state = a.solve(&covolumes, &moles, v, t).expect("solves");
        for (site, want) in [
            (0, 0.101_330_316_299_160),
            (1, 0.101_330_316_299_160),
            (2, 0.101_330_316_299_160),
            (3, 0.101_330_316_299_160),
            (4, 0.031_557_041_194_167_5),
            (5, 0.031_557_041_194_167_5),
        ] {
            close(state.fractions[site], want, &format!("xsite[{site}]"));
        }
        close(state.helmholtz_rt, -6.793_473_163_493_31, "FCPA");
        close(state.unbonded_sites, 2.931_561_607_926_681_7, "hCPA");
        close(state.d_helmholtz_dv, 0.919_769_772_389_381, "dFCPAdV");
        close(state.ln_phi[0], -9.807_081_619_647_33, "dFCPAdN[water]");
        close(state.ln_phi[1], -8.298_303_821_393_14, "dFCPAdN[methanol]");

        // The association energy's volume derivative matches too.
        close(state.d_helmholtz_dv, 0.919_769_772_389_381, "dFCPAdV");

        // The composition derivative does *not* match, and that is the finding.
        //
        // NeqSim's `dFCPAdNdN(i,j)` carries `-0.5 h calc_lngij(j)`, and `calc_lngij` is not
        // the derivative of `calc_lngi`. Measured at this state:
        // `calc_lngij(i,j) = d^2 ln g/dn_i dn_j + (b_i + b_j) f`, where `f` is the
        // distribution function's numerator over its denominator - at `(0,0)` the extra term
        // is 0.443178856, exactly `calc_lngi(0)`, and at `(1,1)` it is 0.945834970, exactly
        // `calc_lngi(1)`. So NeqSim's second derivative is not the derivative of its own
        // first: at this state it returns -0.920056554 where the exact derivative of
        // `dFCPAdN` - which matches NeqSim's to 1e-14 - is -4.097114380.
        //
        // azoth keeps the exact derivative, because a second-order flash solves with it and
        // an inexact Jacobian only slows convergence. This is recorded rather than forced to
        // agree; the finite-difference tests above are what establish the exactness.
        let d = a
            .derivatives(&covolumes, &moles, v, t, &state)
            .expect("differentiable");
        assert!(
            (d.d_ln_phi_dn[0][0] / -4.097_114_380_096_661 - 1.0).abs() < 1.0e-9,
            "azoth's dFCPAdNdN[0][0] is the exact derivative: {}",
            d.d_ln_phi_dn[0][0]
        );
        assert!(
            (d.d_ln_phi_dn[0][0] / -0.920_056_554_036_658 - 1.0).abs() > 1.0,
            "and it is not NeqSim's, whose that gap is the finding recorded above"
        );
    }

    /// The second derivatives, against finite differences of the first.
    ///
    /// `d2X/dV2` first, because it is the piece a partial `J_V` gets wrong and the energy
    /// derivatives cannot localise: the test below would still pass with the Hessian terms
    /// dropped if only the energy were checked.
    #[test]
    fn the_second_derivatives_match_finite_differences_of_the_first() {
        let (a, b, n, v, t) = water_and_methanol();
        let state = a.solve(&b, &n, v, t).expect("solves");
        let d = a.derivatives(&b, &n, v, t, &state).expect("differentiable");

        let x_v = |vol: f64| -> Vec<f64> {
            let s = a.solve(&b, &n, vol, t).expect("solves");
            a.derivatives(&b, &n, vol, t, &s)
                .expect("differentiable")
                .d_fractions_dv
        };
        let phi_v = |vol: f64| -> f64 { a.solve(&b, &n, vol, t).expect("solves").d_helmholtz_dv };

        let h = 1.0e-10;
        let (up, down) = (x_v(v + h), x_v(v - h));
        for site in 0..a.site_count() {
            let numerical = (up[site] - down[site]) / (2.0 * h);
            assert!(
                (d.d2_fractions_dv2[site] / numerical - 1.0).abs() < 1.0e-4,
                "d2X_{site}/dV2: analytic {} vs numerical {numerical}",
                d.d2_fractions_dv2[site]
            );
        }

        let numerical = (phi_v(v + h) - phi_v(v - h)) / (2.0 * h);
        assert!(
            (d.d2_helmholtz_dv2 / numerical - 1.0).abs() < 1.0e-4,
            "d2(A/(RT))/dV2: analytic {} vs numerical {numerical}",
            d.d2_helmholtz_dv2
        );
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
