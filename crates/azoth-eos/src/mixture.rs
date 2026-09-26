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
    Schwartzentruber, Soave, SoreideWhitsonWater, TwuCoon, matcop_kappa, umr_kappa,
};
use crate::association::{
    Association, AssociationCubic, AssociationRecord, NON_ASSOCIATING, R, SiteDerivatives,
};
use crate::cubic::Cubic;
use crate::mixing_rule::{MixingRule, SoreideWhitsonRole, UMR_HWFC};
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

/// One phase's `ln phi` together with its first derivatives at a state.
///
/// See [`Mixture::phase_derivatives`], which computes it, for what each family is
/// differentiated at and why the composition derivative is a mole-number one.
#[derive(Debug, Clone, PartialEq)]
pub struct PhaseDerivatives {
    /// `ln phi_i`, the same values [`PhaseState::ln_phi`] carries at this state.
    pub ln_phi: Vec<f64>,
    /// `d ln phi_i / d n_j` at constant `T`, `P` and the other mole numbers, row `i`.
    ///
    /// At unit total mole number, and carrying the Gibbs-Duhem sum rule
    /// `sum_i x_i d ln phi_i / d n_j = 0`.
    pub d_ln_phi_dn: Vec<Vec<f64>>,
    /// `d ln phi_i / d T` at constant `P` and composition.
    pub d_ln_phi_dt: Vec<f64>,
    /// `d ln phi_i / d P` at constant `T` and composition.
    pub d_ln_phi_dp: Vec<f64>,
}

/// An associating mixture's association terms at one state.
///
/// Carried together because [`Mixture::phase_derivatives`] needs all of them at once and
/// none is meaningful alone: the root sensitivities say how `Z` responds to the state,
/// `alpha` converts one to a molar volume, and the kernel's derivatives are what the
/// association adds at constant volume.
struct AssociationDerivatives {
    /// `R T/P`.
    alpha: f64,
    /// `dZ/dn_j`, at constant `T` and `P`.
    d_z_n: Vec<f64>,
    /// `T dZ/dT`, at constant `P` and composition.
    t_d_z: f64,
    /// `P dZ/dP`, at constant `T` and composition.
    p_d_z: f64,
    /// The kernel's derivatives at the converged state.
    kernel: SiteDerivatives,
    /// The association's `ln phi_i`, which adds to the cubic's.
    ln_phi: Vec<f64>,
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
    /// The four liquid-viscosity parameters `LIQVISC1`-`LIQVISC4`, whose meaning is
    /// The three liquid-conductivity coefficients `LIQCOND1`-`LIQCOND3`, whose polynomial is
    /// `c0 + c1 T + c2 T^2`. Read by `eos.liquid_conductivity_polynom`.
    pub liquid_conductivity: [f64; 3],
    /// [`Self::liqvisc_model`]. Read by `eos.aqueous_viscosity`; a component built from
    /// critical constants alone carries zeros and a model of zero, which that model reads as
    /// NeqSim's own default branch.
    pub liqvisc: [f64; 4],
    /// Which of NeqSim's four liquid-viscosity expressions those four are; zero names none.
    pub liqvisc_model: u32,
    /// The volume-translation parameter in m³/mol, subtracted from the untranslated
    /// molar volume: `v_corr = v - c`. Zero for a component without a translation,
    /// which is the plain PR/SRK/RK forms. The Peneloux shift of
    /// [`crate::pr_peneloux_shift`] or [`crate::srk_peneloux_shift`] is one source; a
    /// fitted constant is another.
    pub volume_shift: f64,
    /// Whether the substance precipitates as wax: the databank's `waxformer`, and the flag
    /// that decides whether it can be in a wax phase at all.
    ///
    /// **NeqSim's wax model reads it where a cubic cannot see it.** `ComponentWax.fugcoef`
    /// returns a `1e50` marker for a substance that is not a former, so the flag is what
    /// makes the exclusion a number the fraction solve can carry rather than a phase list
    /// that has to be assembled per state.
    pub wax_former: bool,
    /// NeqSim's `COMPTYPE` for the substance, lower-cased: `hc`, `inert`, `glycol`, `water`
    /// and the rest.
    ///
    /// **Carried because a phase's *label* is a function of it.** `PhaseEos.init` decides
    /// whether a phase is a gas, an oil or an aqueous one from its volume ratio and from
    /// whether its hydrocarbons outweigh its aqueous components - and "is this a
    /// hydrocarbon" is `COMPTYPE == "HC"` plus the two flags a TBP cut carries, which a
    /// databank row does not need. Empty for a card's substance, which has no table row.
    pub class: String,
    /// Heat of fusion, in J/mol, and the triple-point temperature in K, for a substance the
    /// databank states them for.
    ///
    /// Read by `eos.wax_solid_fugacity`'s fusion term. **Zero means the table carries no
    /// value**, which is every substance the wax model does not treat as a solid; a model
    /// must refuse a zero rather than compute with it, the rule the other zero-means-absent
    /// fields here follow.
    pub heat_of_fusion: f64,
    /// The triple-point temperature, in K.
    pub triple_point_temperature: f64,
    /// The solid's heat-capacity coefficients, in the table's own units.
    pub cp_solid: [f64; 4],
    /// The liquid's heat-capacity coefficients.
    pub cp_liquid: [f64; 5],
    /// The solid's density-correlation coefficients, zero where the table states none.
    pub solid_density_coefs: [f64; 4],
    /// The liquid's density-correlation coefficients, the same table's other half.
    pub liquid_density_coefs: [f64; 4],
    /// The association parameters, or `None` for a component with no site scheme.
    ///
    /// **A component that carries this is not described by its critical constants
    /// alone.** Its attraction and covolume are fitted - water's `b` is 1.4515e-5
    /// m³/mol against `0.08664 R Tc/Pc`'s 2.11e-5 - so an associating equation of state
    /// reads this in place of what the cubic would derive.
    pub association: Option<AssociationRecord>,
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
            liqvisc: [0.0; 4],
            liqvisc_model: 0,
            liquid_conductivity: [0.0; 3],
            volume_shift: 0.0,
            wax_former: false,
            class: String::new(),
            heat_of_fusion: 0.0,
            triple_point_temperature: 0.0,
            cp_solid: [0.0; 4],
            cp_liquid: [0.0; 5],
            solid_density_coefs: [0.0; 4],
            liquid_density_coefs: [0.0; 4],
            association: None,
        })
    }

    /// This component, with its association parameters attached.
    ///
    /// Only an associating equation of state reads them; a cubic ignores the field, as
    /// it ignores the alpha parameters of a correlation it is not using.
    #[must_use]
    pub fn with_association(mut self, association: Option<AssociationRecord>) -> Self {
        self.association = association;
        self
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

    /// This component, with the substance's `COMPTYPE` attached.
    ///
    /// The databank populates it; a card's substance has none, and a model that classifies
    /// by it refuses rather than guessing.
    #[must_use]
    pub fn with_class(mut self, class: String) -> Self {
        self.class = class;
        self
    }

    /// This component, with the wax data attached.
    ///
    /// The databank populates it; a caller-supplied pseudo-component may set it from
    /// `eos.tbp_fraction_properties`'s fits, which is what NeqSim's `addTBPWax` does for a
    /// cut. A model reading it refuses a zero heat of fusion rather than treating it as a
    /// substance that does not melt.
    #[must_use]
    pub fn with_wax_data(
        mut self,
        wax_former: bool,
        heat_of_fusion: f64,
        triple_point_temperature: f64,
    ) -> Self {
        self.wax_former = wax_former;
        self.heat_of_fusion = heat_of_fusion;
        self.triple_point_temperature = triple_point_temperature;
        self
    }

    /// This component, with the solid route's tabulated polynomials attached.
    ///
    /// Only `eos.tp_solid_flash` reads them, and only for a component it is asked to
    /// precipitate; a cubic ignores the field entirely.
    #[must_use]
    pub fn with_solid_tables(
        mut self,
        cp_solid: [f64; 4],
        cp_liquid: [f64; 5],
        solid_density_coefs: [f64; 4],
        liquid_density_coefs: [f64; 4],
    ) -> Self {
        self.cp_solid = cp_solid;
        self.cp_liquid = cp_liquid;
        self.solid_density_coefs = solid_density_coefs;
        self.liquid_density_coefs = liquid_density_coefs;
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

    /// Attach a liquid-viscosity set: the four parameters and the model they belong to.
    /// The three liquid-conductivity coefficients, from the databank's own columns.
    #[must_use]
    pub fn with_liquid_conductivity(mut self, conductivity: [f64; 3]) -> Self {
        self.liquid_conductivity = conductivity;
        self
    }

    pub fn with_liquid_viscosity(mut self, liqvisc: [f64; 4], model: u32) -> Self {
        self.liqvisc = liqvisc;
        self.liqvisc_model = model;
        self
    }
}

/// The UMR rule's state at a composition: the fugacity's three quantities and the
/// departure's one.
///
/// `t_d_alpha_mix` is `T d alpha_mix/dT`, which the enthalpy's excess-Gibbs branch is
/// built from. It travels with `ader` because both are derived from the same
/// `ln gamma` and its temperature derivative, and computing them apart would evaluate
/// the activity coefficients twice.
#[derive(Debug, Clone, PartialEq)]
pub struct UmrState {
    /// The per-component attraction coefficients `qPure_i + hwfc ln gamma_i`.
    pub ader: Vec<f64>,
    /// `sum_i x_i ader_i`, which `A` is `B` times.
    pub alpha_mix: f64,
    /// The co-volume vector the fugacity's volume-derivative term carries.
    pub b_der: Vec<f64>,
    /// `T d alpha_mix/dT` at this composition.
    pub t_d_alpha_mix: f64,
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
    associating: bool,
    /// The Fürst electrolyte term, or `None` for every mixture that is not one.
    ///
    /// Carried whole rather than as a flag because it is not a property of the components:
    /// the short-range table is built from the *whole* composition's names, and the
    /// dielectric coefficients are read per component - so the term is assembled by the
    /// resolver that knows the names, exactly as [`crate::association::Association`] is
    /// assembled from the records the databank supplies.
    furst: Option<crate::furst_electrolyte::FurstElectrolyte>,
    /// The hydrate's own tables, when a mixture was resolved from names for one.
    ///
    /// The `Mixture` holds critical constants and no names, and the hydrate's guest data is
    /// keyed by name - so it travels here, assembled by the resolver that knows them. The
    /// same reason [`Self::furst`] is a field rather than a lookup.
    hydration: Option<crate::hydrate::Hydration>,
    /// The components' names, in the order every vector here is indexed by.
    ///
    /// **A `Component` carries critical constants and no name**, because the databank's key is
    /// not a parameter and a caller building one from constants alone has none to give. A
    /// model that has to look a substance up *by name* therefore has nowhere to read it: the
    /// hydrate's guest tables travel in [`Self::hydration`] for exactly that reason and the
    /// Fürst short-range table in [`Self::furst`], and this is the same thing for a model
    /// whose key is a substance the databank knows and no specialised record carries.
    ///
    /// `None` for a mixture built from constants, which is every caller-supplied one.
    names: Option<Vec<String>>,
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
        let mixture = Self {
            components,
            mixing_rule: MixingRule::Classic { kij },
            cubic: Cubic::default(),
            alpha: Alpha::default(),
            associating: false,
            furst: None,
            hydration: None,
            names: None,
        };
        // Refused here rather than at the first phase evaluation: a mixture whose
        // association cannot be computed is a mistake in what was asked for, and the
        // caller is the one who can fix it.
        Ok(mixture)
    }

    /// This mixture, running the Fürst electrolyte term.
    ///
    /// The term carries its own components and its own short-range table, so this is a
    /// statement that the mixture *is* a Fürst one - the same kind of statement
    /// [`Self::with_association`] makes, and for the same reason: the same substances are
    /// a classical fluid under another model.
    #[must_use]
    pub fn with_furst(mut self, furst: crate::furst_electrolyte::FurstElectrolyte) -> Self {
        self.furst = Some(furst);
        self
    }

    /// The Fürst electrolyte term, if this mixture has one.
    #[must_use]
    pub fn furst(&self) -> Option<&crate::furst_electrolyte::FurstElectrolyte> {
        self.furst.as_ref()
    }

    /// This mixture, with the hydrate's guest tables attached.
    ///
    /// A statement that the mixture *is* one a hydrate can be calculated for - the species
    /// are the same substances under another model.
    #[must_use]
    pub fn with_hydration(mut self, hydration: crate::hydrate::Hydration) -> Self {
        self.hydration = Some(hydration);
        self
    }

    /// The hydrate's tables, if this mixture was resolved for one.
    #[must_use]
    pub fn hydration(&self) -> Option<&crate::hydrate::Hydration> {
        self.hydration.as_ref()
    }

    /// This mixture, carrying the names its components were resolved from.
    #[must_use]
    pub fn with_names(mut self, names: Vec<String>) -> Self {
        self.names = Some(names);
        self
    }

    /// The components' names, in the order every vector here is indexed by, or `None` for a
    /// mixture built from constants.
    #[must_use]
    pub fn names(&self) -> Option<&[String]> {
        self.names.as_deref()
    }

    /// The index of a component by name, or `None` where this mixture carries no names or
    /// none matches.
    ///
    /// The lookup every name-keyed model starts with, and it is here rather than in each of
    /// them so that the case-insensitivity is decided once: the databank's keys are lower
    /// case and a caller's spelling need not be.
    #[must_use]
    pub fn index_of(&self, name: &str) -> Option<usize> {
        let key = name.trim().to_lowercase();
        self.names
            .as_ref()?
            .iter()
            .position(|candidate| candidate.to_lowercase() == key)
    }

    /// This mixture, running the Wertheim association contribution.
    ///
    /// **Whether a mixture associates is the phase model's decision, not the
    /// components'.** A component carries its fitted CPA parameters whenever the table
    /// has them - the databank resolves them by name - but a cubic that is not CPA must
    /// not use them: NeqSim's `SystemSrkCPA` builds `ComponentSrkCPA` and substitutes the
    /// fitted `a` and `b`, while `SystemNRTL`'s `PhaseSrkEos` builds `ComponentSrkEos`
    /// over the *same* substances and uses the cubic's own. So this is opt-in, and a
    /// mixture that does not call it ignores the field.
    ///
    /// # Errors
    pub fn with_association(mut self) -> Result<Self> {
        self.associating = true;
        Ok(self)
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

        for (index, component) in self.components.iter().enumerate() {
            let reduced_temperature = t.value / component.tc.value;
            let reduced_pressure = p.value / component.pc.value;
            reduced_temperatures.push(reduced_temperature);
            // **A Fürst ion's attraction and covolume are not the cubic's.** NeqSim's
            // `ComponentModifiedFurstElectrolyteEos` overwrites both in its constructor -
            // `b = (p0 d^3 + p1)` from the fitted parameters and `a = 1e-35` - so a port
            // that left them at Peng-Robinson's would give every ion a real attraction and
            // a covolume of the wrong size, and the root would be plausible and wrong.
            // Substituted here rather than on the component because it is the *model's*
            // statement, exactly as the association's fitted `a`/`b` substitution is.
            if let Some(furst) = self.furst.as_ref()
                && let Some(covolume) = furst.ion_covolume(index)
            {
                // The reduced forms of the two: `A_i = a_i P/(R T)^2` and `B_i = b_i P/(R T)`,
                // both already dimensionless, so this is the same reduction the cubic's own
                // branch does with `omega_a alpha Pr/Tr^2` and `omega_b Pr/Tr`.
                let r_t = R * t.value;
                a.push(
                    crate::furst_electrolyte::FurstElectrolyte::ion_attraction() * p.value
                        / (r_t * r_t),
                );
                b.push(covolume * p.value / r_t);
                // The ion's alpha is Schwartzentruber's in NeqSim and contributes nothing:
                // `psi` is read only through weights carrying `sqrt(A_i A_j)`, and `A` is
                // `1e-35`.
                psi.push(0.0);
                psi_t.push(0.0);
                continue;
            }
            // The kappa correlation belongs to the alpha term; the Omega constants to
            // the cubic. SRK and PR share the Soave alpha form and differ in both; RK
            // is kappa-free.
            let cubic = self.cubic;
            let non_soave = |term: &dyn AlphaTerm| -> (f64, f64, f64, f64, f64) {
                let alpha = term.alpha(reduced_temperature);
                let a_reduced = cubic.omega_a() * alpha * reduced_pressure
                    / (reduced_temperature * reduced_temperature);
                let b_reduced = cubic.omega_b() * reduced_pressure / reduced_temperature;
                (
                    a_reduced,
                    b_reduced,
                    term.psi(reduced_temperature),
                    term.psi_t(reduced_temperature),
                    alpha,
                )
            };
            // The alpha itself comes back with the reduced pair because a fitted
            // association set *replaces* `a` and `b` and then has to be scaled by the
            // model's alpha - `ComponentUMRCPA.setAttractiveTerm` installs the
            // Mathias-Copeman term over the substituted set, where `ComponentSrkCPA`
            // installs a Soave coefficient. Recomputing it in that branch would be a
            // second expression of the same match.
            let (a_reduced, b_reduced, psi_value, psi_t_value, alpha_value) = match self.cubic {
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
                    Alpha::SoreideWhitson => {
                        // **Water by role, and the salinity from the same rule.** The role
                        // vector is the only per-component identity this layer has, which is
                        // why `Mixture::new` refuses this alpha beside any other rule: there
                        // would be nothing to read the role from.
                        // **Refused rather than resolved against a zero.** This is the one
                        // alpha that is not a function of the component alone - it needs the
                        // salinity and the role vector, and beside another rule there is
                        // nothing to read either from. An error rather than a panic, per this
                        // crate's no-panic rule.
                        let MixingRule::SoreideWhitson {
                            roles, salinity, ..
                        } = &self.mixing_rule
                        else {
                            return Err(AzothError::invalid_input(
                                "mixing_rule",
                                "the Soreide-Whitson alpha takes its salinity from the \
                                 Soreide-Whitson mixing rule and reads which component is \
                                 water from that rule's roles, so it cannot be resolved \
                                 beside any other rule",
                            ));
                        };
                        let salinity = *salinity;
                        if roles[index] == SoreideWhitsonRole::Water {
                            non_soave(&SoreideWhitsonWater {
                                salinity,
                                critical_temperature: component.tc.value,
                            })
                        } else {
                            // Every other component is Peng-Robinson 1978, which is what
                            // `AttractiveTermSoreideWhitson.alpha` delegates to.
                            let kappa = pr78_kappa(component.omega)?;
                            warnings.extend(kappa.warnings);
                            non_soave(&Soave { kappa: kappa.kappa })
                        }
                    }
                    Alpha::Pr | Alpha::Srk | Alpha::Pr78 | Alpha::Twu | Alpha::SrkFitted => {
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
                            // **The coefficient is supplied rather than derived**, which is
                            // what makes this a variant of its own rather than a flag on
                            // `Alpha::Srk`: there is no correlation to fall back to, so a
                            // component without one is refused rather than given Soave's
                            // default for an acentric factor nobody stated.
                            Alpha::SrkFitted => (
                                component.alpha_params.first().copied().ok_or_else(|| {
                                    AzothError::invalid_input(
                                        "alpha_params",
                                        format!(
                                            "component {index} takes the fitted Soave alpha \
                                                 and carries no coefficient; `alpha_params[0]` \
                                                 is where it is read from, and \
                                                 `eos.tbp_fraction_properties` is what gives a \
                                                 cut its own"
                                        ),
                                    )
                                })?,
                                Vec::new(),
                            ),
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
                            term.alpha(reduced_temperature),
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
                        term.alpha(reduced_temperature),
                    )
                }
                Cubic::Tst => {
                    let kappa = twu_kappa(component.omega)?;
                    warnings.extend(kappa.warnings);
                    non_soave(&Soave { kappa: kappa.kappa })
                }
            };
            // An associating component carries its own attraction and covolume, and
            // NeqSim's `ComponentSrkCPA` substitutes them for the cubic's:
            // `if (|aCPA| > 1e-6) { a = aCPA; b = bCPA; }`. This is not a refinement of
            // the critical-constant values - water's fitted covolume is 1.4515e-5
            // m3/mol against `0.08664 R Tc/Pc`'s 2.11e-5 - and its alpha coefficient is
            // fitted too, so the whole `a_i(T)` differs.
            let family = AssociationCubic::of(self.cubic);
            // **UMR-CPA substitutes a third set, and scales it with the model's own
            // alpha.** `ComponentUMRCPA.setAttractiveTerm` installs the five-parameter
            // Mathias-Copeman term (22) over the substituted `a` and `b`, where
            // `ComponentSrkCPA` calls `setm(mCPA)` and gets a Soave coefficient - so
            // which alpha scales the fitted set is the model's choice, and the model
            // states it by asking for `Alpha::MatCop5PrUmr`.
            if self.associating
                && let Some(record) = &component.association
                && let (Alpha::MatCop5PrUmr, Some(umr)) = (self.alpha, &record.umr_cpa)
            {
                let r_t = R * t.value;
                a.push(umr.attraction() * alpha_value * p.value / (r_t * r_t));
                b.push(umr.covolume() * p.value / r_t);
                psi.push(psi_value);
                psi_t.push(psi_t_value);
                continue;
            }
            if self.associating
                && let Some(record) = &component.association
                && record.has_fitted_set(family)
            {
                let term = Soave {
                    kappa: record.alpha_m(family),
                };
                let alpha = term.alpha(reduced_temperature);
                let r_t = R * t.value;
                a.push(record.attraction(family) * alpha * p.value / (r_t * r_t));
                b.push(record.covolume(family) * p.value / r_t);
                psi.push(term.psi(reduced_temperature));
                psi_t.push(term.psi_t(reduced_temperature));
                continue;
            }
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
            pressure: p.value,
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
        let seed = match side {
            RootSide::Liquid => z_min,
            RootSide::Vapour => z_max,
        };
        let _ = seed;
        // An associating mixture's root is **not** a root of the cubic: the association
        // contributes to the pressure, and NeqSim's `PhaseSrkCPA.molarVolume` solves for
        // the volume where the *total* pressure equals the specified one. See
        // [`Self::associating_root`].
        let z = match (self.association(), self.furst.as_ref()) {
            (Some(_), Some(_)) => {
                return Err(AzothError::invalid_input(
                    "mixture",
                    "this mixture carries both the Wertheim association and the Fürst \
                     electrolyte term. NeqSim's `PhaseElectrolyteCPA` is the model that \
                     has both, and it is not ported - so a mixture naming two terms is a \
                     request for a model that does not exist here rather than a phase \
                     whose pressure is the sum of the two",
                ));
            }
            (Some(association), None) => self.associating_root(reduced, x, &association, side)?,
            (None, Some(furst)) => self.furst_root(reduced, x, furst, side)?,
            (None, None) => seed,
        };
        self.phase_state_at(reduced, x, z)
    }

    /// The root of an associating mixture's equation of state.
    ///
    /// **The association carries a pressure, so the root is not the cubic's.**
    /// `PhaseSrkCPA.molarVolume` solves `BonV - (B/n) dFdV() - P B/(n R T) = 0`, and
    /// `PhaseSrkCPA.dFdV()` is `super.dFdV() + dFCPAdV()` - the cubic *plus* the
    /// association. Measured on water/methanol at 300 K and 100 bar that is the whole of
    /// the difference between this library's `Z = 0.15437` and NeqSim's `0.10505`, and it
    /// is not a small correction to the cubic root: NeqSim's own `A` and `B` give a cubic
    /// root of `0.15229`, so the association moves the root by 31%.
    ///
    /// The residual, with `A` and `B` the reduced parameters at the specified pressure and
    /// `P_a = -R T d(A_assoc/(R T))/dV` evaluated at `V = Z R T/P`:
    ///
    /// ```text
    /// f(Z) = Z - Z/(Z - B) + A Z/((Z + d1 B)(Z + d2 B)) - P_a Z/P
    /// ```
    ///
    /// The first three terms are the cubic's own residual, so a non-associating component
    /// of the sum reduces to `f_cubic` exactly.
    ///
    /// **Newton from the cubic root, with a bisection fallback.** NeqSim Newton-solves
    /// with the second and third volume derivatives and a damping rule; a numerical
    /// derivative and a bracket is this library's house style where the analytic form
    /// would be another derivation, and the seed being the cubic root means the fallback
    /// is only reached when the association has moved the root out of the basin.
    ///
    /// # Errors
    /// * Propagates the association solve's errors.
    /// * [`AzothError::OutOfRange`] if no root is found above `B`.
    fn associating_root(
        &self,
        reduced: &ReducedParameters,
        x: &[f64],
        association: &Association,
        side: RootSide,
    ) -> Result<f64> {
        let (a_mix, b_mix) = self.mixture_parameters(reduced, x);
        // The cubic's own root for this side is the seed, and its ordering is what picks
        // the branch: the association moves a root, it does not choose between them.
        let seed = self.cubic_root(a_mix, b_mix, side)?;
        let (r, t, p) = (R, reduced.t_kelvin, reduced.pressure);
        let delta1 = self.cubic.delta1();
        let delta2 = self.cubic.delta2();
        // The covolumes in m3/mol, recovered from their reduced form. The kernel is
        // invariant under a common rescaling of covolume and volume, and `V = Z R T/P`
        // below is SI, so these must be too.
        let covolumes: Vec<f64> = reduced.b.iter().map(|b| b * r * t / p).collect();

        // One association solve per evaluation: `d(A/(RT))/dV` is on the state the solve
        // already returns, so this needs no derivative call of its own.
        let residual = |z: f64| -> Result<f64> {
            let v = z * r * t / p;
            let state = association.solve(&covolumes, x, v, t)?;
            let p_assoc = -r * t * state.d_helmholtz_dv;
            let cubic =
                z - z / (z - b_mix) + a_mix * z / ((z + delta1 * b_mix) * (z + delta2 * b_mix));
            Ok(cubic - p_assoc * z / p)
        };

        Self::root_with_extra_pressure(&residual, seed, b_mix, side, "associating")
    }

    /// The root of a Fürst electrolyte mixture's equation of state.
    ///
    /// **The electrolyte terms carry a pressure, so the root is not the cubic's**, and the
    /// residual has the same shape as the association's:
    ///
    /// ```text
    /// f(Z) = Z - Z/(Z - B) + A Z/((Z + d1 B)(Z + d2 B)) - P_e Z/P
    /// ```
    ///
    /// with `P_e = -R T d(A_e/(R T))/dV_si` evaluated at `V = Z R T/P`. NeqSim gives this
    /// model its own Halley iteration rather than the CPA one, but the *equation* is the
    /// same pressure balance: at the probe's converged root `dFdV = 0.415155577019407` and
    /// `Z = 1 - V 1e5 dFdV` gives `0.00963585` against the printed
    /// `0.00963585200923298`. The root of an equation does not depend on the iteration that
    /// found it, so this reuses the shared solver and the probe's `Z` is the check.
    ///
    /// **The composition is passed as mole numbers summing to one**, which is what makes the
    /// term's *intensive* outputs - the pressure and the `ln phi` - independent of how large
    /// the phase is. Its Helmholtz energy is extensive and is not read here.
    ///
    /// # Errors
    /// * Propagates the term's errors.
    /// * [`AzothError::OutOfRange`] if no root is found above `B`.
    fn furst_root(
        &self,
        reduced: &ReducedParameters,
        x: &[f64],
        furst: &crate::furst_electrolyte::FurstElectrolyte,
        side: RootSide,
    ) -> Result<f64> {
        let (a_mix, b_mix) = self.mixture_parameters(reduced, x);
        let seed = self.cubic_root(a_mix, b_mix, side)?;
        let (r, t, p) = (R, reduced.t_kelvin, reduced.pressure);
        let delta1 = self.cubic.delta1();
        let delta2 = self.cubic.delta2();

        let residual = |z: f64| -> Result<f64> {
            let solution = furst.solve(t, x, z * r * t / p)?;
            let p_term = -r * t * solution.helmholtz_rt_dv;
            let cubic =
                z - z / (z - b_mix) + a_mix * z / ((z + delta1 * b_mix) * (z + delta2 * b_mix));
            Ok(cubic - p_term * z / p)
        };
        Self::root_with_extra_pressure(&residual, seed, b_mix, side, "electrolyte")
    }

    /// The cubic's own root for a side, which both pressure roots seed from.
    ///
    /// # Errors
    /// * [`AzothError::OutOfRange`] if the cubic has no admissible root.
    fn cubic_root(&self, a_mix: f64, b_mix: f64, side: RootSide) -> Result<f64> {
        match self.cubic {
            Cubic::Pr | Cubic::Tst => match side {
                RootSide::Liquid => Ok(pr_z_factor(a_mix, b_mix)?.z_min),
                RootSide::Vapour => Ok(pr_z_factor(a_mix, b_mix)?.z_max),
            },
            Cubic::Srk | Cubic::Rk => match side {
                RootSide::Liquid => Ok(srk_z_factor(a_mix, b_mix)?.z_min),
                RootSide::Vapour => Ok(srk_z_factor(a_mix, b_mix)?.z_max),
            },
        }
    }

    /// The walk, Newton and scan a cubic-plus-pressure root solve shares.
    ///
    /// Both terms that carry a pressure need the same three passes - a geometric walk up
    /// from the covolume for the liquid side, Newton from the cubic's own root, and a
    /// linear scan with a bisection for whatever neither finds - and they differ only in
    /// the residual. Factored rather than copied because the *solver* is the part that
    /// must not drift: a fix to the walk in one is a fix the other needs.
    ///
    /// `seed` is the cubic's own root for the side wanted, which is what makes the walk
    /// and Newton safe to start from; `what` names the term in the failure message, so a
    /// refusal says which pressure carried the equation away.
    ///
    /// # Errors
    /// * Propagates the residual's own errors.
    /// * [`AzothError::OutOfRange`] if no sign change is found above `B`.
    fn root_with_extra_pressure(
        residual: &impl Fn(f64) -> Result<f64>,
        seed: f64,
        b_mix: f64,
        side: RootSide,
        what: &str,
    ) -> Result<f64> {
        /// The zero of `f` in `[lo, hi]`, by bisection. The count is fixed rather than
        /// converged on a residual: near the covolume `f` is a difference of two large
        /// terms, and a bisection that stopped at an absolute tolerance there would stop
        /// early.
        fn bisect(f: &impl Fn(f64) -> Result<f64>, lo: f64, hi: f64) -> Result<f64> {
            let (mut lo, mut hi) = (lo, hi);
            let sign = f(lo)?.signum();
            for _ in 0..200 {
                let mid = 0.5 * (lo + hi);
                if f(mid)?.signum() == sign {
                    lo = mid;
                } else {
                    hi = mid;
                }
            }
            Ok(0.5 * (lo + hi))
        }

        // `RootSide::Liquid` wants the **lowest** zero above `B`, and the seed cannot
        // find it. At low pressure the substituted cubic's `a` is small enough that it
        // has one real root, so `z_min` and `z_max` are the same number - the vapour root
        // - and Newton from it converges there and stops. **The association's pressure is
        // what creates the liquid root**, so it has to be *found* rather than seeded: the
        // residual is `-inf` at `B+`, and the zero is a fraction of a per cent above `B`,
        // which a linear scan steps straight over and a cubic root never points at.
        //
        // The walk is geometric because the root's *ratio* to `B` is what is bounded, not
        // its distance. Newton still serves the vapour side, and both sides reach the
        // fallback for the case neither handles.
        if matches!(side, RootSide::Liquid) {
            let floor = b_mix * 1.000_001;
            let steps = 64;
            let mut previous = floor;
            let mut previous_value = residual(previous)?;
            for step in 1..=steps {
                let z_try = floor * (seed / floor).powf((step as f64) / (steps as f64));
                let value = residual(z_try)?;
                if previous_value.signum() != value.signum() {
                    return bisect(&residual, previous, z_try);
                }
                previous = z_try;
                previous_value = value;
            }
        }

        // Newton, with the step capped to a tenth of the distance to the covolume so it
        // cannot cross the `Z > B` boundary where the cubic term is singular.
        let mut z = seed.max(b_mix * 1.000_001);
        for _ in 0..60 {
            let value = residual(z)?;
            if value.abs() < 1.0e-12 {
                return Ok(z);
            }
            let step = 1.0e-7 * z.abs().max(1.0e-6);
            let slope = (residual(z + step)? - value) / step;
            if !slope.is_finite() || slope == 0.0 {
                break;
            }
            let mut delta = -value / slope;
            let ceiling = 0.1 * (z - b_mix);
            if delta.abs() > ceiling {
                delta = ceiling * delta.signum();
            }
            let next = z + delta;
            if !next.is_finite() || next <= b_mix {
                break;
            }
            z = next;
        }

        // The fallback: scan for a sign change and bisect. NeqSim reaches
        // `molarVolumeChangePhase` on the same failure, so arriving here is upstream's
        // behaviour and not a port artifact.
        let ceiling = match side {
            RootSide::Liquid => seed.max(1.0),
            RootSide::Vapour => (seed * 8.0).max(1.0),
        };
        let floor = b_mix * 1.000_001;
        let steps = 400;
        let mut previous = floor;
        let mut previous_value = residual(previous)?;
        for step in 1..=steps {
            let z_try = floor + (ceiling - floor) * (step as f64) / (steps as f64);
            let value = residual(z_try)?;
            if previous_value == 0.0 {
                return Ok(previous);
            }
            if previous_value.signum() != value.signum() {
                return bisect(&residual, previous, z_try);
            }
            previous = z_try;
            previous_value = value;
        }
        Err(AzothError::OutOfRange {
            field: "z".to_string(),
            value: seed,
            detail: format!(
                "no volume root for this {what} mixture above B = {b_mix}: the residual \
             has no sign change between the covolume and {ceiling}"
            ),
        })
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
        // `z - 1` as the equation of state writes it. Equal to `z - 1` at the cubic's own
        // root and not at an associating mixture's, which is the whole of the difference
        // between this model and NeqSim's at every state where the two were compared.
        let z_minus_one = self.cubic.eos_z_minus_one(z, a_mix, b_mix);

        let mut ln_phi: Vec<f64> = match &self.mixing_rule {
            MixingRule::HuronVidal { .. }
            | MixingRule::WongSandler { .. }
            | MixingRule::Umr { .. } => {
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
                    b_ratio * z_minus_one - ln_z_minus_b - coefficient * factor * i_term
                })
                .collect(),
        };

        // The Wertheim association, when the phase model runs it. Its `ln phi` adds to
        // the cubic's, and its Helmholtz energy to the departure enthalpy.
        //
        // **The two must move together.** `s_dep_r` below is the Gibbs identity
        // `h_dep_rt - sum_i x_i ln phi_i`, not a second formula - so adding the
        // association to `ln phi` alone would fold it silently into the entropy. The
        // enthalpy's association term is `A/(RT) - T d(A/RT)/dT`, which is why the
        // kernel carries that derivative at all.
        let mut assoc_dep_rt = 0.0;
        if let Some(association) = self.association() {
            let r_t = R * reduced.t_kelvin;
            // `reduced.b[i]` is `b_i P/(R T)`, so this is the covolume in m3/mol. The
            // kernel is invariant under a common rescaling of covolume and volume, but
            // its `Delta` carries the covolume and its `S_i` divides by the volume, so
            // the two must be in the same units - which SI gives.
            let covolumes: Vec<f64> = reduced
                .b
                .iter()
                .map(|b| b * r_t / reduced.pressure)
                .collect();
            let v = z * r_t / reduced.pressure;
            let state = association.solve(&covolumes, x, v, reduced.t_kelvin)?;
            for (value, addition) in ln_phi.iter_mut().zip(&state.ln_phi) {
                *value += addition;
            }
            let derivatives =
                association.derivatives(&covolumes, x, v, reduced.t_kelvin, &state)?;
            assoc_dep_rt = -reduced.t_kelvin * derivatives.d_helmholtz_dt;
        }

        // The Fürst electrolyte contribution, at the same root. **`ln phi` and not the
        // Helmholtz energy**: the term's `dFdN` is `ln phi_i = dFdN_i - ln Z`'s own
        // `dFdN`, and this library's cubic `ln phi` is already that identity for its part,
        // so the contributions add. The term takes the composition as mole numbers summing
        // to one, which is what makes its intensive outputs independent of the phase size.
        if let Some(furst) = self.furst.as_ref() {
            let r_t = R * reduced.t_kelvin;
            let solution = furst.solve(reduced.t_kelvin, x, z * r_t / reduced.pressure)?;
            for (value, addition) in ln_phi.iter_mut().zip(&solution.ln_phi) {
                *value += addition;
            }
        }

        // The mixture's departure functions. `psi_bar` is the composition-weighted
        // average of the components' `psi`, and the two lines below are then
        // `pr_departure`'s expressions with `psi_bar` in place of `psi` - which is
        // what makes them reduce to it exactly at one component.
        //
        // **A universal rule's departure is not the weighted one.** `h_dep_rt`'s cubic term
        // is `(T dA/dT + A)/(delta_diff B) * I`, and for the classical rule that is
        // `A (psi_bar - 1)/(delta_diff B)` because `A` carries `(R T)^-2` through every
        // weight. For the UMR rule `A = B alpha_mix` with **both** factors moving, so
        // `T dA/dT + A` collapses to `B * T d alpha_mix/dT` and the term is
        // `T d alpha_mix/dT / delta_diff * I` - no `psi_bar` in it at all. A pure fluid
        // hides the difference: there `alpha_mix` is one component's `qPure`, and the two
        // expressions agree exactly, which is why this was invisible until a mixture was
        // compared against NeqSim.
        let umr_state = match &self.mixing_rule {
            MixingRule::Umr { .. } => Some(self.umr_ader(reduced, x)),
            _ => None,
        };
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
        let excess = match &umr_state {
            Some(state) => state.t_d_alpha_mix / self.cubic.delta_diff(),
            None => coefficient * (psi_bar - 1.0),
        };
        let h_dep_rt = (z - 1.0) + excess * i_term + assoc_dep_rt;
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

    /// An associating mixture's root sensitivities to `n_j`, `T` and `P`.
    ///
    /// **The root is the associating residual's, not the cubic's**, so its sensitivity is
    /// taken from that residual: the cubic's own `dZ/dn_j` is `-0.8` where the associating
    /// one is `-1.7` at the state `tests/mixture.rs` measures, and the difference is the
    /// association's pressure responding to the composition.
    ///
    /// The residual is `Q(Z) - (Z/P) P_a(x, V(Z))` where `Q` is the cubic's z-cubic and
    /// `P_a = -R T d(A_assoc/(R T))/dV` - the polynomial form `associating_root` solves, so
    /// that `dQ/dZ` is the cubic's own `df_dz` and the two sensitivities are consistent.
    /// `Q` vanishes at the same `Z` as the monic cubic's residual because it *is* that
    /// residual times `q(Z) = (Z + delta1 B)(Z + delta2 B)`; multiplying through by `q` is
    /// what lets one derivative serve both.
    ///
    /// # Errors
    /// * Propagates the association solve's errors.
    // Eight arguments, and each is a sensitivity the caller has already computed at this
    // state. Recomputing them here would be a second derivation of `A`'s and `B`'s
    // temperature derivatives, which is the defect the module is organised to avoid.
    #[allow(clippy::too_many_arguments)]
    fn association_derivatives(
        &self,
        reduced: &ReducedParameters,
        x: &[f64],
        association: &Association,
        z: f64,
        a_mix: f64,
        b_mix: f64,
        abar: &[f64],
        t_d_a: f64,
        t_d_b: f64,
    ) -> Result<AssociationDerivatives> {
        let n = self.len();
        let (r, t, p) = (R, reduced.t_kelvin, reduced.pressure);
        // `R T/P`, which carries a `Z` sensitivity to a molar volume and back.
        let alpha = r * t / p;
        let covolumes: Vec<f64> = reduced.b.iter().map(|b| b * alpha).collect();
        let v = z * alpha;
        let state = association.solve(&covolumes, x, v, t)?;
        let kernel = association.derivatives(&covolumes, x, v, t, &state)?;

        // `Q(Z) := (Z - B) q(Z) f_cubic(Z)` is the monic z-cubic whose derivative the
        // cubic's own `df_dz` already is, and `(Z - B) q(Z) Z P_a/P` is the association's
        // term written the same way. So the residual this differentiates is
        //
        //   R(Z) = monic(Z) - S(Z) P_a(V(Z))/P,   S(Z) = (Z - B) q(Z)
        //
        // **`Z q(Z)`, not `S`.** Multiplying the cubic's residual through by `q` alone
        // leaves a residual whose roots are not the equation of state's, and the two
        // differ by `-B q P_a/P` - a difference that vanishes only when the association
        // does. The first version of this used `Z q` and disagreed with the root's own
        // finite difference in both sign and magnitude.
        let b = b_mix;
        let (ds, dp) = (self.cubic.delta_sum(), self.cubic.delta_prod());
        let q = z * z + ds * b * z + dp * b * b;
        let s_factor = (z - b) * q;
        let s_prime = q + (z - b) * (2.0 * z + ds * b);
        let d_s_db = -q + (z - b) * (ds * z + 2.0 * dp * b);

        let p_assoc = -r * t * state.d_helmholtz_dv;
        let a_vv = kernel.d2_helmholtz_dv2;
        let a_vt = kernel.d2_helmholtz_dv_dt;
        let fa = self.cubic.df_da(z, b);
        let fb = self.cubic.df_db(z, a_mix, b);
        let fz = self.cubic.df_dz(z, a_mix, b);

        // `P_a` depends on `V`, and `V = Z R T/P` moves with every one of the three state
        // variables - which is why each sensitivity below carries an `A_vv` term beside
        // the explicit one.
        let g_z = fz - (s_prime / p) * p_assoc + s_factor * alpha * alpha * a_vv;

        let d_a_dn: Vec<f64> = (0..n).map(|j| 2.0 * (abar[j] - a_mix)).collect();
        let d_b_dn: Vec<f64> = (0..n).map(|j| reduced.b[j] - b_mix).collect();
        // `dA/dn_j` and `dB/dn_j` are the *normalised* partials - the ones that hold the
        // total mole number at one, and the ones NeqSim reports - so the association's
        // composition derivative has to be taken the same way. The kernel's is raw, so
        // the chain rule from `x = n/(sum n)` subtracts the composition-weighted sum of
        // all of them, exactly as `2 abar_j` becomes `2 (abar_j - A)`.
        let weighted_a_vn: f64 = (0..n).map(|i| x[i] * kernel.d2_helmholtz_dv_dn[i]).sum();
        let d_z_n: Vec<f64> = (0..n)
            .map(|j| {
                -((fa * d_a_dn[j] + fb * d_b_dn[j]) - d_s_db * d_b_dn[j] * p_assoc / p
                    + s_factor * alpha * (kernel.d2_helmholtz_dv_dn[j] - weighted_a_vn))
                    / g_z
            })
            .collect();

        // `T dB/dT = -B`, so `T dS/dT = -B dS/dB`, and
        // `T dP_a/dT = P_a - R T alpha Z A_vv - R T^2 A_vt` at fixed `Z`.
        let t_d_p_assoc = p_assoc - r * t * alpha * z * a_vv - r * t * t * a_vt;
        let t_d_z = -(self.cubic.t_dfdt(z, a_mix, b, t_d_a, t_d_b) + (b * d_s_db / p) * p_assoc
            - (s_factor / p) * t_d_p_assoc)
            / g_z;

        // `P dB/dP = B`, and `P dP_a/dP = Z alpha^2 A_vv` at fixed `Z`.
        let p_d_z = -(fa * a_mix + fb * b - d_s_db * b * p_assoc / p + s_factor * p_assoc / p
            - s_factor * z * alpha * alpha * a_vv)
            / g_z;

        Ok(AssociationDerivatives {
            alpha,
            d_z_n,
            t_d_z,
            p_d_z,
            ln_phi: state.ln_phi.clone(),
            kernel,
        })
    }

    /// One phase's `ln phi` together with its first derivatives at a state.
    ///
    /// The three families a second-order flash is written in, all analytic:
    ///
    /// * `d_ln_phi_dn[i][j]` - `d ln phi_i / d n_j` at constant `T`, `P` and the other
    ///   mole numbers.
    /// * `d_ln_phi_dt[i]` - `d ln phi_i / d T` at constant `P` and composition.
    /// * `d_ln_phi_dp[i]` - `d ln phi_i / d P` at constant `T` and composition.
    ///
    /// # Mole numbers, not mole fractions
    ///
    /// The composition derivative is taken holding the *other mole numbers* fixed, so
    /// the total moves with `j` and the result is not intensive. The values returned are
    /// for a total mole number of one, which is what a composition summing to one is,
    /// and they carry the Gibbs-Duhem sum rule `sum_i x_i d ln phi_i / d n_j = 0` at
    /// that scale - the property a constant-composition derivative would not have.
    ///
    /// # Scope
    ///
    /// The classic van der Waals one-fluid rules. The activity-coefficient rules write
    /// `ln phi` in a different form - the `ader`/`alpha_mix`/`b_der` branch of
    /// [`Self::phase_state_at`] - and their composition derivative is a second
    /// derivation over the excess Gibbs energy's own Hessian, which is not this one.
    /// Asking for it is an error naming that, not a silent return of the wrong matrix.
    ///
    /// # Errors
    /// * [`AzothError::InvalidInput`] if `x` is not one entry per component, or if the
    ///   mixing rule is one of the three this does not differentiate.
    /// * [`AzothError::OutOfRange`] if the root is not admissible, on the same test
    ///   [`Self::phase_state_at`] makes.
    pub fn phase_derivatives(
        &self,
        reduced: &ReducedParameters,
        x: &[f64],
        z: f64,
    ) -> Result<PhaseDerivatives> {
        let n = self.len();
        if x.len() != n {
            return Err(AzothError::invalid_input(
                "x",
                format!("a composition for {n} components has {} entries", x.len()),
            ));
        }
        match &self.mixing_rule {
            MixingRule::HuronVidal { .. }
            | MixingRule::WongSandler { .. }
            | MixingRule::Umr { .. } => {
                return Err(AzothError::invalid_input(
                    "mixing_rule",
                    "the activity-coefficient rules write ln phi as `ader`, `alpha_mix` \
                     and `b_der`, whose composition derivative is the excess Gibbs \
                     energy's second derivative rather than this classical one. The \
                     derivative surface covers the cubic family, and this rule is \
                     outside it",
                ));
            }
            MixingRule::SoreideWhitson { .. } => {
                return Err(AzothError::invalid_input(
                    "mixing_rule",
                    "the Soreide-Whitson rule resolves its interaction matrix against \
                     the *phase composition*, so its `A_ij` carries a composition \
                     derivative this classical one does not have. The derivative \
                     surface covers the rules whose matrix is fixed by the state",
                ));
            }
            _ => {}
        }

        let kij = self.phase_kij(reduced, x);
        let (a_mix, b_mix) = self.mixture_parameters(reduced, x);
        if !z.is_finite() || z <= b_mix {
            return Err(AzothError::OutOfRange {
                field: "z".to_string(),
                value: z,
                detail: format!(
                    "the root must exceed the mixture's B = {b_mix}, because `ln(z - B)` \
                     is otherwise the logarithm of a negative number"
                ),
            });
        }

        let delta1 = self.cubic.delta1();
        let delta2 = self.cubic.delta2();
        let delta_diff = self.cubic.delta_diff();
        let i_term = self.cubic.i_term(z, b_mix);
        let coefficient = self.cubic.coefficient(a_mix, b_mix);
        let ln_z_minus_b = (z - b_mix).ln();
        let fz = self.cubic.df_dz(z, a_mix, b_mix);
        let fa = self.cubic.df_da(z, b_mix);
        let fb = self.cubic.df_db(z, a_mix, b_mix);
        let (z_minus_one, zmo_dz, zmo_da, zmo_db) =
            self.cubic.eos_z_minus_one_partials(z, a_mix, b_mix);

        // `A_ij` and its two row sums. `abar[i]` is `sum_j x_j A_ij`, the same quantity
        // `phase_state_at` hoists as `cross`.
        let mut a_ij = vec![vec![0.0; n]; n];
        let mut abar = vec![0.0; n];
        for i in 0..n {
            for j in 0..n {
                a_ij[i][j] = (1.0 - kij[i * n + j]) * (reduced.a[i] * reduced.a[j]).sqrt();
                abar[i] += x[j] * a_ij[i][j];
            }
        }
        let b_ratio: Vec<f64> = (0..n).map(|i| reduced.b[i] / b_mix).collect();
        let factor: Vec<f64> = (0..n).map(|i| 2.0 * abar[i] / a_mix - b_ratio[i]).collect();
        let ln_phi: Vec<f64> = (0..n)
            .map(|i| b_ratio[i] * z_minus_one - ln_z_minus_b - coefficient * factor[i] * i_term)
            .collect();

        let t_dkij = self.mixing_rule.t_d_effective_kij(reduced.t_kelvin);

        // --- the state's sensitivities -------------------------------------------
        //
        // `A_i` and `B_i` are functions of `(T, P)` alone, so composition moves `A` and
        // `B` alone: at unit total mole number `dA/dn_j = 2 (abar_j - A)` and
        // `dB/dn_j = B_j - B`. `A_ij` moves with temperature through both the components'
        // alphas and the interaction parameter, and `B_ij` scales as `1/T`, so
        // `T dB/dT = -B`.
        let d_a_dn: Vec<f64> = (0..n).map(|j| 2.0 * (abar[j] - a_mix)).collect();
        let d_b_dn: Vec<f64> = (0..n).map(|j| reduced.b[j] - b_mix).collect();

        let mut t_d_a_ij = vec![vec![0.0; n]; n];
        let mut t_d_abar = vec![0.0; n];
        for i in 0..n {
            for j in 0..n {
                let psi_pair = 0.5 * (reduced.psi[i] + reduced.psi[j]);
                // `sqrt(A_i A_j)` is the part the alpha functions carry; the interaction
                // parameter multiplies it, so its own temperature derivative enters
                // with the opposite sign.
                let root = (reduced.a[i] * reduced.a[j]).sqrt();
                t_d_a_ij[i][j] = a_ij[i][j] * (psi_pair - 2.0) - t_dkij[i * n + j] * root;
                t_d_abar[i] += x[j] * t_d_a_ij[i][j];
            }
        }
        let mut t_d_a = 0.0;
        for i in 0..n {
            t_d_a += x[i] * t_d_abar[i];
        }
        let t_d_b = -b_mix;

        // --- the root's sensitivities --------------------------------------------
        //
        // **An associating mixture's root is not the cubic's, so its sensitivity is not
        // either.** The association carries a pressure that responds to the composition
        // and to the temperature, and the whole of the difference between this surface
        // and the cubic's is here.
        let association = self.association();
        let sensitivities = match &association {
            Some(association) => Some(self.association_derivatives(
                reduced,
                x,
                association,
                z,
                a_mix,
                b_mix,
                &abar,
                t_d_a,
                t_d_b,
            )?),
            None => None,
        };
        let (d_z_dn, t_d_z, p_d_z) = match &sensitivities {
            Some(association) => (
                association.d_z_n.clone(),
                association.t_d_z,
                association.p_d_z,
            ),
            None => (
                (0..n)
                    .map(|j| -(fa * d_a_dn[j] + fb * d_b_dn[j]) / fz)
                    .collect::<Vec<f64>>(),
                -(fa * t_d_a + fb * t_d_b) / fz,
                -(fa * a_mix + fb * b_mix) / fz,
            ),
        };

        // --- composition, at constant temperature and pressure -------------------
        let d_coeff_dn: Vec<f64> = (0..n)
            .map(|j| {
                d_a_dn[j] / (delta_diff * b_mix) - a_mix * d_b_dn[j] / (delta_diff * b_mix * b_mix)
            })
            .collect();
        let d_i_dn: Vec<f64> = (0..n)
            .map(|j| {
                (d_z_dn[j] + delta1 * d_b_dn[j]) / (z + delta1 * b_mix)
                    - (d_z_dn[j] + delta2 * d_b_dn[j]) / (z + delta2 * b_mix)
            })
            .collect();

        // The association's `ln phi`, which `PhaseState::ln_phi` carries too: the surface
        // is the same values, so a caller comparing the two must find them equal.
        let mut ln_phi = ln_phi;
        if let Some(association) = &sensitivities {
            for (value, addition) in ln_phi.iter_mut().zip(&association.ln_phi) {
                *value += addition;
            }
        }

        let mut d_ln_phi_dn = vec![vec![0.0; n]; n];
        for i in 0..n {
            for j in 0..n {
                // `d(abar_i)/dn_j = A_ij - abar_i`, and the `B_i/B` term's derivative is
                // `-B_i/B^2 dB/dn_j`, which is `factor_i`'s `-b_ratio_i` differentiated.
                let d_factor = 2.0 * ((a_ij[i][j] - abar[i]) * a_mix - abar[i] * d_a_dn[j])
                    / (a_mix * a_mix)
                    + b_ratio[i] * d_b_dn[j] / b_mix;
                // `B_i/B` times the *equation of state's* `z - 1`, so the composition moves
                // it through `z`, through `A` and through `B`, and the `B` in the
                // denominator moves it again.
                let d_z_minus_one = zmo_dz * d_z_dn[j] + zmo_da * d_a_dn[j] + zmo_db * d_b_dn[j];
                d_ln_phi_dn[i][j] = b_ratio[i] * d_z_minus_one
                    - (d_z_dn[j] - d_b_dn[j]) / (z - b_mix)
                    // `b_ratio_i` is `B_i/B`, and `B` moves with the composition, so
                    // the term above is not the whole of it. The `B_i/B` inside
                    // `factor_i` carries the same derivative with the opposite sign,
                    // and the two do not cancel.
                    - b_ratio[i] * z_minus_one * d_b_dn[j] / b_mix
                    - d_coeff_dn[j] * factor[i] * i_term
                    - coefficient * d_factor * i_term
                    - coefficient * factor[i] * d_i_dn[j];
            }
        }
        // The association's own contribution: the kernel's `d ln phi_i/dn_j` at constant
        // volume, plus what the volume does - `dV/dn_j = (R T/P) dZ/dn_j`.
        // The same normalisation as the root's: the kernel's composition derivative holds
        // the other mole numbers, not the total, so the chain rule from `x = n/(sum n)`
        // subtracts the composition-weighted sum of the row.
        if let Some(association) = &sensitivities {
            for (i, row) in d_ln_phi_dn.iter_mut().enumerate() {
                let kernel_row = &association.kernel.d_ln_phi_dn[i];
                let weighted: f64 = (0..n).map(|k| x[k] * kernel_row[k]).sum();
                let volume = association.alpha * association.kernel.d_ln_phi_dv[i];
                for (j, value) in row.iter_mut().enumerate() {
                    *value += kernel_row[j] - weighted + volume * association.d_z_n[j];
                }
            }
        }

        // --- temperature, at constant pressure and composition --------------------
        let t_d_coeff = t_d_a / (delta_diff * b_mix) - a_mix * t_d_b / (delta_diff * b_mix * b_mix);
        let t_d_i = (t_d_z + delta1 * t_d_b) / (z + delta1 * b_mix)
            - (t_d_z + delta2 * t_d_b) / (z + delta2 * b_mix);
        let mut d_ln_phi_dt = vec![0.0; n];
        for i in 0..n {
            // `b_ratio_i` is `B_i/B`, two quantities that both scale as `1/T`, so it is
            // temperature-independent and drops out of the sum below.
            let t_d_factor = 2.0 * (t_d_abar[i] * a_mix - abar[i] * t_d_a) / (a_mix * a_mix);
            let t_d_z_minus_one = zmo_dz * t_d_z + zmo_da * t_d_a + zmo_db * t_d_b;
            let t_d_ln_phi = b_ratio[i] * t_d_z_minus_one
                - (t_d_z - t_d_b) / (z - b_mix)
                - t_d_coeff * factor[i] * i_term
                - coefficient * t_d_factor * i_term
                - coefficient * factor[i] * t_d_i;
            d_ln_phi_dt[i] = t_d_ln_phi / reduced.t_kelvin;
        }
        // And the association's, whose volume moves with the temperature as
        // `dV/dT = (R/P)(Z + T dZ/dT)`.
        if let Some(association) = &sensitivities {
            let volume = association.alpha / reduced.t_kelvin * (z + association.t_d_z);
            for (i, value) in d_ln_phi_dt.iter_mut().enumerate() {
                *value +=
                    association.kernel.d_ln_phi_dt[i] + volume * association.kernel.d_ln_phi_dv[i];
            }
        }

        // --- pressure, at constant temperature and composition --------------------
        //
        // Every `A_i` and `B_i` is linear in `P`, so `A`, `B` and their row sums scale
        // as `P` while `coefficient` and `factor_i` - ratios of two such quantities -
        // do not move at all.
        let p_d_i = (p_d_z + delta1 * b_mix) / (z + delta1 * b_mix)
            - (p_d_z + delta2 * b_mix) / (z + delta2 * b_mix);
        let mut d_ln_phi_dp = vec![0.0; n];
        for i in 0..n {
            let p_d_z_minus_one = zmo_dz * p_d_z + zmo_da * a_mix + zmo_db * b_mix;
            let p_d_ln_phi = b_ratio[i] * p_d_z_minus_one
                - (p_d_z - b_mix) / (z - b_mix)
                - coefficient * factor[i] * p_d_i;
            d_ln_phi_dp[i] = p_d_ln_phi / reduced.pressure;
        }
        // The association carries no explicit pressure, so its whole contribution is the
        // volume's: `dV/dP = (R T/P)(dZ/dP - Z/P)`.
        if let Some(association) = &sensitivities {
            let volume = association.alpha * (association.p_d_z - z) / reduced.pressure;
            for (i, value) in d_ln_phi_dp.iter_mut().enumerate() {
                *value += volume * association.kernel.d_ln_phi_dv[i];
            }
        }

        Ok(PhaseDerivatives {
            ln_phi,
            d_ln_phi_dn,
            d_ln_phi_dt,
            d_ln_phi_dp,
        })
    }

    /// The mixture's association parameters at this cubic, or `None` when no component
    /// carries a site scheme.
    ///
    /// Built from the components rather than stored, because the fitted families differ
    /// by cubic and [`Self::with_cubic`] can change the cubic after construction.
    pub fn association(&self) -> Option<Association> {
        if !self.associating {
            return None;
        }
        let cubic = AssociationCubic::of(self.cubic);
        if !self.components.iter().any(|c| c.association.is_some()) {
            return None;
        }
        // **The UMR-CPA set replaces the family's for a component that carries one**, the
        // same substitution `reduced_parameters` makes for `a` and `b`, and gated the same
        // way: the model states that it is UMR-CPA by asking for `Alpha::MatCop5PrUmr`.
        // Water is the check - this set's `kappa_AB` is 0.125 and its `eps` 14177 J/mol
        // against the PR family's, so reading the wrong one moves its `ln phi` in the
        // third decimal while leaving the methane beside it almost untouched.
        let umr_cpa = self.alpha == Alpha::MatCop5PrUmr;
        Association::new(
            self.components
                .iter()
                .map(|c| {
                    c.association.as_ref().map_or(NON_ASSOCIATING, |r| {
                        match (umr_cpa, &r.umr_cpa) {
                            (true, Some(umr)) => umr.at(r.scheme),
                            _ => r.at(cubic),
                        }
                    })
                })
                .collect(),
            // No cross-rule overrides: `INTER.csv`'s `cpaBetaCross`/`cpaEpsCross` are
            // not read yet, and NeqSim's default for every pair is the Elliott rule.
            Vec::new(),
        )
        .ok()
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

    /// The UMR attraction coefficients `A_i/B_i + hwfc ln gamma_i`, and `T d alpha_mix/dT`.
    ///
    /// `EosMixingRuleHandler.init`'s `qPure[i] + hwfc * ln(gamma_i)`, with
    /// `qPure[i] = a_i^T/(b_i R T)` and `hwfc = -1/0.53`. The reduced attraction `A_i`
    /// that [`ReducedParameters::a`] carries is `a_i^T P/(R T)^2` and the co-volume
    /// `B_i` is `b_i P/(R T)`, so their ratio **is** `a_i^T/(b_i R T)` and the division
    /// cancels the state; that is why this needs no pressure.
    ///
    /// **The temperature derivative comes back with it because the departure needs it and
    /// the two must be one expression.** `T d alpha_mix/dT = sum_i x_i [qPure_i (psi_i - 1)
    /// + hwfc T d ln gamma_i/dT]` - the first term from `a_i^T/(b_i R T)` carrying `1/T`,
    /// the second from the residual only, since the Flory-Huggins combinatorial is built
    /// from `r` and `q` and does not move.
    fn umr_ader(&self, reduced: &ReducedParameters, x: &[f64]) -> UmrState {
        let MixingRule::Umr { unifac, .. } = &self.mixing_rule else {
            unreachable!("umr_ader is only called for the UMR rule");
        };
        let (ln_gamma, t_d_ln_gamma) =
            crate::unifac_umrpru_activity_coefficients::umrpru_ln_gamma_and_dt(
                unifac,
                reduced.t_kelvin,
                x,
            );
        let qpure: Vec<f64> = (0..self.len())
            .map(|i| reduced.a[i] / reduced.b[i])
            .collect();
        let ader = crate::mixing_rule::umr_ader(&qpure, &ln_gamma);
        let alpha_mix: f64 = (0..self.len()).map(|i| x[i] * ader[i]).sum();
        let t_d_alpha_mix: f64 = (0..self.len())
            .map(|i| x[i] * (qpure[i] * (reduced.psi[i] - 1.0) + UMR_HWFC * t_d_ln_gamma[i]))
            .sum();
        UmrState {
            ader,
            alpha_mix,
            b_der: reduced.b.clone(),
            t_d_alpha_mix,
        }
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

    /// The Wong-Sandler attraction coefficients.
    ///
    /// The `DijT` is the rule's own, from `WSGIJT`/`WSGJIT`, and **this was measured
    /// rather than reasoned**. NeqSim's `WongSandlerMixingRule` is handed `NRTLDijT` -
    /// loaded from those two columns - and passes it into a `PhaseGENRTLmodifiedHV`,
    /// while the rule's own `getHVDijTParameter` reports the *Huron-Vidal* `HVGIJT`
    /// instead, so the two readings of the same rule disagree about which column is in
    /// force. NeqSim settles it: for `CO2`/`water` at 350 K, 5 bar, whose `WSGIJT` is
    /// 0.96 and `HVGIJT` is -0.842, its Wong-Sandler GE evaluation gives
    /// `ln gamma = [2.9589427251536327, 2.2321445185372717]`, which is this function on
    /// `WSGIJT` to the last digit and is not either of the alternatives.
    ///
    /// This used to force zeros, on the reading that the rule took a DijT-free form. It
    /// did so *and passed its test*, because that test is `water`/`ethanol`, where both
    /// columns are zero and the choice is unobservable.
    fn ws_ader(&self, reduced: &ReducedParameters, x: &[f64]) -> Vec<f64> {
        let MixingRule::WongSandler {
            kij,
            hv_gij,
            hv_gij_t,
            hv_alpha,
            hv_pairs,
        } = &self.mixing_rule
        else {
            unreachable!("ws_ader is only called for the Wong-Sandler rule");
        };
        self.ge_ader(reduced, x, kij, hv_gij, hv_gij_t, hv_alpha, hv_pairs)
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
            MixingRule::Umr { .. } => {
                let state = self.umr_ader(reduced, x);
                (state.ader, state.alpha_mix, state.b_der)
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
            MixingRule::Umr { .. } => {
                let state = self.umr_ader(reduced, x);
                let b_mix: f64 = (0..n).map(|i| x[i] * reduced.b[i]).sum();
                (b_mix * state.alpha_mix, b_mix)
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
        let cubic = -total * (1.0 - b_sum).ln() - (a_sum / (self.cubic.delta_diff() * b_sum)) * g;
        // **The association's own `A^R/RT`, at the same volume.** Without it this surface
        // is a cubic's, and the two things that read it - `eos.stability_test`'s choice of
        // the feed's root, and every RootSide choice behind it - compare two roots by the
        // wrong energy. An associating mixture's two roots at 356 K and 1 bar differ by
        // `0.157 RT`, and the association is most of that.
        let Some(association) = self.association() else {
            return Ok(cubic);
        };
        let r_t = R * reduced.t_kelvin;
        // **Mole numbers, and the volume `V = Z R T/P`.** This function's `n` are mole
        // numbers and its volume is fixed by `Z` alone - `b_i/V` is `B_i/Z` for every
        // component, whatever the moles - so a perturbation at fixed `Z` *is* a
        // perturbation at fixed volume. The kernel's `A_assoc/(RT)` is
        // `sum_i n_i sum_A (ln X_A - X_A/2 + 1/2)`, extensive in the moles it is given;
        // fractions there were wrong by the total, which the finite difference caught, and
        // a volume scaled by the total was wrong again.
        let covolumes: Vec<f64> = reduced
            .b
            .iter()
            .map(|b| b * r_t / reduced.pressure)
            .collect();
        let state = association.solve(
            &covolumes,
            n,
            compressibility * r_t / reduced.pressure,
            reduced.t_kelvin,
        )?;
        Ok(cubic + state.helmholtz_rt)
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
        // **Refused rather than answered from the cubic.** The energy above carries the
        // association and this does not, so for an associating mixture the two would be
        // inconsistent - the Hessian would not be the energy's second derivative, and
        // `criticality_matrix` would place a critical point on a cubic the mixture is not.
        // The same rule the heat-capacity departure follows: report the association's
        // absence rather than the cubic's value.
        if self.association().is_some() {
            return Err(AzothError::invalid_input(
                "mixture",
                "this mixture associates, and the association's second composition \
                 derivative is not assembled - so the Hessian here would be the cubic's \
                 and not the energy's",
            ));
        }
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
pub fn wilson_k(mixture: &Mixture, t: ThermodynamicTemperature, p: Pressure) -> Vec<f64> {
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
    /// The state's absolute pressure, in pascals.
    ///
    /// Carried for the same reason as the temperature: a derivative with respect to
    /// pressure needs the pressure, and a caller holding a `ReducedParameters` should
    /// not have to hand back a second number that could describe a different state.
    pub pressure: f64,
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
