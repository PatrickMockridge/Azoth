//! `KineticReaction` - a configured rate law, and not the one P10 landed.
//!
//! **P10's `Kinetics` is a different class answering a different question.** That one is
//! NeqSim's `thermodynamics.chemicalreaction.Kinetics`, whose two laws sit behind a selector
//! and whose legacy form is the only one a fluid's own reactions reach. This is
//! `process.equipment.reactor.KineticReaction`, an object a **caller configures** and hands to
//! a reactor: pre-exponential factor, activation energy, a stoichiometry map, per-component
//! orders, and optionally adsorption terms. Nothing in a fluid reaches it, so it has no home in
//! `azoth-reactions` - it is the process layer's.
//!
//! **Three rate types, and one of them is the first.** `RateType::Equilibrium` is documented as
//! "use `GibbsReactor` instead for full equilibrium" and its `calculateRate` dispatches to the
//! **power-law** body, so an equilibrium-typed reaction is a power law here. That is the
//! class's behaviour and this reproduces it rather than refusing it.
//!
//! **The gas constant is the class's own third value.** `R_GAS = 8.31446` in
//! `KineticReaction.java:67` - against ISO 6976's `8.314510` and the Sm³ conversion's
//! `8.3144621`, both already in this tree. The three are within `6e-6` of each other and the
//! rate constant's exponential makes that visible, so this keeps the class's number and says so
//! rather than reaching for a shared one.

use azoth_core::{AzothError, Result};

/// The class's own gas constant, J/(mol·K).
pub const R_GAS: f64 = 8.314_46;

/// Which of the class's three rate expressions evaluates the rate.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RateType {
    /// `r = k(T) * Π Cᵢ^αᵢ`, with the reverse term when the reaction is reversible.
    PowerLaw,
    /// Langmuir-Hinshelwood-Hougen-Watson: the power law over a denominator of adsorption
    /// terms.
    Lhhw,
    /// The class's third type, which its own `calculateRate` evaluates as a power law.
    Equilibrium,
}

/// What a rate's units are per, which the reactor converts before it integrates.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RateBasis {
    /// mol/(m³ reactor · s): already volumetric, and the class's default.
    Volume,
    /// mol/(kg catalyst · s).
    CatalystMass,
    /// mol/(m² catalyst · s).
    CatalystArea,
}

/// A rate law a reactor is configured with.
///
/// The fields are the class's own setters and defaults (`KineticReaction.java:98-143`), and the
/// species collections are **ordered vectors rather than maps**, because the class's are Java
/// `LinkedHashMap`s whose insertion order is read: the forward rate multiplies its factors in
/// that order, floating-point multiplication is not associative, and the species a reactor adds
/// to its feed are appended in it too. A `BTreeMap` would be name-ordered and would answer
/// differently in the last bits, so the order is the class's and
/// [`add_reactant`](Self::add_reactant) and [`add_product`](Self::add_product) are how it is set.
#[derive(Debug, Clone, PartialEq)]
pub struct KineticReaction {
    /// The reaction's name, for a spec's `notes` and the class's `getName`.
    pub name: String,
    /// Which expression evaluates the rate.
    pub rate_type: RateType,
    /// What the rate's units are per.
    pub rate_basis: RateBasis,
    /// `A` in `k(T) = A · Tⁿ · exp(-Ea/(R·T))`, the class's default `1.0e10`.
    pub pre_exponential_factor: f64,
    /// `Ea`, J/mol, the class's default `100000.0`.
    pub activation_energy: f64,
    /// `n`, the class's default `0.0`.
    pub temperature_exponent: f64,
    /// ΔH of reaction, J/mol, which the energy balance consumes. The class's default `0.0`.
    pub heat_of_reaction: f64,
    /// Whether the reverse term is subtracted.
    pub reversible: bool,
    /// `[a, b, c, d]` for `ln Keq = a + b/T + c·ln T + d·T`. The class's default is four zeros,
    /// which makes `Keq = 1`.
    pub equilibrium_coefficients: [f64; 4],
    /// Stoichiometric coefficients in insertion order, negative for reactants. The species set a
    /// reactor must carry is the union of this and the feed.
    pub stoichiometry: Vec<(String, f64)>,
    /// Kinetic orders for the forward rate, in insertion order; empty means `k(T)` alone.
    pub reaction_orders: Vec<(String, f64)>,
    /// Orders for the reverse term of a reversible reaction.
    pub product_orders: Vec<(String, f64)>,
    /// LHHW adsorption terms: component → (Kᵢ factor, adsorption order).
    pub adsorption_terms: Vec<(String, (f64, f64))>,
    /// The LHHW denominator's exponent `m`, the class's default `1`.
    pub adsorption_exponent: i32,
    /// The adsorption pre-exponential factor, the class's default `1.0`.
    pub adsorption_pre_exponential_factor: f64,
    /// The adsorption activation energy, J/mol, the class's default `0.0`.
    pub adsorption_activation_energy: f64,
}

impl KineticReaction {
    /// A reaction with the class's defaults and no stoichiometry.
    ///
    /// The defaults are real numbers and not a blank: `A = 1.0e10` with `Ea = 100000 J/mol`
    /// gives `k(500 K) = 0.36`, so an unconfigured reaction has a rate rather than none.
    #[must_use]
    pub fn new(name: impl Into<String>) -> Self {
        Self {
            name: name.into(),
            rate_type: RateType::PowerLaw,
            rate_basis: RateBasis::Volume,
            pre_exponential_factor: 1.0e10,
            activation_energy: 100_000.0,
            temperature_exponent: 0.0,
            heat_of_reaction: 0.0,
            reversible: false,
            equilibrium_coefficients: [0.0; 4],
            stoichiometry: Vec::new(),
            reaction_orders: Vec::new(),
            product_orders: Vec::new(),
            adsorption_terms: Vec::new(),
            adsorption_exponent: 1,
            adsorption_pre_exponential_factor: 1.0,
            adsorption_activation_energy: 0.0,
        }
    }

    /// A reactant: its stoichiometry is stored **negative**, and its order set, exactly as
    /// `addReactant(name, stoichCoeff, order)` does.
    pub fn add_reactant(&mut self, component: impl Into<String>, stoich_coeff: f64, order: f64) {
        let component = component.into();
        self.stoichiometry
            .push((component.clone(), -stoich_coeff.abs()));
        self.reaction_orders.push((component, order));
    }

    /// A product, with its stoichiometry stored positive.
    pub fn add_product(&mut self, component: impl Into<String>, stoich_coeff: f64) {
        self.stoichiometry
            .push((component.into(), stoich_coeff.abs()));
    }

    /// A product with a reverse-order, for a reversible reaction.
    pub fn add_product_with_order(
        &mut self,
        component: impl Into<String>,
        stoich_coeff: f64,
        reverse_order: f64,
    ) {
        let component = component.into();
        self.stoichiometry
            .push((component.clone(), stoich_coeff.abs()));
        self.product_orders.push((component, reverse_order));
    }

    /// An LHHW adsorption term: the component, its `Kᵢ` factor and its adsorption order.
    pub fn add_adsorption_term(&mut self, component: impl Into<String>, factor: f64, order: f64) {
        self.adsorption_terms
            .push((component.into(), (factor, order)));
    }

    /// `k(T) = A · Tⁿ · exp(-Ea/(R·T))`, the class's Arrhenius form.
    ///
    /// **A zero exponent skips the power rather than raising to it.** The class guards with
    /// `|n| > 1e-10` (`KineticReaction.java:227`) and takes `T⁰ = 1` otherwise, so an exponent
    /// of `1e-11` is not a very small correction but *exactly* one - which the port reproduces,
    /// because `powf` on a near-zero exponent does not give one bit-exactly.
    #[must_use]
    pub fn rate_constant(&self, temperature: f64) -> f64 {
        let t_power = if self.temperature_exponent.abs() > 1.0e-10 {
            temperature.powf(self.temperature_exponent)
        } else {
            1.0
        };
        self.pre_exponential_factor
            * t_power
            * (-self.activation_energy / (R_GAS * temperature)).exp()
    }

    /// `Keq(T)`, from `ln Keq = a + b/T + c·ln T + d·T`.
    #[must_use]
    pub fn equilibrium_constant(&self, temperature: f64) -> f64 {
        let [a, b, c, d] = self.equilibrium_coefficients;
        (a + b / temperature + c * temperature.ln() + d * temperature).exp()
    }

    /// The rate at a state, given each component's molar concentration in mol/m³.
    ///
    /// `concentration` returns `0.0` for a component the state does not carry, which is the
    /// class's own `getConcentration` behaviour for an absent component.
    ///
    /// # Errors
    /// [`AzothError::InvalidInput`] if `temperature` is not positive, since the Arrhenius term
    /// divides by it.
    pub fn rate(&self, temperature: f64, concentration: impl Fn(&str) -> f64) -> Result<f64> {
        if temperature <= 0.0 {
            return Err(AzothError::invalid_input(
                "temperature",
                format!("a rate constant needs a positive temperature, not {temperature}"),
            ));
        }
        let k = self.rate_constant(temperature);
        match self.rate_type {
            RateType::PowerLaw | RateType::Equilibrium => {
                self.power_law(temperature, k, &concentration)
            }
            RateType::Lhhw => self.lhhw(temperature, k, &concentration),
        }
    }

    /// The forward power-law product, or `None` where the class returns a hard zero.
    ///
    /// The zero is the class's: a reactant whose concentration is non-positive and whose order
    /// is positive makes the whole rate zero rather than a NaN, which is the `0^n` case stated
    /// as a convention.
    fn forward(&self, k: f64, concentration: &impl Fn(&str) -> f64) -> Option<f64> {
        let mut forward = k;
        for (component, order) in &self.reaction_orders {
            let c = concentration(component);
            if c <= 0.0 && *order > 0.0 {
                return None;
            }
            forward *= c.max(0.0).powf(*order);
        }
        Some(forward)
    }

    fn power_law(
        &self,
        temperature: f64,
        k: f64,
        concentration: &impl Fn(&str) -> f64,
    ) -> Result<f64> {
        let Some(forward) = self.forward(k, concentration) else {
            return Ok(0.0);
        };
        if !self.reversible {
            return Ok(forward);
        }
        let keq = self.equilibrium_constant(temperature);
        if keq <= 0.0 {
            return Ok(forward);
        }
        let mut reverse = 1.0;
        for (component, order) in &self.product_orders {
            // **The reverse term floors at `1e-30` where the forward one floors at zero.** The
            // class writes `max(C, 1e-30)` here (`KineticReaction.java:289`) and `max(C, 0.0)`
            // there, so an absent product makes the forward rate zero and leaves the reverse
            // finite.
            reverse *= concentration(component).max(1.0e-30).powf(*order);
        }
        Ok(forward - k * reverse / keq)
    }

    fn lhhw(&self, temperature: f64, k: f64, concentration: &impl Fn(&str) -> f64) -> Result<f64> {
        let Some(numerator) = self.forward(k, concentration) else {
            return Ok(0.0);
        };
        let k_ads = self.adsorption_pre_exponential_factor
            * (-self.adsorption_activation_energy / (R_GAS * temperature)).exp();
        let mut sum = 1.0;
        for (component, (factor, order)) in &self.adsorption_terms {
            sum += factor * k_ads * concentration(component).max(0.0).powf(*order);
        }
        // The class floors the denominator at `1e-30` rather than refusing.
        let denominator = sum.powi(self.adsorption_exponent).max(1.0e-30);
        Ok(numerator / denominator)
    }
}
