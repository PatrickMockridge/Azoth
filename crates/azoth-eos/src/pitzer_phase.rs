//! Pitzer's activity-coefficient arithmetic, as NeqSim's `PhasePitzer` and
//! `ComponentGePitzer` carry it.
//!
//! The model is Harvie and Weare's (1980) mixed-electrolyte formulation: a single-ion
//! activity coefficient built from a Debye-Huckel term and binary cation-anion
//! interactions, plus same-sign `theta` and ternary `psi` mixing terms, plus an optional
//! neutral-solute layer.
//!
//! # The pieces, and the order they are ported in
//!
//! The whole is large and most of it is private, so `validation/neqsim/PitzerArithmetic.java`
//! prints the intermediates one at a time and this module grows against that oracle:
//!
//! * [`debye_huckel_a_phi`] - the Debye-Huckel parameter `A_phi(T)`.
//! * [`g`] and [`g_prime`] - the two functions every binary term is built from.
//!
//! # The name that does not mean what it says
//!
//! NeqSim's `debyeHuckelAphi(T)` **returns `3 * A_phi`**, and every caller divides by
//! three to get back what the name promises - `getGamma` and `getWaterGamma` both carry a
//! comment saying so. This port returns `A_phi` and drops the round trip, because a
//! function named for one number and returning another is a factor of three waiting to
//! happen, and three times an activity coefficient is a plausible number rather than a
//! failure. The oracle's own two columns are the check that the two agree.

/// The Debye-Huckel parameter `A_phi(T)`, in `(mol/kg)^-1/2`.
///
/// `ComponentGePitzer.debyeHuckelAphi`, divided by the three its callers divide by. Its
/// inputs, in the order NeqSim assembles them:
///
/// * **Water density**, Kell (1975)'s polynomial for `0..150 C`, overwritten above
///   `100 C` by an IAPWS-style approximation because the Kell form under-predicts a
///   pressurised liquid, and floored at `700 kg/m3`.
/// * **The dielectric constant**, Archer and Wang (1990), floored at `20`.
/// * `A_phi = 1.4006e6 sqrt(rho_g_per_cm3) / (eps T)^1.5`.
///
/// **There is no clamp on the temperature below `0 C`** - the Kell polynomial is
/// extrapolated - and that is NeqSim's behaviour rather than an oversight here. The
/// floor on the density is what keeps the bracketed terms from diverging at the cold end.
#[must_use]
pub fn debye_huckel_a_phi(t: f64) -> f64 {
    let celsius = t - 273.15;
    let mut density = 999.83 + 5.0948e-2 * celsius - 7.5722e-3 * celsius * celsius
        + 3.8907e-5 * celsius * celsius * celsius
        - 1.2e-7 * celsius * celsius * celsius * celsius;
    if celsius > 100.0 {
        let excess = celsius - 100.0;
        density = 958.0 - 1.08 * excess - 0.0028 * excess * excess;
    }
    if density < 700.0 {
        density = 700.0;
    }
    let per_cm3 = density / 1000.0;

    let mut dielectric = 87.740 - 0.40008 * celsius + 9.398e-4 * celsius * celsius
        - 1.410e-6 * celsius * celsius * celsius;
    if dielectric < 20.0 {
        dielectric = 20.0;
    }

    1.4006e6 * per_cm3.sqrt() / (dielectric * t).powf(1.5)
}

/// Below this `x` the two `g` functions are zero rather than their limit.
///
/// NeqSim's own guard (`x1 > 1e-12`), kept as a guard rather than replaced by the limit
/// `g(0) = 1`: the two differ when `x` is small and not zero, and reproducing the branch
/// is the point.
const G_FLOOR: f64 = 1.0e-12;

/// `g(x) = 2 (1 - (1 + x) exp(-x)) / x^2`, the Pitzer binary-term function.
///
/// Zero below [`G_FLOOR`], which is NeqSim's branch and not the limit.
#[must_use]
pub fn g(x: f64) -> f64 {
    if x <= G_FLOOR {
        return 0.0;
    }
    2.0 * (1.0 - (1.0 + x) * (-x).exp()) / (x * x)
}

/// `g'(x) = -2 (1 - (1 + x + x^2/2) exp(-x)) / x^2`, the function's derivative.
///
/// Zero below [`G_FLOOR`], as [`g`] is.
#[must_use]
pub fn g_prime(x: f64) -> f64 {
    if x <= G_FLOOR {
        return 0.0;
    }
    -2.0 * (1.0 - (1.0 + x + x * x / 2.0) * (-x).exp()) / (x * x)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The oracle's own column, temperature by temperature.
    ///
    /// `validation/neqsim/PitzerArithmetic.java` prints `Aphi = Agamma / 3` for each of
    /// these, and the two agree to the last digit it prints - so this is a check against
    /// NeqSim rather than against a literature value. The literature's `A_phi` at 298.15 K
    /// is `0.3915` against NeqSim's `0.392034`, and the difference is NeqSim's.
    #[test]
    fn the_debye_huckel_parameter_matches_neqsim() {
        for (t, expected) in [
            (273.15, 0.377466966949275),
            (283.15, 0.382907418517215),
            (298.15, 0.392034451863750),
            (313.15, 0.402349556023220),
            (323.15, 0.409909254205050),
            (348.15, 0.431348257867475),
            (373.15, 0.456800085250619),
            (398.15, 0.485617249619569),
        ] {
            let got = debye_huckel_a_phi(t);
            assert!(
                (got - expected).abs() < 1.0e-14,
                "A_phi({t}) = {got}, and NeqSim gives {expected}"
            );
        }
    }

    /// **The density floor is what keeps the cold end finite**, and it is reachable.
    ///
    /// The Kell polynomial is extrapolated below `0 C` - NeqSim has no clamp there - and at
    /// `250 K` below freezing it has gone **negative**, which is what the floor catches. A
    /// port that dropped the floor would take the square root of a negative number and
    /// return `NaN`, so this is the difference between a number and a failure.
    ///
    /// **The dielectric floor does not apply at that temperature**, which the first draft
    /// of this test assumed. The Archer-Wang polynomial's cubic term dominates once
    /// Celsius is large and negative, so `eps` *grows* - `268.5` at `-250 C`, far above
    /// the `20` floor.
    #[test]
    fn the_density_floor_is_reachable_and_the_dielectric_floor_is_not() {
        let celsius: f64 = -250.0;
        let t = celsius + 273.15;

        let unfloored = 999.83 + 5.0948e-2 * celsius - 7.5722e-3 * celsius * celsius
            + 3.8907e-5 * celsius * celsius * celsius
            - 1.2e-7 * celsius * celsius * celsius * celsius;
        assert!(
            unfloored < 0.0,
            "the Kell polynomial is negative here ({unfloored}), which is what the floor is for"
        );

        let dielectric = 87.740 - 0.40008 * celsius + 9.398e-4 * celsius * celsius
            - 1.410e-6 * celsius * celsius * celsius;
        assert!(
            dielectric > 20.0,
            "the dielectric polynomial grows at the cold end ({dielectric})"
        );

        let expected = 1.4006e6 * 0.7f64.sqrt() / (dielectric * t).powf(1.5);
        let got = debye_huckel_a_phi(t);
        assert!(
            (got - expected).abs() < 1.0e-12,
            "A_phi = {got}, and the floored density gives {expected}"
        );
    }

    /// The two `g` functions against the oracle's table, at the same `x`.
    #[test]
    fn the_binary_functions_match_neqsim() {
        for (x, want_g, want_g_prime) in [
            (0.10, 0.935768032088879, -0.0309306140529486),
            (0.50, 0.721632083448399, -0.115101423735766),
            (1.00, 0.528482235314231, -0.160602794142788),
            (2.00, 0.296997075145081, -0.161661791908468),
            (4.00, 0.113552725694541, -0.0952370868058070),
            (8.00, 0.0311556511359024, -0.0308201885079999),
            (12.00, 0.0138877795172140, -0.0138816353048607),
        ] {
            assert!(
                (g(x) - want_g).abs() < 1.0e-14,
                "g({x}) = {}, and NeqSim gives {want_g}",
                g(x)
            );
            assert!(
                (g_prime(x) - want_g_prime).abs() < 1.0e-14,
                "g'({x}) = {}, and NeqSim gives {want_g_prime}",
                g_prime(x)
            );
        }
    }

    /// **The guard is at `1e-12`, and the formula is already meaningless at `1e-8`.**
    ///
    /// `g`'s numerator is `1 - (1 + x) exp(-x)`, which is `x^2/2` evaluated as a
    /// cancellation of two numbers near one. That costs about sixteen digits, so the
    /// expression returns garbage long before `x` reaches the guard: measured,
    /// `g(1e-8) = 2.22` and `g(1e-9) = 0.0`, where the true value is `1`. The guard does
    /// not protect against this - it is seven orders of magnitude too low.
    ///
    /// **It is not reachable in practice, and that is measured too.** `x = alpha sqrt(I)`
    /// with `alpha >= 1.4`, so `x < 1e-5` needs `I < 2.5e-11` mol/kg. A sodium chloride
    /// mole fraction of `1e-9` gives `I = 5.55e-8` and `x = 4.7e-4` - four orders clear.
    /// Reaching the broken region needs a mole fraction near `1e-13`.
    ///
    /// So this test records the behaviour rather than asserting it is right: the port
    /// reproduces NeqSim, and the defect is NeqSim's.
    #[test]
    fn the_small_x_guard_is_below_where_the_formula_stops_working() {
        assert_eq!(g(0.0), 0.0);
        assert_eq!(g(1.0e-13), 0.0, "NeqSim's own branch, and not the limit");
        assert_eq!(g_prime(0.0), 0.0);

        // The region between the guard and the precision floor, where the value is a
        // number with no relationship to `g`.
        assert!(g(1.0e-8) > 2.0, "cancellation, not a value: {}", g(1.0e-8));
        assert_eq!(g(1.0e-9), 0.0, "and here it is exactly zero");

        // And where it is trustworthy, which is the ordinary range.
        for x in [1.0e-3f64, 1.0e-2, 0.5, 1.0] {
            let expected = 2.0 * (1.0 - (1.0 + x) * (-x).exp()) / (x * x);
            assert!((g(x) - expected).abs() < 1.0e-15);
        }
    }
}

/// The temperature both of NeqSim's Pitzer forms state their value at, in K.
///
/// `PitzerParameterDatasets.PHREEQC_REFERENCE_TEMPERATURE_K` and the `298.15` the
/// Silvester form divides by are the same number, which is why one constant serves both.
pub const REFERENCE_TEMPERATURE_K: f64 = 298.15;

/// Within this many kelvin of the reference, the catalogue form returns its constant term
/// unchanged rather than evaluating the polynomial.
///
/// `PitzerTemperatureFunction.REFERENCE_TOLERANCE_K`. It is a *branch* and not a rounding
/// guard: at `298.1501 K` the six coefficients evaluate to `a0` plus about `1e-8`, and
/// NeqSim returns exactly `a0`. Measured on the built phase, which is where the oracle for
/// it lives.
const REFERENCE_TOLERANCE_K: f64 = 1.0e-3;

/// How a Pitzer parameter varies with temperature.
///
/// **The two datasets use different forms, and this is where that shows.** The PHREEQC
/// catalogue carries six coefficients and evaluates a polynomial; `PitzerParameters.csv`
/// carries two and evaluates Silvester and Pitzer's log form. NeqSim keeps them in one
/// class behind two getters, so a pair knows which it is only by which dataset answered.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum TemperatureForm {
    /// The PHREEQC six-coefficient form `a0`..`a5`, from the vendored catalogue.
    Catalog([f64; 6]),
    /// Silvester and Pitzer's `value_25 + t1 (1/T - 1/Tr) + t2 ln(T/Tr)`, from the CSV.
    ///
    /// **Zero `t1` and `t2` is the flat case**, and NeqSim returns `value_25` for it
    /// rather than evaluating the sum - which differs from the sum only in the last bit,
    /// but differs, so the branch is reproduced.
    Silvester { at_25: f64, t1: f64, t2: f64 },
}

impl TemperatureForm {
    /// The parameter at a temperature.
    #[must_use]
    pub fn value_at(&self, t: f64) -> f64 {
        match *self {
            TemperatureForm::Catalog(a) => catalog_value(&a, t),
            TemperatureForm::Silvester { at_25, t1, t2 } => {
                if t1.abs() < 1.0e-20 && t2.abs() < 1.0e-20 {
                    return at_25;
                }
                at_25
                    + t1 * (1.0 / t - 1.0 / REFERENCE_TEMPERATURE_K)
                    + t2 * (t / REFERENCE_TEMPERATURE_K).ln()
            }
        }
    }
}

/// `PitzerTemperatureFunction.valueAt`, whose form is
/// `a0 + a1(1/T - 1/Tr) + a2 ln(T/Tr) + a3(T - Tr) + a4(T^2 - Tr^2) + a5(1/T^2 - 1/Tr^2)`.
fn catalog_value(a: &[f64; 6], t: f64) -> f64 {
    if (t - REFERENCE_TEMPERATURE_K).abs() < REFERENCE_TOLERANCE_K {
        return a[0];
    }
    let inverse = 1.0 / t;
    let inverse_reference = 1.0 / REFERENCE_TEMPERATURE_K;
    a[0] + a[1] * (inverse - inverse_reference)
        + a[2] * (t / REFERENCE_TEMPERATURE_K).ln()
        + a[3] * (t - REFERENCE_TEMPERATURE_K)
        + a[4] * (t * t - REFERENCE_TEMPERATURE_K * REFERENCE_TEMPERATURE_K)
        + a[5] * (inverse * inverse - inverse_reference * inverse_reference)
}

/// `PhasePitzer.getPitzerAlpha1`'s **bounded** form: `1.4` for a 2:2 pair, else `2.0`.
///
/// **Three call sites in NeqSim use an unbounded form instead** - `getGamma`,
/// `getWaterGamma` and `phreeqcBinaryBprime` all test `|z| >= 1.5` with no upper bound,
/// so a 3-valent ion gets `1.4` there and `2.0` here. The two agree for every charge the
/// vendored data carries - the highest is `2` - so this is a latent inconsistency rather
/// than a defect, and it is recorded rather than resolved. [`alpha1_as_used`] is what the
/// activity coefficient is actually built from.
#[must_use]
pub fn alpha1(first_charge: f64, second_charge: f64) -> f64 {
    let bounded = |z: f64| (1.5..2.5).contains(&z);
    if bounded(first_charge.abs()) && bounded(second_charge.abs()) {
        1.4
    } else {
        2.0
    }
}

/// The `alpha1` the activity coefficient is actually built from, at the three sites that
/// use the unbounded test. See [`alpha1`] for why both exist.
#[must_use]
pub fn alpha1_as_used(first_charge: f64, second_charge: f64) -> f64 {
    if first_charge.abs() >= 1.5 && second_charge.abs() >= 1.5 {
        1.4
    } else {
        2.0
    }
}

/// `PhasePitzer.getPitzerAlpha2`: `12.0` when either ion is monovalent **or** the pair is
/// 2:2, else `50.0`.
///
/// The 2:2 clause is not redundant with the monovalent one and is why `CaCl2`'s qualified
/// `B2` row survives: that row is a 2:1 term with `alpha2 = 12`, and a 2:2-only branch
/// would have discarded it.
#[must_use]
pub fn alpha2(first_charge: f64, second_charge: f64) -> f64 {
    let (one, two) = (first_charge.abs(), second_charge.abs());
    let monovalent = one < 1.5 || two < 1.5;
    let two_two = (1.5..2.5).contains(&one) && (1.5..2.5).contains(&two);
    if monovalent || two_two { 12.0 } else { 50.0 }
}

#[cfg(test)]
mod parameter_tests {
    use super::*;

    /// **Both temperature forms against the built phase.**
    ///
    /// `validation/neqsim/PitzerArithmetic.java` prints each across temperature from a
    /// live `SystemPitzer`, so these are NeqSim's numbers rather than a restatement of
    /// its formula. The catalogue pair is `Na+/Cl-` (the PHREEQC six coefficients) and the
    /// CSV pair is `Na+/HCO3-`, whose two coefficients are zero and which is therefore
    /// flat.
    #[test]
    fn both_temperature_forms_match_neqsim() {
        // Na+/Cl- from the PHREEQC catalogue: a0 = 0.07534, a1 = 9598.4, ... (the row
        // `Cl-|Na+` of B0). The other coefficients come from the compiled table.
        let catalogue =
            crate::pitzer_catalog::find(crate::pitzer_catalog::Family::B0, &["Na+", "Cl-"])
                .expect("the catalogue carries Na+/Cl-");
        let form = TemperatureForm::Catalog(catalogue);

        // At the reference the value is `a0`, and it stays `a0` inside the 1e-3 K window.
        assert!((form.value_at(298.15) - 0.07534).abs() < 1.0e-15);
        assert!(
            (form.value_at(298.1501) - 0.07534).abs() < 1.0e-15,
            "within the reference tolerance the constant term is returned unchanged"
        );

        // And outside it, the oracle's own column.
        for (t, expected) in [
            (273.15, 0.0493895679117),
            (323.15, 0.0892387746618),
            (373.15, 0.100154205617),
        ] {
            let got = form.value_at(t);
            assert!(
                (got - expected).abs() < 1.0e-12,
                "beta0({t}) = {got}, and NeqSim gives {expected}"
            );
        }

        // The CSV's flat case: `Na+/HCO3-` carries beta0_25 = 0.0277 with no temperature
        // coefficients, so every temperature gives the same number.
        let flat = TemperatureForm::Silvester {
            at_25: 0.0277,
            t1: 0.0,
            t2: 0.0,
        };
        for t in [273.15, 298.15, 323.15, 373.15] {
            assert_eq!(flat.value_at(t), 0.0277);
        }
    }

    /// The two `alpha` functions, at the charges the vendored data carries.
    #[test]
    fn the_alpha_coefficients_follow_the_charge() {
        // Monovalent against monovalent: alpha1 = 2.0, alpha2 = 12.0. The oracle prints
        // exactly this for `Na+/Cl-`.
        assert_eq!(alpha1_as_used(1.0, -1.0), 2.0);
        assert_eq!(alpha2(1.0, -1.0), 12.0);

        // A 2:2 pair: alpha1 = 1.4, alpha2 = 12.0 (the redundant-looking second clause).
        assert_eq!(alpha1_as_used(2.0, -2.0), 1.4);
        assert_eq!(alpha2(2.0, -2.0), 12.0);

        // A 2:1 pair: alpha1 = 2.0, alpha2 = 12.0, which is what keeps CaCl2's B2 row.
        assert_eq!(alpha1_as_used(2.0, -1.0), 2.0);
        assert_eq!(alpha2(2.0, -1.0), 12.0);

        // A 3:1 pair: alpha1 = 2.0 and alpha2 = 50.0, since nothing is monovalent.
        assert_eq!(alpha1_as_used(3.0, -1.0), 2.0);
        assert_eq!(alpha2(3.0, -1.0), 12.0, "the chloride is monovalent");

        // **And the one charge where the two alpha1 definitions part.** A 3-valent pair
        // gets 1.4 from the sites that compute the activity coefficient and 2.0 from the
        // getter. Nothing vendored carries a charge above 2, so this is unreachable
        // rather than wrong.
        assert_eq!(alpha1_as_used(3.0, -3.0), 1.4);
        assert_eq!(alpha1(3.0, -3.0), 2.0);
        assert_ne!(alpha1_as_used(3.0, -3.0), alpha1(3.0, -3.0));
    }
}

/// `b`, the Debye-Huckel denominator constant, the same `1.2` in every Pitzer expression
/// NeqSim carries.
pub const B: f64 = 1.2;

/// Below this an ionic strength or an `x` makes a derivative term zero rather than small.
const DERIVATIVE_FLOOR: f64 = 1.0e-12;

/// Below this a fitted `beta2` is treated as absent.
const BETA2_FLOOR: f64 = 1.0e-20;

/// Below this a component's charge makes it a neutral rather than an ion.
const ION_CHARGE: f64 = 0.5;

/// The pair parameters the activity coefficient reads, from whichever dataset is in force.
///
/// **A trait and not a struct, because the two datasets answer differently.** The PHREEQC
/// catalogue evaluates six coefficients through [`TemperatureForm::Catalog`] and
/// `PitzerParameters.csv` evaluates Silvester-Pitzer's two; both are reachable from the
/// selection rule, and which one answers a given pair is a property of the phase's
/// topology rather than of the pair.
pub trait PairParameters {
    /// `PhasePitzer.getBeta0ij(i, j, T)`.
    fn beta0(&self, first: usize, second: usize, temperature: f64) -> f64;
    /// `PhasePitzer.getBeta1ij(i, j, T)`.
    fn beta1(&self, first: usize, second: usize, temperature: f64) -> f64;
    /// `PhasePitzer.getCphiij(i, j, T)`.
    fn cphi(&self, first: usize, second: usize, temperature: f64) -> f64;
    /// `PhasePitzer.getBeta2ij(i, j, T)`.
    fn beta2(&self, first: usize, second: usize, temperature: f64) -> f64;
    /// `PhasePitzer.getThetaij(i, j, T)`.
    fn theta(&self, first: usize, second: usize, temperature: f64) -> f64;
    /// `PhasePitzer.getPsiijk(i, j, k, T)`.
    fn psi(&self, first: usize, second: usize, third: usize, temperature: f64) -> f64;
}

/// The phase state the ion branch is a function of, and the three topology flags that
/// decide which of its terms are live.
///
/// The flags are NeqSim's, read off the phase: `isPhreeqcCommonIonTermsActive`,
/// `isNonTwoTwoBeta2Active` and `hasUnequalChargeSameSignPair`. They are not derivable
/// from the composition - `nonTwoTwoBeta2Active` is a property of which rows the dataset
/// *configured*, not of which ions are present - which is why they are passed in.
#[derive(Debug, Clone, Copy)]
pub struct IonActivityContext<'a> {
    /// Molality per component, from [`crate::electrolyte::composition`].
    pub molality: &'a [f64],
    /// Ionic charge per component, in units of the elementary charge.
    pub charge: &'a [f64],
    /// Temperature, in K.
    pub temperature: f64,
    /// `A_phi` at that temperature, from [`debye_huckel_a_phi`].
    pub a_phi: f64,
    /// `PhasePitzer.isPhreeqcCommonIonTermsActive`.
    pub common_ion_terms: bool,
    /// `PhasePitzer.isNonTwoTwoBeta2Active`.
    pub non_two_two_beta2: bool,
    /// `PhasePitzer.hasUnequalChargeSameSignPair`.
    pub unequal_charge_same_sign: bool,
    /// Whether the neutral Pitzer layer is active, `ComponentGePitzer`'s
    /// `neutralPitzerInteractionsActive`.
    ///
    /// Read by the water node as well as the ion one, and for a second reason there: it
    /// decides whether `sumMolalities` counts the neutral solutes or only the ions.
    pub neutral_interactions_active: bool,
}

impl IonActivityContext<'_> {
    /// `I = 1/2 sum m_i z_i^2`, recomputed here so the context cannot disagree with
    /// [`crate::electrolyte::ionic_strength`] about what it holds.
    #[must_use]
    pub fn ionic_strength(&self) -> f64 {
        crate::electrolyte::ionic_strength(self.molality, self.charge)
    }

    /// `Z = sum m_i |z_i|` over the ions, the `Z` of the `Z C` term.
    #[must_use]
    pub fn charge_sum(&self) -> f64 {
        self.molality
            .iter()
            .zip(self.charge)
            .filter(|&(_, z)| z.abs() > ION_CHARGE)
            .map(|(m, z)| m * z.abs())
            .sum()
    }
}

/// `ln gamma_M` for one ion, as `ComponentGePitzer.getGamma`'s ion branch computes it.
///
/// NeqSim's own comment gives the shape:
///
/// ```text
/// ln(gamma_M) = z_M^2 F + sum_a m_a (2 B_Ma + Z C_Ma)
/// F = -A_phi [ sqrt(I)/(1 + b sqrt(I)) + (2/b) ln(1 + b sqrt(I)) ] + sum_c sum_a m_c m_a B'_ca
/// ```
///
/// and the pieces that are easy to lose are the three the flags gate: `fBprime` is *not*
/// accumulated when the common-ion terms are active, because
/// [`common_ion_contribution`] computes its own; `beta2` is applied only for a 2:2 pair or
/// when the dataset activated it elsewhere; and `E_theta` enters twice, once through `F`
/// and once through the theta sum.
#[must_use]
pub fn ln_gamma(
    context: &IonActivityContext<'_>,
    parameters: &dyn PairParameters,
    component: usize,
) -> f64 {
    let charge = context.charge[component];
    let i = context.ionic_strength();
    let sqrt_i = i.sqrt();
    let a_phi = context.a_phi;

    // f^phi = -A_phi [ sqrt(I)/(1 + b sqrt(I)) + (2/b) ln(1 + b sqrt(I)) ]
    let debye_huckel = -a_phi * (sqrt_i / (1.0 + B * sqrt_i) + (2.0 / B) * (1.0 + B * sqrt_i).ln());

    let z_sum = context.charge_sum();
    let common_ion = if context.common_ion_terms {
        common_ion_contribution(context, parameters, charge)
    } else {
        0.0
    };

    let m_this = context.molality[component];
    let mut sum = 0.0;
    let mut b_prime_sum = 0.0;

    // The binary terms: this ion against every opposite-sign ion.
    for j in 0..context.molality.len() {
        if j == component {
            continue;
        }
        let other_charge = context.charge[j];
        if other_charge * charge >= 0.0 {
            continue;
        }
        let m_j = context.molality[j];
        let beta0 = parameters.beta0(component, j, context.temperature);
        let beta1 = parameters.beta1(component, j, context.temperature);
        let cphi = parameters.cphi(component, j, context.temperature);

        // **The unbounded charge test**, which is what the three sites that compute an
        // activity coefficient use. See [`alpha1`] and [`alpha1_as_used`].
        let alpha_one = alpha1_as_used(charge, other_charge);
        let is_two_two = charge.abs() >= 1.5 && other_charge.abs() >= 1.5;
        let x1 = alpha_one * sqrt_i;
        let mut b_value = beta0 + beta1 * g(x1);
        let mut b_derivative = if i > DERIVATIVE_FLOOR {
            beta1 * g_prime(x1) / i
        } else {
            0.0
        };

        if is_two_two || context.non_two_two_beta2 {
            let beta2 = parameters.beta2(component, j, context.temperature);
            if beta2.abs() > BETA2_FLOOR {
                let x2 = alpha2(charge, other_charge) * sqrt_i;
                b_value += beta2 * g(x2);
                if i > DERIVATIVE_FLOOR {
                    b_derivative += beta2 * g_prime(x2) / i;
                }
            }
        }

        // C_ca = Cphi_ca / (2 sqrt(|z_c z_a|))
        let c_value = cphi / (2.0 * (charge * other_charge).abs().sqrt());
        sum += m_j * (2.0 * b_value + z_sum * c_value);

        if !context.common_ion_terms {
            b_prime_sum += m_this * m_j * b_derivative;
        }
    }

    // E_theta's derivative through F, over the unequal-charge same-sign pairs.
    let mut e_theta_prime = 0.0;
    if context.unequal_charge_same_sign {
        for first in 0..context.molality.len() {
            let first_charge = context.charge[first];
            if first_charge.abs() < ION_CHARGE {
                continue;
            }
            for second in (first + 1)..context.molality.len() {
                let second_charge = context.charge[second];
                if first_charge * second_charge <= 0.0
                    || (first_charge - second_charge).abs() < DERIVATIVE_FLOOR
                {
                    continue;
                }
                let (_, derivative) =
                    crate::pitzer_electrostatic::calculate(first_charge, second_charge, i, a_phi);
                e_theta_prime += context.molality[first] * context.molality[second] * derivative;
            }
        }
    }

    // The same-sign theta and opposite-sign psi terms.
    for j in 0..context.molality.len() {
        if j == component {
            continue;
        }
        let other_charge = context.charge[j];
        if other_charge.abs() < ION_CHARGE || other_charge * charge <= 0.0 {
            continue;
        }
        let m_j = context.molality[j];
        let theta = parameters.theta(component, j, context.temperature);
        let mut e_theta = 0.0;
        if context.unequal_charge_same_sign && (charge - other_charge).abs() >= DERIVATIVE_FLOOR {
            let (value, _) = crate::pitzer_electrostatic::calculate(charge, other_charge, i, a_phi);
            e_theta = value;
        }
        sum += m_j * 2.0 * (theta + e_theta);

        for k in 0..context.molality.len() {
            let third_charge = context.charge[k];
            if third_charge * charge >= 0.0 || third_charge.abs() < ION_CHARGE {
                continue;
            }
            let m_k = context.molality[k];
            let psi = parameters.psi(component, j, k, context.temperature);
            sum += m_j * m_k * psi;
        }
    }

    let f = debye_huckel + b_prime_sum + e_theta_prime;
    charge * charge * f + sum + common_ion
}

/// PHREEQC's common B-prime and C0 contributions to one ion's `ln gamma`.
///
/// `ComponentGePitzer.phreeqcCommonIonContribution`, which replaces the `fBprime` the main
/// loop would otherwise accumulate: it sums over **every** cation-anion pair rather than
/// over this ion's pairs, and adds a `C0` term the main loop has no counterpart for.
fn common_ion_contribution(
    context: &IonActivityContext<'_>,
    parameters: &dyn PairParameters,
    target_charge: f64,
) -> f64 {
    let i = context.ionic_strength();
    let sqrt_i = i.sqrt();
    let mut b_prime = 0.0;
    let mut cphi_sum = 0.0;
    for cation in 0..context.molality.len() {
        let cation_charge = context.charge[cation];
        if cation_charge <= 0.0 {
            continue;
        }
        for anion in 0..context.molality.len() {
            let anion_charge = context.charge[anion];
            if anion_charge >= 0.0 {
                continue;
            }
            let first = context.molality[cation];
            let second = context.molality[anion];
            b_prime +=
                first * second * binary_b_prime(context, parameters, cation, anion, i, sqrt_i);
            cphi_sum += first * second * parameters.cphi(cation, anion, context.temperature)
                / (2.0 * (cation_charge * anion_charge).abs().sqrt());
        }
    }
    target_charge * target_charge * b_prime + target_charge.abs() * cphi_sum
}

/// One pair's `B'` as the common-ion sum computes it: `PitzerParameterDatasets`' shape,
/// which differs from the main loop's in that it applies `beta2` wherever the catalogue
/// carries one rather than gating on the 2:2 test.
fn binary_b_prime(
    context: &IonActivityContext<'_>,
    parameters: &dyn PairParameters,
    first: usize,
    second: usize,
    i: f64,
    sqrt_i: f64,
) -> f64 {
    let first_charge = context.charge[first];
    let second_charge = context.charge[second];
    let x1 = alpha1_as_used(first_charge, second_charge) * sqrt_i;
    let mut derivative = 0.0;
    if x1 > DERIVATIVE_FLOOR && i > DERIVATIVE_FLOOR {
        derivative = parameters.beta1(first, second, context.temperature) * g_prime(x1) / i;
    }
    let beta2 = parameters.beta2(first, second, context.temperature);
    let x2 = alpha2(first_charge, second_charge) * sqrt_i;
    if beta2.abs() > BETA2_FLOOR && x2 > DERIVATIVE_FLOOR && i > DERIVATIVE_FLOOR {
        derivative += beta2 * g_prime(x2) / i;
    }
    derivative
}

#[cfg(test)]
mod ion_tests {
    use super::*;
    use crate::pitzer_catalog::{Family, find};

    /// One pair's parameter as the PHREEQC catalogue states it: six coefficients, or zero
    /// where the catalogue carries no row.
    ///
    /// A zero for an absent row is NeqSim's own behaviour - its parameter arrays are
    /// zero-initialised and only the rows it read are written - so this reproduces the
    /// arithmetic rather than papering over a gap. A gap that *matters* is caught earlier,
    /// by the selection rule's coverage requirement.
    pub(super) fn form(family: Family, names: &[&str]) -> TemperatureForm {
        match find(family, names) {
            Some(a) => TemperatureForm::Catalog(a),
            None => TemperatureForm::Catalog([0.0; 6]),
        }
    }

    /// The catalogue-backed [`PairParameters`], which is what a selected brine reads.
    pub(super) struct CatalogParameters {
        species: Vec<String>,
    }

    impl CatalogParameters {
        pub(super) fn new(names: &[&str]) -> Self {
            Self {
                species: names
                    .iter()
                    .map(|name| crate::pitzer_catalog::canonical_species(name))
                    .collect(),
            }
        }

        fn pair(&self, first: usize, second: usize) -> [&str; 2] {
            [&self.species[first], &self.species[second]]
        }
    }

    impl PairParameters for CatalogParameters {
        fn beta0(&self, first: usize, second: usize, t: f64) -> f64 {
            form(Family::B0, &self.pair(first, second)).value_at(t)
        }
        fn beta1(&self, first: usize, second: usize, t: f64) -> f64 {
            form(Family::B1, &self.pair(first, second)).value_at(t)
        }
        fn cphi(&self, first: usize, second: usize, t: f64) -> f64 {
            form(Family::C0, &self.pair(first, second)).value_at(t)
        }
        fn beta2(&self, first: usize, second: usize, t: f64) -> f64 {
            form(Family::B2, &self.pair(first, second)).value_at(t)
        }
        fn theta(&self, first: usize, second: usize, t: f64) -> f64 {
            form(Family::Theta, &self.pair(first, second)).value_at(t)
        }
        fn psi(&self, first: usize, second: usize, third: usize, t: f64) -> f64 {
            let names = [
                self.species[first].as_str(),
                self.species[second].as_str(),
                self.species[third].as_str(),
            ];
            form(Family::Psi, &names).value_at(t)
        }
    }

    /// The molalities of a composition, per kilogram of water.
    ///
    /// One mole of mixture's worth: `m_i = x_i / (x_water M_water)`, with the water molar
    /// mass the databank carries.
    pub(super) fn molalities(x: &[f64]) -> Vec<f64> {
        let solvent = 0.018015;
        let mass_of_water = x[0] * solvent;
        x.iter().map(|&xi| xi / mass_of_water).collect()
    }

    /// **The end-to-end oracle**, from `validation/neqsim/PitzerArithmetic.java` on a live
    /// `SystemPitzer` at 298.15 K.
    ///
    /// The compositions are the probe's own. The first brine has no same-sign pair and no
    /// non-2:2 `beta2`, so it exercises the plain path; the second has both, and is the
    /// case the two flags exist for.
    #[test]
    fn the_ion_activity_coefficient_matches_neqsim() {
        // water + Na+ + Cl- at mole fractions 0.88 / 0.06 / 0.06, and the molalities
        // **computed from them** rather than copied from the probe's twelve-digit
        // printout - which is what the first draft did, and it failed in the thirteenth
        // digit of `I`.
        let names = ["water", "Na+", "Cl-"];
        let mole_fraction = [0.88, 0.06, 0.06];
        let molality = molalities(&mole_fraction);
        let charge = [0.0, 1.0, -1.0];
        let parameters = CatalogParameters::new(&names);
        let context = IonActivityContext {
            molality: &molality,
            charge: &charge,
            temperature: 298.15,
            a_phi: debye_huckel_a_phi(298.15),
            common_ion_terms: true,
            non_two_two_beta2: false,
            unequal_charge_same_sign: false,
            neutral_interactions_active: false,
        };
        assert!(
            (context.ionic_strength() - 3.784_724_850_503_37).abs() < 1.0e-12,
            "I = {}",
            context.ionic_strength()
        );
        for (index, expected) in [(1, -0.267_512_432_075_052), (2, -0.267_512_432_075_052)] {
            let got = ln_gamma(&context, &parameters, index);
            assert!(
                (got - expected).abs() < 1.0e-10,
                "ln gamma({}) = {got}, and NeqSim gives {expected}",
                names[index]
            );
        }
    }

    /// The brine with a same-sign pair of unequal charge and a non-2:2 `beta2`, which is
    /// the case `CaCl2` in the catalogue creates.
    #[test]
    fn the_same_sign_and_beta2_terms_match_neqsim() {
        // The probe adds `Cl-` twice and NeqSim aggregates the two entries, so the
        // effective composition is water 0.88 / Na+ 0.03 / Ca++ 0.03 / Cl- 0.06.
        let names = ["water", "Na+", "Ca++", "Cl-"];
        let mole_fraction = [0.88, 0.03, 0.03, 0.06];
        let molality = molalities(&mole_fraction);
        let charge = [0.0, 1.0, 2.0, -1.0];
        let parameters = CatalogParameters::new(&names);
        let context = IonActivityContext {
            molality: &molality,
            charge: &charge,
            temperature: 298.15,
            a_phi: debye_huckel_a_phi(298.15),
            common_ion_terms: true,
            non_two_two_beta2: true,
            unequal_charge_same_sign: true,
            neutral_interactions_active: false,
        };
        assert!((context.ionic_strength() - 6.623_268_488_380_89).abs() < 1.0e-12);

        for (index, expected) in [
            (1, -0.504_259_131_023_400),
            (2, -1.841_321_926_770_47),
            (3, 0.726_599_630_955_942),
        ] {
            let got = ln_gamma(&context, &parameters, index);
            assert!(
                (got - expected).abs() < 1.0e-10,
                "ln gamma({}) = {got}, and NeqSim gives {expected}",
                names[index]
            );
        }
    }
}

/// The water molar mass NeqSim hard-codes in the osmotic expression, in kg/mol.
///
/// `getWaterGamma` writes `double Mw = 18.015` in g/mol and divides by `1000`, and the
/// databank's water row is `0.018015 kg/mol` - the two agree exactly, so this is a named
/// restatement rather than a second value, and
/// [`crate::electrolyte::ln_water_activity`] takes it as an argument for that reason.
pub const WATER_MOLAR_MASS: f64 = 0.018015;

/// `ln gamma_w`, as `ComponentGePitzer.getWaterGamma` computes it.
///
/// **This route does not go through the ion expression at all.** The solvent's activity
/// coefficient comes from the Pitzer *osmotic* coefficient, and the osmotic coefficient is
/// built from a different binary function:
///
/// ```text
/// phi - 1 = (2 / sum m) [ -A_phi I^1.5/(1 + b sqrt(I)) + sum_c sum_a m_c m_a (B^phi_ca + Z C_ca)
///                         + theta + psi ]          with   B^phi_ca = beta0 + beta1 exp(-alpha1 sqrt(I))
/// ln a_w = -phi M_w sum m
/// ln gamma_w = ln a_w - ln x_w
/// ```
///
/// where the ion branch's `B` is `beta0 + beta1 g(alpha1 sqrt(I))`. **`exp(-alpha sqrt(I))`
/// and `g(alpha sqrt(I))` are different functions**, not two spellings of one, so a port
/// that shared them would be wrong at every state.
///
/// `neutral_osmotic` is the neutral layer's contribution, zero when the layer is
/// inactive; it is passed in because it is the one term here that is not a function of the
/// ions.
///
/// # The two policies for a vanished solvent
///
/// NeqSim returns `gamma = 1` - the ideal value - when the water mole fraction falls below
/// `1e-10`, and that is reproduced. [`crate::electrolyte::water_activity_coefficient`]
/// *refuses* the same state, and the two differ deliberately: that helper is the shared
/// relation, where a silently-ideal answer is the failure this library is organised
/// against, and this is a port of a model that chose otherwise.
#[must_use]
pub fn ln_gamma_water(
    context: &IonActivityContext<'_>,
    parameters: &dyn PairParameters,
    solvent: usize,
    x_water: f64,
    neutral_osmotic: f64,
) -> f64 {
    let i = context.ionic_strength();
    let sqrt_i = i.sqrt();
    let a_phi = context.a_phi;
    let n = context.molality.len();

    // Charged components, plus the neutral solutes when that layer is active - and
    // **`charge != 0` exactly**, which is a different test from the `|z| >= 0.5` the ion
    // branch uses, so a component with a tiny charge counts here and not there.
    let sum_molalities: f64 = (0..n)
        .filter(|&k| {
            context.charge[k] != 0.0 || (context.neutral_interactions_active && k != solvent)
        })
        .map(|k| context.molality[k])
        .sum();

    // NeqSim's own guard, and its own answer: the ideal value rather than a refusal.
    if sum_molalities < 1.0e-10 {
        return 0.0;
    }

    // f^phi = -A_phi I^1.5 / (1 + b sqrt(I)). **The base is non-negative by construction**
    // - an ionic strength is a sum of non-negative terms - so the fractional power is
    // total here and needs no guard.
    let f_phi = -a_phi * i.powf(1.5) / (1.0 + B * sqrt_i);

    // Z = sum |z_i| m_i over **every** component, water included; its charge is zero.
    let z: f64 = context
        .charge
        .iter()
        .zip(context.molality)
        .map(|(&charge, &m)| charge.abs() * m)
        .sum();

    let mut binary = 0.0;
    for cation in 0..n {
        let cation_charge = context.charge[cation];
        if cation_charge <= 0.0 {
            continue;
        }
        for anion in 0..n {
            let anion_charge = context.charge[anion];
            if anion_charge >= 0.0 {
                continue;
            }
            let beta0 = parameters.beta0(cation, anion, context.temperature);
            let beta1 = parameters.beta1(cation, anion, context.temperature);
            let cphi = parameters.cphi(cation, anion, context.temperature);
            let alpha_one = alpha1_as_used(cation_charge, anion_charge);
            let mut b_phi = beta0 + beta1 * (-alpha_one * sqrt_i).exp();
            if cation_charge.abs() >= 1.5 && anion_charge.abs() >= 1.5 || context.non_two_two_beta2
            {
                let beta2 = parameters.beta2(cation, anion, context.temperature);
                if beta2.abs() > BETA2_FLOOR {
                    b_phi += beta2 * (-alpha2(cation_charge, anion_charge) * sqrt_i).exp();
                }
            }
            let product = (cation_charge * anion_charge).abs();
            let c = if product > 0.0 {
                cphi / (2.0 * product.sqrt())
            } else {
                0.0
            };
            binary += context.molality[cation] * context.molality[anion] * (b_phi + z * c);
        }
    }

    let mut theta_psi = 0.0;
    // Cation-cation, then cation-cation-anion.
    for first in 0..n {
        let first_charge = context.charge[first];
        if first_charge <= 0.0 {
            continue;
        }
        for second in (first + 1)..n {
            let second_charge = context.charge[second];
            if second_charge <= 0.0 {
                continue;
            }
            let m_first = context.molality[first];
            let m_second = context.molality[second];
            let theta = parameters.theta(first, second, context.temperature);
            let electrostatic = electrostatic_phi(context, first_charge, second_charge, i);
            theta_psi += m_first * m_second * (theta + electrostatic);
            for anion in 0..n {
                if context.charge[anion] >= 0.0 {
                    continue;
                }
                let m = context.molality[anion];
                theta_psi += m_first
                    * m_second
                    * m
                    * parameters.psi(first, second, anion, context.temperature);
            }
        }
    }
    // Anion-anion, then anion-anion-cation.
    for first in 0..n {
        let first_charge = context.charge[first];
        if first_charge >= 0.0 {
            continue;
        }
        for second in (first + 1)..n {
            let second_charge = context.charge[second];
            if second_charge >= 0.0 {
                continue;
            }
            let m_first = context.molality[first];
            let m_second = context.molality[second];
            let theta = parameters.theta(first, second, context.temperature);
            let electrostatic = electrostatic_phi(context, first_charge, second_charge, i);
            theta_psi += m_first * m_second * (theta + electrostatic);
            for cation in 0..n {
                if context.charge[cation] <= 0.0 {
                    continue;
                }
                let m = context.molality[cation];
                theta_psi += m_first
                    * m_second
                    * m
                    * parameters.psi(first, second, cation, context.temperature);
            }
        }
    }

    let phi = 1.0 + (2.0 / sum_molalities) * (f_phi + binary + theta_psi + neutral_osmotic);
    let ln_a_w = crate::electrolyte::ln_water_activity(phi, sum_molalities, WATER_MOLAR_MASS);

    if x_water < 1.0e-10 {
        return 0.0;
    }
    ln_a_w - x_water.ln()
}

/// `E_theta + I dE_theta/dI`, the combination the *osmotic* route uses.
///
/// The ion branch takes `E_theta` and its derivative separately - the value through the
/// theta sum, the derivative through `F` - so this sum appears here and nowhere else.
fn electrostatic_phi(
    context: &IonActivityContext<'_>,
    first_charge: f64,
    second_charge: f64,
    ionic_strength: f64,
) -> f64 {
    if !context.unequal_charge_same_sign || (first_charge - second_charge).abs() < DERIVATIVE_FLOOR
    {
        return 0.0;
    }
    let (value, derivative) = crate::pitzer_electrostatic::calculate(
        first_charge,
        second_charge,
        ionic_strength,
        context.a_phi,
    );
    value + ionic_strength * derivative
}

#[cfg(test)]
mod water_tests {
    use super::ion_tests::{CatalogParameters, molalities};
    use super::*;

    /// **The water node against the built phase**, at the two brines the ion tests use.
    ///
    /// The osmotic coefficient and `ln gamma_w` both come from the capture, and they are
    /// the check that this route is not the ion one: the ion branch's `B` is
    /// `beta0 + beta1 g(alpha sqrt(I))` and this is `beta0 + beta1 exp(-alpha sqrt(I))`, so
    /// a shared implementation would agree on neither number.
    #[test]
    fn the_water_activity_coefficient_matches_neqsim() {
        // water + Na+ + Cl-.
        let names = ["water", "Na+", "Cl-"];
        let molality = molalities(&[0.88, 0.06, 0.06]);
        let charge = [0.0, 1.0, -1.0];
        let parameters = CatalogParameters::new(&names);
        let context = IonActivityContext {
            molality: &molality,
            charge: &charge,
            temperature: 298.15,
            a_phi: debye_huckel_a_phi(298.15),
            common_ion_terms: true,
            non_two_two_beta2: false,
            unequal_charge_same_sign: false,
            neutral_interactions_active: false,
        };
        let got = ln_gamma_water(&context, &parameters, 0, 0.88, 0.0);
        assert!(
            (got - -0.022_033_938_564_997_1).abs() < 1.0e-9,
            "ln gamma_w = {got}, and NeqSim gives -0.0220339385649971"
        );

        // water + Na+ + Ca++ + Cl-, where `E_theta` reaches the osmotic route too.
        let names = ["water", "Na+", "Ca++", "Cl-"];
        let molality = molalities(&[0.88, 0.03, 0.03, 0.06]);
        let charge = [0.0, 1.0, 2.0, -1.0];
        let parameters = CatalogParameters::new(&names);
        let context = IonActivityContext {
            molality: &molality,
            charge: &charge,
            temperature: 298.15,
            a_phi: debye_huckel_a_phi(298.15),
            common_ion_terms: true,
            non_two_two_beta2: true,
            unequal_charge_same_sign: true,
            neutral_interactions_active: false,
        };
        let got = ln_gamma_water(&context, &parameters, 0, 0.88, 0.0);
        assert!(
            (got - -0.056_315_124_549_584_6).abs() < 1.0e-9,
            "ln gamma_w = {got}, and NeqSim gives -0.0563151245495846"
        );
    }

    /// **The osmotic coefficient itself**, which is the number both `gamma_w` and
    /// `getOsmoticCoefficientOfWater` report.
    ///
    /// Inverted from `ln a_w = -phi M_w sum m` so the number checked is the one the model
    /// reports. The first draft of this sliced one name list for both cases and paired
    /// `Cl-` with `Ca++`'s parameters, which is why the two are written out per case.
    #[test]
    fn the_osmotic_coefficient_matches_neqsim() {
        for (names, mole_fraction, charges, expected) in [
            (
                vec!["water", "Na+", "Cl-"],
                vec![0.88, 0.06, 0.06],
                vec![0.0, 1.0, -1.0],
                1.099_026_940_549_14,
            ),
            (
                vec!["water", "Na+", "Ca++", "Cl-"],
                vec![0.88, 0.03, 0.03, 0.06],
                vec![0.0, 1.0, 2.0, -1.0],
                1.350_422_304_436_11,
            ),
        ] {
            let molality = molalities(&mole_fraction);
            let parameters = CatalogParameters::new(&names);
            let context = IonActivityContext {
                molality: &molality,
                charge: &charges,
                temperature: 298.15,
                a_phi: debye_huckel_a_phi(298.15),
                common_ion_terms: true,
                non_two_two_beta2: charges.len() == 4,
                unequal_charge_same_sign: charges.len() == 4,
                neutral_interactions_active: false,
            };
            let ln_a_w = ln_gamma_water(&context, &parameters, 0, 0.88, 0.0) + 0.88f64.ln();
            let sum_m: f64 = molality
                .iter()
                .zip(&charges)
                .filter(|&(_, z)| *z != 0.0)
                .map(|(&m, _)| m)
                .sum();
            let phi = -ln_a_w / (WATER_MOLAR_MASS * sum_m);
            assert!(
                (phi - expected).abs() < 1.0e-9,
                "phi = {phi}, and NeqSim gives {expected}"
            );
        }
    }
}
