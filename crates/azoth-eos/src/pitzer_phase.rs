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
