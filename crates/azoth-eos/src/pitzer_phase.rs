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
