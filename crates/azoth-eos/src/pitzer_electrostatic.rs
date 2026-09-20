//! The `E_theta` integral of Pitzer's same-sign electrostatic mixing.
//!
//! `PitzerElectrostaticMixing`, and it is the one piece of the Pitzer kernel that is
//! **numerical rather than algebraic**: a Clenshaw recurrence over a 42-term Chebyshev
//! expansion, evaluating
//!
//! ```text
//! E_theta(z_j, z_k, I) = z_j z_k ( J(x_jk) - J(x_jj)/2 - J(x_kk)/2 ) / (4 I)
//! ```
//!
//! with `x_ab = 6 A_phi sqrt(I) z_a z_b` and `J` the integral the recurrence computes. Its
//! second entry is `dE_theta/dI`, which NeqSim computes analytically alongside by carrying
//! the recurrence's derivative in the same loop.
//!
//! # The two branches of the recurrence
//!
//! `J` is evaluated on `x in (0, 1]` through the substitution `x^0.2` and on `x > 1`
//! through `x^-0.1`, each against its own half of the coefficient table - which is what
//! makes the expansion converge over the whole range rather than only near zero.
//!
//! # The fractional exponents here are total
//!
//! `x^0.2` and `x^-0.1` are the rule's case of a fractional power on a base that is
//! **positive by construction**: `x = 6 A_phi sqrt(I) z_a z_b`, and the function returns
//! before reaching them unless `A_phi > 0`, `I > 0` and the two charges are same-sign and
//! non-zero. So `powf` is the right call and there is no domain to state - which is why
//! this note is here rather than a guard.

/// The Chebyshev coefficients, in NeqSim's order: `[0..21]` for `x <= 1`, `[21..42]` for
/// `x > 1`.
const COEFFICIENTS: [f64; 42] = [
    1.925_154_014_814_667,
    -0.060_076_477_753_119,
    -0.029_779_077_456_514,
    -0.007_299_499_690_937,
    0.000_388_260_636_404,
    0.000_636_874_599_598,
    0.000_036_583_601_823,
    -0.000_045_036_975_204,
    -0.000_004_537_895_710,
    0.000_002_937_706_971,
    0.000_000_396_566_462,
    -0.000_000_202_099_617,
    -0.000_000_025_267_769,
    0.000_000_013_522_610,
    0.000_000_001_229_405,
    -0.000_000_000_821_969,
    -0.000_000_000_050_847,
    0.000_000_000_046_333,
    0.000_000_000_001_943,
    -0.000_000_000_002_563,
    -0.000_000_000_010_991,
    0.628_023_320_520_852,
    0.462_762_985_338_493,
    0.150_044_637_187_895,
    -0.028_796_057_604_906,
    -0.036_552_745_910_311,
    -0.001_668_087_945_272,
    0.006_519_840_398_744,
    0.001_130_378_079_086,
    -0.000_887_171_310_131,
    -0.000_242_107_641_309,
    0.000_087_294_451_594,
    0.000_034_682_122_751,
    -0.000_004_583_768_938,
    -0.000_003_548_684_306,
    -0.000_000_250_453_880,
    0.000_000_216_991_779,
    0.000_000_080_779_570,
    0.000_000_004_558_555,
    -0.000_000_006_944_757,
    -0.000_000_002_849_257,
    0.000_000_000_237_816,
];

/// Below this the two charges are the same ion, and `E_theta` is zero.
const CHARGE_TOLERANCE: f64 = 1.0e-12;

/// Where the recurrence's `b` and `d` arrays sit in the workspace, and how long each is.
///
/// NeqSim's own offsets, kept rather than renumbered: the recurrence writes
/// `b[i+1]`, `b[i+2]`, `d[i+1]` and `d[i+2]` at both ends, so the arrays need two guard
/// slots each and the offsets are what keep them from colliding. The two results go at
/// `0` and `2` and the `b` array starts at `6`, which is why the workspace is 50 long.
const B_OFFSET: usize = 6;
const D_OFFSET: usize = 28;
const WORKSPACE: usize = 50;

/// `E_theta` and `dE_theta/dI` for a same-sign pair, the two entries NeqSim's `result`
/// array carries.
///
/// # Panics
/// Never: the guards are the caller's, expressed as an early `(0.0, 0.0)`. NeqSim throws
/// for a non-finite charge, an opposite-sign pair or a negative ionic strength; this
/// returns the zero pair, because a caller that reached here with those has a topology
/// error that the selection rule already refuses upstream.
#[must_use]
pub fn calculate(charge_j: f64, charge_k: f64, ionic_strength: f64, a_phi: f64) -> (f64, f64) {
    if (charge_j - charge_k).abs() < CHARGE_TOLERANCE || ionic_strength == 0.0 || a_phi == 0.0 {
        return (0.0, 0.0);
    }

    let x_constant = 6.0 * a_phi * ionic_strength.sqrt();
    let product = charge_j * charge_k;

    let mut workspace = [0.0; WORKSPACE];
    evaluate_integral(x_constant * product, &mut workspace, 0);
    evaluate_integral(x_constant * charge_j * charge_j, &mut workspace, 2);
    evaluate_integral(x_constant * charge_k * charge_k, &mut workspace, 4);

    let first =
        product * (workspace[0] - 0.5 * workspace[2] - 0.5 * workspace[4]) / (4.0 * ionic_strength);
    let second = product * (workspace[1] - 0.5 * workspace[3] - 0.5 * workspace[5])
        / (8.0 * ionic_strength * ionic_strength)
        - first / ionic_strength;
    (first, second)
}

/// One `J` evaluation, filling `workspace[offset]` and `[offset + 1]` with the integral
/// and its `I`-derivative.
fn evaluate_integral(x: f64, workspace: &mut [f64; WORKSPACE], offset: usize) {
    // The two substitutions, and the derivative each contributes. `x` is positive by
    // construction, so both powers are total - see the module doc.
    let (transformed, derivative_scale, coefficients) = if x <= 1.0 {
        let power = x.powf(0.2);
        (4.0 * power - 2.0, 0.4 * power, 0)
    } else {
        let power = x.powf(-0.1);
        ((40.0 * power - 22.0) / 9.0, -2.0 * power / 9.0, 21)
    };

    let b = B_OFFSET;
    let d = D_OFFSET;
    workspace[b + 21] = 0.0;
    workspace[d + 20] = 0.0;
    workspace[d + 21] = 0.0;
    workspace[b + 20] = COEFFICIENTS[coefficients + 20];
    workspace[b + 19] =
        transformed * COEFFICIENTS[coefficients + 20] + COEFFICIENTS[coefficients + 19];
    workspace[d + 19] = COEFFICIENTS[coefficients + 20];
    for i in (0..=18).rev() {
        workspace[b + i] = transformed * workspace[b + i + 1] - workspace[b + i + 2]
            + COEFFICIENTS[coefficients + i];
        workspace[d + i] =
            workspace[b + i + 1] + transformed * workspace[d + i + 1] - workspace[d + i + 2];
    }

    workspace[offset] = x / 4.0 - 1.0 + 0.5 * (workspace[b] - workspace[b + 2]);
    workspace[offset + 1] = x / 4.0 + derivative_scale * (workspace[d] - workspace[d + 2]);
}

#[cfg(test)]
mod tests {
    use super::*;

    /// **The oracle's own table**, from `PitzerElectrostaticMixing.calculate` at
    /// `A_phi = 0.392034451863750` (298.15 K).
    ///
    /// The third row is the case that matters for the guard: `z_j = z_k` returns exactly
    /// zero rather than the integral of a pair with itself.
    #[test]
    fn the_integral_matches_neqsim() {
        let a_phi = 0.392_034_451_863_750;
        for (zj, zk, i, want_e, want_de) in [
            (
                1.0,
                2.0,
                6.0,
                -0.060_026_090_218_998_0,
                0.004_957_269_079_589_56,
            ),
            (1.0, 2.0, 0.5, -0.200_307_694_629_454, 0.187_941_360_122_380),
            (2.0, 2.0, 6.0, 0.0, 0.0),
            (1.0, 3.0, 1.0, -0.876_198_487_983_861, 0.429_161_203_530_903),
            (
                1.0,
                2.0,
                100.0,
                -0.014_736_237_150_057_6,
                0.000_073_788_398_797_269_3,
            ),
        ] {
            let (got_e, got_de) = calculate(zj, zk, i, a_phi);
            assert!(
                (got_e - want_e).abs() < 1.0e-14,
                "E_theta({zj}, {zk}, {i}) = {got_e}, and NeqSim gives {want_e}"
            );
            assert!(
                (got_de - want_de).abs() < 1.0e-14 * want_de.abs().max(1.0),
                "dE_theta/dI({zj}, {zk}, {i}) = {got_de}, and NeqSim gives {want_de}"
            );
        }
    }

    /// **The substitution branches, checked at the seam.**
    ///
    /// `x <= 1` uses `x^0.2` and the first 21 coefficients, `x > 1` uses `x^-0.1` and the
    /// second 21. The two must agree at `x = 1`, or a state whose `x` is near one would
    /// step across a discontinuity that is not in the physics.
    #[test]
    fn the_two_substitutions_agree_at_the_seam() {
        // x = 6 A_phi sqrt(I) z_j z_k = 1 for these, with the charges chosen so the
        // product is exact.
        let a_phi = 0.392_034_451_863_750;
        let build = |product: f64| 1.0 / (6.0 * a_phi * product);
        // z_j z_k = 1 gives a very small I; z_j z_k = 4 gives one four times smaller.
        for (zj, zk) in [(1.0, 1.0), (2.0, 2.0), (1.0, 2.0)] {
            let i = build(zj * zk);
            let just_below = calculate(zj, zk, i * (1.0 - 1.0e-9), a_phi);
            let just_above = calculate(zj, zk, i * (1.0 + 1.0e-9), a_phi);
            assert!(
                (just_below.0 - just_above.0).abs() < 1.0e-6 * just_below.0.abs().max(1.0),
                "the substitution seams disagree at x = 1 for ({zj}, {zk}): {just_below:?} \
                 against {just_above:?}"
            );
        }
    }

    /// The same ion twice, and a state with no ionic strength, are the zero pair.
    #[test]
    fn the_degenerate_cases_are_zero() {
        assert_eq!(calculate(2.0, 2.0, 6.0, 0.4), (0.0, 0.0));
        assert_eq!(calculate(1.0, 2.0, 0.0, 0.4), (0.0, 0.0));
        assert_eq!(calculate(1.0, 2.0, 6.0, 0.0), (0.0, 0.0));
    }
}
