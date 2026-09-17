//! Spec-driven tests for the `eos.ge_nrtl_flash` model.

use azoth_core::AzothError;
use azoth_core::units::{kelvins, pascals};
use azoth_eos::Phase;
use azoth_eos::cubic::Cubic;
use azoth_eos::databank::{ge_nrtl_phase_parameters, mixture_of};
use azoth_eos::ge_nrtl_flash::ge_nrtl_flash;
use azoth_test_support as common;

const MODEL_ID: &str = "eos.ge_nrtl_flash";

/// The cubic a case declares, from the spec's optional `eos` input.
///
/// Read from the case rather than fixed here, because the cubic is half of what a
/// gamma-phi flash *is*: a case that says SRK is NeqSim's `SystemNRTL` pairing, and a
/// case that says nothing is `databank::mixture_of`'s default.
fn cubic_of(case: &azoth_core::spec::TestCase) -> Cubic {
    match case.string("eos").unwrap_or("pr") {
        "srk" => Cubic::Srk,
        "rk" => Cubic::Rk,
        "tst" => Cubic::Tst,
        _ => Cubic::Pr,
    }
}

/// The mixture the hand-written tests use: NeqSim's `SystemNRTL` pairs the NRTL liquid
/// with `PhaseSrkEos`, so the cubic is SRK and not `databank::mixture_of`'s default.
fn srk(names: &[&str]) -> azoth_eos::mixture::Mixture {
    let (mixture, _) = mixture_of(names, None).expect("the components resolve");
    mixture.with_cubic(Cubic::Srk)
}

fn call(case: &azoth_core::spec::TestCase) -> azoth_eos::GeNrtlFlashResult {
    let names = case.list("components").expect("components");
    let params = ge_nrtl_phase_parameters(names, None)
        .unwrap_or_else(|e| panic!("case `{}` should resolve but failed: {e}", case.id));
    let (mixture, _) = mixture_of(names, None).expect("the components resolve");
    let mixture = mixture.with_cubic(cubic_of(case));
    ge_nrtl_flash(
        &params,
        &mixture,
        kelvins(common::input(case, "T")),
        pascals(common::input(case, "P")),
        case.vector("z").expect("z"),
    )
    .unwrap_or_else(|e| panic!("case `{}` should compute but failed: {e}", case.id))
}

#[test]
fn every_case_in_the_spec() {
    let spec = azoth_eos::model_gen::model(MODEL_ID).expect("the model should be in its own table");
    assert!(!spec.cases.is_empty(), "the model should have cases");

    for case in spec.cases {
        let result = call(case);
        let context = &format!("{}::{}", spec.id, case.id);
        for name in ["x", "y", "k", "ln_phi_liquid", "ln_phi_vapour"] {
            let actual = match name {
                "x" => &result.x,
                "y" => &result.y,
                "k" => &result.k,
                "ln_phi_liquid" => &result.ln_phi_liquid,
                _ => &result.ln_phi_vapour,
            };
            let expected = case
                .expected_vector(name)
                .unwrap_or_else(|| panic!("the case declares {name}"));
            for (i, (&got, &want)) in actual.iter().zip(expected).enumerate() {
                common::assert_close(
                    got,
                    want,
                    case.tolerance,
                    &format!("{context} ({name}[{i}])"),
                );
            }
        }
        for name in ["beta", "z_vapour"] {
            let expected = case
                .expected_value(name)
                .unwrap_or_else(|| panic!("the case declares {name}"));
            let got = if name == "beta" {
                result.beta.expect("a two-phase case has a vapour fraction")
            } else {
                result.z_vapour
            };
            common::assert_close(
                got,
                expected,
                case.tolerance,
                &format!("{context} ({name})"),
            );
        }
        common::assert_consistent(&result, context);
    }
}

/// The three flash equations hold at the state the result reports.
///
/// This is what makes the result a flash rather than a plausible-looking pair of
/// vectors: `x` and `y` are *defined* by `beta` and `K`, and `K` is *defined* by the
/// two fugacity coefficients. A model that returned a converged-looking `beta` with
/// compositions from a different iterate would fail here and nowhere else.
#[test]
fn the_flash_equations_hold_at_the_reported_state() {
    let names = ["methanol", "water"];
    let params = ge_nrtl_phase_parameters(&names, None).unwrap();
    let mixture = srk(&names);
    let (t, p) = (350.0, 1.0e5);
    let z = [0.5, 0.5];
    let r = ge_nrtl_flash(&params, &mixture, kelvins(t), pascals(p), &z).unwrap();
    let beta = r.beta.expect("a two-phase case has a vapour fraction");

    // `x_i = z_i / (1 + beta (K_i - 1))` and `y_i = K_i x_i`, and the material balance
    // those two imply.
    for (i, &zi) in z.iter().enumerate() {
        let x = zi / (1.0 + beta * (r.k[i] - 1.0));
        assert!(
            (r.x[i] - x).abs() < 1e-12,
            "component {i}: x is {} but the Rachford-Rice composition is {x}",
            r.x[i]
        );
        assert!(
            (r.y[i] - r.k[i] * r.x[i]).abs() < 1e-15,
            "component {i}: y is {} but K x is {}",
            r.y[i],
            r.k[i] * r.x[i]
        );
        let balance = (1.0 - beta) * r.x[i] + beta * r.y[i];
        assert!(
            (balance - z[i]).abs() < 1e-12,
            "component {i}: the material balance gives {balance}, the feed is {zi}"
        );
        // `K_i = phi_i^L / phi_i^V` at the reported compositions.
        let ln_k = r.ln_phi_liquid[i] - r.ln_phi_vapour[i];
        assert!(
            (ln_k - r.k[i].ln()).abs() < 1e-9,
            "component {i}: ln K is {} but ln(phi_L / phi_V) is {ln_k}",
            r.k[i].ln()
        );
    }
}

/// The liquid half is `eos.ge_nrtl_phase`'s, not a second copy of it.
///
/// The phase model is the one that carries the NeqSim-checked NRTL matrices and
/// Antoine correlations, so a flash that built its own liquid could drift from it.
#[test]
fn the_liquid_half_is_the_nrtl_phase_model() {
    let names = ["methanol", "water"];
    let params = ge_nrtl_phase_parameters(&names, None).unwrap();
    let mixture = srk(&names);
    let r = ge_nrtl_flash(
        &params,
        &mixture,
        kelvins(350.0),
        pascals(1.0e5),
        &[0.5, 0.5],
    )
    .unwrap();

    let phase = azoth_eos::ge_nrtl_phase::ge_nrtl_phase(&params, 350.0, 1.0e5, &r.x).unwrap();
    assert_eq!(r.ln_phi_liquid, phase.ln_phi);
}

/// A subcooled feed is single-phase liquid, and its `beta` is absent rather than zero.
///
/// At 298.15 K every K-value is below one, so the Rachford-Rice equation has no root
/// at all - the feed is single phase *by proof*, and there is no vapour fraction to
/// report. This is the route a `beta = 0.0` would misdescribe: zero would say the
/// flash found a saturated liquid, which it did not.
#[test]
fn a_subcooled_feed_has_no_vapour_fraction_at_all() {
    let names = ["methanol", "water"];
    let params = ge_nrtl_phase_parameters(&names, None).unwrap();
    let mixture = srk(&names);
    let r = ge_nrtl_flash(
        &params,
        &mixture,
        kelvins(298.15),
        pascals(1.0e5),
        &[0.5, 0.5],
    )
    .unwrap();
    assert_eq!(r.phase, Phase::AllLiquid);
    assert_eq!(r.beta, None, "no root exists, so `beta` is absent");
    assert_eq!(r.x, vec![0.5, 0.5], "a single-phase feed is the feed");
    assert_eq!(
        r.iterations, 1,
        "the bracket is refused before any evaluation"
    );
}

/// The other single-phase route: a root that exists and is negative.
///
/// At 340 K the K-values straddle one, so Rachford-Rice has a root and the flash
/// converges to it - the negative flash, the amount of vapour that would have to be
/// added to bring the feed to saturation. `beta` is present and the result carries
/// `OUT_OF_VALID_RANGE`, which is a different statement from the case above.
#[test]
fn a_negative_flash_reports_its_vapour_fraction() {
    let names = ["methanol", "water"];
    let params = ge_nrtl_phase_parameters(&names, None).unwrap();
    let mixture = srk(&names);
    let r = ge_nrtl_flash(
        &params,
        &mixture,
        kelvins(340.0),
        pascals(1.0e5),
        &[0.5, 0.5],
    )
    .unwrap();
    assert_eq!(r.phase, Phase::AllLiquid);
    let beta = r
        .beta
        .expect("a negative flash reports its vapour fraction");
    assert!(
        beta < 0.0,
        "the flash is subcooled, so `beta` is negative: {beta}"
    );
    assert!(
        r.warnings
            .iter()
            .any(|w| w.code == azoth_core::WarningCode::OutOfValidRange),
        "a negative flash carries OUT_OF_VALID_RANGE: {:?}",
        r.warnings
    );
}

/// The split grows with temperature, and crosses into all-vapour above the dew point.
///
/// A flash that converged to a stationary point of the wrong sign, or that lost the
/// root, would show up here as a `beta` that is out of order or a phase that does not
/// change.
#[test]
fn the_vapour_fraction_grows_with_temperature() {
    let names = ["methanol", "water"];
    let params = ge_nrtl_phase_parameters(&names, None).unwrap();
    let mixture = srk(&names);
    let mut previous = f64::NEG_INFINITY;
    for t in [345.0_f64, 348.0, 350.0, 352.0, 355.0] {
        let r = ge_nrtl_flash(&params, &mixture, kelvins(t), pascals(1.0e5), &[0.5, 0.5]).unwrap();
        assert_eq!(r.phase, Phase::TwoPhase, "at {t} K");
        let beta = r.beta.expect("a split has a vapour fraction");
        assert!(
            beta > previous,
            "beta at {t} K is {beta}, not above {previous}"
        );
        previous = beta;
    }
    assert!(
        (0.409..0.956).contains(&previous),
        "the isotherm from 345 K to 355 K crosses the two-phase region: {previous}"
    );
}

#[test]
fn a_composition_that_does_not_sum_to_one_is_refused() {
    let names = ["methanol", "water"];
    let params = ge_nrtl_phase_parameters(&names, None).unwrap();
    let mixture = srk(&names);
    let err = ge_nrtl_flash(
        &params,
        &mixture,
        kelvins(350.0),
        pascals(1.0e5),
        &[0.6, 0.6],
    )
    .unwrap_err();
    assert!(matches!(err, AzothError::InvalidInput { .. }), "{err:?}");
    assert_eq!(err.field(), Some("z"));
}

#[test]
fn a_component_count_mismatch_is_refused() {
    let names = ["methanol", "water"];
    let params = ge_nrtl_phase_parameters(&names, None).unwrap();
    let mixture = srk(&names);
    let err = ge_nrtl_flash(&params, &mixture, kelvins(350.0), pascals(1.0e5), &[1.0]).unwrap_err();
    assert!(matches!(err, AzothError::InvalidInput { .. }), "{err:?}");
    assert_eq!(err.field(), Some("z"));
}
