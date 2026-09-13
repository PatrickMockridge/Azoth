//! Spec-driven tests for the `eos.molar_enthalpy_entropy` model.

use azoth_core::units::{kelvins, pascals};
use azoth_core::{AzothError, CalcResult};
use azoth_eos::mixture::{Component, Mixture};
use azoth_eos::molar_enthalpy_entropy::{IdealGasModel, molar_enthalpy_entropy};
use azoth_eos::{model_gen, pr_molar_volume::MOLAR_GAS_CONSTANT};

const MODEL_ID: &str = "eos.molar_enthalpy_entropy";

fn mixture_of(tc: &[f64], pc: &[f64], omega: &[f64], kij: Vec<f64>) -> Mixture {
    let components = (0..tc.len())
        .map(|i| Component::new(kelvins(tc[i]), pascals(pc[i]), omega[i]).expect("valid"))
        .collect();
    Mixture::new(components, kij).expect("a valid mixture")
}

fn methane_butane() -> Mixture {
    mixture_of(
        &[190.56, 425.12],
        &[4_599_200.0, 3_796_000.0],
        &[0.01142, 0.2002],
        vec![0.0, 0.05, 0.05, 0.0],
    )
}

fn from_case(case: &azoth_core::spec::TestCase) -> (Mixture, IdealGasModel, f64) {
    let mixture = mixture_of(
        case.vector("Tc").expect("Tc"),
        case.vector("Pc").expect("Pc"),
        case.vector("omega").expect("omega"),
        case.matrix("kij").expect("kij").to_vec(),
    );
    let ideal_gas = IdealGasModel {
        cp_a: case.vector("cp_a").expect("cp_a").to_vec(),
        cp_b: case.vector("cp_b").expect("cp_b").to_vec(),
        cp_c: case.vector("cp_c").expect("cp_c").to_vec(),
        cp_d: case.vector("cp_d").expect("cp_d").to_vec(),
        h_ref: case.vector("h_ref").expect("h_ref").to_vec(),
        s_ref: case.vector("s_ref").expect("s_ref").to_vec(),
        t_ref: kelvins(common_input(case, "T_ref")),
        p_ref: pascals(common_input(case, "P_ref")),
    };
    (mixture, ideal_gas, common_input(case, "compressibility"))
}

fn common_input(case: &azoth_core::spec::TestCase, name: &str) -> f64 {
    case.input(name)
        .unwrap_or_else(|| panic!("test `{}` should supply `{name}`", case.id))
}

#[test]
fn every_case_in_the_spec() {
    let spec = model_gen::model(MODEL_ID).expect("the model should be in its own table");
    assert!(!spec.cases.is_empty(), "the model should have cases");

    for case in spec.cases {
        let (mixture, ideal_gas, compressibility) = from_case(case);
        let result = molar_enthalpy_entropy(
            &mixture,
            &ideal_gas,
            kelvins(common_input(case, "T")),
            pascals(common_input(case, "P")),
            case.vector("z").expect("z"),
            compressibility,
        )
        .unwrap_or_else(|e| panic!("case `{}` should compute but failed: {e}", case.id));

        // Absolute rather than relative: these are enthalpies on a zero datum, so a
        // relative comparison would divide by a number that is near a sign change.
        for (name, actual) in [
            ("h", result.h.value),
            ("s", result.s.value),
            ("h_ideal", result.h_ideal.value),
            ("s_ideal", result.s_ideal.value),
            ("h_departure", result.h_departure.value),
            ("s_departure", result.s_departure.value),
        ] {
            let expected = case
                .expected_value(name)
                .unwrap_or_else(|| panic!("case `{}` should declare `{name}`", case.id));
            assert!(
                (actual - expected).abs() <= 1e-6,
                "{}::{} ({name}): got {actual}, expected {expected}",
                spec.id,
                case.id
            );
        }
        assert!(
            (result.psi_bar - case.expected_value("psi_bar").expect("psi_bar")).abs() < 1e-12,
            "{}::{}: psi_bar",
            spec.id,
            case.id
        );
        assert!(
            result.is_clean(),
            "unexpected warnings: {:?}",
            result.warnings
        );
    }
}

/// With the ideal-gas terms off, the enthalpy **is** the departure.
///
/// The identity that ties this model to `eos.pr_departure`: `h = R*T*h_dep_rt`. Exact
/// rather than approximate, because both sides are the same product - so a defect in
/// the departure arithmetic or an extra term in the enthalpy shows up as a disagreement
/// rather than as a tolerance.
#[test]
fn the_enthalpy_is_the_departure_when_the_ideal_gas_terms_are_off() {
    let mixture = methane_butane();
    let ideal_gas = IdealGasModel {
        cp_a: vec![0.0; 2],
        cp_b: vec![0.0; 2],
        cp_c: vec![0.0; 2],
        cp_d: vec![0.0; 2],
        h_ref: vec![0.0; 2],
        s_ref: vec![0.0; 2],
        t_ref: kelvins(298.15),
        p_ref: pascals(101_325.0),
    };
    let result = molar_enthalpy_entropy(
        &mixture,
        &ideal_gas,
        kelvins(330.0),
        pascals(2_500_000.0),
        &[0.6, 0.4],
        0.8274482588400789,
    )
    .expect("a state");

    assert_eq!(
        result.h_ideal.value, 0.0,
        "the ideal-gas enthalpy should be zero"
    );
    assert_eq!(result.h.value, result.h_departure.value);
    assert!(
        (result.h.value / (MOLAR_GAS_CONSTANT * 330.0) - -0.5399247170906822).abs() < 1e-12,
        "h/(R T) is {}, not the departure the mixture module computes",
        result.h.value / (MOLAR_GAS_CONSTANT * 330.0)
    );
    // And the entropy is the departure *plus* the two ideal-gas terms that no
    // coefficient switches off - which is the pair a reader might mistake for a defect.
    let pressure_and_mixing = -MOLAR_GAS_CONSTANT * (2_500_000.0f64 / 101_325.0).ln()
        - MOLAR_GAS_CONSTANT * (0.6 * 0.6f64.ln() + 0.4 * 0.4f64.ln());
    assert!(
        (result.s_ideal.value - pressure_and_mixing).abs() < 1e-9,
        "s_ideal is {}, but the pressure and mixing terms are {pressure_and_mixing}",
        result.s_ideal.value
    );
    assert!((result.s.value - (pressure_and_mixing + result.s_departure.value)).abs() < 1e-9);
}

/// At one component the departures are `eos.pr_departure`'s, exactly.
///
/// The reduction that holds the mixture form, run through this model rather than
/// through the mixture module directly: `h_departure / (R*T)` must equal the registered
/// calc's `h_dep_rt` bit for bit, because the mixture form is the pure form with `psi`
/// replaced by a weighted average that at `N = 1` is that component's own `psi`.
#[test]
fn the_departures_reduce_to_pr_departure_at_one_component() {
    for (tc, pc, omega, t, p, z) in [
        (
            369.83,
            4_248_000.0,
            0.1523,
            300.0,
            1_000_000.0,
            0.964_968_034_8,
        ),
        (
            425.12,
            3_796_000.0,
            0.2002,
            350.0,
            1_000_000.0,
            0.988_866_701_4,
        ),
    ] {
        let mixture = mixture_of(&[tc], &[pc], &[omega], vec![0.0]);
        let ideal_gas = IdealGasModel {
            cp_a: vec![4.0],
            cp_b: vec![1.0],
            cp_c: vec![-0.5],
            cp_d: vec![0.1],
            h_ref: vec![0.0],
            s_ref: vec![0.0],
            t_ref: kelvins(298.15),
            p_ref: pascals(101_325.0),
        };
        let reduced = mixture
            .reduced_parameters(kelvins(t), pascals(p))
            .expect("a state");
        let roots = azoth_eos::pr_z_factor(
            mixture.mixture_parameters(&reduced, &[1.0]).0,
            mixture.mixture_parameters(&reduced, &[1.0]).1,
        )
        .expect("roots");
        let compressibility = roots.z_max;

        let result = molar_enthalpy_entropy(
            &mixture,
            &ideal_gas,
            kelvins(t),
            pascals(p),
            &[1.0],
            compressibility,
        )
        .expect("a state");

        let kappa = azoth_eos::pr_kappa(omega).expect("kappa").kappa;
        let pure =
            azoth_eos::pr_departure(reduced.a[0], reduced.b[0], compressibility, kappa, t / tc)
                .expect("a departure");

        assert_eq!(
            result.psi_bar, reduced.psi[0],
            "Tc={tc}: psi_bar should be psi"
        );
        assert_eq!(
            result.h_departure.value / (MOLAR_GAS_CONSTANT * t),
            pure.h_dep_rt,
            "Tc={tc}, T={t}: the departure enthalpy should be bit-identical"
        );
        assert!(
            (result.s_departure.value / MOLAR_GAS_CONSTANT - pure.s_dep_r).abs() < 1e-14,
            "Tc={tc}, T={t}: the entropy departures differ by more than one ulp"
        );
        let _ = z;
    }
}

/// The datum shifts the enthalpy, and only by what the caller put in.
///
/// A model whose absolute enthalpy is meaningful only relative to a datum must move
/// with it exactly: raising every `h_ref` by a constant raises `h` by that constant,
/// leaving both the ideal-gas integrals and the departures untouched. It is what would
/// fail if the reference were mixed into the integral rather than added to it, which is
/// invisible at a zero datum - the only datum a spec case can state.
#[test]
fn the_datum_shifts_the_enthalpy_exactly() {
    let mixture = methane_butane();
    let build = |offset: f64| IdealGasModel {
        cp_a: vec![4.0, 4.0],
        cp_b: vec![1.0, 1.0],
        cp_c: vec![-0.5, -0.5],
        cp_d: vec![0.1, 0.1],
        h_ref: vec![offset, offset],
        s_ref: vec![0.0, 0.0],
        t_ref: kelvins(298.15),
        p_ref: pascals(101_325.0),
    };
    let at = |model: &IdealGasModel| {
        molar_enthalpy_entropy(
            &mixture,
            model,
            kelvins(330.0),
            pascals(2_500_000.0),
            &[0.6, 0.4],
            0.8274482588400789,
        )
        .expect("a state")
    };

    let plain = at(&build(0.0));
    let shifted = at(&build(-285_830.0));
    assert!(
        (shifted.h.value - (plain.h.value - 285_830.0)).abs() < 1e-9,
        "shifting the datum by -285830 moved h by {}, not by that",
        shifted.h.value - plain.h.value
    );
    // Relative rather than exact: adding 285830 back to a value near -284700 is a
    // cancellation, and the relative error of the recovery is the arithmetic's rather
    // than the model's. `s` and `h_departure` below need no such allowance, because
    // the datum does not enter either of them.
    assert!(
        (shifted.h_ideal.value + 285_830.0 - plain.h_ideal.value).abs() < 1e-9,
        "recovering the ideal-gas enthalpy from the shifted one lost more than the \
         cancellation explains"
    );
    assert_eq!(shifted.h_departure.value, plain.h_departure.value);
    assert_eq!(
        shifted.s.value, plain.s.value,
        "the datum is an enthalpy only"
    );
}

#[test]
fn a_malformed_input_is_refused() {
    let mixture = methane_butane();
    let good = IdealGasModel {
        cp_a: vec![4.0, 4.0],
        cp_b: vec![1.0, 1.0],
        cp_c: vec![-0.5, -0.5],
        cp_d: vec![0.1, 0.1],
        h_ref: vec![0.0, 0.0],
        s_ref: vec![0.0, 0.0],
        t_ref: kelvins(298.15),
        p_ref: pascals(101_325.0),
    };
    let call = |model: &IdealGasModel, z: &[f64], compressibility: f64| {
        molar_enthalpy_entropy(
            &mixture,
            model,
            kelvins(330.0),
            pascals(2_500_000.0),
            z,
            compressibility,
        )
    };

    // A coefficient vector of the wrong length.
    let mut short = good.clone();
    short.cp_b = vec![1.0];
    assert!(matches!(
        call(&short, &[0.6, 0.4], 0.827),
        Err(AzothError::InvalidInput { .. })
    ));

    // A composition that is not one.
    for z in [vec![0.6, 0.5], vec![0.6, -0.6], vec![0.6, 0.4, 0.0]] {
        assert!(call(&good, &z, 0.827).is_err(), "z = {z:?}");
    }

    // A root at or below the mixture's B, where ln(z - B) is undefined.
    let err = call(&good, &[0.6, 0.4], 0.01).unwrap_err();
    assert!(matches!(err, AzothError::OutOfRange { .. }));
    assert_eq!(err.field(), Some("z"));

    // A reference state that is not a state.
    let mut bad_ref = good.clone();
    bad_ref.t_ref = kelvins(0.0);
    assert!(matches!(
        call(&bad_ref, &[0.6, 0.4], 0.827),
        Err(AzothError::OutOfRange { .. })
    ));
}
