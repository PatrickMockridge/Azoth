//! Spec-driven tests for the `eos.molar_enthalpy_entropy` model.

use azoth_core::units::{kelvins, pascals};
use azoth_core::{AzothError, CalcResult};
use azoth_eos::mixture::Mixture;
use azoth_eos::molar_enthalpy_entropy::{IdealGasModel, molar_enthalpy_entropy};
use azoth_eos::{databank, model_gen, pr_molar_volume::MOLAR_GAS_CONSTANT};

const MODEL_ID: &str = "eos.molar_enthalpy_entropy";

/// The methane/n-butane pair the spec's cases use, resolved through the databank.
fn methane_butane() -> Mixture {
    databank::mixture_of(&["methane", "n-butane"], None)
        .expect("the pair resolves")
        .0
}

fn from_case(case: &azoth_core::spec::TestCase) -> (Mixture, IdealGasModel, f64) {
    let names = case
        .list("components")
        .expect("the case declares components");
    let (mixture, ideal_gas) =
        databank::mixture_of(names, None).expect("the case's components resolve");
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
        cp_e: vec![0.0; 2],
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
        (result.h.value / (MOLAR_GAS_CONSTANT * 330.0) - -0.5461443674029838).abs() < 1e-12,
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
    for (name, t, p, z) in [
        ("propane", 300.0, 1_000_000.0, 0.964_968_034_8),
        ("n-butane", 350.0, 1_000_000.0, 0.988_866_701_4),
    ] {
        let entry = databank::entry(name, None).expect("a databank entry");
        let (tc, omega) = (entry.tc, entry.omega);
        let mixture = databank::mixture_of(&[name], None).expect("resolves").0;
        let ideal_gas = IdealGasModel {
            cp_a: vec![4.0],
            cp_b: vec![1.0],
            cp_c: vec![-0.5],
            cp_d: vec![0.1],
            cp_e: vec![0.0],
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

/// `integral Cp dT` from `T_ref` to `T_ref` is zero, whatever the coefficients.
///
/// `h_ideal` is measured from NeqSim's fixed `referenceTemperature` of 273.15 K, so at
/// that temperature the ideal-gas enthalpy is exactly nothing, and the entropy is then
/// only the pressure and mixing terms. This is the check that would fail if a reference
/// were mixed into the integral rather than being its lower limit, and it pins the
/// constant: 298.15 K, a rounder number, gives a non-zero answer.
#[test]
fn the_ideal_gas_enthalpy_is_zero_at_the_reference_temperature() {
    let mixture = methane_butane();
    let ideal_gas = IdealGasModel {
        cp_a: vec![4.0, 4.0],
        cp_b: vec![1.0, 1.0],
        cp_c: vec![-0.5, -0.5],
        cp_d: vec![0.1, 0.1],
        cp_e: vec![0.0, 0.0],
    };
    let at = |t: f64| {
        molar_enthalpy_entropy(
            &mixture,
            &ideal_gas,
            kelvins(t),
            pascals(2_500_000.0),
            &[0.6, 0.4],
            0.8274482588400789,
        )
        .expect("a state")
    };

    let reference = at(273.15);
    assert_eq!(
        reference.h_ideal.value, 0.0,
        "the integral from T_ref to T_ref is zero"
    );

    let other = at(298.15);
    assert!(
        other.h_ideal.value.abs() > 1.0,
        "298.15 K is not the reference"
    );
    assert_ne!(other.h_departure.value, reference.h_departure.value);
}

#[test]
fn a_malformed_input_is_refused() {
    let mixture = methane_butane();
    let good = IdealGasModel {
        cp_a: vec![4.0, 4.0],
        cp_b: vec![1.0, 1.0],
        cp_c: vec![-0.5, -0.5],
        cp_d: vec![0.1, 0.1],
        cp_e: vec![0.0, 0.0],
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
}
