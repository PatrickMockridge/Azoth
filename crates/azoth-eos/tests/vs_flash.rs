//! Spec-driven tests for the `eos.vs_flash` model.

use azoth_core::units::{cubic_meters_per_mole, joules_per_mole_kelvin};
use azoth_eos::Cubic;
use azoth_eos::databank;
use azoth_eos::vs_flash::vs_flash;
use azoth_test_support as common;

const MODEL_ID: &str = "eos.vs_flash";

#[test]
fn every_case_in_the_spec() {
    let spec = azoth_eos::model_gen::model(MODEL_ID).expect("the model");
    assert!(!spec.cases.is_empty(), "the model should have cases");
    for case in spec.cases {
        let (mixture, ideal_gas) = databank::mixture_of(
            case.list("components").expect("components"),
            Cubic::Pr,
            None,
        )
        .expect("the case's fluid resolves");
        let context = &format!("{}::{}", spec.id, case.id);
        let result = vs_flash(
            &mixture,
            &ideal_gas,
            cubic_meters_per_mole(common::input(case, "V")),
            joules_per_mole_kelvin(common::input(case, "S")),
            case.vector("z").expect("z"),
        )
        .unwrap_or_else(|e| panic!("{context} should compute but failed: {e}"));
        common::assert_close(
            result.pressure.value,
            common::expected(case, "P"),
            case.tolerance,
            &format!("{context} (P)"),
        );
        common::assert_close(
            result.temperature.value,
            common::expected(case, "T"),
            case.tolerance,
            &format!("{context} (T)"),
        );
        common::assert_consistent(&result, context);
    }
}

/// The answer is the state whose volume and entropy are the ones asked for.
///
/// Two more states than the spec carries, and the check is the specification itself
/// rather than a round trip: the volume at the answer must be the volume asked for, and
/// the entropy must be too. A model that satisfied one and not the other would pass
/// every round-trip case whose two happened to be consistent.
#[test]
fn both_specifications_hold_at_the_answer() {
    use azoth_eos::flash_property::{Property, property_at};
    let (mixture, ideal_gas) =
        databank::mixture_of(&["methane", "propane", "n-butane"], Cubic::Pr, None)
            .expect("the ternary");
    let z = [0.5, 0.3, 0.2];
    for (t, pbar) in [(380.0, 1.5e6), (320.0, 4.0e6), (430.0, 8.0e6)] {
        let (v, _) = property_at(
            &mixture,
            &ideal_gas,
            azoth_core::units::kelvins(t),
            azoth_core::units::pascals(pbar),
            &z,
            Property::Volume,
        )
        .expect("a volume");
        let (s, _) = property_at(
            &mixture,
            &ideal_gas,
            azoth_core::units::kelvins(t),
            azoth_core::units::pascals(pbar),
            &z,
            Property::Entropy,
        )
        .expect("an entropy");
        let result = vs_flash(
            &mixture,
            &ideal_gas,
            cubic_meters_per_mole(v),
            joules_per_mole_kelvin(s),
            &z,
        )
        .expect("the state has an answer");
        let (found_v, _) = property_at(
            &mixture,
            &ideal_gas,
            result.temperature,
            result.pressure,
            &z,
            Property::Volume,
        )
        .expect("a volume");
        let (found_s, _) = property_at(
            &mixture,
            &ideal_gas,
            result.temperature,
            result.pressure,
            &z,
            Property::Entropy,
        )
        .expect("an entropy");
        // Both are looser than the arithmetic's last bit because the iteration *accepts*
        // 1e-3 relative - NeqSim's own acceptance, recorded in the spec - and stops
        // there rather than at the answer's last digit. Measured at ~5e-6 relative, so
        // the tolerance is that with margin and not the model's own bound: a test that
        // asserted 1e-3 would pass on a state that had barely moved.
        assert!(
            (found_v - v).abs() <= 1.0e-5 * v.abs(),
            "T={t}, P={pbar}: the volume at the answer is {found_v} against the {v} asked for"
        );
        // The entropy is the looser of the two specifications, and by more than the
        // volume: the iteration's temperature derivative is the flashed heat capacity,
        // which is the frozen-composition slope rather than the equilibrium one, so the
        // entropy converges slowly. Measured at ~2e-6 relative, and the model *accepts*
        // 1e-3 - see the spec's assumption, which records both.
        assert!(
            (found_s - s).abs() <= 1.0e-5 * s.abs().max(1.0),
            "T={t}, P={pbar}: the entropy at the answer is {found_s} against the {s} asked for"
        );
    }
}
