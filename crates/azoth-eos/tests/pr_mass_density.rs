//! Spec-driven tests for `eos.pr_mass_density`.

use azoth_core::AzothError;
use azoth_core::units::{MolarVolume, cubic_meters_per_mole, kilograms_per_mole};
use azoth_eos::pr_mass_density;
use azoth_eos::spec_gen;
use azoth_test_support as common;

const CALC_ID: &str = "eos.pr_mass_density";

fn call(case: &azoth_core::spec::TestCase) -> azoth_eos::PrMassDensityResult {
    pr_mass_density(
        kilograms_per_mole(common::input(case, "M")),
        cubic_meters_per_mole(common::input(case, "v")),
    )
    .unwrap_or_else(|e| panic!("test `{}` should compute but failed: {e}", case.id))
}

#[test]
fn every_case_in_the_spec() {
    let spec = common::spec(spec_gen::specs(), CALC_ID);
    common::assert_skips_are_explained(spec);

    let mut executed = 0;
    for case in spec.all_tests() {
        if !case.is_active() {
            continue;
        }
        match case.kind {
            "worked_example" | "reference" => {
                let result = call(case);
                common::assert_close(
                    result.rho.value,
                    common::expected(case, "rho"),
                    case.tolerance,
                    &format!("{}::{} (rho)", spec.id, case.id),
                );
                common::assert_consistent(&result, &format!("{}::{}", spec.id, case.id));
                common::assert_warnings_agree_with_spec(
                    spec,
                    &result.warnings,
                    |quantity| match quantity {
                        "M" => Some(common::input(case, "M")),
                        "v" => Some(common::input(case, "v")),
                        "rho" => Some(result.rho.value),
                        _ => None,
                    },
                    &format!("{}::{}", spec.id, case.id),
                );
            }
            "property" => match case.property {
                Some("consistency_with") | Some("unit_round_trip") => {}
                other => panic!("{}::{}: unhandled property {other:?}", spec.id, case.id),
            },
            other => panic!("{}::{}: unknown test kind {other:?}", spec.id, case.id),
        }
        executed += 1;
    }
    assert!(
        executed >= 3,
        "expected several active cases, ran {executed}"
    );
}

/// Doubling `M` doubles `rho`, and doubling `v` halves it.
///
/// The two relations that pin the formula itself, using no reference value at all. A
/// point-value case cannot tell `M/v` from `M*v` or from `v/M`; these can, and they
/// are the reason the spec declares this as its `consistency_with` property rather
/// than another number.
#[test]
fn the_relation_is_linear_in_molar_mass_and_inverse_in_volume() {
    let (m, v) = (0.0440956, 0.0018317107825229842);
    let base = pr_mass_density(kilograms_per_mole(m), cubic_meters_per_mole(v))
        .unwrap()
        .rho
        .value;

    let double_m = pr_mass_density(kilograms_per_mole(2.0 * m), cubic_meters_per_mole(v))
        .unwrap()
        .rho
        .value;
    assert!(
        (double_m - 2.0 * base).abs() < 1e-12,
        "rho should double with M"
    );

    let double_v = pr_mass_density(kilograms_per_mole(m), cubic_meters_per_mole(2.0 * v))
        .unwrap()
        .rho
        .value;
    assert!(
        (double_v - base / 2.0).abs() < 1e-12,
        "rho should halve with v"
    );
}

/// The same state in different units gives the same density.
///
/// The conversion factor between `m**3/mol` and `cm**3/mol` is 1e-6, so an
/// implementation that got the direction wrong fails by a factor of a million rather
/// than marginally - which is what makes this worth asserting rather than assuming.
#[test]
fn the_same_state_in_other_units() {
    let si = pr_mass_density(kilograms_per_mole(0.0440956), cubic_meters_per_mole(0.002))
        .unwrap()
        .rho
        .value;

    let from_cm3 = pr_mass_density(
        kilograms_per_mole(0.0440956),
        MolarVolume::new::<uom::si::molar_volume::cubic_centimeter_per_mole>(2000.0),
    )
    .unwrap()
    .rho
    .value;

    // **Relative, not bit-for-bit, and the reason is `uom`'s rather than this
    // crate's**: its `cubic_centimeter_per_mole` factor is computed as `0.01 * 0.01
    // * 0.01`, which in binary is `1.0000000000000002e-6` and not `1e-6`. So 2000
    // cm**3/mol stores as `0.0020000000000000005` - one ulp above `0.002` - and a
    // bit-exact assertion fails on the *input*, not on the arithmetic. Measured
    // rather than guessed; the first version of this test asserted bit-equality and
    // that is what it tripped over.
    //
    // The tolerance is still far tighter than the mistake the test is for: a
    // conversion applied in the wrong direction is a factor of a million, and a
    // molar mass in g/mol is a factor of a thousand.
    let relative = (from_cm3 - si).abs() / si;
    assert!(
        relative < 1e-15,
        "expected agreement to 1e-15, got {relative:e}"
    );
    assert!(
        relative < 1e-6 / 1e9,
        "and certainly not a conversion error"
    );

    // The mistake the input's description warns about, asserted so the warning is
    // attached to something: a molar mass in g/mol gives a density 1000x out.
    let in_grams = pr_mass_density(kilograms_per_mole(44.0956), cubic_meters_per_mole(0.002))
        .unwrap()
        .rho
        .value;
    assert!((in_grams - 1000.0 * si).abs() < 1e-6);
}

#[test]
fn a_non_positive_input_is_an_error() {
    for (m, v) in [(0.0, 1e-3), (-0.044, 1e-3), (0.044, 0.0), (0.044, -1e-3)] {
        let err = pr_mass_density(kilograms_per_mole(m), cubic_meters_per_mole(v)).unwrap_err();
        assert!(matches!(err, AzothError::OutOfRange { .. }), "M={m}, v={v}");
    }
}
