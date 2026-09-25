//! Spec-driven tests for `eos.fuller_schettler_giddings_diffusivity`.

use azoth_core::units::{cubic_meters_per_mole, kelvins, kilograms_per_mole, pascals};
use azoth_eos::databank::{FullerVolumeSource, fuller_diffusion_volume};
use azoth_eos::spec_gen;
use azoth_eos::{FullerSchettlerGiddingsDiffusivityResult, fuller_schettler_giddings_diffusivity};
use azoth_test_support as common;

const CALC_ID: &str = "eos.fuller_schettler_giddings_diffusivity";

fn call(case: &azoth_core::spec::TestCase) -> FullerSchettlerGiddingsDiffusivityResult {
    fuller_schettler_giddings_diffusivity(
        kilograms_per_mole(common::input(case, "MA")),
        kilograms_per_mole(common::input(case, "MB")),
        cubic_meters_per_mole(common::input(case, "VA")),
        cubic_meters_per_mole(common::input(case, "VB")),
        kelvins(common::input(case, "T")),
        pascals(common::input(case, "P")),
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
                    result.d.value,
                    common::expected(case, "d"),
                    case.tolerance,
                    &format!("{}::{} (d)", spec.id, case.id),
                );
                common::assert_consistent(&result, &format!("{}::{}", spec.id, case.id));
                common::assert_warnings_agree_with_spec(
                    spec,
                    &result.warnings,
                    |quantity| match quantity {
                        "MA" => Some(common::input(case, "MA")),
                        "MB" => Some(common::input(case, "MB")),
                        "VA" => Some(common::input(case, "VA")),
                        "VB" => Some(common::input(case, "VB")),
                        "T" => Some(common::input(case, "T")),
                        "P" => Some(common::input(case, "P")),
                        "d" => Some(result.d.value),
                        _ => None,
                    },
                    &format!("{}::{}", spec.id, case.id),
                );
            }
            "property" => match case.property {
                Some("unit_round_trip") => {}
                other => panic!("{}::{}: unhandled property {other:?}", spec.id, case.id),
            },
            other => panic!("{}::{}: unknown test kind {other:?}", spec.id, case.id),
        }
        executed += 1;
    }
    assert!(
        executed >= 2,
        "expected several active cases, ran {executed}"
    );
}

/// **The diffusion-volume ladder, against the class's own three rungs.**
///
/// The special table answers first (methane's 25.14, nitrogen's 18.5 - the numbers the
/// capture's rows are reproducible from), `0.285*Vc` second (ammonia, whose table key is the
/// formula `NH3` and whose name therefore misses it), and the molar-mass estimate last, for a
/// component with neither.
#[test]
fn the_diffusion_volume_ladder_is_the_classs_own() {
    // The first rung: the table, matched case-insensitively, so `MEG` and `TEG` are found.
    assert_eq!(
        fuller_diffusion_volume("methane", Some(9.9e-5), 0.016043),
        (25.14, FullerVolumeSource::SpecialTable)
    );
    assert_eq!(
        fuller_diffusion_volume("nitrogen", Some(8.98e-5), 0.0280135),
        (18.5, FullerVolumeSource::SpecialTable)
    );
    assert_eq!(
        fuller_diffusion_volume("TEG", None, 0.15017),
        (129.18, FullerVolumeSource::SpecialTable)
    );
    // The second rung: `0.285*Vc`, for a name the table does not carry.
    let (ammonia, source) = fuller_diffusion_volume("ammonia", Some(9.9e-5), 0.017031);
    assert_eq!(source, FullerVolumeSource::CriticalVolume);
    assert!((ammonia - 0.285 * 99.0).abs() < 1e-12);
    // The third: a molar-mass estimate, which the capture's probe rows never reach because
    // every component it names carries a critical volume. `max(10, 0.95*M)` with M in g/mol.
    assert_eq!(
        fuller_diffusion_volume("unlisted", None, 0.016043),
        (15.24085, FullerVolumeSource::MolarMass)
    );
    assert_eq!(
        fuller_diffusion_volume("unlisted", Some(0.0), 0.005),
        (10.0, FullerVolumeSource::MolarMass)
    );
}
