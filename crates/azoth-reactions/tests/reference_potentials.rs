//! Spec-driven tests for `reactions.reference_potentials`.
//!
//! The oracle is `validation/neqsim/ReactionProbe.java`'s sibling,
//! `ReferencePotentialProbe`, whose capture is `captures/reference_potential_probe.tsv`.
//! Its three fluids are the reason the case file has two of them: the first exercises a
//! three-reaction basis and the second a five-reaction one, so a port that stopped at the
//! first fluid's reaction count is caught by the second.
//!
//! **The basis is checked as well as the numbers**, because a wrong rank rule would
//! still produce plausible potentials - it would solve over a different set of components
//! and propagate the rest. `independent` is that check.

use azoth_core::spec::TestCase;
use azoth_core::units::kelvins;
use azoth_reactions::databank::ReactionDataSource;
use azoth_reactions::model_gen;
use azoth_reactions::reference_potentials::reference_potentials;
use azoth_test_support as common;

const MODEL_ID: &str = "reactions.reference_potentials";

fn call(case: &TestCase) -> azoth_reactions::ReferencePotentialsResult {
    let names: Vec<String> = case
        .list("components")
        .expect("the case declares components")
        .iter()
        .map(|name| (*name).to_string())
        .collect();
    let source: ReactionDataSource = common::input_str(case, "source")
        .parse()
        .unwrap_or_else(|e| panic!("test `{}` names a source: {e}", case.id));
    reference_potentials(&names, source, kelvins(common::input(case, "T")))
        .unwrap_or_else(|e| panic!("test `{}` should compute but failed: {e}", case.id))
}

#[test]
fn every_case_in_the_spec() {
    let spec = model_gen::model(MODEL_ID).expect("the model should be in its own table");
    assert!(!spec.cases.is_empty());

    for case in spec.cases {
        let result = call(case);
        let context = format!("{}::{}", spec.id, case.id);

        let want = case
            .expected_vector("potentials")
            .expect("the case states potentials");
        assert_eq!(
            result.potentials.len(),
            want.len(),
            "{context}: one potential per component"
        );
        for (i, expected) in want.iter().enumerate() {
            common::assert_close(
                result.potentials[i],
                *expected,
                case.tolerance,
                &format!("{context} potentials[{i}]"),
            );
        }

        // **The basis**, which the numbers alone do not pin: a port that solved over a
        // different independent set and propagated the rest lands on the same answer for
        // a consistent system, and differs only where it does not.
        for name in ["independent", "survivors"] {
            let want = case
                .expected_vector(name)
                .unwrap_or_else(|| panic!("{context}: the case states {name}"));
            let got = if name == "independent" {
                &result.independent
            } else {
                &result.survivors
            };
            assert_eq!(got.len(), want.len(), "{context}: the {name} mask's length");
            for (i, expected) in want.iter().enumerate() {
                assert_eq!(got[i], *expected, "{context}: {name}[{i}]");
            }
        }

        assert_eq!(
            result.rank as f64,
            common::expected(case, "rank"),
            "{context}: the basis rank"
        );
    }
}

/// The capture read directly, keyed by component name rather than by position.
///
/// The case file asserts positions, because a case has to name the order; this asserts
/// **which component got which potential**, which is the claim that survives NeqSim's own
/// component order being a `HashSet`'s. The two are different tests: a port that permuted
/// its output would pass one and fail the other, and the capture is what says which.
#[test]
fn the_potentials_are_what_the_capture_attributes_to_each_component() {
    let names: Vec<String> = ["CO2", "water", "OH-", "H3O+", "HCO3-", "CO3--"]
        .iter()
        .map(|name| (*name).to_string())
        .collect();
    let result = reference_potentials(&names, ReactionDataSource::Standard, kelvins(298.15))
        .expect("the CO2-water fluid resolves");

    // From `captures/reference_potential_probe.tsv`, `fluid=co2-water`.
    let expected = [
        ("CO2", -135_846.032_628_876_45),
        ("water", 86_061.160_036_685_31),
        ("OH-", 72_306.228_492_481_24),
        ("H3O+", 199_632.183_161_778_75),
        ("HCO3-", -127_079.608_272_790_43),
        ("CO3--", -154_589.471_361_198_58),
    ];
    for (i, (name, want)) in expected.iter().enumerate() {
        assert_eq!(&names[i], name, "the case's component order");
        common::assert_close(
            result.potentials[i],
            *want,
            1e-12,
            &format!("the capture's potential for {name}"),
        );
    }

    // The capture's `independent_columns=0 1 2` and `dependent_columns=3 4 5`.
    assert_eq!(result.independent, [1.0, 1.0, 1.0, 0.0, 0.0, 0.0]);
    assert_eq!(result.rank, 3);
}

/// The three loaded reactions all survive for this fluid, and the four the fluid cannot
/// run are not candidates at all - so the mask is over the table's loaded rows, not over
/// the reactions that happened to survive.
#[test]
fn the_survivor_mask_is_over_the_tables_loaded_rows() {
    let names: Vec<String> = ["CO2", "water", "OH-", "H3O+", "HCO3-", "CO3--"]
        .iter()
        .map(|name| (*name).to_string())
        .collect();
    let result = reference_potentials(&names, ReactionDataSource::Standard, kelvins(298.15))
        .expect("the CO2-water fluid resolves");

    // Seven rows of `REACTIONDATA.csv` carry `USEREACTION=1`: CO2water, waterreac,
    // carbonate, water-H2S, water-HS, MDEAprot and DEAprot. The first three are the ones
    // whose reactants this fluid has.
    assert_eq!(
        result.survivors,
        [1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0],
        "the mask is the loaded rows in table order"
    );
}

/// A fluid carrying none of a reaction's reactants is not a candidate for it, so a
/// methane-oxygen fluid - which the standard source has no loaded reaction for - has an
/// empty basis rather than a failure.
#[test]
fn a_fluid_the_table_cannot_react_has_an_empty_basis() {
    let names: Vec<String> = ["methane", "oxygen"]
        .iter()
        .map(|name| (*name).to_string())
        .collect();
    let result = reference_potentials(&names, ReactionDataSource::Standard, kelvins(303.3))
        .expect("no reaction is not an error");
    assert_eq!(result.rank, 0);
    assert_eq!(result.potentials, vec![0.0, 0.0]);
    assert!(
        result.survivors.iter().all(|kept| *kept == 0.0),
        "the one combustion row is dormant, so nothing survives: {:?}",
        result.survivors
    );
}
