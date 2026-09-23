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
use azoth_reactions::databank::{
    ReactionDataSource, formation_properties, reactions, stoichiometry,
};
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

/// The deadlock fallback is reproduced, and this is the case that exercises it.
///
/// `water/H3O+/OH-` reduces to one independent column and two dependent ones, and no
/// surviving reaction has exactly one unknown - so NeqSim seeds the first uncomputed
/// dependent component with its Gibbs energy of formation and carries on from there.
/// `water` is solved for, `H3O+` takes its Gibbs energy of formation raw and not negated,
/// and `OH-` follows from it by ordinary propagation.
#[test]
fn the_deadlock_fallback_seeds_one_component_and_propagates_from_it() {
    let names: Vec<String> = ["water", "H3O+", "OH-"]
        .iter()
        .map(|name| (*name).to_string())
        .collect();
    let result = reference_potentials(&names, ReactionDataSource::Standard, kelvins(298.15))
        .expect("the deadlock branch answers rather than refusing");

    assert_eq!(result.rank, 1);
    assert_eq!(
        result.independent,
        vec![1.0, 0.0, 0.0],
        "only the first column is solved for"
    );

    let seeded = formation_properties("H3O+")
        .expect("parses")
        .expect("the component databank carries it");
    assert_eq!(
        result.potentials[1], seeded.gibbs_energy_of_formation,
        "the seed is the component's Gibbs energy of formation, not its negative"
    );

    let propagated = formation_properties("OH-")
        .expect("parses")
        .expect("the component databank carries it");
    assert_ne!(
        result.potentials[2], propagated.gibbs_energy_of_formation,
        "OH- is propagated from the seeded pair, not seeded itself"
    );
}

/// **Every subset of the species the three tables' loaded reactions name answers**, and
/// the deadlock fallback is what makes that true: 3,320 of these 14,333 subsets seed a
/// component the propagation cannot reach, measured.
///
/// The counts are asserted rather than described, in both directions. A change that left
/// the propagation unable to finish would stop the sweep from answering at all, and a
/// change that made the reduction complete on its own would leave the fallback
/// unexercised and drop the second count to zero - so the branch cannot quietly stop
/// being reached, and a pin move that adds a reaction to any of the three tables has to
/// re-measure both numbers.
#[test]
fn every_reaction_set_the_vendored_tables_admit_answers_and_the_seed_is_exercised() {
    let mut subsets = 0usize;
    let mut refused = 0usize;
    let mut seeded = 0usize;
    for source in [
        ReactionDataSource::Standard,
        ReactionDataSource::Pitzer,
        ReactionDataSource::KentEisenberg,
    ] {
        let mut species: Vec<String> = Vec::new();
        for row in reactions(source).expect("the table parses") {
            if !row.use_reaction {
                continue;
            }
            for (component, _) in stoichiometry(&row.name).expect("the table parses") {
                if !species.contains(&component) {
                    species.push(component);
                }
            }
        }
        let width = species.len();
        assert!(
            (1..31).contains(&width),
            "{source:?} names {width} species, which the sweep's bitmask does not cover"
        );
        for mask in 1u32..(1u32 << width) {
            let set: Vec<String> = (0..width)
                .filter(|index| mask & (1 << index) != 0)
                .map(|index| species[index].clone())
                .collect();
            let result = match reference_potentials(&set, source, kelvins(298.15)) {
                Ok(result) => {
                    subsets += 1;
                    result
                }
                Err(error) => {
                    // **Only the pitzer source refuses**, because it is the only one whose
                    // rows carry evidence. A refusal anywhere else is a defect, not the
                    // gate.
                    assert_eq!(
                        source,
                        ReactionDataSource::Pitzer,
                        "{source:?} with {set:?} is refused: {error}"
                    );
                    refused += 1;
                    continue;
                }
            };
            // A seed is visible as a potential that *is* a component's Gibbs energy of
            // formation. A zero-valued one is not counted: 36 of the databank's rows hold
            // a zero there, so equality would credit a solved-for zero to the seed.
            let did_seed = set.iter().enumerate().any(|(index, name)| {
                formation_properties(name)
                    .expect("the component table parses")
                    .is_some_and(|row| {
                        row.gibbs_energy_of_formation != 0.0
                            && result.potentials[index] == row.gibbs_energy_of_formation
                    })
            });
            if did_seed {
                seeded += 1;
            }
        }
    }
    // The three counts, each of which a silent change would hide:
    assert_eq!(
        subsets + refused,
        14_333,
        "the three tables' loaded species, subset by subset"
    );
    assert_eq!(
        refused, 2_496,
        "subsets the pitzer source's evidence gate refuses - every one of them a set whose \
         survivors include a row that is not VALIDATED"
    );
    assert_eq!(
        seeded, 3_320,
        "answered subsets where the deadlock fallback seeds a component"
    );
}

/// **The gate: the `pitzer` source refuses when any *survivor* is not `VALIDATED`.**
///
/// `ChemicalReactionList.requireValidatedEvidenceForActiveReactions` runs *after* both
/// removals, so what it judges is the survivor set and not the table. The oracle is
/// `PitzerStrictnessProbe.java` / `captures/pitzer_strictness_probe.tsv`: six `SystemPitzer`
/// fluids, of which exactly one is refused, naming `MDEAprot`.
///
/// **It fires on the constructor's second pass.** `MDEAprot` names `MDEA+` among its
/// reactants, and `MDEA+` is not in the feed - the operations constructor's own loop adds
/// the ions its chemistry needs and re-reads the names, so the refusal is raised on a
/// component list the fluid did not start with. That is why the case below supplies the
/// augmented names: they are what `system.getComponentNames()` reports after the attempt.
#[test]
fn the_pitzer_source_refuses_an_unvalidated_active_row() {
    let augmented: Vec<String> = ["MDEA", "water", "CO2", "OH-", "H3O+", "HCO3-"]
        .iter()
        .map(|name| (*name).to_string())
        .collect();
    let error = reference_potentials(&augmented, ReactionDataSource::Pitzer, kelvins(313.15))
        .expect_err("MDEAprot is an active row and is not VALIDATED");
    let message = error.to_string();
    assert!(
        message.contains("MDEAprot"),
        "the refusal names the row it rejected: {message}"
    );
    assert!(
        !message.contains("CO2water") && !message.contains("carbonate"),
        "the four VALIDATED rows are not named: {message}"
    );

    // The same fluid on the standard source is not refused: only the pitzer table carries
    // the evidence column, and `requiresValidatedActiveReactions` is set on that source
    // alone.
    reference_potentials(&augmented, ReactionDataSource::Standard, kelvins(313.15))
        .expect("the standard source requires no evidence");
}

/// **The fluid the gate does not fire on**, from the same capture: CO2/water survives on
/// the four `VALIDATED` rows alone - `CO2water`, `waterreac` and `carbonate` here, the
/// fourth needing H2S the fluid does not carry. So the gate is not "the pitzer source
/// always refuses", which is what a port that judged the *table* rather than the survivors
/// would answer.
#[test]
fn the_four_validated_rows_carry_a_co2_water_fluid() {
    let names: Vec<String> = ["CO2", "water", "OH-", "H3O+", "HCO3-", "CO3--"]
        .iter()
        .map(|name| (*name).to_string())
        .collect();
    let result = reference_potentials(&names, ReactionDataSource::Pitzer, kelvins(298.15))
        .expect("every active row is VALIDATED");
    assert_eq!(
        result.survivors.iter().filter(|mask| **mask == 1.0).count(),
        3,
        "CO2water, waterreac and carbonate survive; the capture is the oracle"
    );
    assert_eq!(result.rank, 3);
}
