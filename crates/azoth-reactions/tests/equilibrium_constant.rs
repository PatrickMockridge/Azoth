//! Spec-driven tests for `reactions.equilibrium_constant`.
//!
//! Two oracles, both from `validation/neqsim/ReactionProbe.java`: the arithmetic of the
//! `CO2water` rows, from `captures/reaction_probe.tsv`.
//!
//! **The cases are placed away from the reference temperature on purpose.** `CO2water`'s
//! `K3` is `-39.440767`, so the `ln T` term contributes about `-225` to a `ln K` of
//! `-14.6`: a port that dropped it would be wrong by a factor of `e**225` and would still
//! miss the worked example by only a few percent, because at 298.15 K the term is close to
//! what the fit's other constants absorb. The two reference tests are at 423.15 K and on
//! the Pitzer source, which is where a wrong coefficient shows.

use azoth_core::spec::TestCase;
use azoth_core::units::kelvins;
use azoth_reactions::databank::{ReactionDataSource, reactions};
use azoth_reactions::equilibrium_constant::{equilibrium_constant, ln_k_of};
use azoth_reactions::spec_gen;
use azoth_test_support as common;

const CALC_ID: &str = "reactions.equilibrium_constant";

fn source_of(case: &TestCase) -> ReactionDataSource {
    common::input_str(case, "source")
        .parse()
        .unwrap_or_else(|e| panic!("test `{}` names a source: {e}", case.id))
}

fn call(case: &TestCase) -> azoth_reactions::EquilibriumConstantResult {
    equilibrium_constant(
        source_of(case),
        common::input_str(case, "reaction"),
        kelvins(common::input(case, "T")),
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
                    result.ln_k,
                    common::expected(case, "ln_k"),
                    case.tolerance,
                    &format!("{}::{} (ln K)", spec.id, case.id),
                );
                common::assert_close(
                    result.k,
                    common::expected(case, "k"),
                    case.tolerance,
                    &format!("{}::{} (K)", spec.id, case.id),
                );
                common::assert_close(
                    result.reaction_heat.value,
                    common::expected(case, "reaction_heat"),
                    case.tolerance,
                    &format!("{}::{} (heat of reaction)", spec.id, case.id),
                );
                common::assert_consistent(&result, &format!("{}::{}", spec.id, case.id));
                common::assert_warnings_agree_with_spec(
                    spec,
                    &result.warnings,
                    |quantity| match quantity {
                        "T" => Some(common::input(case, "T")),
                        "ln_k" => Some(result.ln_k),
                        "k" => Some(result.k),
                        _ => None,
                    },
                    &format!("{}::{}", spec.id, case.id),
                );
            }
            "property" => match case.property {
                // `ln K` carries no dimension to round-trip: its inputs are a name and a
                // temperature and its outputs are dimensionless but for the heat. The
                // declaration is kept so the spec states that rather than omitting it.
                Some("unit_round_trip") => {}
                other => panic!("{}::{}: unhandled property {other:?}", spec.id, case.id),
            },
            other => panic!("{}::{}: unknown test kind {other:?}", spec.id, case.id),
        }
        executed += 1;
    }
    assert!(
        executed >= 4,
        "expected several active cases, ran {executed}"
    );
}

/// **The `ln T` term is what a wrong port drops**, and this is the state that shows it.
///
/// Dropping `K3 ln T` changes `ln K` by about `+225` here, not by a rounding: the term is
/// the largest single contribution to the logarithm of a constant that is itself `1e-7`.
/// So the assertion is on the size of the term rather than on the answer, because it is
/// the term that a reader checking this file wants to see is real.
#[test]
fn the_logarithmic_term_is_not_negligible() {
    let result = equilibrium_constant(ReactionDataSource::Standard, "CO2water", kelvins(423.15))
        .expect("the row exists");
    let row = azoth_reactions::databank::reaction(ReactionDataSource::Standard, "CO2water")
        .expect("the table parses")
        .expect("the row exists");
    let term = row.coefficients[2] * 423.15_f64.ln();

    assert!(
        term.abs() > 200.0,
        "the `K3 ln T` term is {term}, and the spec's assumptions say a port that \
         dropped it fails only away from the reference temperature"
    );
    assert!(
        result.ln_k < 0.0,
        "the constant is far below one: {result:?}"
    );
}

/// The heat and the derivative are one relation, so a sign error in either shows here.
#[test]
fn the_heat_follows_the_derivative() {
    let result = equilibrium_constant(ReactionDataSource::Standard, "CO2water", kelvins(273.15))
        .expect("the row exists");
    let from_identity = result.ln_k_derivative
        * 273.15
        * 273.15
        * azoth_reactions::equilibrium_constant::GAS_CONSTANT;
    assert!(
        (result.reaction_heat.value - from_identity).abs() / from_identity.abs() < 1.0e-14,
        "the heat is `d(ln K)/dT * R T**2` and nothing else"
    );
    // **Exothermic in the direction the table writes it** at this temperature, and the
    // sign is a statement about that direction rather than about the chemistry.
    assert!(
        result.reaction_heat.value > 0.0,
        "CO2water's forward direction is endothermic at 273.15 K: {result:?}"
    );
}

/// A name no source carries, and a source this library does not know, are different
/// failures - one is the data, the other is the caller.
#[test]
fn the_two_refusals_are_distinct() {
    let unknown = equilibrium_constant(ReactionDataSource::Standard, "no-such", kelvins(298.15))
        .expect_err("no row");
    assert!(
        matches!(unknown, azoth_core::AzothError::PropertyUnavailable { .. }),
        "{unknown:?}"
    );

    assert!(
        "phreeqc".parse::<ReactionDataSource>().is_err(),
        "a source this library does not carry is refused rather than defaulted"
    );
}

/// A dormant row is answered rather than hidden, matching `getChemicalReaction`, which
/// selects on the name alone and never reads the flag.
#[test]
fn a_dormant_row_is_answered() {
    let result = equilibrium_constant(
        ReactionDataSource::Standard,
        "methanecombustion",
        kelvins(298.15),
    )
    .expect("the factory answers for it, so this does");
    assert!(
        result.ln_k.is_finite(),
        "a finite answer for a row the reaction set will not load: {result:?}"
    );
}

/// **The table carries `MDEAprot` twice, and the two fits are about 2.76 apart at 298 K.**
///
/// This is the premise behind `ln_k_of`. `equilibrium_constant` selects on the name alone and
/// answers the *first* row - which for this name is the disabled Austgen fit - while the
/// reaction set loads the enabled Huttenhuis row. A caller that has already selected a row has
/// to evaluate *that* row: `reactive_phase_equilibrium`'s certification filters on `usereaction`
/// and then measures the kept row's species against a constant, so a by-name relookup measures
/// one fit against another's `ln K` - 2.76, against a `2e-6` tolerance.
#[test]
fn a_duplicated_reaction_name_carries_two_different_fits() {
    let rows: Vec<_> = reactions(ReactionDataSource::Standard)
        .expect("the table loads")
        .into_iter()
        .filter(|row| row.name == "MDEAprot")
        .collect();
    assert_eq!(rows.len(), 2, "the table carries MDEAprot twice");
    let first = &rows[0];
    let enabled = rows
        .iter()
        .find(|row| row.use_reaction)
        .expect("one of the two rows is enabled");
    assert!(!first.use_reaction, "the first row is the disabled one");

    let t = kelvins(298.15);
    let by_name = equilibrium_constant(ReactionDataSource::Standard, "MDEAprot", t)
        .expect("the name resolves")
        .ln_k;
    assert!(
        (by_name - ln_k_of(first, t)).abs() < 1e-12,
        "the by-name lookup answers the first row"
    );
    assert!(
        (ln_k_of(enabled, t) - by_name).abs() > 2.5,
        "the enabled fit is a different constant: {} against {by_name}",
        ln_k_of(enabled, t)
    );
}
