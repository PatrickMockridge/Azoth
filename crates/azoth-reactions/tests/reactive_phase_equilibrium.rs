//! Spec-driven tests for `reactions.reactive_phase_equilibrium`.
//!
//! The oracle is `validation/neqsim/ReactionOperationsProbe.java`, whose capture is
//! `captures/reaction_operations_probe.tsv`.
//!
//! **The three cases are three different arguments.** The first is the ordinary one, an
//! aqueous phase whose charge row happens to be zero. The second is a brine whose
//! chloride is not modelled, so the row carries its sodium and a port that assumed a
//! zero there passes the first case and fails this one. The third never solves at all:
//! one gas phase, no aqueous and no liquid, so `getReactivePhaseIndex` returns `-1` and
//! NeqSim's facade returns without solving - and the matrix, the element amounts and the
//! reference potentials are still built, as its constructor builds them.
//!
//! Whether a case was skipped is a boolean, and the generator carries numbers, vectors
//! and strings, so it is asserted here - which is where a flag belongs anyway.

use azoth_core::spec::TestCase;
use azoth_reactions::databank::ReactionDataSource;
use azoth_reactions::model_gen;
use azoth_reactions::reactive_phase::reactive_phase_index;
use azoth_reactions::chemical_equilibrium::ConcentrationBasis;
use azoth_reactions::reactive_phase_equilibrium::{
    ELEMENT_BALANCE_RESIDUAL_TOLERANCE_MOLES, REACTION_LOG_RESIDUAL_TOLERANCE,
    REACTIVE_PHASE_CHARGE_TOLERANCE_MOLES, ReactionSeed, ReactivePhaseEquilibriumResult,
    reactive_phase_equilibrium,
};
use azoth_test_support as common;

const MODEL_ID: &str = "reactions.reactive_phase_equilibrium";

/// The cases that take the solve, and the one that does not, by id.
fn takes_the_solve(id: &str) -> bool {
    id != "a_single_gas_phase_is_skipped"
}

fn call(case: &TestCase) -> ReactivePhaseEquilibriumResult {
    let source: ReactionDataSource = common::input_str(case, "source")
        .parse()
        .expect("the case names one of the three sources");
    let seed: ReactionSeed = common::input_str(case, "seed")
        .parse()
        .expect("the case says where the solve starts");
    let components: Vec<String> = case
        .list("components")
        .expect("the case states its components")
        .iter()
        .map(|name| (*name).to_string())
        .collect();
    let moles = case
        .vector("moles")
        .expect("the case states the phase composition");
    let log_activity = case
        .vector("log_activity")
        .expect("the case states the activity coefficients");

    reactive_phase_equilibrium(
        &components,
        source,
        common::input_str(case, "phase"),
        moles,
        common::input(case, "phase_charge"),
        common::input(case, "phase_moles"),
        case.flag("whole_system")
            .expect("the case says whether the phase is the whole system"),
        log_activity,
        common::input(case, "T"),
        common::input(case, "max_iterations") as u32,
        common::input(case, "tolerance"),
        seed,
        common::input_str(case, "concentration_basis")
            .parse()
            .expect("the case states the basis"),
    )
    .unwrap_or_else(|e| panic!("test `{}` should compute but failed: {e}", case.id))
}

#[test]
fn every_case_in_the_spec() {
    let spec = model_gen::model(MODEL_ID).expect("the model should be in its own table");
    assert!(!spec.cases.is_empty(), "the model has no cases");

    for case in spec.cases {
        let tolerance = case.tolerance;
        let result = call(case);

        assert_eq!(
            result.skipped,
            !takes_the_solve(case.id),
            "{}: skipped is {}",
            case.id,
            result.skipped
        );

        let expected_moles = case
            .expected_vector("moles")
            .expect("the case states expected moles");
        assert_eq!(result.moles.len(), expected_moles.len(), "{}", case.id);
        for (i, (actual, expected)) in result.moles.iter().zip(expected_moles).enumerate() {
            common::assert_close(
                *actual,
                *expected,
                tolerance,
                &format!("{}: moles[{i}]", case.id),
            );
        }

        // **The seed, which is where the case and the code can disagree about the start.**
        // `seed_applied` is a flag on the result and `TestCase` carries flags for inputs
        // and not for expectations, so it is read from the raw spec here - the entry that
        // says so is `SELF_ASSERTED_EXPECTATIONS` in `tools/gen_registry.py`.
        if let Some(flag) = case.flag("seed_applied") {
            assert_eq!(result.seed_applied, flag, "{}: seed_applied", case.id);
        }
        if let Some(expected_seed) = case.expected_vector("seed_moles") {
            assert_eq!(result.seed_moles.len(), expected_seed.len(), "{}", case.id);
            for (i, (actual, expected)) in result.seed_moles.iter().zip(expected_seed).enumerate() {
                common::assert_close(
                    *actual,
                    *expected,
                    tolerance,
                    &format!("{}: seed_moles[{i}]", case.id),
                );
            }
        }

        let expected_b = case.expected_vector("b").expect("the case states b");
        assert_eq!(result.b.len(), expected_b.len(), "{}", case.id);
        for (i, (actual, expected)) in result.b.iter().zip(expected_b).enumerate() {
            common::assert_close(
                *actual,
                *expected,
                tolerance,
                &format!("{}: b[{i}]", case.id),
            );
        }

        let expected_potentials = case
            .expected_vector("chem_ref")
            .expect("the case states the reference potentials");
        for (i, (actual, expected)) in result.chem_ref.iter().zip(expected_potentials).enumerate() {
            common::assert_close(
                *actual,
                *expected,
                tolerance,
                &format!("{}: chem_ref[{i}]", case.id),
            );
        }

        // **The matrix's own values are asserted by the Python case runner**, which reads
        // `a_matrix` from the spec: the generator that turns a spec into this table carries
        // scalars, vectors and strings and has no field for a matrix expectation, so the
        // case's block is read there and registered in `SELF_ASSERTED_EXPECTATIONS` rather
        // than here. Its *shape* is checked here, because it is what the operation builds.
        let species = expected_moles.len();
        assert_eq!(
            result.a_matrix.len(),
            result.b.len(),
            "{}: one row per element amount, and the charge row is one of them",
            case.id
        );
        for (row, values) in result.a_matrix.iter().enumerate() {
            assert_eq!(
                values.len(),
                species,
                "{}: a_matrix row {row} has one column per component",
                case.id
            );
        }

        common::assert_close(
            f64::from(result.iterations),
            common::expected(case, "iterations"),
            1.0,
            &format!("{}: iterations", case.id),
        );
    }
}

/// The phase search itself, against the capture's own listings.
///
/// The capture prints each fluid's phase type names in system order beside the index
/// NeqSim's private method returned, both before and after a flash. That is ten
/// (listing, answer) pairs, and four of them take the fallback branch rather than the
/// aqueous one - so this is where the two-pass order is checked rather than assumed.
#[test]
fn the_phase_search_reproduces_the_capture() {
    // (labels, the index the capture recorded, which fluid and when)
    let captured: &[(&[&str], Option<usize>, &str)] = &[
        (
            &["gas", "liquid"],
            Some(1),
            "co2-water-aqueous before its flash",
        ),
        (
            &["gas", "aqueous"],
            Some(1),
            "co2-water-aqueous after its flash",
        ),
        (
            &["gas", "liquid"],
            Some(1),
            "co2-h2s-water before its flash",
        ),
        (
            &["gas", "aqueous"],
            Some(1),
            "co2-h2s-water after its flash",
        ),
        (
            &["gas", "liquid"],
            Some(1),
            "co2-water-warm before its flash",
        ),
        (
            &["gas", "aqueous"],
            Some(1),
            "co2-water-warm after its flash",
        ),
        (
            &["gas", "liquid"],
            Some(1),
            "gas-condensate-no-aqueous before its flash",
        ),
        // The skip: one gas phase, so no aqueous and nothing for the fallback.
        (&["gas"], None, "gas-condensate-no-aqueous after its flash"),
        (
            &["gas", "liquid"],
            Some(1),
            "bicarbonate-brine before its flash",
        ),
        (
            &["gas", "aqueous"],
            Some(1),
            "bicarbonate-brine after its flash",
        ),
    ];

    for (labels, expected, what) in captured {
        assert_eq!(
            reactive_phase_index(labels),
            *expected,
            "{what}: {labels:?} should search to {expected:?}"
        );
    }

    // And the ordering rule on its own: a liquid listed first does not pre-empt an
    // aqueous phase later, because the first pass is over every phase.
    assert_eq!(reactive_phase_index(&["liquid", "aqueous"]), Some(1));
    assert_eq!(reactive_phase_index(&["oil", "liquid"]), Some(0));
}

/// A skipped phase still comes back with its matrix, its amounts and its potentials.
///
/// NeqSim's `solveChemEq` returns before rebuilding them, but its *constructor* built
/// them, and the capture prints them for the skipped fluid. This is the assertion that
/// the skip is a conditional on the solve rather than on the operation.
#[test]
fn a_skip_still_carries_what_the_facade_built() {
    let spec = model_gen::model(MODEL_ID).expect("the model is in its table");
    let case = spec
        .cases
        .iter()
        .find(|case| case.id == "a_single_gas_phase_is_skipped")
        .expect("the spec carries the skip case");

    let result = call(case);
    assert!(result.skipped);
    assert!(!result.converged, "a skip is not a convergence");
    assert_eq!(result.iterations, 0);
    assert_eq!(result.a_matrix.len(), 4, "elements plus the charge row");
    assert_eq!(result.b.len(), 4);
    assert_eq!(result.chem_ref.len(), 6);
    assert!(
        result.chem_ref.iter().any(|potential| *potential != 0.0),
        "the potentials are built at the skip, not left at zero"
    );
}

/// A component with no element row is refused rather than quietly dropped.
#[test]
fn a_substance_outside_the_element_table_is_refused() {
    // MEG is a P9 inhibitor staple and has no row, so a fluid carrying it cannot enter
    // the balance. NeqSim builds the row from whatever components do have one; this
    // refuses, because a smaller balance still looks like an answer.
    let error = reactive_phase_equilibrium(
        &["CO2".to_string(), "MEG".to_string()],
        ReactionDataSource::Standard,
        "aqueous",
        &[0.01, 1.0],
        0.0,
        1.01,
        false,
        &[0.0, 0.0],
        298.15,
        100,
        1e-8,
        ReactionSeed::None,
        ConcentrationBasis::MoleFraction,
    )
    .expect_err("MEG has no formula row");
    assert!(error.to_string().contains("MEG"), "{error}");
}

/// **The certificate is NeqSim's gate, and it is a different question from convergence.**
///
/// `solveChemEq` returns `converged && r1 && r2 && r3` - a converged solve whose reaction log
/// residual is over `2e-6` is reported as a failure - and on the five captured fluids that is
/// what happens. These cases model the *direct* solve, which passes its own gate: the first
/// converges and is certified, the second converges nowhere and is not, and the skipped phase
/// certifies nothing because there is no phase to certify on.
///
/// The four numbers themselves are not pinned here. The captures' residuals come from
/// `solveChemEq` and this path does not, so a value taken from them would be a claim about a
/// different call; what is pinned is what the flags mean.
#[test]
fn the_certificate_is_the_gate_and_not_the_solver_flag() {
    let spec = model_gen::model(MODEL_ID).expect("the model should be in its own table");
    for case in spec.cases {
        let result = call(case);
        if result.skipped {
            assert_eq!(result.refinements, 0, "{}: nothing ran", case.id);
            assert!(!result.certified, "{}: nothing to certify", case.id);
            assert!(
                result.net_charge_moles.is_nan(),
                "{}: NeqSim's net charge is NaN without a reactive phase",
                case.id
            );
            continue;
        }
        assert_eq!(result.refinements, 1, "{}: the first refinement", case.id);

        let under = result.max_reaction_log_residual <= REACTION_LOG_RESIDUAL_TOLERANCE
            && result.net_charge_moles.abs() <= REACTIVE_PHASE_CHARGE_TOLERANCE_MOLES
            && result.max_element_residual <= ELEMENT_BALANCE_RESIDUAL_TOLERANCE_MOLES;
        assert_eq!(
            result.certified,
            under,
            "{}: certified={} against residuals {} / {} / {}",
            case.id,
            result.certified,
            result.max_reaction_log_residual,
            result.net_charge_moles,
            result.max_element_residual
        );
        assert_eq!(
            result.converged, result.certified,
            "{}: the returned flag is the gate's",
            case.id
        );
        // The element balance is conserved by construction, so its residual is rounding.
        assert!(
            result.max_element_residual < 1e-9,
            "{}: the elements came out at {}",
            case.id,
            result.max_element_residual
        );
    }
}
