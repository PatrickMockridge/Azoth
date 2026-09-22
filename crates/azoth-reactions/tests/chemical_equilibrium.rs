//! Spec-driven tests for `reactions.chemical_equilibrium`.
//!
//! The oracle is `validation/neqsim/ChemicalEquilibriumProbe.java`, whose capture is
//! `captures/chemical_equilibrium_probe.tsv`.
//!
//! **One of the two cases is a failure, and that is the point of it.** NeqSim's solver
//! does not converge on the CO2-H2S-water fluid - the error settles at 117 against a
//! tolerance of `1e-8` - and it returns its last iterate with `converged` false. A port
//! that converged there would be answering a question NeqSim refuses to, so the case pins
//! the iteration count, the error and the composition the failure leaves behind.
//!
//! The `converged` flag cannot live in the case file - the generator carries numbers,
//! vectors and strings - so it is asserted here, which is where a boolean belongs anyway.

use azoth_core::spec::TestCase;
use azoth_reactions::chemical_equilibrium::{ConcentrationBasis, chemical_equilibrium};
use azoth_reactions::model_gen;
use azoth_test_support as common;

const MODEL_ID: &str = "reactions.chemical_equilibrium";

fn call(case: &TestCase) -> azoth_reactions::ChemicalEquilibriumResult {
    let b = case.vector("b").expect("the case states b");
    let moles = case
        .vector("moles")
        .expect("the case states the starting moles");
    let flat = case
        .matrix("a_matrix")
        .expect("the case states the element matrix");

    let n_elements = b.len();
    let species = moles.len();
    assert_eq!(
        flat.len(),
        n_elements * species,
        "{}: the element matrix is elements x components",
        case.id
    );
    let a_matrix: Vec<Vec<f64>> = (0..n_elements)
        .map(|e| flat[e * species..(e + 1) * species].to_vec())
        .collect();

    chemical_equilibrium(
        &a_matrix,
        b,
        case.flag("whole_system")
            .expect("the case says whether the phase is the whole system"),
        moles,
        case.vector("chem_ref")
            .expect("the case states the reference potentials"),
        case.vector("log_activity")
            .expect("the case states the activity coefficients"),
        common::input(case, "T"),
        common::input(case, "max_iterations") as u32,
        common::input(case, "tolerance"),
        common::input_str(case, "concentration_basis")
            .parse()
            .expect("the case names a concentration basis"),
        common::input(case, "solvent_weight"),
        case.vector("solvent_mask")
            .expect("the case states the reference-state split"),
        common::input(case, "phase_moles"),
    )
    .unwrap_or_else(|e| panic!("test `{}` should compute but failed: {e}", case.id))
}

#[test]
fn every_case_in_the_spec() {
    let spec = model_gen::model(MODEL_ID).expect("the model should be in its own table");
    assert!(!spec.cases.is_empty());

    for case in spec.cases {
        let result = call(case);
        let context = format!("{}::{}", spec.id, case.id);

        // **A case may decline to pin the composition, and the one that does is the one
        // that never converged.** An iterate the solver gave up on is where two rounding
        // paths diverge rather than an answer, so the failure is asserted - the flag and
        // the large error - and the point it stopped at is not.
        match case.expected_vector("moles") {
            Some(want) => {
                assert_eq!(
                    result.moles.len(),
                    want.len(),
                    "{context}: one mole number per entry"
                );
                for (i, expected) in want.iter().enumerate() {
                    common::assert_close(
                        result.moles[i],
                        *expected,
                        case.tolerance,
                        &format!("{context} moles[{i}]"),
                    );
                }
            }
            None => assert!(
                !result.converged,
                "{context}: a case that pins no composition must be one that failed",
            ),
        }

        // **The pass count is asserted as a band, and the reason is measured.** The
        // solve is a Newton iteration whose two `M`-matrix solves are diagonal and exact,
        // but whose bordered `(NELE+1)` system goes through LU - and this port's LU takes
        // a different rounding path from Jama's. On the captured fluid that costs one
        // extra pass (16 against 15) and lands `1.5e-06` relative away on the trace ions,
        // `6e-11` on water. A fixed point reached through a different rounding path is a
        // different point, and the band says so rather than pretending to bit equality.
        let iterations = common::expected(case, "iterations");
        assert!(
            (f64::from(result.iterations) - iterations).abs() <= 2.0,
            "{context}: {} passes against NeqSim's {iterations}",
            result.iterations
        );
        let error = common::expected(case, "error");
        assert!(
            (result.error - error).abs() <= 1e-6 * error.abs().max(1.0),
            "{context}: error {} against NeqSim's {error}",
            result.error
        );
    }
}

/// The flag the case file cannot carry, asserted per case by name.
///
/// **Both directions matter.** A solve that reports success where the error never came in
/// is worse than one that fails, because a caller reads `moles` on the strength of it.
#[test]
fn the_convergence_flags_are_what_neqsim_reports() {
    let spec = model_gen::model(MODEL_ID).expect("the model should be in its own table");
    for case in spec.cases {
        let result = call(case);
        let expected = match case.id {
            "co2_water_from_the_phase_composition" => true,
            "pitzer_solute_molality_basis" => true,
            "co2_h2s_water_does_not_converge" => false,
            other => panic!("{other}: a case with no stated convergence"),
        };
        assert_eq!(
            result.converged, expected,
            "{}::{}: converged={} with error {} against tolerance {}",
            spec.id, case.id, result.converged, result.error, case.tolerance
        );
        // The two flags must also agree with the error itself, which is what NeqSim's
        // `converged = error < maxError` says - so a case whose numbers moved must move
        // this with them.
        if expected {
            assert!(
                result.error < common::input(case, "tolerance"),
                "{}::{}: reported convergence at {}",
                spec.id,
                case.id,
                result.error
            );
        } else {
            assert!(
                result.error > common::input(case, "tolerance"),
                "{}::{}: reported failure at {}",
                spec.id,
                case.id,
                result.error
            );
        }
    }
}

/// The element matrix's last row is the charge balance, and its right-hand side is zero.
///
/// Asserted rather than assumed, because a port that dropped the row would still solve -
/// it would conserve the elements and not the charge, and every captured number would
/// move.
#[test]
fn the_captured_element_matrices_carry_a_charge_row() {
    let spec = model_gen::model(MODEL_ID).expect("the model should be in its own table");
    for case in spec.cases {
        let b = case.vector("b").expect("the case states b");
        assert_eq!(
            b.last(),
            Some(&0.0),
            "{}: the last element amount is the electroneutrality target",
            case.id
        );
        let flat = case.matrix("a_matrix").expect("the element matrix");
        let species = case.vector("moles").expect("moles").len();
        let last = &flat[flat.len() - species..];
        assert!(
            last.iter().any(|value| *value < 0.0),
            "{}: an electroneutrality row has anions in it: {last:?}",
            case.id
        );
    }
}

/// A shape the solver cannot read is refused rather than indexed out of bounds.
#[test]
fn a_matrix_with_the_wrong_width_is_refused() {
    let error = chemical_equilibrium(
        &[vec![1.0, 0.0]],
        &[1.0],
        false,
        &[1.0, 1.0, 1.0],
        &[0.0, 0.0, 0.0],
        &[0.0, 0.0, 0.0],
        298.15,
        100,
        1e-8,
        ConcentrationBasis::MoleFraction,
        0.0,
        &[],
        3.0,
    )
    .expect_err("one column against three components");
    assert!(
        matches!(error, azoth_core::AzothError::InvalidInput { .. }),
        "{error:?}"
    );
}
