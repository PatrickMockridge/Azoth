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
    solve(
        case,
        case.vector("log_activity")
            .expect("the case states the activity coefficients"),
    )
}

/// The same call with the activity vector **substituted**, which is what localises the
/// residue below: the vector is an input, so a run with a different one is a run of the
/// same equations under the state NeqSim had.
fn solve(case: &TestCase, log_activity: &[f64]) -> azoth_reactions::ChemicalEquilibriumResult {
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
        log_activity,
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

/// **The residue is localised, and it is the activity vector the port holds fixed.**
///
/// The case's own reason for its loose tolerance used to be that the bordered matrix is
/// ill-conditioned. That is not the explanation, and this measures three things:
///
/// * the conditioning's floor is **~1e-10**: the Rust and Python kernels agree to `5.2e-10`
///   on this state, and a one-ulp change to an input moves the answer `6.9e-10`, so an
///   amplification of about `4e6` is what the matrix does - three to five orders of
///   magnitude below the `4.2e-5` the port is away from NeqSim by;
/// * the **largest single contributor is the fixed activity vector**, which
///   `specs/models/reactions/chemical_equilibrium.toml` already records as a divergence:
///   `ChemicalEquilibrium.solve` refreshes `logactivityVec` on every accepted step and this
///   port holds the caller's for the whole solve. Substituting the vector NeqSim has at its
///   answer - which `PitzerReactionBasisProbe` now prints beside the initial one - cuts the
///   gap `8.3x` on CO2 and `7.5x` worst case;
/// * **a residue of `5.0e-6 .. 7.7e-3` is left and conditioning does not explain it.** The
///   next measurement that would split it is the activity vector at *this port's* answer,
///   which is one probe round trip away.
#[test]
fn the_residue_is_the_fixed_activity_vector_and_not_the_conditioning() {
    let spec = model_gen::model(MODEL_ID).expect("the model is in its own table");
    let case = spec
        .cases
        .iter()
        .find(|case| case.id == "pitzer_solute_molality_basis")
        .expect("the molality case is in the spec");
    let captured = |key: &str| -> Vec<f64> {
        let text = std::fs::read_to_string(
            std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
                .join("../../validation/neqsim/captures/pitzer_reaction_basis_probe.tsv"),
        )
        .expect("the capture is committed");
        let block = &text[text
            .find("fluid=co2-water-pitzer")
            .expect("the block is in the capture")..];
        block
            .lines()
            .find_map(|line| line.trim().strip_prefix(&format!("{key}=")))
            .unwrap_or_else(|| panic!("no `{key}` in the capture"))
            .split_whitespace()
            .map(|value| value.parse().expect("a number"))
            .collect()
    };

    let converged = captured("log_activity_after");
    let neqsim = case
        .expected_vector("moles")
        .expect("the case states the answer it pins");

    let worst = |moles: &[f64]| -> Vec<f64> {
        moles
            .iter()
            .zip(neqsim)
            .map(|(got, want)| (got - want).abs() / want.abs().max(f64::MIN_POSITIVE))
            .collect()
    };

    let initial = worst(&solve(case, case.vector("log_activity").expect("the case")).moles);
    let fixed = worst(&solve(case, &converged).moles);

    // The vector NeqSim has at its answer is not the one it started with: the trace ions'
    // activity coefficients move by a factor of about 400 between the two, which is why
    // holding the initial one is not a small approximation.
    let start = case.vector("log_activity").expect("the case");
    let drift = (converged[2] - start[2]).abs() / start[2].abs().max(f64::MIN_POSITIVE);
    assert!(
        drift > 100.0,
        "the ion's activity coefficient moves: {drift}"
    );

    assert!(
        initial[1] > 4.0e-5,
        "CO2 starts 4.2e-5 from NeqSim: {}",
        initial[1]
    );
    assert!(
        fixed[1] < 6.0e-6,
        "the converged vector brings CO2 to 5.0e-6: {}",
        fixed[1]
    );
    let factor = initial[1] / fixed[1];
    assert!(factor > 5.0, "and that is an improvement of {factor}");

    // **And the residue does not go to the kernels' own agreement**, so the vector is not
    // the whole of it.
    assert!(
        fixed[1] > 1.0e-9,
        "a residue remains: {} - if this ever reaches the conditioning's floor, the reason \
         in the case file has been superseded and this test is the one that says so",
        fixed[1]
    );
}
