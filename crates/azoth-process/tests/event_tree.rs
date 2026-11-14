//! `process.event_tree` against `validation/neqsim/captures/process_event_tree.tsv`.

use azoth_process::models::event_tree::{EventTreeBranch, EventTreeResult, event_tree};

fn close(actual: f64, expected: f64, what: &str) {
    let error = (actual - expected).abs() / expected.abs().max(f64::MIN_POSITIVE);
    assert!(
        error <= 1.0e-9,
        "{what}: {actual} against {expected}, relative error {error:.3e}"
    );
}

fn branch(name: &str, probability_of_yes: f64) -> EventTreeBranch {
    EventTreeBranch {
        name: name.to_string(),
        probability_of_yes,
    }
}

fn tree(frequency: f64, branches: &[EventTreeBranch]) -> EventTreeResult {
    event_tree(frequency, branches).expect("the captured state is valid")
}

#[test]
fn the_two_branch_tree_reproduces_the_capture() {
    let result = tree(1.0e-3, &[branch("Ignition", 0.1)]);
    assert_eq!(result.outcomes.len(), 2);
    assert_eq!(result.outcomes[0].path, "Ignition=N");
    close(result.outcomes[0].probability, 0.9, "the no branch");
    close(
        result.outcomes[0].frequency_per_year,
        9.000_000_000_000_001e-4,
        "the no branch's frequency",
    );
    assert_eq!(result.outcomes[1].path, "Ignition=Y");
    close(result.outcomes[1].probability, 0.1, "the yes branch");
    close(
        result.outcomes[1].frequency_per_year,
        1.0e-4,
        "the yes frequency",
    );
    let probability: f64 = result.outcomes.iter().map(|row| row.probability).sum();
    close(probability, 1.0, "the probabilities sum to one");
    let frequency: f64 = result
        .outcomes
        .iter()
        .map(|row| row.frequency_per_year)
        .sum();
    close(frequency, 1.0e-3, "the frequencies sum to the initiating");
}

#[test]
fn the_three_branch_tree_reproduces_the_capture() {
    let result = tree(
        1.0e-2,
        &[
            branch("Immediate ignition", 0.05),
            branch("Delayed ignition", 0.10),
            branch("ESD success", 0.99),
        ],
    );
    assert_eq!(result.outcomes.len(), 8);
    assert_eq!(
        result.outcomes[0].path,
        "Immediate ignition=N; Delayed ignition=N; ESD success=N"
    );
    close(
        result.outcomes[0].probability,
        0.008_550_000_000_000_007,
        "the first leaf",
    );
    close(
        result.outcomes[4].probability,
        0.846_449_999_999_999_9,
        "the most likely leaf",
    );
    close(
        result.outcomes[7].frequency_per_year,
        4.950_000_000_000_001e-5,
        "the last leaf's frequency",
    );
    let probability: f64 = result.outcomes.iter().map(|row| row.probability).sum();
    close(probability, 1.0, "the probabilities sum to one");
    let frequency: f64 = result
        .outcomes
        .iter()
        .map(|row| row.frequency_per_year)
        .sum();
    close(frequency, 1.0e-2, "the frequencies sum to the initiating");
}

#[test]
fn a_tree_with_no_branches_is_one_leaf() {
    let result = tree(1.0e-3, &[]);
    assert_eq!(result.outcomes.len(), 1);
    assert_eq!(result.outcomes[0].path, "");
    close(result.outcomes[0].probability, 1.0, "the certain leaf");
    close(
        result.outcomes[0].frequency_per_year,
        1.0e-3,
        "the initiating frequency",
    );
}

#[test]
fn certain_and_impossible_branches_leave_one_path() {
    let result = tree(1.0e-3, &[branch("Certain", 1.0), branch("Impossible", 0.0)]);
    assert_eq!(result.outcomes.len(), 4);
    close(result.outcomes[0].probability, 0.0, "Certain=N");
    close(
        result.outcomes[1].probability,
        1.0,
        "Certain=Y; Impossible=N",
    );
    close(result.outcomes[2].probability, 0.0, "Impossible=Y");
    close(result.outcomes[3].probability, 0.0, "both yes");
}

#[test]
fn the_guards_are_the_class_and_the_cap() {
    assert!(
        event_tree(-1.0e-3, &[]).is_err(),
        "the class constructor refuses a negative frequency"
    );
    let above = [branch("Bad", 1.5)];
    assert!(
        event_tree(1.0e-3, &above).is_err(),
        "addBranch refuses a probability above one"
    );
    let nan = [branch("Bad", f64::NAN)];
    assert!(
        event_tree(1.0e-3, &nan).is_err(),
        "a NaN probability is refused"
    );
    let too_many: Vec<EventTreeBranch> = (0..21)
        .map(|index| branch(&format!("B{index}"), 0.5))
        .collect();
    assert!(
        event_tree(1.0e-3, &too_many).is_err(),
        "21 branches exceed the procedure's 20-branch cap"
    );
}
