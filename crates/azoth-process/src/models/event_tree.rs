//! `process.event_tree` - the event-tree enumeration as a registered id.
//!
//! Spec: `specs/models/process/event_tree.toml`. The arithmetic is
//! `process/safety/risk/eta/EventTreeAnalyzer.java` at the pin: the full `2^n`
//! enumeration in bit order, the path spelling and the products.

use azoth_core::{AzothError, Result, Warning, apply_checks};
use serde::Serialize;

use crate::executor::json::warnings as wire_warnings;
use crate::model_gen;

/// One pivotal-event branch, as the caller states it.
#[derive(Debug, Clone, PartialEq)]
pub struct EventTreeBranch {
    /// The pivotal event's name, used only in the path.
    pub name: String,
    /// The probability the event occurs; the no branch carries its complement.
    pub probability_of_yes: f64,
}

/// One leaf of the tree.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct EventTreeOutcome {
    /// The branch states joined by `"; "`, each `name=Y` or `name=N`.
    pub path: String,
    /// The product of the branch probabilities and complements on this path.
    pub probability: f64,
    /// The probability times the initiating frequency, in 1/year.
    pub frequency_per_year: f64,
}

/// Result of `process.event_tree`.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct EventTreeResult {
    /// One row per leaf, in the class's bit order.
    pub outcomes: Vec<EventTreeOutcome>,
    /// Caveats.
    #[serde(serialize_with = "wire_warnings")]
    pub warnings: Vec<Warning>,
}

/// The full `2^n` event tree in bit order.
///
/// # Errors
/// `InvalidInput` for a negative initiating frequency, a branch probability outside
/// `[0, 1]` or a NaN, or more than the algorithm's `max_iterations` branches.
pub fn event_tree(
    initiating_frequency_per_year: f64,
    branches: &[EventTreeBranch],
) -> Result<EventTreeResult> {
    let spec = &model_gen::EVENT_TREE_SPEC;
    let mut warnings = Vec::new();
    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "initiating_frequency_per_year" => Some(initiating_frequency_per_year),
            _ => None,
        },
        &mut warnings,
    )?;
    if !initiating_frequency_per_year.is_finite() {
        return Err(AzothError::invalid_input(
            "initiating_frequency_per_year",
            "a non-finite frequency is not a number this calculation publishes",
        ));
    }

    let cap = spec
        .algorithm
        .map_or(0, |algorithm| algorithm.max_iterations) as usize;
    if branches.len() > cap {
        return Err(AzothError::invalid_input(
            "branches",
            format!(
                "the tree enumerates 2^n leaves, and {} branches exceeds the {cap} the \
                 procedure carries",
                branches.len()
            ),
        ));
    }
    for (index, branch) in branches.iter().enumerate() {
        let probability = branch.probability_of_yes;
        if probability.is_nan() || !(0.0..=1.0).contains(&probability) {
            return Err(AzothError::invalid_input(
                format!("branches[{index}].probability_of_yes"),
                format!("a branch probability is in [0, 1], and {probability} is not"),
            ));
        }
    }

    let leaves = 1usize << branches.len();
    let mut outcomes = Vec::with_capacity(leaves);
    for index in 0..leaves {
        let mut probability = 1.0;
        let mut path = String::new();
        for (bit, branch) in branches.iter().enumerate() {
            let yes = (index >> bit) & 1 == 1;
            probability *= if yes {
                branch.probability_of_yes
            } else {
                1.0 - branch.probability_of_yes
            };
            if bit > 0 {
                path.push_str("; ");
            }
            path.push_str(&branch.name);
            path.push_str(if yes { "=Y" } else { "=N" });
        }
        outcomes.push(EventTreeOutcome {
            path,
            probability,
            frequency_per_year: probability * initiating_frequency_per_year,
        });
    }

    Ok(EventTreeResult { outcomes, warnings })
}
