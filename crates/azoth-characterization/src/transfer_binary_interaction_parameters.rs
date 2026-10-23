//! `characterization.transfer_binary_interaction_parameters` - a table moved between slates.
//!
//! Spec: `specs/models/characterization/transfer_binary_interaction_parameters.toml`. Oracle:
//! `validation/neqsim/captures/characterization_fluid_rebuild_probe.tsv`.
//!
//! # What moves, and what stays
//!
//! The reference's **nonzero** pairs are written onto the target: base components by name,
//! pseudo-components by position - their names are their own fluids' and do not match - and
//! base-to-pseudo pairs alongside the base-to-base ones. Everything the reference does not state
//! is left at the target's value, which is why the answer starts from the target's own matrix
//! rather than from zeros: a target built from the databank has fitted defaults, and dropping
//! them silently would be the failure this shape exists to prevent.
//!
//! # A missing name is skipped once and refused once
//!
//! The class looks the reference's first name up in the target and `continue`s when it is not
//! there, but its setter throws when the *second* name is missing. That asymmetry is reproduced:
//! a reference base component the target does not carry contributes nothing, and a pair whose
//! other end is missing is a mistake in what was asked for.

use azoth_core::{AzothError, Result};

use crate::component_row::InteractionComponent;
use crate::results::TransferBinaryInteractionParametersResult;

/// The reference's nonzero interaction pairs, written onto the target's own matrix.
///
/// # Errors
/// * [`azoth_core::AzothError::InvalidInput`] if either matrix is not `N x N` for its component
///   list, or if a pair the reference states names a base component the target does not carry.
///
/// # Example
/// ```
/// use azoth_characterization::{
///     InteractionComponent, transfer_binary_interaction_parameters,
/// };
///
/// let reference = [
///     InteractionComponent { name: "methane".to_owned(), pseudo_component: false },
///     InteractionComponent { name: "CO2".to_owned(), pseudo_component: false },
/// ];
/// let target = reference.clone();
/// let transferred = transfer_binary_interaction_parameters(
///     &reference,
///     &target,
///     &[vec![0.0, 0.09], vec![0.09, 0.0]],
///     &[vec![0.0, 0.05], vec![0.05, 0.0]],
/// )?;
/// assert_eq!(transferred.target_kij, vec![vec![0.0, 0.09], vec![0.09, 0.0]]);
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
pub fn transfer_binary_interaction_parameters(
    reference_components: &[InteractionComponent],
    target_components: &[InteractionComponent],
    reference_kij: &[Vec<f64>],
    target_kij: &[Vec<f64>],
) -> Result<TransferBinaryInteractionParametersResult> {
    check_matrix("reference_kij", reference_kij, reference_components.len())?;
    check_matrix("target_kij", target_kij, target_components.len())?;

    let mut transferred = target_kij.to_vec();
    let reference_pseudo: Vec<usize> = pseudo_indices(reference_components);
    let target_pseudo: Vec<usize> = pseudo_indices(target_components);

    for (first_index, first) in reference_components.iter().enumerate() {
        if first.pseudo_component {
            continue;
        }
        let Some(first_in_target) = index_of(target_components, &first.name) else {
            continue;
        };
        for (second_index, second) in reference_components.iter().enumerate() {
            if first.name == second.name {
                continue;
            }
            let kij = reference_kij[first_index][second_index];
            if kij == 0.0 {
                continue;
            }
            let second_in_target = if second.pseudo_component {
                let position = reference_pseudo
                    .iter()
                    .position(|index| *index == second_index)
                    .expect("a pseudo-component's index is in its own list");
                match target_pseudo.get(position) {
                    Some(index) => *index,
                    None => continue,
                }
            } else {
                match index_of(target_components, &second.name) {
                    Some(index) => index,
                    None => {
                        return Err(AzothError::invalid_input(
                            "target_components",
                            format!(
                                "no component with name: {} in the target fluid, and the \
                                 reference states a nonzero interaction with it",
                                second.name
                            ),
                        ));
                    }
                }
            };
            write_pair(&mut transferred, first_in_target, second_in_target, kij);
        }
    }

    let common = reference_pseudo.len().min(target_pseudo.len());
    for first in 0..common {
        for second in (first + 1)..common {
            let kij = reference_kij[reference_pseudo[first]][reference_pseudo[second]];
            if kij != 0.0 {
                write_pair(
                    &mut transferred,
                    target_pseudo[first],
                    target_pseudo[second],
                    kij,
                );
            }
        }
    }

    Ok(TransferBinaryInteractionParametersResult {
        target_kij: transferred,
        warnings: Vec::new(),
    })
}

/// The indices of the pseudo-components, in component order.
fn pseudo_indices(components: &[InteractionComponent]) -> Vec<usize> {
    components
        .iter()
        .enumerate()
        .filter(|(_, component)| component.pseudo_component)
        .map(|(index, _)| index)
        .collect()
}

/// The index of a component by name, or `None`. A repeated name resolves to the first entry,
/// as `getComponent(name)` does.
fn index_of(components: &[InteractionComponent], name: &str) -> Option<usize> {
    components
        .iter()
        .position(|component| component.name == name)
}

/// One pair, both triangles, the diagonal never touched.
fn write_pair(matrix: &mut [Vec<f64>], first: usize, second: usize, kij: f64) {
    matrix[first][second] = kij;
    matrix[second][first] = kij;
}

/// A matrix is `N x N` for its component list, or the refusal naming it.
fn check_matrix(field: &str, matrix: &[Vec<f64>], count: usize) -> Result<()> {
    let square = matrix.len() == count && matrix.iter().all(|row| row.len() == count);
    if square {
        return Ok(());
    }
    Err(AzothError::invalid_input(
        field,
        format!("an interaction matrix is {count} x {count}, one row and column per component"),
    ))
}
