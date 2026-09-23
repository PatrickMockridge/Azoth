//! The formula matrix a reactive flash balances on, from
//! `flashops/reactiveflash/FormulaMatrix.java`.
//!
//! **This is not `ChemicalReactionOperations`' element matrix**, and the differences are
//! deliberate. `calcAmatrix` sorts the element rows by symbol and appends an
//! electroneutrality row whose right-hand side is a corrected charge; this collects the
//! elements in **the order they are first encountered** over the phase's components, adds a
//! `Charge` row only when some component is an ion, and gives a component with no element row
//! a pseudo-element of its own so that its conservation is still enforced. The two are the
//! same matrix up to a row permutation wherever both exist, and the row order is the port's
//! business because the flash's own numbers were computed with this one.
//!
//! # The rank decides whether anything reacts at all
//!
//! `NR = NC - rank(A)` independent reactions: with four components over three elements the
//! water-gas shift has one, and a fluid whose rank reaches its component count has none - the
//! flash then falls back to an ordinary VLE flash. The rank is a **Gaussian elimination with
//! a pivot floor of `1e-12`**, reproduced rather than replaced by this crate's exact integer
//! rank: the matrix is not integer (the charge row and the pseudo-elements are not), and a
//! different rank rule would change which fluids take the reactive branch at all.

use azoth_core::{AzothError, Result};

use crate::databank::{element_composition, ionic_charge};

/// The pivot floor `FormulaMatrix.getRank` uses, from `Math.abs(M[row][col]) > maxVal` with
/// `maxVal` starting at `1e-12`.
pub const RANK_PIVOT_FLOOR: f64 = 1e-12;

/// The element names a component with no element row is given, one per component.
///
/// NeqSim's `"Pseudo_" + componentNames[i]`: the component still has to conserve something,
/// so it gets an element of its own rather than a row of zeros - which would leave it
/// unconstrained and let the flash create or destroy it.
pub const PSEUDO_ELEMENT_PREFIX: &str = "Pseudo_";

/// The formula matrix `A[element][component]`, its element names, and whether an ion is
/// present.
#[derive(Debug, Clone, PartialEq)]
pub struct FormulaMatrix {
    /// Element names in the order they are first encountered, **with `Charge` last** when any
    /// component carries a non-zero charge. Not sorted: this is the order NeqSim's
    /// `LinkedHashSet` produced and the order the captured matrices are printed in.
    pub element_names: Vec<String>,
    /// The components, in the caller's order.
    pub component_names: Vec<String>,
    /// `A[element][component]`.
    pub matrix: Vec<Vec<f64>>,
    /// Whether the last row is `Charge` rather than an element.
    pub has_ionic_species: bool,
}

impl FormulaMatrix {
    /// The matrix a phase's components build.
    ///
    /// # Errors
    /// [`AzothError::InvalidInput`] where a component is named and the component databank has
    /// no row for it, because its charge is what the electroneutrality row is built from.
    pub fn build(components: &[String]) -> Result<Self> {
        let mut element_names: Vec<String> = Vec::new();
        let mut per_component: Vec<Vec<(String, f64)>> = Vec::new();
        let mut charges: Vec<f64> = Vec::new();
        let mut has_ionic_species = false;

        for name in components {
            let charge = ionic_charge(name)?.ok_or_else(|| AzothError::InvalidInput {
                field: "components".to_string(),
                reason: format!(
                    "the component databank has no row for {name}, so its ionic charge - \
                     which the electroneutrality row is built from - is not stated anywhere"
                ),
            })?;
            if charge != 0.0 {
                has_ionic_species = true;
            }
            charges.push(charge);

            let composition = element_composition(name)?;
            let pairs = match composition {
                Some(pairs) => pairs,
                // A component with no element row is given an element of its own.
                None => vec![(format!("{PSEUDO_ELEMENT_PREFIX}{name}"), 1.0)],
            };
            for (element, _) in &pairs {
                if !element_names.contains(element) {
                    element_names.push(element.clone());
                }
            }
            per_component.push(pairs);
        }

        // The charge row is added only when something is charged, and it is the last row.
        if has_ionic_species {
            element_names.push("Charge".to_string());
        }

        let base_elements = element_names.len() - usize::from(has_ionic_species);
        let mut matrix = vec![vec![0.0; components.len()]; element_names.len()];
        for (i, pairs) in per_component.iter().enumerate() {
            for (element, coefficient) in pairs {
                if let Some(row) = element_names.iter().position(|held| held == element) {
                    matrix[row][i] = *coefficient;
                }
            }
        }
        if has_ionic_species {
            for (i, charge) in charges.iter().enumerate() {
                matrix[base_elements][i] = *charge;
            }
        }

        Ok(Self {
            element_names,
            component_names: components.to_vec(),
            matrix,
            has_ionic_species,
        })
    }

    /// The matrix's rank, by the elimination `FormulaMatrix.getRank` runs.
    #[must_use]
    pub fn rank(&self) -> usize {
        matrix_rank(&self.matrix, self.component_names.len())
    }

    /// How many independent reactions the fluid can run: `NC - rank(A)`.
    #[must_use]
    pub fn independent_reactions(&self) -> usize {
        self.component_names.len().saturating_sub(self.rank())
    }
}

/// The rank of a bare matrix, which is the elimination `FormulaMatrix.getRank` runs - the same
/// one [`FormulaMatrix::rank`] delegates to, and the one the reactive solve asks for so that it
/// can short-circuit at `NR = 0` where the class does.
///
/// Indexed rather than iterated: a column's elimination mutates every row but the pivot's, and
/// the pivot is chosen *within* the pass, so the loop bounds are the matrix's own and an
/// iterator form would have to rebuild them.
#[must_use]
#[allow(clippy::needless_range_loop)]
pub fn matrix_rank(matrix: &[Vec<f64>], columns: usize) -> usize {
    let rows = matrix.len();
    let mut working = matrix.to_vec();
    let mut row_used = vec![false; rows];
    let mut rank = 0;

    for column in 0..columns {
        if rank == rows {
            break;
        }
        // The pivot is the largest entry in the column, accepted only above the floor.
        let mut pivot_row = None;
        let mut largest = RANK_PIVOT_FLOOR;
        for row in 0..rows {
            if !row_used[row] && working[row][column].abs() > largest {
                largest = working[row][column].abs();
                pivot_row = Some(row);
            }
        }
        let Some(pivot_row) = pivot_row else {
            continue;
        };
        row_used[pivot_row] = true;
        rank += 1;
        let pivot = working[pivot_row][column];
        let pivot_values = working[pivot_row][column..columns].to_vec();
        for row in 0..rows {
            if row != pivot_row && working[row][column] != 0.0 {
                let factor = working[row][column] / pivot;
                for (offset, pivot_value) in pivot_values.iter().enumerate() {
                    working[row][column + offset] -= factor * pivot_value;
                }
            }
        }
    }
    rank
}

#[cfg(test)]
mod tests {
    use super::*;

    fn names(list: &[&str]) -> Vec<String> {
        list.iter().map(|name| (*name).to_string()).collect()
    }

    /// The two states `ReactiveFlashProbe` drives, from
    /// `captures/reactive_flash_probe.tsv`.
    ///
    /// **The element order is the finding**: the water-gas shift's rows are `C`, `O`, `H` -
    /// first encountered over `CO`, `water`, `CO2`, `hydrogen` - and a sorted order would put
    /// `H` before `O`. The matrices are the capture's, row for row.
    #[test]
    fn the_captured_matrices_are_reproduced_in_first_encountered_order() {
        let wgs = FormulaMatrix::build(&names(&["CO", "water", "CO2", "hydrogen"]))
            .expect("the components are in the databank");
        assert_eq!(wgs.element_names, ["C", "O", "H"]);
        assert_eq!(
            wgs.matrix,
            vec![
                vec![1.0, 0.0, 1.0, 0.0],
                vec![1.0, 1.0, 2.0, 0.0],
                vec![0.0, 2.0, 0.0, 2.0],
            ]
        );
        assert!(!wgs.has_ionic_species, "no ion is present");
        assert_eq!(wgs.rank(), 3);
        assert_eq!(wgs.independent_reactions(), 1);

        // The same shape again, on the fluid whose element order is `C`, `H`, `O`.
        let methane = FormulaMatrix::build(&names(&["methane", "water", "CO2", "hydrogen"]))
            .expect("the components are in the databank");
        assert_eq!(methane.element_names, ["C", "H", "O"]);
        assert_eq!(
            methane.matrix,
            vec![
                vec![1.0, 0.0, 1.0, 0.0],
                vec![4.0, 2.0, 0.0, 2.0],
                vec![0.0, 1.0, 2.0, 0.0],
            ]
        );
        assert_eq!(methane.rank(), 3);
        assert_eq!(methane.independent_reactions(), 1);
    }

    /// An ion puts a `Charge` row last, and it is the component databank's charge rather than
    /// anything the caller states.
    ///
    /// From the capture's `co2-water-ions` block, which is the matrix alone - the flash's
    /// ionic branch is one this port refuses, and the matrix is what decides how many
    /// reactions the fluid has. **The rank is 3 and not 4**: `H`, `O`, `C` and `Charge` over
    /// water/CO2/OH-/H3O+ leave exactly one independent reaction, which is water's own
    /// autoionisation.
    #[test]
    fn an_ion_puts_a_charge_row_last() {
        let ionic = FormulaMatrix::build(&names(&["water", "CO2", "OH-", "H3O+"]))
            .expect("the components are in the databank");
        assert!(ionic.has_ionic_species);
        assert_eq!(ionic.element_names, ["H", "O", "C", "Charge"]);
        assert_eq!(
            ionic.matrix,
            vec![
                vec![2.0, 0.0, 1.0, 3.0],
                vec![1.0, 2.0, 1.0, 1.0],
                vec![0.0, 1.0, 0.0, 0.0],
                vec![0.0, 0.0, -1.0, 1.0],
            ]
        );
        assert_eq!(ionic.rank(), 3);
        assert_eq!(ionic.independent_reactions(), 1);
    }

    #[test]
    fn a_component_the_element_table_does_not_carry_conserves_its_own_pseudo_element() {
        // MEG has no element row - the P9 inhibitor staple - and NeqSim gives it one of its
        // own rather than a row of zeros, which would leave it free to appear or vanish.
        let matrix = FormulaMatrix::build(&names(&["water", "MEG"])).expect("both have rows");
        assert_eq!(
            matrix.element_names,
            ["H", "O", "Pseudo_MEG"],
            "the pseudo-element joins the linked set in the order it was reached"
        );
        assert_eq!(
            matrix.matrix,
            vec![vec![2.0, 0.0], vec![1.0, 0.0], vec![0.0, 1.0]]
        );
        // **The rank is 2, not 3**: water's `H` and `O` rows are the same direction scaled -
        // `[2, 0]` against `[1, 0]` - so they contribute one row between them, and the
        // pseudo-element contributes the second. No reaction, which is what a glycol in water
        // is: the capture's `water-meg` block is this matrix.
        assert_eq!(matrix.rank(), 2);
        assert_eq!(matrix.independent_reactions(), 0);
    }
}
