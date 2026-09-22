//! The two linear-algebra routines the reference-potential solve needs.
//!
//! NeqSim calls `Jama.Matrix.rank()` and `Jama.Matrix.solve()` here, and JAMA's rank is
//! a singular-value count against `max(m, n) * s[0] * 2**-52`. **This does not port
//! JAMA**, and the reason is measured rather than stylistic.
//!
//! # The rank is exact, and the tolerance is not load-bearing
//!
//! Both rank tests run on matrices built from **stoichiometric coefficients alone** -
//! `reacGMatrix`'s last column holds `-R T ln K` and is excluded, because the column
//! count is one less than the row width. Every coefficient in the whole of NeqSim's
//! `stoccoefdata` is an integer: the table's 151 rows carry only `0`, `1`, `-1`, `2` and
//! `-2`.
//!
//! So the rank question is asked of an integer matrix, where the answer is exact.
//! `validation/neqsim/captures/reference_potential_probe.tsv` records what JAMA's own
//! test does with three real fluids: the smallest singular value it keeps is `0.2461`,
//! `0.2536` and `0.3178`, while the cutoff it compares against is `2.1e-15`, `3.8e-15`
//! and `4.8e-15`. **Fourteen orders of magnitude**, on every matrix, and in all three
//! cases the test finds a full-rank basis and rejects nothing.
//!
//! A port of the Golub-Reinsch SVD would therefore be a large amount of arithmetic whose
//! tolerance could not change a single decision here - and, being floating point, could
//! only be *less* reliable than the exact test below. What would overturn this is a
//! coefficient table that stopped being integral, which is why
//! [`rank_of_integer_matrix`] refuses one rather than rounding it.

use azoth_core::{AzothError, Result};

/// The exact rank of a matrix whose entries are all integers.
///
/// **Fraction-free (Bareiss) elimination**: every intermediate stays an integer, so no
/// pivot is ever compared against a tolerance and the answer is exact rather than
/// numerical. Intermediate products can exceed `i64` for a large matrix, so the
/// elimination carries `i128`.
///
/// # Errors
/// [`AzothError::InvalidInput`] if any entry is not integral. A rounded coefficient
/// would make this an approximation of a rank rather than the rank, and the caller
/// cannot tell the difference from the answer.
pub fn rank_of_integer_matrix(matrix: &[Vec<f64>]) -> Result<usize> {
    let rows = matrix.len();
    if rows == 0 {
        return Ok(0);
    }
    let cols = matrix[0].len();
    if cols == 0 {
        return Ok(0);
    }

    let mut a = Vec::with_capacity(rows);
    for (i, row) in matrix.iter().enumerate() {
        if row.len() != cols {
            return Err(AzothError::InvalidInput {
                field: "reaction_matrix".to_string(),
                reason: format!("row {i} has {} entries against {cols}", row.len()),
            });
        }
        let mut converted = Vec::with_capacity(cols);
        for (j, value) in row.iter().enumerate() {
            if value.fract() != 0.0 || !value.is_finite() {
                return Err(AzothError::InvalidInput {
                    field: "reaction_matrix".to_string(),
                    reason: format!(
                        "entry [{i}][{j}] is {value}, and the rank test is exact only over \
                         integers. NeqSim's rank is a singular-value count and would accept \
                         it; this refuses rather than approximating a rank"
                    ),
                });
            }
            converted.push(*value as i128);
        }
        a.push(converted);
    }

    let mut rank = 0usize;
    let mut previous_pivot: i128 = 1;

    for col in 0..cols {
        if rank == rows {
            break;
        }
        // Partial pivoting on the largest magnitude, so the elimination does not depend
        // on row order. Rank is invariant under row permutation, so any pivot choice
        // gives the same answer; this one keeps the intermediates small.
        let mut pivot_row = None;
        let mut best = 0i128;
        for (i, row) in a.iter().enumerate().skip(rank) {
            if row[col].abs() > best {
                best = row[col].abs();
                pivot_row = Some(i);
            }
        }
        let Some(pivot_row) = pivot_row else { continue };
        if best == 0 {
            continue;
        }
        a.swap(rank, pivot_row);

        for i in (rank + 1)..rows {
            // **Read the multiplier before the loop writes over it.** The update below
            // zeroes `a[i][col]` on its first step, and Bareiss's formula needs the
            // original value in every step after that.
            let multiplier = a[i][col];
            for j in col..cols {
                // Bareiss: a[i][j] = (a[i][j]*pivot - multiplier*a[rank][j]) / previous
                a[i][j] = (a[i][j] * a[rank][col] - multiplier * a[rank][j]) / previous_pivot;
            }
        }
        previous_pivot = a[rank][col];
        rank += 1;
    }

    Ok(rank)
}

/// Solve `a x = b` for a square `a`, by LU with partial pivoting.
///
/// The same shape JAMA's `LUDecomposition.solve` takes, and the reason a port of
/// `Matrix.solve` is not needed either: the matrix it is handed is the independent-
/// column submatrix, which is square by construction and integer-valued, so the only
/// floating point is the right-hand side - `-R T ln K` - carried through the
/// elimination. There is no tolerance decision in it at all.
///
/// # Errors
/// [`AzothError::InvalidInput`] if `a` is not square or its shape disagrees with `b`, or
/// if a pivot is zero - which cannot happen for a matrix whose rank was just tested, and
/// is refused rather than divided by.
pub fn solve_lu(a: &[Vec<f64>], b: &[f64]) -> Result<Vec<f64>> {
    let n = a.len();
    if b.len() != n {
        return Err(AzothError::InvalidInput {
            field: "reaction_matrix".to_string(),
            reason: format!("a is {n}x? and b has {} entries", b.len()),
        });
    }
    if n == 0 {
        return Ok(Vec::new());
    }

    let mut lu: Vec<Vec<f64>> = Vec::with_capacity(n);
    for row in a {
        if row.len() != n {
            return Err(AzothError::InvalidInput {
                field: "reaction_matrix".to_string(),
                reason: format!(
                    "a must be square: a row has {} entries against {n}",
                    row.len()
                ),
            });
        }
        lu.push(row.clone());
    }
    let mut x = b.to_vec();

    for k in 0..n {
        // Partial pivoting on the largest magnitude: the largest available pivot is the
        // one that keeps the multipliers below one.
        let mut pivot = k;
        let mut best = lu[k][k].abs();
        for (i, row) in lu.iter().enumerate().skip(k + 1) {
            if row[k].abs() > best {
                best = row[k].abs();
                pivot = i;
            }
        }
        if best == 0.0 {
            // Unreachable in this crate: the caller solves with the submatrix its rank
            // test has just accepted as full rank. Refused rather than divided by.
            return Err(AzothError::InvalidInput {
                field: "reaction_matrix".to_string(),
                reason: format!(
                    "the matrix is singular at column {k}, and a full-rank matrix is the \
                     only thing the rank test hands to this solve"
                ),
            });
        }
        if pivot != k {
            lu.swap(k, pivot);
            x.swap(k, pivot);
        }

        // `split_at_mut` rather than indexing twice: the pivot row is read while the
        // rows below it are written, and the split is what lets the borrow checker see
        // that those are disjoint.
        let pivot_value = x[k];
        let (above, below) = lu.split_at_mut(k + 1);
        let pivot_row = &above[k];
        for (offset, row) in below.iter_mut().enumerate() {
            let factor = row[k] / pivot_row[k];
            row[k] = factor;
            for (value, reference) in row[(k + 1)..].iter_mut().zip(&pivot_row[(k + 1)..]) {
                *value -= factor * reference;
            }
            x[k + 1 + offset] -= factor * pivot_value;
        }
    }

    // Back substitution.
    let mut solution = vec![0.0; n];
    for k in (0..n).rev() {
        let mut sum = x[k];
        for j in (k + 1)..n {
            sum -= lu[k][j] * solution[j];
        }
        solution[k] = sum / lu[k][k];
    }
    Ok(solution)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rank_counts_independent_rows() {
        let full = vec![
            vec![-1.0, -2.0, 0.0, 1.0],
            vec![0.0, -2.0, 1.0, 1.0],
            vec![0.0, -1.0, -1.0, 1.0],
        ];
        assert_eq!(rank_of_integer_matrix(&full).unwrap(), 3);

        // The third row is the sum of the first two, so it contributes nothing.
        let dependent = vec![
            vec![-1.0, -2.0, 0.0, 1.0],
            vec![0.0, -2.0, 1.0, 1.0],
            vec![-1.0, -4.0, 1.0, 2.0],
        ];
        assert_eq!(rank_of_integer_matrix(&dependent).unwrap(), 2);
    }

    #[test]
    fn rank_is_invariant_under_column_permutation() {
        // Which is why NeqSim's `HashSet` component order cannot change the answer to
        // `removeDependentReactions`, whose only use of the matrix is this call.
        let a = vec![vec![1.0, 2.0, 3.0], vec![4.0, 5.0, 6.0]];
        let b = vec![vec![3.0, 1.0, 2.0], vec![6.0, 4.0, 5.0]];
        assert_eq!(
            rank_of_integer_matrix(&a).unwrap(),
            rank_of_integer_matrix(&b).unwrap()
        );
    }

    #[test]
    fn a_non_integral_matrix_is_refused_rather_than_rounded() {
        let matrix = vec![vec![1.0, 0.5], vec![0.0, 1.0]];
        let error = rank_of_integer_matrix(&matrix).expect_err("0.5 is not a coefficient");
        assert!(
            matches!(error, AzothError::InvalidInput { .. }),
            "{error:?}"
        );
    }

    #[test]
    fn the_solve_recovers_a_known_solution() {
        // 2x + y = 5, x + 3y = 10  ->  x = 1, y = 3
        let a = vec![vec![2.0, 1.0], vec![1.0, 3.0]];
        let x = solve_lu(&a, &[5.0, 10.0]).unwrap();
        assert!((x[0] - 1.0).abs() < 1e-14);
        assert!((x[1] - 3.0).abs() < 1e-14);
    }

    #[test]
    fn a_singular_system_is_refused() {
        let a = vec![vec![1.0, 2.0], vec![2.0, 4.0]];
        let error = solve_lu(&a, &[1.0, 2.0]).expect_err("the rows are dependent");
        assert!(
            matches!(error, AzothError::InvalidInput { .. }),
            "{error:?}"
        );
    }
}
