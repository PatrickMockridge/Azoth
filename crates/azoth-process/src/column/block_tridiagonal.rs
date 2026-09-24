//! A block-tridiagonal solve, for a staged column's linear system.
//!
//! NeqSim's `BlockTridiagonalMatrix` (100 lines) is storage - one dense lower, diagonal and
//! upper block per row - and `TDMAsolve` (51) is the *scalar* Thomas algorithm; the solver that
//! consumes the blocks is [`NaphtaliSandholmSolver`](super::naphtali_sandholm)'s own
//! `solveBlockTridiagonal`, which is ported here. azoth carried no banded or block-banded
//! linear algebra before this: `azoth-reactions`' `linalg` is LU, rank and null-space work.
//!
//! **The arithmetic is the class's, not a better one.** The forward sweep *inverts* each
//! reduced diagonal block and multiplies, where an LU of the block would be cheaper and
//! better conditioned; the port keeps the inverse because the guard against a singular block
//! (`max |pivot| < 1e-30` inside the Gauss-Jordan) is the class's own test, and an LU would
//! refuse a different set of systems.

// The equations and the elimination are the class's own index arithmetic: the MESH residual is
// written over tray and component indices and the linear algebra over rows and columns, so the
// loops index rather than iterate. Renaming them to iterators would rewrite the port rather
// than clean it up - the precedent `azoth-eos`'s `gerg2008` sets.
#![allow(clippy::needless_range_loop)]

use azoth_core::{AzothError, Result};

/// One row's blocks: the coupling to the row below, the row's own block, and the coupling to
/// the row above.
type Block = Vec<Vec<f64>>;

/// A block-tridiagonal system, held as one dense block per row and column band.
///
/// The system is `N` block-rows of `m` equations and `m` variables each, and the row `j`'s
/// equations may involve the variables of rows `j - 1`, `j` and `j + 1` - which is what a
/// staged column's MESH Jacobian is, one block-row per tray.
#[derive(Debug, Clone)]
pub struct BlockTridiagonal {
    lower: Vec<Block>,
    diagonal: Vec<Block>,
    upper: Vec<Block>,
    m: usize,
}

impl BlockTridiagonal {
    /// Storage for `block_count` rows of `block_size` equations each.
    ///
    /// # Errors
    /// [`AzothError::InvalidInput`] if either count is zero.
    pub fn new(block_count: usize, block_size: usize) -> Result<Self> {
        if block_count == 0 || block_size == 0 {
            return Err(AzothError::invalid_input(
                "block_tridiagonal",
                format!(
                    "a block-tridiagonal system needs at least one block of at least one \
                     equation, and {block_count} blocks of {block_size} is not one"
                ),
            ));
        }
        let band = || vec![vec![0.0; block_size]; block_size];
        Ok(Self {
            lower: vec![band(); block_count],
            diagonal: vec![band(); block_count],
            upper: vec![band(); block_count],
            m: block_size,
        })
    }

    /// The number of rows.
    #[must_use]
    pub fn block_count(&self) -> usize {
        self.diagonal.len()
    }

    /// The number of equations per row.
    #[must_use]
    pub fn block_size(&self) -> usize {
        self.m
    }

    /// The mutable block of the band that couples row `row` to column row `column`.
    ///
    /// # Errors
    /// [`AzothError::InvalidInput`] if the coupling is outside the tridiagonal band, or if
    /// either index is past the end - the class throws for the first and would throw an
    /// array bounds exception for the second.
    pub fn block_mut(&mut self, row: usize, column: usize) -> Result<&mut Block> {
        let rows = self.diagonal.len();
        let band = match row.checked_sub(1) {
            Some(below) if below == column => &mut self.lower,
            _ if row == column => &mut self.diagonal,
            _ if row + 1 == column => &mut self.upper,
            _ => {
                return Err(AzothError::invalid_input(
                    "block_tridiagonal",
                    format!(
                        "the coupling of row {row} to column {column} is outside the \
                         tridiagonal band"
                    ),
                ));
            }
        };
        band.get_mut(row).ok_or_else(|| {
            AzothError::invalid_input(
                "block_tridiagonal",
                format!("row {row} is past the end of {rows} rows"),
            )
        })
    }

    /// The whole matrix, one dense row per equation.
    ///
    /// `BlockTridiagonalMatrix.toDense`: the bands copied into place, which is what a reader
    /// checking a solution by hand - or a solve that has to fall back to a dense elimination -
    /// needs.
    #[must_use]
    pub fn to_dense(&self) -> Vec<Vec<f64>> {
        let n = self.diagonal.len();
        let size = n * self.m;
        let mut dense = vec![vec![0.0; size]; size];
        for block in 0..n {
            copy(&mut dense, &self.diagonal[block], block, block);
            if block > 0 {
                copy(&mut dense, &self.lower[block], block, block - 1);
            }
            if block + 1 < n {
                copy(&mut dense, &self.upper[block], block, block + 1);
            }
        }
        dense
    }

    /// Solve `A x = rhs` by block elimination, or say the system is singular.
    ///
    /// # Errors
    /// [`AzothError::InvalidInput`] if the right-hand side is not one block per row.
    pub fn solve(&self, rhs: &[Vec<f64>]) -> Result<Option<Vec<Vec<f64>>>> {
        let n = self.block_count();
        if rhs.len() != n {
            return Err(AzothError::invalid_input(
                "rhs",
                format!(
                    "a system of {n} block-rows needs {n} right-hand blocks and got {}",
                    rhs.len()
                ),
            ));
        }

        let mut prime = self.diagonal.clone();
        let mut prime_rhs = rhs.to_vec();

        for j in 1..n {
            let Some(inverse) = inverse(&prime[j - 1]) else {
                return Ok(None);
            };
            let factor = multiply(&self.lower[j], &inverse);
            let correction = multiply(&factor, &self.upper[j - 1]);
            for r in 0..self.m {
                for c in 0..self.m {
                    prime[j][r][c] = self.diagonal[j][r][c] - correction[r][c];
                }
            }
            let applied = multiply_vector(&factor, &prime_rhs[j - 1]);
            for r in 0..self.m {
                prime_rhs[j][r] -= applied[r];
            }
        }

        let mut solution = vec![vec![0.0; self.m]; n];
        let Some(last) = inverse(&prime[n - 1]) else {
            return Ok(None);
        };
        solution[n - 1] = multiply_vector(&last, &prime_rhs[n - 1]);
        for j in (0..n - 1).rev() {
            let carried = multiply_vector(&self.upper[j], &solution[j + 1]);
            let adjusted: Vec<f64> = (0..self.m).map(|r| prime_rhs[j][r] - carried[r]).collect();
            let Some(inverted) = inverse(&prime[j]) else {
                return Ok(None);
            };
            solution[j] = multiply_vector(&inverted, &adjusted);
        }
        Ok(Some(solution))
    }
}

/// Copy one block into its place in a dense matrix.
fn copy(dense: &mut [Vec<f64>], block: &Block, row_block: usize, column_block: usize) {
    let m = block.len();
    for row in 0..m {
        let target = row_block * m + row;
        let start = column_block * m;
        dense[target][start..start + m].copy_from_slice(&block[row][..m]);
    }
}

/// The inverse of a small dense matrix by Gauss-Jordan, or `None` where it is singular.
///
/// The singularity test is the class's: the largest absolute entry in the column being
/// eliminated must exceed `1e-30`, and a column that does not is a system this returns `None`
/// for rather than a division by a pivot near zero.
fn inverse(a: &Block) -> Option<Block> {
    let m = a.len();
    let mut augmented = vec![vec![0.0; 2 * m]; m];
    for i in 0..m {
        augmented[i][..m].copy_from_slice(&a[i][..m]);
        augmented[i][m + i] = 1.0;
    }
    for column in 0..m {
        let mut pivot_row = column;
        let mut largest = augmented[column][column].abs();
        for row in column + 1..m {
            if augmented[row][column].abs() > largest {
                largest = augmented[row][column].abs();
                pivot_row = row;
            }
        }
        if largest < 1e-30 {
            return None;
        }
        augmented.swap(column, pivot_row);
        let pivot = augmented[column][column];
        for k in 0..2 * m {
            augmented[column][k] /= pivot;
        }
        for row in 0..m {
            if row == column {
                continue;
            }
            let factor = augmented[row][column];
            for k in 0..2 * m {
                augmented[row][k] -= factor * augmented[column][k];
            }
        }
    }
    Some((0..m).map(|i| augmented[i][m..].to_vec()).collect())
}

/// The product of two square blocks.
fn multiply(a: &Block, b: &Block) -> Block {
    let m = a.len();
    let mut result = vec![vec![0.0; m]; m];
    for (i, row) in result.iter_mut().enumerate() {
        for (j, value) in row.iter_mut().enumerate() {
            let mut sum = 0.0;
            for k in 0..m {
                sum += a[i][k] * b[k][j];
            }
            *value = sum;
        }
    }
    result
}

/// The product of a square block and a vector.
fn multiply_vector(a: &Block, v: &[f64]) -> Vec<f64> {
    let m = a.len();
    let mut result = vec![0.0; m];
    for (i, value) in result.iter_mut().enumerate() {
        let mut sum = 0.0;
        for k in 0..m {
            sum += a[i][k] * v[k];
        }
        *value = sum;
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The block-tridiagonal route and a dense elimination agree on a system with a known
    /// answer: a constant-coefficient tridiagonal built so that `x = 1` solves it.
    #[test]
    fn a_block_system_agrees_with_its_own_residual() {
        let (n, m) = (4_usize, 3_usize);
        let mut system = BlockTridiagonal::new(n, m).expect("a system");
        for j in 0..n {
            let diagonal = system.block_mut(j, j).expect("the diagonal band");
            for (r, row) in diagonal.iter_mut().enumerate() {
                row[r] = 4.0;
                if r > 0 {
                    row[r - 1] = -1.0;
                }
            }
            if j + 1 < n {
                let upper = system.block_mut(j, j + 1).expect("the upper band");
                for row in upper.iter_mut() {
                    row[0] = -1.0;
                }
            }
            if j > 0 {
                let lower = system.block_mut(j, j - 1).expect("the lower band");
                for row in lower.iter_mut() {
                    row[m - 1] = -1.0;
                }
            }
        }
        // rhs = A * 1
        let mut rhs = vec![vec![0.0; m]; n];
        for j in 0..n {
            for r in 0..m {
                rhs[j][r] = system.block_mut(j, j).expect("the diagonal")[r]
                    .iter()
                    .sum::<f64>();
                if j + 1 < n {
                    rhs[j][r] += system.block_mut(j, j + 1).expect("the upper band")[r]
                        .iter()
                        .sum::<f64>();
                }
                if j > 0 {
                    rhs[j][r] += system.block_mut(j, j - 1).expect("the lower band")[r]
                        .iter()
                        .sum::<f64>();
                }
            }
        }
        let solution = system
            .solve(&rhs)
            .expect("the system is well formed")
            .expect("the system is not singular");
        for block in &solution {
            for value in block {
                assert!((value - 1.0).abs() < 1.0e-12, "{value}");
            }
        }
    }

    /// A coupling outside the band is refused rather than silently written into a neighbour.
    #[test]
    fn a_coupling_outside_the_band_is_refused() {
        let mut system = BlockTridiagonal::new(3, 2).expect("a system");
        assert!(system.block_mut(0, 2).is_err());
        assert!(system.block_mut(2, 0).is_err());
        assert!(system.block_mut(1, 1).is_ok());
    }
}
