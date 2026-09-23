//! The DIIS accelerator, from `flashops/reactiveflash/DIISAccelerator.java`.
//!
//! Pulay's direct inversion in the iterative subspace: a rolling history of iterate-residual
//! pairs, and the coefficients `c_i` with `sum(c_i) = 1` that minimise `||sum(c_i r_i)||`,
//! read off the augmented system
//!
//! ```text
//! [ B  -1 ] [c ]   [ 0]
//! [-1   0 ] [mu] = [-1]      B_ij = r_i . r_j
//! ```
//!
//! and applied as `x_new = sum(c_i x_i)`. NeqSim drives it on the Lagrange multipliers, with
//! the element-balance residual vector as the error, and accepts an extrapolated step only
//! when the full residual does not rise.
//!
//! **It is load-bearing in the multiphase solve, and the capture says so**: the forced-one-phase
//! 300 K water-gas shift accepts 19 extrapolated steps in 35 iterations, and the single-phase
//! 600 K state accepts none - so a port that dropped DIIS would follow a different trajectory
//! on one branch and the same one on the other.
//!
//! # Where this diverges from the class
//!
//! `addEntry` copies `vectorLength` entries with `System.arraycopy` and throws on a short
//! array; this returns [`AzothError::InvalidInput`] instead. The refusal is the same and the
//! arithmetic is not touched - a wrapped slot keeps the values it had rather than being
//! half-filled with a stale tail.

use azoth_core::{AzothError, Result};

/// The pivot floor `extrapolate` refuses a Pulay system at, from its own `1.0e-30`.
pub const DIIS_PIVOT_FLOOR: f64 = 1.0e-30;

/// The Pulay accelerator, with the class's circular buffer.
#[derive(Debug, Clone, PartialEq)]
pub struct DiisAccelerator {
    vector_length: usize,
    max_history: usize,
    iterate_history: Vec<Vec<f64>>,
    residual_history: Vec<Vec<f64>>,
    count: usize,
    next_slot: usize,
}

impl DiisAccelerator {
    /// A buffer of `max_history` slots, each a vector of `vector_length` - the class's
    /// `DIIS_DEPTH` of 6 in the solver, and its `ne` for the width.
    #[must_use]
    pub fn new(vector_length: usize, max_history: usize) -> Self {
        Self {
            vector_length,
            max_history,
            iterate_history: vec![vec![0.0; vector_length]; max_history],
            residual_history: vec![vec![0.0; vector_length]; max_history],
            count: 0,
            next_slot: 0,
        }
    }

    /// `addEntry`: store a pair, overwriting the oldest slot once the buffer is full.
    ///
    /// # Errors
    /// [`AzothError::InvalidInput`] where either vector is not [`Self::vector_length`] long -
    /// the class throws where this refuses.
    pub fn add_entry(&mut self, iterate: &[f64], residual: &[f64]) -> Result<()> {
        if iterate.len() != self.vector_length || residual.len() != self.vector_length {
            return Err(AzothError::InvalidInput {
                field: "iterate".to_string(),
                reason: format!(
                    "{} and {} entries against a vector length of {}",
                    iterate.len(),
                    residual.len(),
                    self.vector_length
                ),
            });
        }
        self.iterate_history[self.next_slot].copy_from_slice(iterate);
        self.residual_history[self.next_slot].copy_from_slice(residual);
        self.next_slot = (self.next_slot + 1) % self.max_history;
        if self.count < self.max_history {
            self.count += 1;
        }
        Ok(())
    }

    /// `canExtrapolate`: two pairs are the minimum a combination needs.
    #[must_use]
    pub fn can_extrapolate(&self) -> bool {
        self.count >= 2
    }

    /// `getCount`: how many pairs are stored, up to the history's length.
    #[must_use]
    pub fn count(&self) -> usize {
        self.count
    }

    /// `reset`: discard the history.
    pub fn reset(&mut self) {
        self.count = 0;
        self.next_slot = 0;
    }

    /// `extrapolate`: the combination, or `None` where the class returns `null` - fewer than
    /// two pairs stored, or a Pulay system whose pivot falls under [`DIIS_PIVOT_FLOOR`].
    ///
    /// The elimination and the back-substitution are the class's own, partial pivoting and
    /// pivot floor included, and the inner sweep of the elimination starts at the current
    /// column rather than the next one, as it does there.
    #[must_use]
    #[allow(clippy::needless_range_loop)] // the class's own index loops over the Pulay system
    pub fn extrapolate(&self) -> Option<Vec<f64>> {
        if self.count < 2 {
            return None;
        }

        let m = self.count;
        let dim = m + 1;
        let mut aug = vec![vec![0.0_f64; dim + 1]; dim];
        for i in 0..m {
            let ii = self.buffer_index(i);
            for j in 0..m {
                let jj = self.buffer_index(j);
                let mut dot = 0.0_f64;
                for k in 0..self.vector_length {
                    dot += self.residual_history[ii][k] * self.residual_history[jj][k];
                }
                aug[i][j] = dot;
            }
            aug[i][m] = -1.0;
            aug[m][i] = -1.0;
            aug[i][dim] = 0.0;
        }
        aug[m][m] = 0.0;
        aug[m][dim] = -1.0;

        for column in 0..dim {
            let mut pivot_row = column;
            let mut largest = aug[column][column].abs();
            for row in column + 1..dim {
                if aug[row][column].abs() > largest {
                    largest = aug[row][column].abs();
                    pivot_row = row;
                }
            }
            if largest < DIIS_PIVOT_FLOOR {
                return None;
            }
            if pivot_row != column {
                aug.swap(column, pivot_row);
            }
            for row in column + 1..dim {
                let factor = aug[row][column] / aug[column][column];
                for entry in column..=dim {
                    aug[row][entry] -= factor * aug[column][entry];
                }
            }
        }

        let mut solution = vec![0.0_f64; dim];
        for i in (0..dim).rev() {
            solution[i] = aug[i][dim];
            for j in i + 1..dim {
                solution[i] -= aug[i][j] * solution[j];
            }
            if aug[i][i].abs() < DIIS_PIVOT_FLOOR {
                return None;
            }
            solution[i] /= aug[i][i];
        }

        let mut result = vec![0.0_f64; self.vector_length];
        for i in 0..m {
            let ii = self.buffer_index(i);
            let coefficient = solution[i];
            for k in 0..self.vector_length {
                result[k] += coefficient * self.iterate_history[ii][k];
            }
        }
        Some(result)
    }

    /// `bufferIndex`: the oldest stored pair is logical index zero, which is `nextSlot` once
    /// the buffer has wrapped.
    fn buffer_index(&self, logical: usize) -> usize {
        if self.count < self.max_history {
            logical
        } else {
            (self.next_slot + logical) % self.max_history
        }
    }
}
