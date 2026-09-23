//! The linear-program initial estimate `LinearProgrammingChemicalEquilibrium` seeds the
//! reactive solve with.
//!
//! `ChemicalReactionOperations.solveChemEq` hands its Newton solve a starting composition
//! from a linear program: **minimise `sum(mu_i n_i / R T)` subject to `A n = b` and
//! `n >= 0`**, where `A` is the element matrix with the electroneutrality row last. The
//! answer is written back into the phase through `updateMoles` before the solve starts, so
//! it is a step of the operation and not a suggestion.
//!
//! # Why this is enumeration and not a simplex
//!
//! Commons Math's `SimplexSolver` is epsilon-tolerant Dantzig-with-Bland pivoting, so where
//! an optimal face is degenerate it returns *some* optimal vertex. The captured states have
//! 4 rows against 6 columns and 5 against 9, so the basic feasible solutions are `C(6,4) =
//! 15` and `C(9,5) = 126` - small enough to enumerate exactly. **What this module pins is
//! the optimum**: its value and the set of components it puts moles on are properties of
//! the program. Which vertex of a degenerate optimal face is returned is a property of the
//! pivoting, and where the two differ the port's choice is the lowest basis mask, which is
//! deterministic and is the divergence the spec records.
//!
//! # What is not reproduced
//!
//! NeqSim's objective vector and every constraint carry a leading entry for a **slack
//! variable that appears nowhere else** - coefficient `0.0` in the objective and in every
//! row (`LinearProgrammingChemicalEquilibrium.java:307-320`). A variable with no
//! coefficient cannot affect feasibility or the optimum, so it has no column here. Its
//! `inertMoles` parameter is read by nothing in the method either.

use crate::linalg::solve_lu;

/// The floor every seed entry is lifted to, from NeqSim's `MIN_MOLES`
/// (`LinearProgrammingChemicalEquilibrium.java:41`).
pub const MIN_MOLES: f64 = 1e-60;

/// The most columns a basis sweep will enumerate, so a caller cannot ask for `2**40` subset
/// tests by passing a large fluid.
const MAX_SWEPT_COLUMNS: usize = 20;

/// The LP's answer: the feasible vertex minimising `sum(reduced_potentials[i] * n[i])`, or
/// `None` where the program has no feasible solution.
///
/// `None` is NeqSim's `null` and is **a state and not an error**: `solveChemEq` reads it as
/// "keep the composition the phase already has". Every returned entry is floored at
/// [`MIN_MOLES`], which is what NeqSim does to the solver's point before returning it.
#[must_use]
pub fn initial_estimate(
    a_matrix: &[Vec<f64>],
    b: &[f64],
    reduced_potentials: &[f64],
) -> Option<Vec<f64>> {
    let rows = b.len();
    let columns = reduced_potentials.len();
    if rows == 0 || columns == 0 || rows > columns || a_matrix.len() < rows {
        return None;
    }
    if columns > MAX_SWEPT_COLUMNS {
        return None;
    }

    let mut best: Option<(f64, u32, Vec<f64>)> = None;
    for mask in 1u32..(1u32 << columns) {
        if mask.count_ones() as usize != rows {
            continue;
        }
        let basis: Vec<usize> = (0..columns).filter(|c| mask & (1 << c) != 0).collect();
        let mut square = vec![vec![0.0; rows]; rows];
        for (row, values) in a_matrix.iter().take(rows).enumerate() {
            for (column, &source) in basis.iter().enumerate() {
                square[row][column] = values[source];
            }
        }
        let Ok(solved) = solve_lu(&square, b) else {
            continue;
        };
        if solved.iter().any(|value| *value < 0.0) {
            continue;
        }
        let mut moles = vec![0.0; columns];
        for (position, &column) in basis.iter().enumerate() {
            moles[column] = solved[position];
        }
        let objective: f64 = moles
            .iter()
            .zip(reduced_potentials)
            .map(|(moles, potential)| moles * potential)
            .sum();
        // Ties are broken by the lowest basis mask, which the ascending sweep gives for
        // free. A degenerate optimal face is the case that needs it.
        let better = match &best {
            None => true,
            Some((held, best_mask, _)) => {
                objective < held - OBJECTIVE_TIE
                    || (objective - held).abs() <= OBJECTIVE_TIE && mask < *best_mask
            }
        };
        if better {
            best = Some((objective, mask, moles));
        }
    }

    best.map(|(_, _, moles)| {
        moles
            .into_iter()
            .map(|value| value.max(MIN_MOLES))
            .collect()
    })
}

/// Two objectives within this of each other are the same optimum: the comparison is between
/// sums of products of a potential and a mole number, so it carries the arithmetic's own
/// rounding.
const OBJECTIVE_TIE: f64 = 1e-12;

#[cfg(test)]
mod tests {
    use super::*;

    /// The element matrix of both captured fluids, with the electroneutrality row last.
    fn captured_matrix() -> Vec<Vec<f64>> {
        vec![
            vec![1.0, 0.0, 0.0, 0.0, 1.0, 1.0],
            vec![0.0, 2.0, 1.0, 3.0, 1.0, 0.0],
            vec![2.0, 1.0, 1.0, 1.0, 3.0, 3.0],
            vec![0.0, 0.0, -1.0, 1.0, -1.0, -2.0],
        ]
    }

    /// The reduced reference potentials, from the capture. The same for both fluids: they
    /// are a function of the temperature and the standard state, not of the composition.
    fn captured_objective() -> [f64; 6] {
        [
            -54.799_675_618_630_715,
            34.716_683_013_169_76,
            29.168_005_792_406_905,
            80.530_720_467_865_22,
            -51.263_339_652_447_62,
            -62.360_694_093_973_336,
        ]
    }

    /// **The program NeqSim states has no solution, and that is the finding.**
    ///
    /// On the unflashed CO2-water fluid Commons Math returns a vertex -
    /// `[0.0100000002, 10.00000000025, ~0, ~0, 3.92e-16, ~0]` - and that point does not
    /// satisfy `A n = b`: it misses the oxygen row by `1.5e-10` and the charge row by
    /// `3.92e-16`. The exact program has no non-negative solution at all, which an
    /// exhaustive search over every support of size up to six confirms, so the answer
    /// Commons Math returns is *epsilon-feasible* - its own `DEFAULT_EPSILON` is `1e-6` -
    /// and not a solution of its own constraints.
    ///
    /// The cause is visible in the tree: `calcBVector` builds `b` from a cached `nVector`
    /// that only a solve getting past the phase search refreshes, so `b` is not `A n` of
    /// the composition the phase holds. The port answers the program as stated, which is
    /// `None` here; the divergence is recorded in the spec rather than smoothed over.
    #[test]
    fn the_captured_programs_have_no_exact_solution() {
        let a = captured_matrix();
        let objective = captured_objective();

        // The unflashed fluid: NeqSim seeds, and there is no exact vertex to seed with.
        let unflashed = [
            0.010_000_000_2,
            20.000_000_000_499_995,
            10.020_000_000_8,
            0.0,
        ];
        assert_eq!(initial_estimate(&a, &unflashed, &objective), None);

        // And the point NeqSim returned violates the rows it was solved against.
        let reported = [
            0.010_000_000_199_999_607,
            10.000_000_000_249_997,
            MIN_MOLES,
            MIN_MOLES,
            3.920_475_055_707_584e-16,
            MIN_MOLES,
        ];
        let worst = (0..a.len())
            .map(|row| {
                (0..reported.len())
                    .map(|column| a[row][column] * reported[column])
                    .sum::<f64>()
                    - unflashed[row]
            })
            .fold(0.0_f64, |worst, residual| worst.max(residual.abs()));
        assert!(
            worst > 1e-11,
            "NeqSim's own seed is {} from satisfying the program, so this case no longer \
             records a divergence",
            worst
        );
        assert!(worst < 1e-6, "{worst} is above the epsilon it passes on");

        // The flashed aqueous phase: NeqSim throws `NoFeasibleSolutionException` here, and
        // so does this - the one place the two agree.
        let flashed = [
            1.405_182_302_147_056_8e-4,
            19.999_593_938_541_302,
            10.000_066_636_279_326,
            0.0,
        ];
        assert_eq!(initial_estimate(&a, &flashed, &objective), None);
    }

    /// The enumeration itself, on a program with an answer that can be checked by hand:
    /// `n_0 + n_1 = 2`, `n_0 - n_1 = 0`, minimising `n_0 + 2 n_1`. The only feasible point
    /// is `(1, 1)`, so there is nothing for the tie-break to choose between.
    #[test]
    fn a_program_with_one_feasible_point_returns_it() {
        let a = vec![vec![1.0, 1.0], vec![1.0, -1.0]];
        let seed = initial_estimate(&a, &[2.0, 0.0], &[1.0, 2.0]).expect("feasible");
        assert!((seed[0] - 1.0).abs() < 1e-15, "{}", seed[0]);
        assert!((seed[1] - 1.0).abs() < 1e-15, "{}", seed[1]);
    }

    /// A degenerate optimal face: every point on `n_0 + n_1 = 1` with the objective
    /// `n_0 + n_1` is optimal, so both vertices cost the same and the answer is decided by
    /// the tie-break - the lowest basis mask, which is `n_0 = 1`.
    #[test]
    fn a_degenerate_optimum_is_broken_towards_the_lowest_basis() {
        let a = vec![vec![1.0, 1.0]];
        let seed = initial_estimate(&a, &[1.0], &[1.0, 1.0]).expect("feasible");
        assert!((seed[0] - 1.0).abs() < 1e-15, "{}", seed[0]);
        assert!(seed[1] <= MIN_MOLES, "{}", seed[1]);
    }

    #[test]
    fn a_shape_the_sweep_cannot_cover_is_answered_with_none() {
        let a = vec![vec![1.0, 1.0]];
        assert_eq!(initial_estimate(&[], &[], &[]), None);
        // More rows than columns has no basic solution, and neither has a program whose
        // only basis is singular.
        assert_eq!(initial_estimate(&a, &[1.0, 1.0], &[1.0]), None);
        assert_eq!(
            initial_estimate(&[vec![1.0, 2.0], vec![2.0, 4.0]], &[1.0, 2.0], &[1.0, 1.0]),
            None,
            "a singular basis is skipped rather than solved"
        );
    }

    /// **The component order does not move the optimum, and NeqSim's own order is a sort.**
    ///
    /// `LinearProgrammingChemicalEquilibrium`'s fallback branch sorts its component array
    /// with `ReferencePotComparator` before building the matrix and the potentials - a
    /// comparator that returns `1` or `0` and never `-1`, so it is not transitive and the
    /// permutation it produces is arbitrary. Whatever permutation it produces, the sort
    /// moves the potentials array and the matrix's columns *together* (both are read off
    /// the same sorted array), so the program is the same one written in another column
    /// order.
    ///
    /// This sweeps a permutation of both kinds of program and asserts what that implies:
    /// **the objective value at the optimum is invariant**, and where the optimum is
    /// unique the vertex permutes with the columns. Where it is *not* unique, the two runs
    /// land on different vertices of the same optimal face - which is what makes the
    /// tie-break a convention rather than an answer, and why this module pins the objective
    /// and the mask and records the choice.
    #[test]
    fn a_column_permutation_moves_the_vertex_and_not_the_optimum() {
        // A unique optimum: basis {0, 2} costs 3, and every other basis costs more or is
        // singular.
        let a = vec![vec![1.0, 0.0, 0.0], vec![0.0, 2.0, 1.0]];
        let b = [1.0, 2.0];
        let objective = [1.0, 5.0, 1.0];
        let seed = initial_estimate(&a, &b, &objective).expect("feasible");
        assert!(
            seed[0] > 0.5 && seed[1] <= MIN_MOLES && seed[2] > 1.5,
            "{seed:?}"
        );

        let swap = [2usize, 1, 0];
        let moved: Vec<Vec<f64>> = a
            .iter()
            .map(|row| swap.iter().map(|&column| row[column]).collect())
            .collect();
        let moved_objective: Vec<f64> = swap.iter().map(|&column| objective[column]).collect();
        let moved_seed =
            initial_estimate(&moved, &b, &moved_objective).expect("the permutation is feasible");

        let cost = |seed: &[f64], objective: &[f64]| -> f64 {
            seed.iter().zip(objective).map(|(n, c)| n * c).sum()
        };
        assert!(
            (cost(&seed, &objective) - cost(&moved_seed, &moved_objective)).abs() < 1e-12,
            "the optimum's value is a property of the program and not of the column order"
        );
        for (i, &column) in swap.iter().enumerate() {
            assert!(
                (moved_seed[i] - seed[column]).abs() < 1e-15,
                "column {i} of the permuted program is column {column} of the original"
            );
        }

        // A degenerate face: `n_0 + n_1 + n_2 = 1` with every coefficient costing one, so
        // every vertex is optimal and the tie-break decides which one comes back.
        let a = vec![vec![1.0, 1.0, 1.0]];
        let b = [1.0];
        let objective = [1.0, 1.0, 1.0];
        let seed = initial_estimate(&a, &b, &objective).expect("feasible");
        assert!(seed[0] > 0.5, "the lowest basis mask wins: {seed:?}");

        let swap = [2usize, 0, 1];
        let moved: Vec<Vec<f64>> = a
            .iter()
            .map(|row| swap.iter().map(|&column| row[column]).collect())
            .collect();
        let moved_objective: Vec<f64> = swap.iter().map(|&column| objective[column]).collect();
        let moved_seed =
            initial_estimate(&moved, &b, &moved_objective).expect("the permutation is feasible");
        assert!(
            (cost(&seed, &objective) - cost(&moved_seed, &moved_objective)).abs() < 1e-12,
            "still an optimum, whichever vertex the tie-break lands on"
        );
    }
}
