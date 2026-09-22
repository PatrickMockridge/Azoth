"""The linear-program initial estimate `LinearProgrammingChemicalEquilibrium` seeds with.

The Python twin of ``crates/azoth-reactions/src/lp_seed.rs``, written to mirror it
statement for statement: minimise ``sum(mu_i n_i / R T)`` subject to ``A n = b`` and
``n >= 0``, where ``A`` is the element matrix with the electroneutrality row last.

**Enumeration, not a simplex.** Commons Math's ``SimplexSolver`` is epsilon-tolerant
pivoting, so where an optimal face is degenerate it returns *some* optimal vertex. The
captured states have 4 rows against 6 columns and 5 against 9, so the basic feasible
solutions are ``C(6,4) = 15`` and ``C(9,5) = 126`` - small enough to enumerate exactly.
What this pins is the optimum: its value and the components it puts moles on. Which vertex
of a degenerate optimal face is returned is the pivoting's business, and this takes the
lowest basis mask, which is deterministic and is the divergence the spec records.

**NeqSim's answer on the captured fluids is not a solution of this program.** Commons Math
returns a vertex that misses the oxygen row by ``1.5e-10`` where the exact program has no
non-negative solution at all; this answers the program as stated, so it returns ``None``
there. See the spec's assumptions.
"""

from __future__ import annotations

from typing import Final

from azoth.reactions.reference import _linalg

#: The floor every seed entry is lifted to, from NeqSim's ``MIN_MOLES``.
MIN_MOLES: Final[float] = 1e-60

#: The most columns a basis sweep will enumerate, so a caller cannot ask for ``2**40``
#: subset tests by passing a large fluid.
_MAX_SWEPT_COLUMNS: Final[int] = 20

#: Two objectives within this of each other are the same optimum: the comparison is a sum
#: of products of a potential and a mole number, so it carries the arithmetic's rounding.
_OBJECTIVE_TIE: Final[float] = 1e-12


def initial_estimate(
    a_matrix: list[list[float]],
    b: list[float],
    reduced_potentials: list[float],
) -> list[float] | None:
    """The feasible vertex minimising ``sum(reduced_potentials[i] * n[i])``.

    ``None`` is NeqSim's ``null`` and **is a state and not an error**: ``solveChemEq``
    reads it as "keep the composition the phase already has". Every returned entry is
    floored at :data:`MIN_MOLES`, which is what NeqSim does to the solver's point before
    returning it.
    """
    rows = len(b)
    columns = len(reduced_potentials)
    if rows == 0 or columns == 0 or rows > columns or len(a_matrix) < rows:
        return None
    if columns > _MAX_SWEPT_COLUMNS:
        return None

    best: tuple[float, int, list[float]] | None = None
    for mask in range(1, 1 << columns):
        basis = [column for column in range(columns) if mask & (1 << column)]
        if len(basis) != rows:
            continue
        square = [[a_matrix[row][column] for column in basis] for row in range(rows)]
        try:
            solved = _linalg.solve_lu(square, list(b))
        except Exception:  # a singular basis is skipped, not a failure
            continue
        if any(value < 0.0 for value in solved):
            continue
        moles = [0.0] * columns
        for position, column in enumerate(basis):
            moles[column] = solved[position]
        objective = sum(moles[i] * reduced_potentials[i] for i in range(columns))
        # Ties go to the lowest basis mask, which the ascending sweep gives for free. A
        # degenerate optimal face is the case that needs it.
        if best is None:
            better = True
        else:
            held, best_mask, _ = best
            better = objective < held - _OBJECTIVE_TIE or (
                abs(objective - held) <= _OBJECTIVE_TIE and mask < best_mask
            )
        if better:
            best = (objective, mask, moles)

    if best is None:
        return None
    return [max(value, MIN_MOLES) for value in best[2]]
