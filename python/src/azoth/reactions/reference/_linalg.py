"""The two linear-algebra routines the reference-potential solve needs.

The Python twin of ``crates/azoth-reactions/src/linalg.rs``, and the same argument
applies: **JAMA is not ported**, because the rank it computes here is exact.

NeqSim calls ``Jama.Matrix.rank()``, a singular-value count against
``max(m, n) * s[0] * 2**-52``, and ``Jama.Matrix.solve()``, which is LU with partial
pivoting. Both rank tests run on matrices built from **stoichiometric coefficients alone**
- the matrix is one column wider than the coefficients and the extra column holds
``-R T ln K`` - and every coefficient in NeqSim's ``stoccoefdata`` is an integer. The
table's 151 rows carry only ``0``, ``1``, ``-1``, ``2`` and ``-2``.

``validation/neqsim/captures/reference_potential_probe.tsv`` records what JAMA's own test
does with three real fluids: the smallest singular value it keeps is ``0.2461``, while the
cutoff it compares against is ``4.8e-15``. Fourteen orders of magnitude, on every matrix,
and in all three cases it finds a full-rank basis and rejects nothing.

So the rank is computed exactly instead. Python's integers are arbitrary precision, so the
fraction-free elimination below cannot overflow and needs no scaling.
"""

from __future__ import annotations

from azoth.core.errors import InvalidInputError


def _as_integers(matrix: list[list[float]]) -> list[list[int]]:
    """The matrix as integers, refusing anything that is not one.

    A rounded coefficient would make this an approximation of a rank rather than the
    rank, and the caller cannot tell the difference from the answer.
    """
    converted: list[list[int]] = []
    for i, row in enumerate(matrix):
        converted_row: list[int] = []
        for j, value in enumerate(row):
            if not float(value).is_integer():
                raise InvalidInputError(
                    "reaction_matrix",
                    f"entry [{i}][{j}] is {value}, and the rank test is exact only over "
                    f"integers. NeqSim's rank is a singular-value count and would accept "
                    f"it; this refuses rather than approximating a rank",
                )
            converted_row.append(int(value))
        converted.append(converted_row)
    return converted


def rank_of_integer_matrix(matrix: list[list[float]]) -> int:
    """The exact rank of a matrix whose entries are all integers.

    Fraction-free (Bareiss) elimination: every intermediate stays an integer, so no pivot
    is ever compared against a tolerance and the answer is exact rather than numerical.

    Raises:
        InvalidInputError: if any entry is not integral.
    """
    if not matrix or not matrix[0]:
        return 0
    a = _as_integers(matrix)
    rows = len(a)
    cols = len(a[0])

    rank = 0
    previous_pivot = 1

    for col in range(cols):
        if rank == rows:
            break

        # Partial pivoting on the largest magnitude. Rank is invariant under row
        # permutation, so any pivot choice gives the same answer; this one keeps the
        # intermediates small.
        pivot_row = None
        best = 0
        for i in range(rank, rows):
            if abs(a[i][col]) > best:
                best = abs(a[i][col])
                pivot_row = i
        if pivot_row is None or best == 0:
            continue

        a[rank], a[pivot_row] = a[pivot_row], a[rank]

        for i in range(rank + 1, rows):
            # Read the multiplier before the loop writes over it: the update below zeroes
            # `a[i][col]` on its first step, and Bareiss's formula needs the original
            # value in every step after that.
            multiplier = a[i][col]
            for j in range(col, cols):
                a[i][j] = (a[i][j] * a[rank][col] - multiplier * a[rank][j]) // previous_pivot

        previous_pivot = a[rank][col]
        rank += 1

    return rank


def solve_lu(a: list[list[float]], b: list[float]) -> list[float]:
    """Solve ``a x = b`` for a square ``a``, by LU with partial pivoting.

    The same shape ``Jama.Matrix.solve`` takes for a square matrix, and the reason a port
    of it is not needed either: the matrix it is handed is the independent-column
    submatrix, which is square by construction and integer-valued, so the only floating
    point is the right-hand side - ``-R T ln K`` - carried through the elimination. There
    is no tolerance decision in it at all.

    Raises:
        InvalidInputError: if ``a`` is not square or its shape disagrees with ``b``, or if
            a pivot is zero.
    """
    n = len(a)
    if len(b) != n:
        raise InvalidInputError("reaction_matrix", f"a is {n}x? and b has {len(b)} entries")
    if n == 0:
        return []
    for row in a:
        if len(row) != n:
            raise InvalidInputError(
                "reaction_matrix", f"a must be square: a row has {len(row)} entries against {n}"
            )

    lu = [list(row) for row in a]
    x = list(b)

    for k in range(n):
        pivot = k
        best = abs(lu[k][k])
        for i in range(k + 1, n):
            if abs(lu[i][k]) > best:
                best = abs(lu[i][k])
                pivot = i
        if best == 0.0:
            # Unreachable in this library: the caller solves with the submatrix its rank
            # test has just accepted as full rank.
            raise InvalidInputError(
                "reaction_matrix",
                f"the matrix is singular at column {k}, and a full-rank matrix is the "
                f"only thing the rank test hands to this solve",
            )
        if pivot != k:
            lu[k], lu[pivot] = lu[pivot], lu[k]
            x[k], x[pivot] = x[pivot], x[k]

        for i in range(k + 1, n):
            factor = lu[i][k] / lu[k][k]
            lu[i][k] = factor
            for j in range(k + 1, n):
                lu[i][j] -= factor * lu[k][j]
            x[i] -= factor * x[k]

    solution = [0.0] * n
    for k in range(n - 1, -1, -1):
        total = x[k]
        for j in range(k + 1, n):
            total -= lu[k][j] * solution[j]
        solution[k] = total / lu[k][k]
    return solution
