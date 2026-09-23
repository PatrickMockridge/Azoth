"""The Pulay accelerator, the Python twin of ``crates/azoth-reactions/src/diis.rs``.

Pulay's direct inversion in the iterative subspace: a rolling history of
iterate-residual pairs, and the coefficients ``c_i`` with ``sum(c_i) = 1`` that
minimise ``||sum(c_i r_i)||``, read off the augmented system

.. code-block:: text

    [ B  -1 ] [c ]   [ 0]
    [-1   0 ] [mu] = [-1]      B_ij = r_i . r_j

and applied as ``x_new = sum(c_i x_i)``. NeqSim drives it on the Lagrange
multipliers, with the element-balance residual vector as the error.

**It is load-bearing, and the capture says so**: the reactive flash's
forced-one-phase 300 K water-gas shift accepts 19 extrapolated steps in 35
iterations, and the 600 K state accepts none. A reference implementation that
dropped DIIS would follow a different trajectory on one branch and the same one
on the other, so it is ported rather than approximated.

The one divergence from the class is stated in the Rust twin and repeated here:
``addEntry`` copies with ``System.arraycopy`` and throws on a short array, where
this refuses with an ``InvalidInputError`` - a wrapped slot half-filled from a
short vector would otherwise keep a stale tail.
"""

from __future__ import annotations

from azoth.core.errors import InvalidInputError

#: The pivot floor ``extrapolate`` refuses a Pulay system at, from the class's ``1e-30``.
PIVOT_FLOOR = 1.0e-30


class DiisAccelerator:
    """The Pulay accelerator, with the class's circular buffer.

    ``max_history`` is the class's ``DIIS_DEPTH`` of 6 in the solve, and the vector
    length is its ``ne``.
    """

    def __init__(self, vector_length: int, max_history: int) -> None:
        self.vector_length = vector_length
        self.max_history = max_history
        self._iterate_history = [[0.0] * vector_length for _ in range(max_history)]
        self._residual_history = [[0.0] * vector_length for _ in range(max_history)]
        self._count = 0
        self._next_slot = 0

    def add_entry(self, iterate: list[float], residual: list[float]) -> None:
        """Store a pair, overwriting the oldest slot once the buffer is full.

        Raises:
            InvalidInputError: if either vector is not ``vector_length`` long - the class
                throws where this refuses.
        """
        if len(iterate) != self.vector_length or len(residual) != self.vector_length:
            raise InvalidInputError(
                "iterate",
                f"{len(iterate)} and {len(residual)} entries against a vector length of "
                f"{self.vector_length}",
            )
        self._iterate_history[self._next_slot] = list(iterate)
        self._residual_history[self._next_slot] = list(residual)
        self._next_slot = (self._next_slot + 1) % self.max_history
        if self._count < self.max_history:
            self._count += 1

    def can_extrapolate(self) -> bool:
        """``canExtrapolate``: two pairs are the minimum a combination needs."""
        return self._count >= 2

    def count(self) -> int:
        """``getCount``: how many pairs are stored, up to the history's length."""
        return self._count

    def reset(self) -> None:
        """``reset``: discard the history."""
        self._count = 0
        self._next_slot = 0

    def _buffer_index(self, logical: int) -> int:
        """The oldest stored pair is logical index zero, which is ``nextSlot`` once wrapped."""
        if self._count < self.max_history:
            return logical
        return (self._next_slot + logical) % self.max_history

    def extrapolate(self) -> list[float] | None:
        """The combination, or ``None`` where the class returns ``null``.

        That is fewer than two pairs stored, or a Pulay system whose pivot falls under
        :data:`PIVOT_FLOOR`. The elimination and the back-substitution are the class's own,
        partial pivoting and pivot floor included, and the inner sweep of the elimination
        starts at the current column rather than the next one, as it does there.
        """
        if self._count < 2:
            return None

        m = self._count
        dim = m + 1
        aug = [[0.0] * (dim + 1) for _ in range(dim)]
        for i in range(m):
            ii = self._buffer_index(i)
            for j in range(m):
                jj = self._buffer_index(j)
                dot = 0.0
                for k in range(self.vector_length):
                    dot += self._residual_history[ii][k] * self._residual_history[jj][k]
                aug[i][j] = dot
            aug[i][m] = -1.0
            aug[m][i] = -1.0
            aug[i][dim] = 0.0
        aug[m][m] = 0.0
        aug[m][dim] = -1.0

        for column in range(dim):
            pivot_row = column
            largest = abs(aug[column][column])
            for row in range(column + 1, dim):
                if abs(aug[row][column]) > largest:
                    largest = abs(aug[row][column])
                    pivot_row = row
            if largest < PIVOT_FLOOR:
                return None
            if pivot_row != column:
                aug[column], aug[pivot_row] = aug[pivot_row], aug[column]
            for row in range(column + 1, dim):
                factor = aug[row][column] / aug[column][column]
                for entry in range(column, dim + 1):
                    aug[row][entry] -= factor * aug[column][entry]

        solution = [0.0] * dim
        for i in range(dim - 1, -1, -1):
            value = aug[i][dim]
            for j in range(i + 1, dim):
                value -= aug[i][j] * solution[j]
            if abs(aug[i][i]) < PIVOT_FLOOR:
                return None
            solution[i] = value / aug[i][i]

        result = [0.0] * self.vector_length
        for i in range(m):
            ii = self._buffer_index(i)
            for k in range(self.vector_length):
                result[k] += solution[i] * self._iterate_history[ii][k]
        return result
