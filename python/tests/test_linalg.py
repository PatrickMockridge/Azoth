"""The reactions tier's linear algebra: the LU solve and its pseudo-inverse fallback.

The Python twin of the Rust unit tests in ``crates/azoth-reactions/src/linalg.rs``. Both
assert the same analytic values, because the two implementations are mirrors and the
values - not either implementation - are the oracle.
"""

from __future__ import annotations

import pytest

from azoth.core.errors import InvalidInputError
from azoth.reactions.reference import _linalg
from azoth.reactions.reference.chemical_equilibrium import _solve_bordered


def test_a_singular_system_is_refused() -> None:
    """The branch the rank-tested callers cannot reach and the bordered solve can."""
    with pytest.raises(InvalidInputError, match="singular at column 1"):
        _linalg.solve_lu([[1.0, 2.0], [2.0, 4.0]], [1.0, 2.0])


def test_the_pseudo_inverse_answers_the_singular_system() -> None:
    """``A = u v^T`` with ``u = v = [1, 2]``, so ``A+ = [[1, 2], [2, 4]] / 25``."""
    x = _linalg.solve_least_squares([[1.0, 2.0], [2.0, 4.0]], [1.0, 2.0])
    assert x == pytest.approx([0.2, 0.4], abs=1e-14)
    # The solution solves the system and carries no null-space component.
    assert x[0] * 1.0 + x[1] * 2.0 == pytest.approx(1.0, abs=1e-14)
    assert x[0] * 2.0 - x[1] == pytest.approx(0.0, abs=1e-14)


def test_the_pseudo_inverse_keeps_the_minimum_norm_solution() -> None:
    """Rows 1 and 2 sum to row 3, and ``b`` does too; ``t = 5/3`` is the minimum norm."""
    a = [
        [1.0, 0.0, 1.0],
        [0.0, 1.0, 1.0],
        [1.0, 1.0, 2.0],
    ]
    x = _linalg.solve_least_squares(a, [2.0, 3.0, 5.0])
    assert x == pytest.approx([1.0 / 3.0, 4.0 / 3.0, 5.0 / 3.0], abs=1e-14)


def test_the_bordered_solve_falls_back_instead_of_refusing() -> None:
    """The state the reactive brine reaches: the LU arm fails and the answer must exist."""
    solved = _solve_bordered([[1.0, 2.0], [2.0, 4.0]], [1.0, 2.0])
    assert solved == pytest.approx([0.2, 0.4], abs=1e-14)


def test_a_non_finite_entry_is_refused_rather_than_propagated() -> None:
    with pytest.raises(InvalidInputError, match="not a finite number"):
        _linalg.solve_least_squares([[1.0, float("nan")], [0.0, 1.0]], [1.0, 1.0])
