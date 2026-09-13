"""The iterative solver used by implicit equations.

Mirrors ``azoth_hydraulics::solver`` in Rust, statement for statement.

This module is small and deliberately not general. It exists because
``hydraulics.friction_factor_colebrook`` is implicit in ``f``, and the two
implementations have to agree numerically. Two implementations running the *same*
algorithm converge to the same floating-point value; two running
different-but-equally-valid algorithms - Newton against fixed point, say - agree
only to within their difference, which is far larger than a tight tolerance
allows.

So the scheme is fixed by the spec (kind, tolerance, iteration cap, initial guess,
convergence rule) and both languages implement exactly it. Nothing here should
grow into a general-purpose numerical library.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass
from enum import StrEnum

from azoth.core.errors import InvalidInputError, SolverNotConvergedError


class Convergence(StrEnum):
    """How successive iterates are compared to decide convergence.

    Named in the spec, because the same tolerance means different things under
    the two rules and both implementations must agree on which is meant.
    """

    #: ``|x_{k+1} - x_k| <= tolerance``
    ABSOLUTE = "absolute"
    #: ``|x_{k+1} - x_k| <= tolerance * |x_{k+1}|``
    RELATIVE = "relative"

    @classmethod
    def parse(cls, name: str) -> Convergence:
        """Parse the spec's spelling.

        Unknown names are an error rather than a silent default: guessing which
        convergence rule was meant would change the answer.
        """
        try:
            return cls(name)
        except ValueError:
            raise InvalidInputError(
                "convergence",
                f"unknown convergence rule `{name}`; expected `absolute` or `relative`",
            ) from None


@dataclass(frozen=True, slots=True)
class SolverOutcome:
    """What the iteration produced."""

    #: The final iterate. When `converged` is false this is the last value
    #: computed, and it is not a solution to anything.
    x: float
    #: Iterations performed.
    iterations: int
    #: Whether the stopping rule was met.
    converged: bool
    #: ``|x_k - x_{k-1}|`` at the final step.
    residual: float


def fixed_point(
    x0: float,
    tolerance: float,
    max_iterations: int,
    convergence: Convergence,
    f: Callable[[float], float],
) -> SolverOutcome:
    """Iterate ``x -> f(x)`` from ``x0`` until the stopping rule is met or the cap
    is reached.

    Iteration counts are deterministic: the same inputs always perform the same
    number of steps, which is what makes the Python and Rust results comparable
    rather than merely close.

    A ``NaN`` iterate is not converged - and cannot be, since every comparison
    with ``NaN`` is false - so it falls through to the iteration cap rather than
    masquerading as a solution.
    """
    x = x0
    residual = float("inf")

    for iteration in range(1, max_iterations + 1):
        next_x = f(x)
        residual = abs(next_x - x)
        if convergence is Convergence.ABSOLUTE:
            close_enough = residual <= tolerance
        else:
            close_enough = residual <= tolerance * abs(next_x)
        x = next_x
        if close_enough:
            return SolverOutcome(x=x, iterations=iteration, converged=True, residual=residual)

    return SolverOutcome(x=x, iterations=max_iterations, converged=False, residual=residual)


def require_converged(outcome: SolverOutcome, tolerance: float) -> SolverOutcome:
    """Turn a non-converged outcome into an error.

    A caller that gets a value from a run that never converged has a number that
    solves nothing, so this is an error rather than a warning - unlike an
    out-of-range input, where the number is real and merely untrustworthy.
    """
    if outcome.converged:
        return outcome
    raise SolverNotConvergedError(outcome.iterations, outcome.residual, tolerance)
