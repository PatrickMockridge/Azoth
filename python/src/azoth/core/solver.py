"""The iterative solver used by implicit equations.

Mirrors ``azoth_core::solver`` in Rust, statement for statement.

This module is small and deliberately not general. It exists because
``hydraulics.friction_factor_colebrook`` is implicit in ``f``, and the two
implementations have to agree numerically. Two implementations running the *same*
algorithm converge to the same floating-point value; two running
different-but-equally-valid algorithms - Newton against fixed point, say - agree
only to within their difference, which is far larger than a tight tolerance
allows.

So the scheme is fixed by the spec (kind, tolerance, iteration cap, initial guess,
convergence rule) and both languages implement exactly it.

It lives in :mod:`azoth.core` rather than in a namespace package because every
namespace may need it and no namespace may depend on another. It was in
``azoth.hydraulics.reference`` while hydraulics was the only domain with an
implicit equation; the first calc outside it made that untenable.

The rule that governs what may be added here: **it holds the named schemes that
``specs/schema/calc.schema.json``'s ``solver.kind`` permits, and it grows only
when that enum grows** - each addition bringing the scheme in both languages and
the contract test that holds the enum to them, in one commit. It is not a
general-purpose numerical library and must not become one.
"""

from __future__ import annotations

import math
from collections.abc import Callable
from dataclasses import dataclass
from enum import StrEnum

from azoth.core.errors import InvalidInputError, SolverNotConvergedError


class SolverKind(StrEnum):
    """The named solution schemes the spec schema permits.

    A cross-language contract in the same way :class:`Convergence` and
    :class:`~azoth.core.warnings.WarningCode` are: the schema's ``solver.kind``
    enum, this class, and ``azoth_core::solver::SolverKind::ALL`` in Rust must
    name the same set, and ``python/tests/test_solver_contract.py`` asserts that
    rather than trusting three hand-edited lists to stay in step. The Rust list
    reaches Python through ``azoth._core.solver_kinds()``.

    Adding a name here is a claim that *both* implementations run it. The
    schema's own description of ``solver.kind`` says why that matters: a spec
    naming a kind neither can run would "describe a calculation neither
    implementation can run, which is the 'looks like validation, does nothing'
    failure this schema exists to catch".
    """

    #: ``x -> f(x)``, iterated from the spec's declared initial guess.
    FIXED_POINT = "fixed_point"
    #: The real roots of a cubic, formed analytically and then polished.
    CUBIC_ROOTS = "cubic_roots"

    @classmethod
    def parse(cls, name: str) -> SolverKind:
        """Parse the spec's spelling.

        Unknown names are an error rather than a silent default: guessing which
        scheme was meant would change the answer, and falling back to the only
        implemented kind would run a scheme the spec did not ask for.
        """
        try:
            return cls(name)
        except ValueError:
            known = [kind.value for kind in cls]
            raise InvalidInputError(
                "solver.kind",
                f"unknown solver kind `{name}`; expected one of {known}",
            ) from None


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


@dataclass(frozen=True, slots=True)
class CubicRootsOutcome:
    """What :func:`cubic_roots` produced.

    Not a :class:`SolverOutcome`, because a cubic has up to three answers rather
    than one and the caller - not this module - decides which of them is physical.
    ``roots`` is returned whole so that decision stays with the calc that knows
    what its roots mean.
    """

    #: The real roots, ascending, after polishing.
    #:
    #: Length 1 when the discriminant is positive and 3 when it is not, which
    #: includes the degenerate cases: a double or triple root is returned as
    #: repeated values rather than collapsed, so a caller that counts roots gets the
    #: algebraic count and not a de-duplicated one.
    roots: tuple[float, ...]
    #: Newton steps taken, summed over the roots.
    iterations: int
    #: Whether every root met the stopping rule.
    converged: bool
    #: The largest ``|x_k - x_{k-1}|`` seen on the final step of any root.
    residual: float


def cubic_roots(
    c2: float,
    c1: float,
    c0: float,
    tolerance: float,
    max_iterations: int,
    convergence: Convergence,
) -> CubicRootsOutcome:
    """The real roots of ``x**3 + c2*x**2 + c1*x + c0``, formed analytically then polished.

    # Why analytic first

    An iterative root-finder started from a fixed guess does not reliably find the
    root a caller wants. Measured on the Peng-Robinson cubic for propane at
    ``Tr = 0.8, Pr = 0.25`` - roots 0.0368, 0.1481, 0.7908 - Newton started at 0.30
    converges to the *middle* root, at 0.50 to the liquid root and at 0.70 to the
    vapour one. Which root comes back depends on where the search began, so an
    implementation that guesses cannot be told what it found. Forming the roots
    analytically and ordering them removes the guess.

    # Why the polish is not optional

    Three reasons, all measured rather than anticipated:

    * **Cardano cancels catastrophically when the roots are close.** The subtraction
      ``cbrt(-q/2 + sqrt(d)) - cbrt(-q/2 - sqrt(d))`` loses precision as the two
      terms approach each other, which is exactly the near-triple-root case.
    * **The discriminant's sign is not reliably computable near zero.** With the
      rounded Peng-Robinson constants at the critical point it is ``1.8e-12`` -
      positive, but small enough that a different rounding puts it on the other
      side, selecting the three-real-root branch for a cubic whose roots are mostly
      complex.
    * **``cbrt``, ``acos`` and ``cos`` are not correctly rounded** in any libm, so
      the analytic value is only close to begin with. The polish is what makes the
      two languages agree to a tolerance worth asserting.

    A Newton step on the polynomial is cheap, is exactly representable arithmetic,
    and converges quadratically once it is near - and the analytic formation has
    already put it near.
    """
    # Depress x**3 + c2*x**2 + c1*x + c0 to w**3 + p*w + q by x = w - c2/3.
    p = c1 - c2 * c2 / 3.0
    q = 2.0 * c2 * c2 * c2 / 27.0 - c2 * c1 / 3.0 + c0
    half_q = q / 2.0
    third_p = p / 3.0
    discriminant = half_q * half_q + third_p * third_p * third_p

    if discriminant > 0.0:
        # One real root and two complex conjugates. Cardano, with `cbrt` carrying
        # the sign so a negative argument does not become NaN.
        root_d = math.sqrt(discriminant)
        u = math.cbrt(-half_q + root_d)
        v = math.cbrt(-half_q - root_d)
        roots = [u + v - c2 / 3.0]
    else:
        radius = math.sqrt(-third_p * third_p * third_p)
        if radius == 0.0:
            # p = q = 0: the cubic is (w)**3, a triple root. The general form would
            # divide by `radius` here.
            roots = [-c2 / 3.0, -c2 / 3.0, -c2 / 3.0]
        else:
            # Clamped because rounding can push the ratio a hair outside [-1, 1],
            # and `acos` of that is NaN rather than a slightly wrong angle.
            cosine = max(-1.0, min(1.0, -half_q / radius))
            angle = math.acos(cosine)
            scale = 2.0 * math.sqrt(-third_p)
            roots = sorted(
                scale * math.cos((angle + 2.0 * math.pi * k) / 3.0) - c2 / 3.0 for k in range(3)
            )

    # Newton polish, on the polynomial rather than the depressed form: that is the
    # equation the caller wrote, and a step on it needs no back-substitution.
    iterations = 0
    residual = 0.0
    converged = True
    polished: list[float] = []
    for root in roots:
        x = root
        step_residual = math.inf
        root_converged = False
        for _ in range(max_iterations):
            f = ((x + c2) * x + c1) * x + c0
            df = (3.0 * x + 2.0 * c2) * x + c1
            if df == 0.0:
                # A stationary point: Newton cannot step, and the analytic value is
                # the best available. Leave it and let the stopping rule decide.
                break
            next_x = x - f / df
            delta = abs(next_x - x)
            step_residual = delta
            x = next_x
            iterations += 1
            if convergence is Convergence.ABSOLUTE:
                close_enough = delta <= tolerance
            else:
                close_enough = delta <= tolerance * abs(x)
            if close_enough:
                root_converged = True
                break
        if not root_converged:
            converged = False
        residual = max(residual, step_residual)
        polished.append(x)

    return CubicRootsOutcome(
        roots=tuple(polished),
        iterations=iterations,
        converged=converged,
        residual=residual,
    )


def require_cubic_converged(outcome: CubicRootsOutcome, tolerance: float) -> CubicRootsOutcome:
    """Turn a cubic outcome that did not converge into an error.

    The same rule :func:`require_converged` applies and for the same reason: a root
    that did not meet tolerance is the last iterate of an iteration that did not
    finish, and is not a root of anything.
    """
    if outcome.converged:
        return outcome
    raise SolverNotConvergedError(outcome.iterations, outcome.residual, tolerance)


def require_converged(outcome: SolverOutcome, tolerance: float) -> SolverOutcome:
    """Turn a non-converged outcome into an error.

    A caller that gets a value from a run that never converged has a number that
    solves nothing, so this is an error rather than a warning - unlike an
    out-of-range input, where the number is real and merely untrustworthy.
    """
    if outcome.converged:
        return outcome
    raise SolverNotConvergedError(outcome.iterations, outcome.residual, tolerance)
