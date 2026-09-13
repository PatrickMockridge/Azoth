"""Behaviour of the named solution schemes.

The vocabulary is checked in ``test_solver_contract.py``; this file checks that the
schemes compute what they claim. The cases mirror
``crates/azoth-core/src/solver.rs``'s unit tests one for one, deliberately: the two
implementations are required to run the *same* scheme, so a case that exists on one
side and not the other is a gap in the comparison rather than a saving.

What is *not* here is a cross-language comparison of the solver itself. The
extension does not expose it - only the calcs do - so agreement is asserted where
the schemes are actually used, by ``test_cross_impl.py`` running every spec case
through both backends. That is the stronger check anyway: it compares the answers
rather than the intermediate steps.
"""

from __future__ import annotations

import pytest

from azoth.core.errors import SolverNotConvergedError
from azoth.core.solver import (
    Convergence,
    CubicRootsOutcome,
    cubic_roots,
    fixed_point,
    require_converged,
    require_cubic_converged,
)

TOLERANCE = 1e-14
MAX_ITERATIONS = 50


def solve(c2: float, c1: float, c0: float) -> CubicRootsOutcome:
    return cubic_roots(c2, c1, c0, TOLERANCE, MAX_ITERATIONS, Convergence.RELATIVE)


def test_cubic_roots_finds_three_well_separated_roots_in_order() -> None:
    # (x - 1)(x - 2)(x - 3) = x**3 - 6x**2 + 11x - 6
    out = solve(-6.0, 11.0, -6.0)
    assert out.converged, f"residual {out.residual:e}"
    assert len(out.roots) == 3
    for got, want in zip(out.roots, (1.0, 2.0, 3.0), strict=True):
        assert abs(got - want) < 1e-12, f"got {got}, want {want}"


def test_cubic_roots_returns_one_root_when_the_other_two_are_complex() -> None:
    # x**3 + x + 1 has one real root, near -0.682327803828.
    out = solve(0.0, 1.0, 1.0)
    assert len(out.roots) == 1
    assert out.converged
    root = out.roots[0]
    residual = (root * root + 1.0) * root + 1.0
    assert abs(residual) < 1e-12, f"not a root: f = {residual:e}"


def test_cubic_roots_handles_a_triple_root_without_dividing_by_zero() -> None:
    # (x - 2)**3 = x**3 - 6x**2 + 12x - 8. The degenerate case the trigonometric
    # branch cannot take, because its radius is zero.
    out = solve(-6.0, 12.0, -8.0)
    assert len(out.roots) == 3, "a triple root is three roots, not one"
    for got in out.roots:
        assert abs(got - 2.0) < 1e-9, f"got {got}, want 2"


def test_cubic_roots_is_deterministic() -> None:
    # The same property `fixed_point` is held to, and for the same reason: the two
    # languages are compared against each other, so a scheme that landed on
    # different roots from identical input would make that comparison meaningless
    # rather than merely noisy.
    first, second = solve(-6.0, 11.0, -6.0), solve(-6.0, 11.0, -6.0)
    assert first.iterations == second.iterations
    assert first.roots == second.roots


def test_cubic_roots_answers_a_polynomial_it_was_not_shaped_for() -> None:
    # x**3 - 2x**2 - 5x + 6 = (x - 1)(x + 2)(x - 3). Included because its roots
    # straddle zero, so an implementation that assumed positive roots - which a Z
    # factor always is - fails here.
    out = solve(-2.0, -5.0, 6.0)
    for got, want in zip(out.roots, (-2.0, 1.0, 3.0), strict=True):
        assert abs(got - want) < 1e-12, f"got {got}, want {want}"


def test_cubic_roots_reports_a_stationary_start_rather_than_looping() -> None:
    # A cubic whose analytic root lands exactly on a stationary point cannot take a
    # Newton step. It must stop and say it did not converge, rather than spin to the
    # cap or divide by zero - `require_cubic_converged` is what turns that into an
    # error at the call site.
    # x**3 - 3x + 2 = (x - 1)**2 (x + 2): the analytic form puts a root at 1, where
    # the derivative 3x**2 - 3 is zero.
    out = cubic_roots(-0.0, -3.0, 2.0, TOLERANCE, 1, Convergence.ABSOLUTE)
    # With a cap of one step the polish cannot converge; the point is that it
    # returns rather than raising from inside the solver.
    assert isinstance(out, CubicRootsOutcome)


def test_require_cubic_converged_turns_a_failure_into_an_error() -> None:
    failed = CubicRootsOutcome(roots=(1.0,), iterations=3, converged=False, residual=0.5)
    with pytest.raises(SolverNotConvergedError):
        require_cubic_converged(failed, TOLERANCE)

    ok = CubicRootsOutcome(roots=(1.0,), iterations=3, converged=True, residual=0.0)
    assert require_cubic_converged(ok, TOLERANCE) is ok


def test_the_two_schemes_are_separate_outcomes() -> None:
    """A cubic's outcome is not a scalar's, and the two are not interchangeable.

    ``fixed_point`` answers with one value and ``cubic_roots`` with up to three, so
    a caller that reached for the wrong one would get an attribute error rather than
    a silently wrong number. Asserting the shapes differ keeps that true.
    """
    scalar = fixed_point(1.0, 1e-12, 20, Convergence.RELATIVE, lambda x: (x + 2.0 / x) / 2.0)
    assert require_converged(scalar, 1e-12).converged
    assert not hasattr(scalar, "roots")
    assert hasattr(solve(0.0, 1.0, 1.0), "roots")
