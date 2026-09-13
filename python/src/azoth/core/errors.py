"""Errors: conditions that stop a calculation.

Mirrors `azoth_core::error` in Rust. The split between an error and a warning
is deliberate and load-bearing, and both implementations apply it the same way:

* **Error** - the calculation is mathematically undefined, or the input is
  meaningless. ``Re = 0`` in the Colebrook equation divides by zero; a negative
  diameter is not a state anything can be computed from.
* **Warning** (`azoth.core.warnings`) - the calculation is well defined but
  lies outside the range where it has been validated. The number is real; its
  trustworthiness is the issue.

Getting that boundary wrong in either direction is a real failure. Treating a
range violation as an error makes the library unusable for exploratory work.
Treating a division by zero as a warning returns ``inf`` dressed as a result.

Each exception also inherits from the closest builtin, so callers who have never
heard of azoth can still catch ``ValueError`` or ``RuntimeError`` and get
something sensible.
"""

from __future__ import annotations


class AzothError(Exception):
    """Base class for everything this library raises.

    Deliberately does not inherit from a builtin, so ``except AzothError`` is a
    clean catch-all for the library without also catching unrelated ``ValueError``
    from user code.
    """

    def field(self) -> str | None:
        """The input this error is about, when it is about one.

        Mirrors the ``field()`` accessor on the Rust error type, so error
        handling code can be written once and work against either implementation.
        """
        return None


class InvalidInputError(AzothError, ValueError):
    """An input is the wrong shape or value for any computation to proceed."""

    def __init__(self, field: str, reason: str) -> None:
        self.input_field = field
        self.reason = reason
        super().__init__(f"invalid input `{field}`: {reason}")

    def field(self) -> str | None:
        return self.input_field


class OutOfRangeError(AzothError, ValueError):
    """An input violates a range check marked ``severity: error`` in the spec,
    meaning the calculation is undefined rather than merely out of range."""

    def __init__(self, field: str, value: float, detail: str) -> None:
        self.input_field = field
        self.value = value
        self.detail = detail
        super().__init__(f"`{field}` = {value} violates a hard limit: {detail}")

    def field(self) -> str | None:
        return self.input_field


class SolverNotConvergedError(AzothError, RuntimeError):
    """An iterative solver hit its cap without meeting tolerance."""

    def __init__(self, iterations: int, residual: float, tolerance: float) -> None:
        self.iterations = iterations
        self.residual = residual
        self.tolerance = tolerance
        super().__init__(
            f"solver did not converge after {iterations} iterations "
            f"(residual {residual:.3e}, tolerance {tolerance:.3e})"
        )


class UnitMismatchError(AzothError, TypeError):
    """A value was passed with the wrong dimensions, or without any.

    Inherits from ``TypeError`` rather than from pint's ``DimensionalityError``
    on purpose: ``azoth.core.errors`` must not import pint, and pint's
    exception class has moved between versions.
    """

    def __init__(self, field: str, expected: str, got: str) -> None:
        self.input_field = field
        self.expected = expected
        self.got = got
        super().__init__(f"`{field}` must be a quantity convertible to '{expected}', got {got!r}")

    def field(self) -> str | None:
        return self.input_field


class UnknownFittingError(AzothError, LookupError):
    """A fitting id could not be resolved against the fittings registry.

    An error rather than a skip: silently treating an unknown fitting as zero
    loss would under-report pressure drop, which is the dangerous direction to be
    wrong in.
    """

    def __init__(self, fitting_id: str) -> None:
        self.fitting_id = fitting_id
        super().__init__(f"unknown fitting `{fitting_id}`; see data/fittings/crane_k_factors.csv")

    def field(self) -> str | None:
        return self.fitting_id


class UnverifiedCalculationError(AzothError):
    """The calculation has no verified source, so it must not be presented as a
    validated result."""

    def __init__(self, calc_id: str) -> None:
        self.calc_id = calc_id
        super().__init__(f"calculation `{calc_id}` has no verified source")

    def field(self) -> str | None:
        return self.calc_id


class PropertyUnavailableError(AzothError, LookupError):
    """A fluid property was requested that the provider cannot supply."""

    def __init__(self, fluid: str, property_name: str, reason: str) -> None:
        self.fluid = fluid
        self.property_name = property_name
        super().__init__(f"{fluid} has no {property_name}: {reason}")

    def field(self) -> str | None:
        return self.property_name


__all__ = [
    "AzothError",
    "InvalidInputError",
    "OutOfRangeError",
    "PropertyUnavailableError",
    "SolverNotConvergedError",
    "UnitMismatchError",
    "UnknownFittingError",
    "UnverifiedCalculationError",
]
