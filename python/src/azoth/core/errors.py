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
    """The calculation has no verified source, and the caller asked for that to be
    a hard failure rather than a caveat.

    **Nothing raises this today.** It is defined, exported, and mapped from the
    Rust variant of the same name, but no code path constructs it. It is kept for
    one specific job: an opt-in strictness gate, in the shape of the existing
    ``AZOTH_REQUIRE_RUST``, through which a caller running design work can say "fail
    rather than hand me a number from an unconfirmed source".

    Whether a calc's inputs are placeholders or cited-but-unconfirmed values is
    recorded in the data files and in the specs, not on the result: no warning
    announces it, so a strictness gate would have to read the record rather than
    look for a warning code.

    One constraint on any future use: this must never be raised *from inside a
    calculation* as a substitute for a caveat. Warnings are not errors here, and a
    calc that raised on an unconfirmed source would break that rule and the
    cross-language warning-parity test that enforces it.
    """

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


class KeycardError(AzothError, ValueError):
    """A keycard cannot be loaded, or names something this build does not implement.

    Python-only, and necessarily so: **nothing else parses a keycard file**. A keycard
    is YAML, the workspace takes no YAML dependency, and `azoth.keycard` is the only
    implementation - so there is no second one for this class to be kept identical to.

    That is a narrower claim than it reads as. The Rust core holds a card when it is
    handed one, and refuses a name it cannot resolve with a different class; what is
    unreachable from Rust is *this* one, because what Rust never does is read the
    file.
    """

    def __init__(self, where: str, reason: str) -> None:
        self.where = where
        self.reason = reason
        super().__init__(f"keycard {where}: {reason}")

    def field(self) -> str | None:
        return self.where


__all__ = [
    "AzothError",
    "InvalidInputError",
    "KeycardError",
    "OutOfRangeError",
    "PropertyUnavailableError",
    "SolverNotConvergedError",
    "UnitMismatchError",
    "UnknownFittingError",
    "UnverifiedCalculationError",
]
