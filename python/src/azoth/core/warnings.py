"""Warnings: conditions that do not stop a calculation but do change how much
its answer should be trusted.

Mirrors `azoth_core::warning` in Rust. The string values are a cross-language
contract - a caller must be able to compare a warning from either implementation
without knowing which one produced it - so `test_registry_contract.py` asserts
this enum is exactly the set the spec schema allows.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import StrEnum


class WarningCode(StrEnum):
    """Machine-readable warning identifiers.

    Values are SCREAMING_SNAKE_CASE, identical to
    `azoth_core::WarningCode::as_str()` on the Rust side. Never change an
    existing value; add a new member instead.
    """

    #: The value is outside the range in which the calculation is validated.
    OUT_OF_VALID_RANGE = "OUT_OF_VALID_RANGE"

    #: A range check could not be evaluated because an optional input was not
    #: supplied. Distinct from OUT_OF_VALID_RANGE: the value was never checked,
    #: which is not the same as having passed.
    RANGE_CHECK_SKIPPED = "RANGE_CHECK_SKIPPED"

    #: Flow is in the transitional band, where the friction factor is
    #: indeterminate rather than merely uncertain.
    TRANSITIONAL_FLOW = "TRANSITIONAL_FLOW"

    #: An iterative solver hit its iteration cap without meeting tolerance. The
    #: returned value is the last iterate, not a converged result.
    SOLVER_NOT_CONVERGED = "SOLVER_NOT_CONVERGED"

    # `UNVERIFIED_SOURCE` and `ESTIMATED_DATA` used to sit here. They existed to
    # stop data that had no source looking like data that had one, which was the
    # right worry while every number this library used was either the caller's or a
    # placeholder somebody typed. Once the data is vendored from a licence that
    # permits it, they fire on every result and mean nothing - and a warning that
    # fires on everything is how a reader learns to skim warnings, which is the
    # failure they were built to prevent.

    #: An iterative phase-equilibrium calculation converged to the trivial
    #: solution, ``x = y = z``. Distinct from ``SOLVER_NOT_CONVERGED``: this one
    #: *did* converge, and what it converged to is not a phase split.
    TRIVIAL_SOLUTION = "TRIVIAL_SOLUTION"


@dataclass(frozen=True, slots=True)
class Warning:
    """A single warning attached to a result.

    Frozen so a result's warnings cannot be edited in place after the fact, and
    comparable so the cross-language tests can assert structural parity rather
    than just comparing messages.
    """

    code: WarningCode
    message: str
    field: str | None = None

    def __str__(self) -> str:
        if self.field is None:
            return f"[{self.code}] {self.message}"
        return f"[{self.code}] {self.field}: {self.message}"
