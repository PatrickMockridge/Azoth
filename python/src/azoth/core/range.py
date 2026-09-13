"""Range checking, driven by the ``valid_range`` block of a calc spec.

Mirrors ``azoth_core::range`` in Rust. The two implementations evaluate the
same typed bounds with the same semantics, which is only possible because the
specs describe bounds as numbers rather than as expressions like ``"> 0"`` - a
free-form expression would need a parser in each language, and the two parsers
would eventually disagree.

The distinction this module exists to make is between a value that was *checked
and passed* and a value that was *never checked*. Those look identical to a
caller unless the second case says so, so :func:`apply_checks` emits a
``RANGE_CHECK_SKIPPED`` warning whenever a check cannot be evaluated. That
happens in practice: ``hydraulics.darcy_weisbach`` takes viscosity as an optional
input, and without it there is no Reynolds number to check against.
"""

from __future__ import annotations

import math
from collections.abc import Callable, Iterable, Mapping
from dataclasses import dataclass
from enum import StrEnum
from typing import Any

from azoth.core.errors import InvalidInputError, OutOfRangeError
from azoth.core.warnings import Warning, WarningCode


class Severity(StrEnum):
    """What happens when a range check is violated. Mirrors ``severity``."""

    WARNING = "warning"
    ERROR = "error"


class Band(StrEnum):
    """Which side of the interval counts as a violation. Mirrors ``when``."""

    #: The interval is the *allowed* range; values beyond it violate. The
    #: default, and what almost every check means.
    OUTSIDE = "outside"

    #: The interval is the *forbidden* band; values inside it violate. Needed for
    #: transitional flow, which is a problem between 2000 and 4000 rather than
    #: outside it.
    INSIDE = "inside"


def flatten(text: str) -> str:
    """Collapse whitespace in a piece of spec prose.

    YAML folded scalars arrive with embedded newlines and a trailing one. Those
    are a formatting artefact of the spec file, not meaning, and both languages
    must render the same rationale the same way - the Rust codegen flattens at
    generation time, so Python flattens at load time.
    """
    return " ".join(text.split())


def _format_number(value: float) -> str:
    """Render a bound for a human.

    Not byte-identical to Rust's ``Display`` for ``f64`` - Python's ``:g`` and
    Rust's default differ for very small or very large values. Nothing depends on
    them matching: the generated docs come from the spec via this module, and the
    Rust text only ever reaches an error message.
    """
    return f"{value:g}"


@dataclass(frozen=True, slots=True)
class RangeCheck:
    """One bound on one quantity."""

    quantity: str
    min: float | None
    min_inclusive: bool
    max: float | None
    max_inclusive: bool
    #: A single forbidden value, if this bound is an exclusion rather than an
    #: interval. Mirrors ``RangeCheck::equals`` in Rust, including the reason it
    #: exists: the schema permits it, nothing implemented it, and a bound that
    #: silently does nothing is worse than an absent one.
    equals: float | None
    band: Band
    severity: Severity
    code: WarningCode
    rationale: str

    @classmethod
    def from_spec(cls, raw: Mapping[str, Any]) -> RangeCheck:
        """Build from one entry of a spec's ``valid_range`` list."""
        quantity = str(raw["quantity"])
        try:
            code = WarningCode(raw.get("code", WarningCode.OUT_OF_VALID_RANGE))
            band = Band(raw.get("when", "outside"))
            severity = Severity(raw["severity"])
        except ValueError as exc:
            raise InvalidInputError(
                quantity, f"spec declares an unreadable range check: {exc}"
            ) from exc
        return cls(
            quantity=quantity,
            min=None if raw.get("min") is None else float(raw["min"]),
            min_inclusive=bool(raw.get("min_inclusive", True)),
            max=None if raw.get("max") is None else float(raw["max"]),
            max_inclusive=bool(raw.get("max_inclusive", True)),
            equals=None if raw.get("equals") is None else float(raw["equals"]),
            band=band,
            severity=severity,
            code=code,
            rationale=flatten(str(raw.get("rationale", ""))),
        )

    def violated(self, value: float) -> bool:
        """Whether `value` violates this bound. ``NaN`` violates everything.

        The ``NaN`` case is not defensive padding. Every comparison with ``NaN``
        is false, so without an explicit check a ``NaN`` input would satisfy
        every bound and produce a confident-looking wrong answer - exactly the
        failure this library exists to prevent.
        """
        if math.isnan(value):
            return True
        # An exclusion rather than an interval. Checked first and returned
        # unconditionally: `band` describes which side of an interval violates, and
        # has no meaning for a single forbidden value.
        if self.equals is not None:
            return value == self.equals
        below = False
        if self.min is not None:
            below = value < self.min if self.min_inclusive else value <= self.min
        above = False
        if self.max is not None:
            above = value > self.max if self.max_inclusive else value >= self.max
        outside = below or above
        return outside if self.band is Band.OUTSIDE else not outside

    def describe(self) -> str:
        """Human-readable form of the bound, e.g. ``'Re >= 4000'``."""
        quantity = self.quantity
        if self.equals is not None:
            return f"{quantity} == {_format_number(self.equals)}"
        parts: list[str] = []
        if self.min is not None:
            operator = ">=" if self.min_inclusive else ">"
            parts.append(f"{operator} {_format_number(self.min)}")
        if self.max is not None:
            operator = "<=" if self.max_inclusive else "<"
            parts.append(f"{operator} {_format_number(self.max)}")

        if not parts:
            bounds = f"{quantity} (unbounded)"
        elif len(parts) == 2:
            bounds = f"{quantity} {parts[0]} and {quantity} {parts[1]}"
        else:
            bounds = f"{quantity} {parts[0]}"

        if self.band is Band.INSIDE:
            return f"{bounds} (violation lies inside this band)"
        return bounds

    def apply(self, value: float, warnings: list[Warning]) -> None:
        """Apply this bound to a value.

        Raises:
            OutOfRangeError: when violated with ``Severity.ERROR``.
        """
        if not self.violated(value):
            return
        if self.severity is Severity.ERROR:
            raise OutOfRangeError(
                self.quantity,
                value,
                f"must satisfy {self.describe()}. {self.rationale}",
            )
        warnings.append(
            Warning(
                code=self.code,
                message=f"value {value} violates {self.describe()}. {self.rationale}",
                field=self.quantity,
            )
        )


@dataclass(frozen=True, slots=True)
class SpecChecks:
    """A spec's bounds, split by when they can be evaluated.

    Input checks run before the calculation and guard the arithmetic, because
    their severity is usually ``error``. Derived checks need the calculation to
    have run - and, when an optional input was omitted, may not be evaluable at
    all.
    """

    on_input: tuple[RangeCheck, ...]
    derived: tuple[RangeCheck, ...]

    @property
    def all(self) -> tuple[RangeCheck, ...]:
        """Every bound, in spec order."""
        return self.on_input + self.derived


def checks_for(spec: Mapping[str, Any]) -> SpecChecks:
    """Split a spec's ``valid_range`` into input and derived checks."""
    inputs = spec["inputs"]
    on_input: list[RangeCheck] = []
    derived: list[RangeCheck] = []
    for raw in spec["valid_range"]:
        check = RangeCheck.from_spec(raw)
        (on_input if check.quantity in inputs else derived).append(check)
    return SpecChecks(on_input=tuple(on_input), derived=tuple(derived))


def apply_checks(
    checks: Iterable[RangeCheck],
    resolve: Callable[[str], float | None],
    warnings: list[Warning],
) -> None:
    """Apply a set of bounds, resolving each quantity through `resolve`.

    `resolve` returns ``None`` for a quantity that cannot be computed from the
    inputs actually supplied - typically because an optional input such as
    viscosity was omitted. That is not a pass: it produces a
    ``RANGE_CHECK_SKIPPED`` warning, so the caller can tell "checked and fine"
    from "never checked".
    """
    for check in checks:
        value = resolve(check.quantity)
        if value is None:
            warnings.append(
                Warning(
                    code=WarningCode.RANGE_CHECK_SKIPPED,
                    message=(
                        f"could not check `{check.quantity}` (must satisfy "
                        f"{check.describe()}) because an input it depends on was not "
                        f"supplied; this range is UNVERIFIED for this result"
                    ),
                    field=check.quantity,
                )
            )
        else:
            check.apply(value, warnings)
