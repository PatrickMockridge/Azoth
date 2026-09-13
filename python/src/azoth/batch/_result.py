"""What a batch call returns.

Struct-of-arrays: one sequence per output field, and a positional sequence of warning
tuples. ``warnings[i]`` is element ``i``'s warnings, empty for a clean element.

# The shape, and what it is not

An *array of result objects* would allocate N dataclasses and N quantities, giving back
most of what crossing the language boundary once saved - the boundary crossing being the
cost the batch API exists to remove. A *table of rows* would leave nowhere to put
per-element warnings.

An output that carries dimension is an `array.array` of **SI base magnitudes**, and
:attr:`BatchResult.units` says what each field is in. That is deliberately not the same
kind of thing the scalar API returns - a `pint` quantity - and the difference is the
price of the shape. :meth:`BatchResult.quantities` converts a field to quantities for a
caller who wants them, which allocates, and says so.

An output that is an *enum* - a flow regime - is a tuple of labels, `None` where the
value is absent. There is no number a regime should be, and an index into a table would
make the caller look the mapping up to read a value and need a sentinel for "not
available" that is not one of the real answers.
"""

from __future__ import annotations

from array import array
from dataclasses import dataclass
from typing import TYPE_CHECKING

from azoth.core.warnings import Warning, WarningCode

if TYPE_CHECKING:
    from collections.abc import Mapping

    from azoth.core.units import Q


@dataclass(frozen=True, slots=True, eq=False)
class BatchResult:
    """Common behaviour for every batch result.

    Not constructed directly. Each calculation has its own subclass carrying that
    calc's fields as named attributes, for the reason the scalar results do: a field
    name is a contract, and `result.dp` should be a field the type checker can see
    rather than a dictionary lookup the reader has to find.
    """

    #: Per element, in order. Empty tuples for clean elements.
    warnings: tuple[tuple[Warning, ...], ...]
    #: Field name -> the spec's canonical unit for it, e.g. ``{"dp": "Pa"}``.
    units: Mapping[str, str]

    def __len__(self) -> int:
        """How many elements the batch covered."""
        return len(self.warnings)

    @property
    def is_clean(self) -> bool:
        """True when no element carries a warning.

        A property rather than a method, where the scalar results use a property and the
        shape differs: a scalar result is clean or not, and a batch is clean only if
        every element is. Reading ``result.is_clean`` on a batch is asking about all of
        it, which is the question the name asks.
        """
        return all(not element for element in self.warnings)

    def has_warning(self, code: WarningCode) -> bool:
        """Whether any element carries a warning with this code."""
        return any(any(w.code == code for w in element) for element in self.warnings)

    def warning_indices(self, code: WarningCode) -> tuple[int, ...]:
        """Which elements carry a warning with this code, in order.

        The method that makes per-element warnings usable: a caller who sees
        ``has_warning`` return True needs to know *which* elements, and scanning
        ``warnings`` by hand for the ones that matter is where a batch API stops being
        worth having.
        """
        return tuple(
            index
            for index, element in enumerate(self.warnings)
            if any(w.code == code for w in element)
        )

    def quantities(self, field: str) -> tuple[Q, ...]:
        """One field as `pint` quantities, for a caller who wants them.

        Named for what it does rather than for convenience, and deliberately not the
        default: this allocates one quantity per element, which is most of what the
        batch shape avoids. Calling it is a deliberate choice to spend that.

        Raises:
            AttributeError: if the field is not one of this calculation's outputs, or if
                it is an enum field - there is no unit to give a label.
        """
        from azoth.core.units import ureg

        unit = self.units.get(field)
        if unit is None:
            raise AttributeError(
                f"{type(self).__name__} has no dimensioned field {field!r}; "
                f"fields with units: {sorted(self.units)}"
            )
        values: array[float] = getattr(self, field)
        return tuple(ureg.Quantity(value, unit) for value in values)

    def __repr__(self) -> str:
        """A summary, not a dump.

        The dataclass-generated repr would print every element of every column, which for
        a real batch is thousands of numbers in a traceback or a log line. A subclass is
        declared ``repr=False`` for exactly this reason: without it the generated method
        shadows this one, silently, and the only symptom is an unreadable log.
        """
        return f"{type(self).__name__}(n={len(self)}, clean={self.is_clean})"


__all__ = ["BatchResult"]
