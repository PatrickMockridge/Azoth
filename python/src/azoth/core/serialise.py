"""A result, as data.

Every calculation in this library returns a frozen dataclass, and until now none of them could
leave Python: a front-end, a notebook cell or an agent had to know each result's fields and read
them one at a time. This module writes any of them, and **there is no per-result code** — the walk
is over ``dataclasses.fields``, so the one thing a new result has to be is a dataclass.

That is the whole answer to the objection the middleware's gap table recorded: "a codec per model
is the second writer `executor::json` exists to avoid". There is no codec per model. There is one
writer, and 187 result classes inherit it through :class:`azoth.core.result._HasWarnings`.

# The shape

A result is an object of its fields; a field that is itself a dataclass (a ``Warning``, a
``KComponent``) is an object; a list or tuple is an array; a ``StrEnum`` is its value; and a
quantity is a **number and a unit name**:

    {"re": 149401.2, "regime": "turbulent", "warnings": [{"code": "RANGE_CHECK_SKIPPED", ...}]}
    {"density": {"magnitude": 998.0, "unit": "kg/m**3"}}

**The magnitude is in the unit named, and the name is the vocabulary's own.** That is what the
generated docs print and what ``azoth.core.units.from_si`` hands back, so a reader of a page and a
reader of this document see the same number. The unit is the *spec's* name — ``from_si(x, "mm")``
answers in millimetres, and this says ``"mm"`` rather than ``"millimeter"`` — which makes
``azoth.core.units.unit_for`` the way back, the same one every input takes.

**Why the key is not `magnitude_si`**, which is what ``executor::json`` writes for a stream
record's five fields: those fields are SI by construction (``Pa``, ``K``, ``mol/s``, ``J/mol``), so
that codec's name is true of them. A result's unit is whatever its spec declares, and ``mm`` is a
unit a spec may declare — so this key does not claim SI. Both layers name the *unit* in the same
key, which is the half a consumer reads the same way in either.

# A magnitude that is not a number

It is written as ``null``, and that is a deliberate difference from ``executor::json``, which
*refuses* a non-finite magnitude for a stream. Both are right about their own case: a session that
ran cannot produce one, so a ``NaN`` reaching that codec would be a computation gone wrong and a
``null`` would hide it; while a result carries one **on purpose** — running every registered case
is what found it, `reactions.reactive_phase_equilibrium`'s `a_single_gas_phase_is_skipped` has no
residual because there was nothing to solve, and it says so in its own case name. ``null`` is what
that means: no value there. Refusing would mean a legitimate state had no JSON form at all.

# What is refused

A field whose type the walk does not know. That is refused rather than stringified, because
``str()`` of an object is a value nothing downstream can read back, and the refusal names the path
to the field so a new result that carries something unserialisable says where.
"""

from __future__ import annotations

import array
import dataclasses
import json
from collections.abc import Mapping
from enum import Enum
from typing import Any, Final

from azoth.core._units_gen import CANONICAL_UNITS
from azoth.core.units import ureg

__all__ = ["pint_name", "to_dict", "to_json", "unit_name"]


def _comparable(unit: str) -> str:
    """A unit name in `pint`'s own spelling, which is the only form two names can be compared in.

    **Comparing the strings as written does not work, and the first run of the test beside this is
    what found it.** The vocabulary stores `joule/(mole*kelvin)` and a quantity made from the same
    spec unit prints `joule / kelvin / mole`: same unit, different arrangement, because `pint`
    normalises a denominator's factors into separate divisions. Passing *both* sides through
    `pint`'s formatter is what removes the difference - and it is the same formatter a result's
    quantity went through, so the two are comparable by construction rather than by a rule about
    spaces that would have to be maintained.
    """
    return str(ureg.Unit(unit))


#: The unit's vocabulary name, from its name in `pint`. The inverse of
#: ``_units_gen.CANONICAL_UNITS``, whose entries are a bijection — a test asserts that, so a future
#: unit whose vocabulary name collides with another's is caught here rather than losing one of the
#: two names silently.
_BY_PINT_NAME: Final[dict[str, str]] = {
    _comparable(pint): spec for spec, pint in CANONICAL_UNITS.items()
}


def unit_name(unit: Any) -> str:
    """The name to write for a unit: the vocabulary's, or `pint`'s where it has none.

    **A fallback and not a guess.** A calc can produce a unit no spec names — a metres-per-second
    over a kilogram, say — and inventing a vocabulary name for it would be worse than writing the
    one `pint` parses. The vocabulary covers every unit a *spec* declares, which is every unit a
    result field is documented in.
    """
    return _BY_PINT_NAME.get(_comparable(str(unit)), str(unit))


def pint_name(name: str) -> str:
    """The `pint` unit a vocabulary name stands for, or the name itself.

    The way back for :func:`unit_name`, and the same mapping `azoth.core.units.unit_for` applies to
    a spec's own unit string.
    """
    return CANONICAL_UNITS.get(name) or _BY_PINT_NAME.get(_comparable(name)) or name


def is_finite(value: float) -> bool:
    """Whether a magnitude has a JSON representation at all.

    ``json.dumps`` writes ``NaN`` and ``Infinity`` by default, neither of which is JSON, and a
    reader that gets one takes it for a number.
    """
    return value == value and value not in (float("inf"), float("-inf"))


def to_dict(value: Any, *, path: str = "result") -> Any:
    """One value, as JSON-shaped data.

    ``path`` is the dotted route to where the value sits, and it is only ever used to say where a
    refusal happened — a new result with an unserialisable field fails naming the field rather than
    the result.
    """
    # **The enum branch comes before the scalar one**, because a `StrEnum` *is* a `str`: checked
    # the other way round, the member itself fell through as a string and `json.dumps` happened to
    # write the right thing while `to_dict` handed back an enum member - which is not JSON-shaped
    # data, and which a caller iterating the dict would find is not a `str` either.
    if isinstance(value, Enum):
        # A `StrEnum`'s value is the string a spec declares, which is what a document should carry
        # rather than the member's name.
        return value.value
    if value is None or isinstance(value, (bool, int, str)):
        return value
    if isinstance(value, float):
        return value if is_finite(value) else None
    if isinstance(value, ureg.Quantity):
        magnitude = float(value.magnitude)
        return {
            "magnitude": magnitude if is_finite(magnitude) else None,
            "unit": unit_name(value.units),
        }
    if dataclasses.is_dataclass(value) and not isinstance(value, type):
        return {
            field.name: to_dict(getattr(value, field.name), path=f"{path}.{field.name}")
            for field in dataclasses.fields(value)
        }
    if isinstance(value, Mapping):
        return {str(key): to_dict(item, path=f"{path}.{key}") for key, item in value.items()}
    if isinstance(value, (tuple, list, array.array)):
        # **`array.array` is here because a batch result's columns are one**, and the walk's refusal
        # is what said so - naming the field (`result.re`) rather than writing a stringified buffer.
        # A JSON array is the same sequence of numbers, and the elements go through the same
        # branches, so a non-finite one inside a column is written as no-value like any other.
        return [to_dict(item, path=f"{path}[{index}]") for index, item in enumerate(value)]
    raise TypeError(
        f"{path}: {type(value).__name__} is not something this writer knows, so it is refused "
        f"rather than stringified - a result field of a new type needs a branch here"
    )


def to_json(value: Any) -> str:
    """One value, as the JSON document.

    Sorted keys and no whitespace, so two runs of the same calculation write byte-identical text
    and a diff between two results is a diff between two numbers.
    """
    return json.dumps(to_dict(value), sort_keys=True, separators=(",", ":"), allow_nan=False)
