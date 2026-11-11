"""A calculation's arguments as a JSON document, and back again.

The declaration in ``_registry_gen.spec(calc_id)["inputs"]`` is the authority, and this
mirrors what ``python/tests/_helpers.py``'s ``kwargs_for`` already does: an untyped input
is a quantity whose magnitude is stated in the spec's unit, ``dimensionless`` is a plain
number, an enum or a string passes through, and a fitting list is a list of names.

The encoded form is what a document stores, so a unitful scalar becomes an SI magnitude
and its vocabulary unit name; decoding rebuilds the pint quantity the reference function
takes. Settling on SI is what lets a validator re-run the calculation without knowing
which unit spelling the claimant used.
"""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any

import pint

from azoth._registry_gen import CALCS
from azoth.core.units import from_si, quantity, to_si

__all__ = ["InputError", "UnknownCalc", "UnsupportedCalc", "decode", "encode"]

_BY_ID: dict[str, dict[str, Any]] = {entry["id"]: entry for entry in CALCS}

_QUANTITY = "quantity"
_LIST_OF_NAMES = frozenset({"fitting_list", "vector_string"})
_REFUSED = frozenset({"record_list", "components", "vector", "matrix"})


class InputError(ValueError):
    """An input the bundle cannot encode or decode, named by its field."""


class UnknownCalc(InputError):
    """An id that is not one of the registered calculations."""


class UnsupportedCalc(InputError):
    """A declared input kind this bundle has no codec for, named rather than guessed."""


def _declared(calc_id: str) -> dict[str, Any]:
    entry = _BY_ID.get(calc_id)
    if entry is None:
        raise UnknownCalc(f"{calc_id!r} is not a registered calculation")
    inputs: dict[str, Any] = entry["inputs"]
    return inputs


def _magnitude(value: Any, unit: str, field: str) -> float:
    if isinstance(value, pint.Quantity):
        return to_si(value, unit, field)
    return to_si(quantity(float(value), unit), unit, field)


def encode(calc_id: str, kwargs: Mapping[str, Any]) -> dict[str, Any]:
    """The canonical document form of a call's keyword arguments.

    A number is read in the spec's declared unit (the form a spec's own case uses); a pint
    quantity is converted from whatever unit it carries. A required input that is absent
    is named here rather than left for the call to report as a ``TypeError``.
    """
    declared = _declared(calc_id)
    unknown = sorted(set(kwargs) - set(declared))
    if unknown:
        raise InputError(f"{calc_id}: inputs the spec does not declare: {unknown}")

    encoded: dict[str, Any] = {}
    for name, value in kwargs.items():
        declaration = declared[name]
        kind = declaration.get("type", _QUANTITY)
        if kind == _QUANTITY:
            unit = declaration.get("unit", "dimensionless")
            magnitude = _magnitude(value, unit, name)
            encoded[name] = (
                magnitude if unit == "dimensionless" else {"value": magnitude, "unit": unit}
            )
        elif kind in ("enum", "string"):
            text = str(value)
            allowed = declaration.get("values")
            if allowed is not None and text not in allowed:
                raise InputError(f"{calc_id}.{name}: {text!r} is not one of {allowed}")
            encoded[name] = text
        elif kind == "boolean":
            encoded[name] = bool(value)
        elif kind in _LIST_OF_NAMES:
            encoded[name] = [str(item) for item in value]
        else:
            raise UnsupportedCalc(f"{calc_id}.{name}: the {kind!r} input kind has no bundle codec")

    missing = sorted(
        name
        for name, declaration in declared.items()
        if not declaration.get("optional", False) and name not in encoded
    )
    if missing:
        raise InputError(f"{calc_id}: required inputs absent: {missing}")
    return encoded


def decode(calc_id: str, document: Mapping[str, Any]) -> dict[str, Any]:
    """The keyword arguments a reference function takes, rebuilt from a document."""
    declared = _declared(calc_id)
    unknown = sorted(set(document) - set(declared))
    if unknown:
        raise InputError(f"{calc_id}: document carries undeclared inputs: {unknown}")

    kwargs: dict[str, Any] = {}
    for name, value in document.items():
        declaration = declared[name]
        kind = declaration.get("type", _QUANTITY)
        if kind == _QUANTITY:
            unit = declaration.get("unit", "dimensionless")
            if unit == "dimensionless":
                kwargs[name] = float(value)
                continue
            if not isinstance(value, Mapping) or set(value) != {"value", "unit"}:
                raise InputError(f"{calc_id}.{name}: expected {{value, unit}}, got {value!r}")
            if value["unit"] != unit:
                raise InputError(
                    f"{calc_id}.{name}: the document says {value['unit']!r}, the spec {unit!r}"
                )
            kwargs[name] = from_si(float(value["value"]), unit)
        elif kind in ("enum", "string"):
            kwargs[name] = str(value)
        elif kind == "boolean":
            kwargs[name] = bool(value)
        elif kind in _LIST_OF_NAMES:
            kwargs[name] = [str(item) for item in value]
        else:
            raise UnsupportedCalc(f"{calc_id}.{name}: the {kind!r} input kind has no bundle codec")
    return kwargs
