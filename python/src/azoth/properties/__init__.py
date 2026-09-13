"""Fluid properties for the calcs that need them.

Two built-ins, water and air, each backed by a small table under ``data/fluids/``
with per-row citations. That is deliberately the whole of it: the specs say
property data beyond what the worked examples need is out of scope, and a large
half-verified property database is a bigger liability than a small cited one.

# Provenance

The values are real published figures read from widely used engineering tables
that trace back to IAPWS-IF97 / NIST and the CRC Handbook. They are marked
``unverified`` because they have not been checked against a primary formulation -
the same discipline the fittings registry uses, though for a different reason:
these are not placeholders, they simply have not been through a human's hands.

# Using a different provider

:class:`~azoth.properties.provider.PropertyProvider` is a protocol, so anything
with ``density`` and ``dynamic_viscosity`` methods satisfies it. The calcs take
quantities, not providers, so a provider is only ever needed by code that has to
turn a fluid name and a temperature into a velocity - the CLI, for instance.
"""

from __future__ import annotations

from azoth.core.errors import PropertyUnavailableError
from azoth.properties.provider import (
    STANDARD_PRESSURE_PA,
    PropertyPoint,
    PropertyProvider,
    TableProvider,
)

__all__ = [
    "AIR",
    "STANDARD_PRESSURE_PA",
    "WATER",
    "PropertyPoint",
    "PropertyProvider",
    "TableProvider",
    "available_fluids",
    "provider_for",
]

#: Water at 1 atm, tabulated 0-100 C.
WATER: TableProvider = TableProvider("water", "data/fluids/water.csv")

#: Dry air at 1 atm, tabulated 0-100 C.
AIR: TableProvider = TableProvider("air", "data/fluids/air.csv")

_BUILTINS: dict[str, TableProvider] = {"water": WATER, "air": AIR}


def available_fluids() -> tuple[str, ...]:
    """Names of the built-in fluids."""
    return tuple(sorted(_BUILTINS))


def provider_for(name: str) -> TableProvider:
    """Look up a built-in provider by fluid name.

    Raises:
        PropertyUnavailableError: if the fluid is not one of the built-ins. An
            error rather than a default: silently substituting water for an
            unknown fluid would produce a plausible pressure drop for the wrong
            substance.
    """
    try:
        return _BUILTINS[name.strip().lower()]
    except KeyError:
        raise PropertyUnavailableError(
            name,
            "fluid",
            f"not a built-in fluid; available: {', '.join(available_fluids())}",
        ) from None
