"""``eos.scale_remediation_advisor`` - NeqSim's ``ScaleRemediationAdvisor``.

A lookup over a vendored pipe-delimited knowledge base: a deposit name is normalised
through the class's alias table, its rows are selected (or the ``UNKNOWN`` fallback
is), and they are ordered by descending effectiveness.

Spec: ``specs/models/eos/scale_remediation_advisor.toml``.
"""

from __future__ import annotations

from functools import cache
from typing import Any

from azoth._data import find
from azoth.core.record_outputs_gen import RemediationOptionRow
from azoth.core.result import ScaleRemediationAdvisorResult
from azoth.core.units import from_si

CALC_ID = "eos.scale_remediation_advisor"

REGISTRY_PATH = "data/scale/scale_remediation.csv"

#: The class's alias table, lower-case alias to canonical key.
ALIASES = {
    "calcite": "CaCO3",
    "aragonite": "CaCO3",
    "calcium carbonate": "CaCO3",
    "siderite": "FeCO3",
    "iron carbonate": "FeCO3",
    "mackinawite": "FeS",
    "pyrrhotite": "FeS",
    "troilite": "FeS",
    "iron sulfide": "FeS",
    "iron sulphide": "FeS",
    "barite": "BaSO4",
    "barytes": "BaSO4",
    "barium sulfate": "BaSO4",
    "barium sulphate": "BaSO4",
    "celestite": "SrSO4",
    "strontium sulfate": "SrSO4",
    "strontium sulphate": "SrSO4",
    "gypsum": "CaSO4",
    "anhydrite": "CaSO4",
    "calcium sulfate": "CaSO4",
    "calcium sulphate": "CaSO4",
    "halite": "NaCl",
    "salt": "NaCl",
    "dithiazine": "Dithiazine",
    "triazine solids": "Dithiazine",
    "scavenger solids": "Dithiazine",
    "asphaltene": "Asphaltene",
    "asphaltenes": "Asphaltene",
    "wax": "Wax",
    "paraffin": "Wax",
    "naphthenate": "Naphthenate",
    "calcium naphthenate": "Naphthenate",
}


def canonical_key(deposit: str) -> str:
    """The class's ``canonicalKey``: trim, map the empty string, then the alias table.

    An unmatched name is returned trimmed and unchanged, not case-folded - the
    class's own behaviour, which is why ``canonicalKey("Barite")`` is ``BaSO4`` but
    ``canonicalKey("Gypsumx")`` keeps its spelling.
    """
    trimmed = deposit.strip()
    if not trimmed:
        return "UNKNOWN"
    return ALIASES.get(trimmed.lower(), trimmed)


@cache
def _registry() -> tuple[dict[str, Any], ...]:
    """The embedded knowledge base, parsed the class's way."""
    rows: list[dict[str, Any]] = []
    with find(REGISTRY_PATH).open(encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("|")
            if len(parts) < 10:
                continue
            rows.append(
                {
                    "scale_type": parts[0].strip(),
                    "category": parts[1].strip(),
                    "dissolver": parts[2].strip(),
                    "concentration": parts[3].strip(),
                    "method": parts[4].strip(),
                    "effectiveness": parts[5].strip(),
                    "min_temperature_c": _parse_double(parts[6].strip(), 0.0),
                    "max_temperature_c": _parse_double(parts[7].strip(), 200.0),
                    "cautions": parts[8].strip(),
                    "standard": parts[9].strip(),
                }
            )
    return tuple(rows)


def _parse_double(value: str, fallback: float) -> float:
    try:
        return float(value)
    except ValueError:
        return fallback


def _effectiveness_rank(effectiveness: str) -> int:
    lowered = effectiveness.lower()
    return {"high": 3, "medium": 2, "low": 1}.get(lowered, 0)


def recommend_for(deposit: str) -> list[dict[str, Any]]:
    """The class's ``recommendFor``: the key's rows, or the ``UNKNOWN`` fallback."""
    key = canonical_key(deposit)
    matches = [row for row in _registry() if row["scale_type"].lower() == key.lower()]
    if not matches:
        matches = [row for row in _registry() if row["scale_type"].lower() == "unknown"]
    # `sorted` is stable, as `List.sort` is: equal ranks keep the CSV's order.
    return sorted(matches, key=lambda row: _effectiveness_rank(row["effectiveness"]), reverse=True)


def scale_remediation_advisor(deposit: str) -> ScaleRemediationAdvisorResult:
    """The candidate remedies for a deposit.

    Args:
        deposit: the deposit's name or formula; aliases are mapped, and an unmatched
            name keeps its trimmed spelling.

    Returns:
        The canonical key and the matching options, ordered by descending effectiveness.

    Example:
        >>> r = scale_remediation_advisor("calcite")
        >>> r.canonical_key
        'CaCO3'
        >>> r.options[0].effectiveness
        'HIGH'
    """
    key = canonical_key(deposit)
    rows = recommend_for(key)
    return ScaleRemediationAdvisorResult(
        canonical_key=key,
        options=tuple(
            RemediationOptionRow(
                scale_type=row["scale_type"],
                category=row["category"],
                dissolver=row["dissolver"],
                concentration=row["concentration"],
                method=row["method"],
                effectiveness=row["effectiveness"],
                min_temperature=from_si(row["min_temperature_c"] + 273.15, "K"),
                max_temperature=from_si(row["max_temperature_c"] + 273.15, "K"),
                cautions=row["cautions"],
                standard=row["standard"],
            )
            for row in rows
        ),
        warnings=(),
    )
