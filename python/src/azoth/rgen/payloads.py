"""The two on-chain payloads and the anchor list, built from a bundle's hashes.

r-gen's wire is integer-only: its renderer refuses a float, and the contracts store what
they are given. The builders therefore write integers and strings only, omit optional keys
rather than sending ``null``, and check that every value is representable before returning.
``deadlineHeight`` and ``credentialUri`` are deliberate placeholders the deployer fills in;
a payload deployed with them still present is refused by the contract.
"""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any, Final

from .canonical import PayloadError, hash_value

__all__ = [
    "QUALITY_KEYS",
    "anchors",
    "claim_request",
    "comment_text",
    "validation_request",
]

QUALITY_KEYS: Final[tuple[str, ...]] = (
    "completeness",
    "correctness",
    "traceability",
    "fitness",
)

_QUANTITY_UNIT: Final[str] = "validated-cases"


def integers_only(value: Any, path: str = "payload") -> None:
    """Refuse anything r-gen's renderer cannot carry, naming the path."""
    if isinstance(value, (bool, int, str)):
        return
    if isinstance(value, Mapping):
        for key, item in value.items():
            if not isinstance(key, str):
                raise PayloadError(f"{path}: key {key!r} is not a string")
            integers_only(item, f"{path}.{key}")
        return
    if isinstance(value, list):
        for index, item in enumerate(value):
            integers_only(item, f"{path}[{index}]")
        return
    raise PayloadError(f"{path}: {type(value).__name__} has no r-gen wire form")


def _check_quality(quality: Mapping[str, Any]) -> dict[str, int]:
    if set(quality) != set(QUALITY_KEYS):
        raise PayloadError(f"quality must be exactly {list(QUALITY_KEYS)}, got {sorted(quality)}")
    checked: dict[str, int] = {}
    for key in QUALITY_KEYS:
        value = quality[key]
        if isinstance(value, bool) or not isinstance(value, int):
            raise PayloadError(f"quality.{key} must be an integer")
        if not 1 <= value <= 5:
            raise PayloadError(f"quality.{key} must be 1..5, got {value}")
        checked[key] = value
    if not 4 <= sum(checked.values()) <= 20:
        raise PayloadError("quality must sum to 4..20")
    return checked


def comment_text(calc_id: str, verdict: str, quantity: int) -> str:
    """The validation's comment, deterministic and carrying the count."""
    return f"azoth: {calc_id}; {verdict}; {quantity} validated case(s)"


def claim_request(
    *,
    calc_id: str,
    definition_hash: str,
    build_hash: str,
    statement_hash: str,
    quantity: int,
    project: str,
    jurisdiction: str,
    discipline: str = "process",
    ctype: str = "Validator",
    role: str = "validator",
    level: int = 1,
    deadline_height: int = 0,
    comment_until_height: int | None = None,
) -> dict[str, Any]:
    """The Claim payload. The quantity counts validated cases, never a float result."""
    if quantity < 1:
        raise PayloadError(
            "quantity must be at least 1; a claim with no validated case has nothing to count"
        )
    if deadline_height < 0 or (comment_until_height is not None and comment_until_height < 0):
        raise PayloadError("heights must be non-negative")
    if comment_until_height is not None and comment_until_height > deadline_height:
        raise PayloadError("commentUntilHeight must not exceed deadlineHeight")

    work = {
        "titleHash": hash_value(f"azoth:{calc_id}"),
        "discipline": discipline,
        "jurisdiction": jurisdiction,
        "projectRef": project,
    }
    criteria = {
        "ctype": ctype,
        "discipline": discipline,
        "jurisdiction": jurisdiction,
        "role": role,
        "minLevel": level,
    }
    request: dict[str, Any] = {
        "work": work,
        "selectionCriteria": criteria,
        "selectionCriteriaHash": hash_value(criteria),
        "statementHash": statement_hash,
        "definitionHash": definition_hash,
        "quantity": {"unit": _QUANTITY_UNIT, "value": quantity},
        "evidence": [
            {"kind": "azoth:statement", "hash": statement_hash},
            {"kind": "azoth:definition", "hash": definition_hash},
            {"kind": "azoth:build", "hash": build_hash},
        ],
        "deadlineHeight": deadline_height,
    }
    if comment_until_height is not None:
        request["commentUntilHeight"] = comment_until_height
    integers_only(request)
    return request


def validation_request(
    *,
    verdict: str,
    quality: Mapping[str, Any],
    characteristics_hash: str,
    agreement_hash: str,
    quality_hash: str,
    cases_hash: str,
    comment: str,
    credential_uri: str = "",
) -> dict[str, Any]:
    """The Validation payload; ``credentialUri`` is the deployer's to fill."""
    if verdict not in ("valid", "partial", "invalid"):
        raise PayloadError(f"verdict must be valid, partial or invalid, got {verdict!r}")
    request: dict[str, Any] = {
        "verdict": verdict,
        "quality": _check_quality(quality),
        "credentialUri": credential_uri,
        "characteristicsHash": characteristics_hash,
        "commentHash": hash_value(comment),
        "evidence": [
            {"kind": "azoth:validation", "hash": characteristics_hash},
            {"kind": "azoth:agreement", "hash": agreement_hash},
            {"kind": "azoth:quality", "hash": quality_hash},
            {"kind": "azoth:cases", "hash": cases_hash},
        ],
    }
    integers_only(request)
    return request


def anchors(hashes: Mapping[str, str]) -> list[dict[str, str]]:
    """The seven anchors, parent before child."""
    chain = (
        ("azoth:definition", "definition", ""),
        ("azoth:build", "build", "definition"),
        ("azoth:statement", "statement", "build"),
        ("azoth:cases", "cases", "statement"),
        ("azoth:agreement", "agreement", "cases"),
        ("azoth:quality", "quality", "agreement"),
        ("azoth:validation", "validation-evidence", "quality"),
    )
    out = [
        {
            "kind": kind,
            "hash": hashes[key],
            "parentHash": "" if parent == "" else hashes[parent],
        }
        for kind, key, parent in chain
    ]
    integers_only(out)
    return out
