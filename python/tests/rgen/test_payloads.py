"""The on-chain payloads: integer-only, the four quality keys, and the anchor chain."""

from __future__ import annotations

from typing import Any

import pytest

from azoth.rgen import canonical, payloads

HASHES = {
    "definition": "a" * 64,
    "build": "b" * 64,
    "statement": "c" * 64,
    "cases": "d" * 64,
    "agreement": "e" * 64,
    "quality": "f" * 64,
    "validation-evidence": "0" * 64,
}


def _claim(
    *,
    quantity: int = 5,
    deadline_height: int = 0,
    comment_until_height: int | None = None,
) -> dict[str, Any]:
    return payloads.claim_request(
        calc_id="eos.antoine_vapor_pressure",
        definition_hash=HASHES["definition"],
        build_hash=HASHES["build"],
        statement_hash=HASHES["statement"],
        quantity=quantity,
        project="azoth-demo",
        jurisdiction="US-TX",
        deadline_height=deadline_height,
        comment_until_height=comment_until_height,
    )


def test_a_claim_carries_the_quantity_as_an_integer_object() -> None:
    request = _claim()
    assert request["quantity"] == {"unit": "validated-cases", "value": 5}
    assert isinstance(request["quantity"]["value"], int)
    assert request["deadlineHeight"] == 0
    assert "commentUntilHeight" not in request
    payloads.integers_only(request)
    assert request["selectionCriteriaHash"] == canonical.hash_value(request["selectionCriteria"])


def test_a_zero_quantity_is_refused() -> None:
    with pytest.raises(canonical.PayloadError):
        _claim(quantity=0)


def test_a_comment_window_beyond_the_deadline_is_refused() -> None:
    with pytest.raises(canonical.PayloadError):
        _claim(deadline_height=10, comment_until_height=11)


def test_a_float_anywhere_in_a_payload_is_refused() -> None:
    with pytest.raises(canonical.PayloadError):
        payloads.integers_only({"quantity": {"value": 1.5}})


def test_none_is_refused_because_optional_keys_are_omitted() -> None:
    with pytest.raises(canonical.PayloadError):
        payloads.integers_only({"commentUntilHeight": None})


def test_quality_must_be_exactly_the_four_keys() -> None:
    with pytest.raises(canonical.PayloadError):
        payloads.validation_request(
            verdict="valid",
            quality={"completeness": 5, "correctness": 5, "traceability": 5},
            characteristics_hash="x",
            agreement_hash="y",
            quality_hash="z",
            cases_hash="w",
            comment="c",
        )


def test_quality_values_must_be_in_range() -> None:
    with pytest.raises(canonical.PayloadError):
        payloads.validation_request(
            verdict="valid",
            quality={
                "completeness": 5,
                "correctness": 6,
                "traceability": 5,
                "fitness": 5,
            },
            characteristics_hash="x",
            agreement_hash="y",
            quality_hash="z",
            cases_hash="w",
            comment="c",
        )


def test_the_verdict_vocabulary_is_the_contracts() -> None:
    with pytest.raises(canonical.PayloadError):
        payloads.validation_request(
            verdict="maybe",
            quality={"completeness": 5, "correctness": 5, "traceability": 5, "fitness": 5},
            characteristics_hash="x",
            agreement_hash="y",
            quality_hash="z",
            cases_hash="w",
            comment="c",
        )


def test_the_anchor_chain_is_seven_and_parent_before_child() -> None:
    chain = payloads.anchors(HASHES)
    assert len(chain) == 7
    assert chain[0]["parentHash"] == ""
    seen = {""}
    for anchor in chain:
        assert anchor["parentHash"] in seen
        seen.add(anchor["hash"])
    assert [anchor["kind"] for anchor in chain] == [
        "azoth:definition",
        "azoth:build",
        "azoth:statement",
        "azoth:cases",
        "azoth:agreement",
        "azoth:quality",
        "azoth:validation",
    ]
