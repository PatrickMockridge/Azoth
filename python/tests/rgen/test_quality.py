"""The four scores and the verdict, over the evidence record they are derived from."""

from __future__ import annotations

from typing import Any

from azoth.rgen import quality


def _evidence(**overrides: Any) -> dict[str, Any]:
    base: dict[str, Any] = {
        "declared": 5,
        "replayed": 5,
        "passed": 5,
        "warnings": 0,
        "run_matched": True,
        "run_disagree": False,
        "traceability_items": 5,
        "tracked_dirty": False,
    }
    base.update(overrides)
    return base


def test_a_clean_matched_run_scores_five_everywhere() -> None:
    scored = quality.score(_evidence())
    assert scored["quality"] == {
        "completeness": 5,
        "correctness": 5,
        "traceability": 5,
        "fitness": 5,
    }
    assert scored["verdict"] == "valid"


def test_a_partial_case_set_lowers_both_counts() -> None:
    scored = quality.score(_evidence(passed=4, replayed=4))
    assert scored["quality"]["correctness"] == 4
    assert scored["quality"]["completeness"] == 4
    assert scored["verdict"] == "invalid"


def test_a_dirty_tree_caps_traceability() -> None:
    scored = quality.score(_evidence(tracked_dirty=True))
    assert scored["quality"]["traceability"] == 3


def test_a_missing_provenance_item_lowers_traceability() -> None:
    scored = quality.score(_evidence(traceability_items=3))
    assert scored["quality"]["traceability"] == 3


def test_an_exploratory_run_cannot_score_above_three_on_correctness() -> None:
    scored = quality.score(_evidence(run_matched=False))
    assert scored["quality"]["correctness"] == 3
    assert scored["verdict"] == "partial"


def test_a_disagreeing_run_is_wrong_and_invalid() -> None:
    scored = quality.score(_evidence(run_disagree=True))
    assert scored["quality"]["correctness"] == 1
    assert scored["verdict"] == "invalid"


def test_warnings_lower_fitness_and_make_the_verdict_partial() -> None:
    one = quality.score(_evidence(warnings=1))
    assert one["quality"]["fitness"] == 4
    assert one["verdict"] == "partial"
    two = quality.score(_evidence(warnings=2))
    assert two["quality"]["fitness"] == 3


def test_an_empty_case_set_is_one_not_a_division() -> None:
    scored = quality.score(_evidence(declared=0, replayed=0, passed=0))
    assert scored["quality"]["completeness"] == 1
    assert scored["quality"]["correctness"] == 1


def test_the_payload_quality_has_exactly_the_four_keys() -> None:
    scored = quality.score(_evidence())
    assert set(scored["quality"]) == {"completeness", "correctness", "traceability", "fitness"}
