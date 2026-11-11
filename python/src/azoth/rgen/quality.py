"""The four scores r-gen's validation carries, derived from recorded evidence.

The scores are the validator's stated judgement; every count they were derived from is
written beside them, so a consumer can recompute the judgement rather than take it. The
payload carries only the four keys ``kasset.conventions/v1`` names; the derivation lives
in the quality document.

Correctness is measured against the declared cases' expected values - the external oracle
runs where the calc has them - not against the other implementation. Two same-spec
implementations agreeing is completeness, not correctness.
"""

from __future__ import annotations

from typing import Any, Final

__all__ = ["VERDICTS", "score", "verdict"]

VERDICTS: Final[tuple[str, ...]] = ("valid", "partial", "invalid")


def _band(numerator: int, denominator: int) -> int:
    """A 1..5 score for a count against its total, guarding the empty total."""
    if denominator <= 0:
        return 1
    ratio = numerator / denominator
    if ratio >= 1.0:
        return 5
    if ratio >= 0.75:
        return 4
    if ratio >= 0.5:
        return 3
    if ratio >= 0.25:
        return 2
    return 1


def _traceability(items: int, tracked_dirty: bool | None) -> int:
    score = max(1, min(5, items))
    if tracked_dirty:
        score = min(score, 3)
    return score


def score(evidence: dict[str, Any]) -> dict[str, Any]:
    """The four scores and the evidence they came from."""
    declared = int(evidence["declared"])
    replayed = int(evidence["replayed"])
    passed = int(evidence["passed"])
    warnings = int(evidence["warnings"])

    completeness = _band(replayed, declared)
    correctness = _band(passed, declared)
    if not bool(evidence["run_matched"]):
        correctness = min(correctness, 3)
    if bool(evidence["run_disagree"]):
        correctness = 1
    traceability = _traceability(int(evidence["traceability_items"]), evidence.get("tracked_dirty"))
    fitness = 5 if warnings == 0 else (4 if warnings == 1 else 3)

    record = {
        "completeness": completeness,
        "correctness": correctness,
        "traceability": traceability,
        "fitness": fitness,
    }
    return {
        "quality": record,
        "evidence": {
            "declared": declared,
            "replayed": replayed,
            "passed": passed,
            "warnings": warnings,
            "run_matched": bool(evidence["run_matched"]),
            "run_disagree": bool(evidence["run_disagree"]),
            "traceability_items": int(evidence["traceability_items"]),
            "tracked_dirty": evidence.get("tracked_dirty"),
        },
        "verdict": verdict(
            run_disagree=bool(evidence["run_disagree"]),
            case_failures=declared - passed,
            warnings=warnings,
            run_matched=bool(evidence["run_matched"]),
        ),
    }


def verdict(*, run_disagree: bool, case_failures: int, warnings: int, run_matched: bool) -> str:
    """The validation's verdict from what the evidence says."""
    if run_disagree or case_failures > 0:
        return "invalid"
    if warnings > 0 or not run_matched:
        return "partial"
    return "valid"
