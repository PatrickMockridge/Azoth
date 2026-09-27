"""The provenance block's own arithmetic.

The generated half is data and is checked where it is generated. What is tested here is
the merge: which warnings become ``skipped_checks``, in what order, and what the block
says when a call was clean. Those are the parts a reader relies on to tell "checked and
fine" from "never checked", so they are worth pinning directly rather than only through
a calculation that happens to exercise them.
"""

from __future__ import annotations

from typing import Any

from azoth.core.provenance import FRAMEWORK_DOCUMENT_KEYS, PROVENANCE_KEY, Provenance
from azoth.core.warnings import Warning, WarningCode

STATIC: dict[str, Any] = {
    "calc_id": "hydraulics.darcy_weisbach",
    "name": "Darcy-Weisbach pressure drop",
    "spec": {"path": "specs/calcs/hydraulics/darcy_weisbach.toml", "sha256": "a" * 64},
    "code": [
        {"path": "python/src/azoth/hydraulics/reference/darcy_weisbach.py", "sha256": "b" * 64},
        {"path": "crates/azoth-hydraulics/src/darcy_weisbach.rs", "sha256": "c" * 64},
    ],
    "source": "Crane TP-410",
    "verification": "source_needed",
    "validation_cases": 1,
    "tests_active": 2,
    "tests_skipped": 2,
}


def warning(code: WarningCode, field: str | None = None) -> Warning:
    return Warning(code=code, message="a message", field=field)


def test_a_clean_call_says_so_and_names_no_skipped_check() -> None:
    block = Provenance.of(STATIC, [])

    assert block.clean
    assert block.skipped_checks == ()
    assert block.warning_codes == ()


def test_the_static_half_passes_through_unchanged() -> None:
    block = Provenance.of(STATIC, [])

    assert block.calc_id == STATIC["calc_id"]
    assert block.spec_sha256 == "a" * 64
    assert block.rust_path == STATIC["code"][1]["path"]
    assert block.verification == "source_needed"
    assert (block.validation_cases, block.tests_active, block.tests_skipped) == (1, 2, 2)


def test_only_a_skipped_check_naming_a_field_becomes_a_skipped_check() -> None:
    """Three exclusions in one case, because they are the same rule read three ways.

    A warning whose code does not mean "a check did not run" is not a skipped check; a
    warning that names no field cannot go in a tuple of field names; and two warnings
    about the same field are one skipped check, not two.
    """
    block = Provenance.of(
        STATIC,
        [
            warning(WarningCode.RANGE_CHECK_SKIPPED, "re"),
            warning(WarningCode.RANGE_CHECK_SKIPPED, "re"),
            warning(WarningCode.RANGE_CHECK_SKIPPED, None),
            warning(WarningCode.OUT_OF_VALID_RANGE, "L"),
        ],
    )

    assert block.skipped_checks == ("re",)
    assert not block.clean


def test_skipped_checks_are_ordered_by_name_and_warning_codes_by_the_vocabulary() -> None:
    """Both are ordered by something other than the order the call raised them.

    Otherwise two runs that raised the same warnings in a different order would produce
    different blocks, and the cross-language test that compares them key for key would
    be comparing the order of a warning list rather than the facts in it.
    """
    block = Provenance.of(
        STATIC,
        [
            warning(WarningCode.OUT_OF_VALID_RANGE, "v"),
            warning(WarningCode.RANGE_CHECK_SKIPPED, "v"),
            warning(WarningCode.RANGE_CHECK_SKIPPED, "D"),
        ],
    )

    assert block.skipped_checks == ("D", "v")
    assert block.warning_codes == ("OUT_OF_VALID_RANGE", "RANGE_CHECK_SKIPPED")


def test_the_document_key_is_exported_rather_than_repeated() -> None:
    """The key a result document carries the block under, and the list of such keys.

    The test that asserts a result's JSON shape needs both, and it must not hold its own
    copy: a second list is a second thing to keep in step.
    """
    assert PROVENANCE_KEY == "provenance"
    assert FRAMEWORK_DOCUMENT_KEYS == (PROVENANCE_KEY,)
