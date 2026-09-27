"""The provenance block's own arithmetic.

The generated half is data and is checked where it is generated. What is tested here is
the merge: which warnings become ``skipped_checks``, in what order, and what the block
says when a call was clean. Those are the parts a reader relies on to tell "checked and
fine" from "never checked", so they are worth pinning directly rather than only through
a calculation that happens to exercise them.
"""

from __future__ import annotations

import importlib
import re
import sys
import tomllib
from pathlib import Path
from types import ModuleType
from typing import Any

import pytest

from azoth._provenance_gen import PROVENANCE, provenance
from azoth.core.provenance import FRAMEWORK_DOCUMENT_KEYS, PROVENANCE_KEY, Provenance
from azoth.core.warnings import Warning, WarningCode

REPO_ROOT = Path(__file__).resolve().parents[2]
TOOLS = REPO_ROOT / "tools"

#: A SHA-256, lowercase hex. A hash that is not this is not a hash.
_SHA256 = re.compile(r"^[0-9a-f]{64}$")


def provenance_tool() -> ModuleType:
    """`tools/provenance.py`, imported by name - the tools are not a package."""
    sys.path.insert(0, str(TOOLS))
    try:
        return importlib.import_module("provenance")
    finally:
        sys.path.pop(0)


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


def test_the_block_exists_for_every_registered_id() -> None:
    """Every calc and every model, with no id left able to produce a bare number."""
    from azoth._models_gen import MODELS
    from azoth._registry_gen import CALCS

    registered = {entry["id"] for entry in [*CALCS, *MODELS]}

    assert set(PROVENANCE) == registered
    assert len(registered) == 192


def test_the_rust_table_carries_exactly_the_same_ids() -> None:
    """Two generated files, one list - checked, because they are two files.

    Both are emitted from the same blocks, so they agree unless a language's emitter drops
    one. Nothing else would notice: the Rust side compiles with a short table and the
    Python side passes with a long one, and the id that fell out would be found only by a
    caller who happened to ask for it.
    """
    text = (REPO_ROOT / "crates" / "azoth-core" / "src" / "provenance_gen.rs").read_text(
        encoding="utf-8"
    )
    rust_ids = re.findall(r'^\s*calc_id: "([^"]+)",$', text, flags=re.MULTILINE)

    assert len(rust_ids) == len(PROVENANCE), f"{len(rust_ids)} vs {len(PROVENANCE)}"
    assert set(rust_ids) == set(PROVENANCE)


def test_an_unknown_id_is_refused_rather_than_answered_with_an_empty_block() -> None:
    """The same rule the databank follows for a substance it cannot name.

    A block that existed but said nothing would be worse than an error: it would let a
    number be presented as if it carried its provenance when it carried a placeholder.
    """
    with pytest.raises(KeyError, match="no provenance for"):
        provenance("hydraulics.not_a_calculation")


def test_every_block_is_complete_and_well_formed() -> None:
    """The shape, for all 192, so the contract is the block's rather than one sample's."""
    for calc_id, block in PROVENANCE.items():
        assert set(block) == {
            "calc_id",
            "name",
            "spec",
            "code",
            "source",
            "verification",
            "validation_cases",
            "tests_active",
            "tests_skipped",
        }, calc_id
        assert block["calc_id"] == calc_id, calc_id
        assert block["name"], calc_id
        assert block["source"], f"{calc_id} names no source"
        assert block["verification"] in provenance_tool().VERIFICATION_STATUSES, calc_id
        assert len(block["code"]) == 2, f"{calc_id} does not name both implementations"


def test_every_embedded_hash_is_the_hash_of_the_file_it_names() -> None:
    """**The claim the block makes, checked rather than asserted.**

    A block whose hash was of something else would be worse than no hash at all: it
    would let a reader believe they had identified the code that produced a number when
    they had identified nothing. The hash is only worth carrying because it can be
    compared, so this is the test that makes it a fact about the tree.
    """
    provenance = provenance_tool()
    for calc_id, block in PROVENANCE.items():
        named = [block["spec"], *block["code"]]
        for item in named:
            assert _SHA256.match(item["sha256"]), f"{calc_id}: {item['sha256']!r} is not a hash"
            actual = provenance.sha256_of(REPO_ROOT / item["path"])
            assert actual == item["sha256"], (
                f"{calc_id}: {item['path']} is recorded as {item['sha256'][:12]} "
                f"but hashes to {actual[:12]} - regenerate with `tools/gen_registry.py`"
            )


def test_the_recorded_status_is_the_one_the_tree_derives() -> None:
    """The status is recomputed from the same evidence rather than trusted.

    Called through the tool's own `derive_verification` rather than by restating the
    cascade here, so a change to the rule cannot leave this test passing against a table
    built by the previous one.
    """
    tool = provenance_tool()
    index = tool.validation_cases()
    for tree, entry_of in (
        ("specs/calcs", tool.calc_entry),
        ("specs/models", tool.model_entry),
    ):
        for path in sorted((REPO_ROOT / tree).rglob("*.toml")):
            spec = tomllib.loads(path.read_text(encoding="utf-8"))
            record = entry_of(spec, index)
            block = PROVENANCE[spec["id"]]
            assert block["verification"] == record["verification"], spec["id"]
            assert block["validation_cases"] == record["validation_cases"], spec["id"]
            assert block["tests_active"] == record["tests_active"], spec["id"]
            assert block["tests_skipped"] == record["tests_skipped"], spec["id"]
