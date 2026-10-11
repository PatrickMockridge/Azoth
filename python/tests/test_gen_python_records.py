"""The generated record rows are current, and their shape is the one pyo3 can extract.

`tools/gen_python_records.py` emits `crates/azoth-python/src/records_gen.rs` - one
`#[derive(FromPyObject)]` struct per `record_list` input's `row_type`, from the spec's own `fields`
block. `docs-drift` regenerates and diffs it; the `--check` below runs the same gate in the suite.

**The two shape rules are the reason this file exists rather than only the drift check.** pyo3's
default extraction is positional, so a struct without `#[pyo3(item)]` per field would be found by
the derive and yield *no* fields - and
`test_registration_completeness.py::test_the_record_row_is_one_shape_in_three_declarations` reads
the source with `#[pyo3(item)]\\s*pub (\\w+):`, so the failure it produces is a field-name mismatch
rather than a missing attribute. The other is that the derive stays on one line: the same test's
`#[derive(...FromPyObject...)]` pattern does not match a split attribute.
"""

from __future__ import annotations

import importlib
import re
import subprocess
import sys
from pathlib import Path
from types import ModuleType
from typing import Any

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
OUT = REPO_ROOT / "crates" / "azoth-python" / "src" / "records_gen.rs"


def gen_records() -> ModuleType:
    sys.path.insert(0, str(REPO_ROOT / "tools"))
    try:
        return importlib.import_module("gen_python_records")
    finally:
        sys.path.pop(0)


def test_the_committed_rows_are_current() -> None:
    """The drift gate `docs-drift` runs, run here so it fails in the suite too."""
    result = subprocess.run(
        [sys.executable, str(REPO_ROOT / "tools" / "gen_python_records.py"), "--check"],
        capture_output=True,
        text=True,
        check=False,
    )

    assert result.returncode == 0, (
        f"the record rows are out of date - run `python tools/gen_python_records.py`\n"
        f"{result.stdout}{result.stderr}"
    )


@pytest.mark.parametrize(
    ("declaration", "expected"),
    [
        ({"unit": "K"}, "f64"),
        ({"type": "quantity", "unit": "kg/mol"}, "f64"),
        ({"unit": "dimensionless", "optional": True}, "Option<f64>"),
        ({"type": "string"}, "String"),
        ({"type": "enum", "optional": True}, "Option<String>"),
        ({"type": "boolean"}, "bool"),
    ],
)
def test_a_field_extracts_into_the_magnitude_it_crosses_as(
    declaration: dict[str, Any], expected: str
) -> None:
    """**A dimensioned field is an `f64` here whatever unit the spec declares.**

    Python converts each field to its SI magnitude before the boundary - `input_records_to_si`,
    the one conversion site - so the struct carries numbers and the spec carries the units.
    """
    tool: Any = gen_records()
    assert tool.field_type(declaration) == expected


def test_every_row_field_is_extracted_by_key_and_the_derive_is_on_one_line() -> None:
    """The two rules pyo3 and the agreement test between them require.

    A bare `pub` field is found by the derive and extracts *positionally*, which is a different
    row rather than a compile error - so the guard has to be on the attribute and not on the
    struct's existence.
    """
    text = OUT.read_text(encoding="utf-8")
    structs = re.findall(r"#\[derive\(([^)]*)\)\]\npub struct (\w+) \{(.*?)^\}", text, re.M | re.S)
    assert structs, "no row struct was generated, so there is nothing to hold"

    for derives, name, body in structs:
        assert "FromPyObject" in derives, f"`{name}` does not derive `FromPyObject`"
        fields = re.findall(r"^\s*pub (\w+):", body, re.M)
        keyed = re.findall(r"#\[pyo3\(item\)\]\s*pub (\w+):", body)
        assert keyed == fields, (
            f"`{name}` extracts {len(keyed)} of its {len(fields)} fields by key; the rest would "
            f"be extracted positionally"
        )
