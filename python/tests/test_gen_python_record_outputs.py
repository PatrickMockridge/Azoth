"""The generated record-output rows are current, and their shape is the one pyo3 hands back.

`tools/gen_python_record_outputs.py` emits `crates/azoth-python/src/record_outputs_gen.rs` and
`python/src/azoth/core/record_outputs_gen.py` - a `#[pyclass]` and a frozen dataclass per
`record_list` output, from the spec's own `fields` block. `docs-drift` regenerates and diffs them;
the `--check` below runs the same gate in the suite.

**This is the mirror of `test_gen_python_records.py`, and the direction is why the file exists.**
An input row derives `FromPyObject` and is built in Python; an output row is a `#[pyclass]` built in
Rust, so the rules that matter are the opposite ones: the field *names* have to be the spec's, the
conversion has to be the one a `PyQty` carries, and a renamed kernel field has to fail rather than
quietly transport a number.
"""

from __future__ import annotations

import importlib
import re
import subprocess
import sys
from pathlib import Path
from types import ModuleType

REPO_ROOT = Path(__file__).resolve().parents[2]
RUST_OUT = REPO_ROOT / "crates" / "azoth-python" / "src" / "record_outputs_gen.rs"
PYTHON_OUT = REPO_ROOT / "python" / "src" / "azoth" / "core" / "record_outputs_gen.py"


def gen_outputs() -> ModuleType:
    sys.path.insert(0, str(REPO_ROOT / "tools"))
    try:
        return importlib.import_module("gen_python_record_outputs")
    finally:
        sys.path.pop(0)


def test_the_committed_rows_are_current() -> None:
    """The drift gate `docs-drift` runs, run here so it fails in the suite too."""
    result = subprocess.run(
        [sys.executable, str(REPO_ROOT / "tools" / "gen_python_record_outputs.py"), "--check"],
        capture_output=True,
        text=True,
        check=False,
    )

    assert result.returncode == 0, (
        f"the record output rows are out of date - run "
        f"`python tools/gen_python_record_outputs.py`\n{result.stdout}{result.stderr}"
    )


def test_every_row_is_a_frozen_pyclass_built_only_in_rust() -> None:
    """`frozen` and `skip_from_py_object`, and no `#[new]`.

    A row is produced by a calculation and read by a caller, and the same reasoning that makes a
    result frozen makes it non-constructible: a hand-built row would be a value with nothing behind
    it. `skip_from_py_object` is what stops one crossing back in as an argument.
    """
    text = RUST_OUT.read_text(encoding="utf-8")
    structs = re.findall(r"#\[pyclass\((.*?)\)\]\npub struct (\w+) \{", text, re.S)
    assert structs, "no record output pyclass was generated, so there is nothing to hold"

    for attributes, name in structs:
        assert "frozen" in attributes, f"`{name}` is not frozen"
        assert "skip_from_py_object" in attributes, f"`{name}` can cross back into Rust"
        assert "#[pyo3(get)]" in text, "no field is readable from Python"
    assert "#[new]" not in text, "a record output row is not built from Python"


def test_every_pyclass_field_is_readable_and_the_conversion_is_a_function() -> None:
    """A `#[pyo3(get)]` per field, and a `From` impl whose keys are the fields.

    **A path rather than a closure is load-bearing**: `clippy -D warnings` refuses
    `redundant_closure`, and the transport calls `PyRow::from` over the kernel vector.
    """
    text = RUST_OUT.read_text(encoding="utf-8")
    for match in re.finditer(
        r"#\[pyclass\((.*?)\)\]\npub struct (Py\w+) \{(.*?)^\}", text, re.M | re.S
    ):
        _, name, body = match.groups()
        fields = re.findall(r"^\s*pub (\w+):", body, re.M)
        readable = re.findall(r"#\[pyo3\(get\)\]\s*pub (\w+):", body)
        assert readable == fields, (
            f"`{name}` exposes {len(readable)} of its {len(fields)} fields; the rest would be "
            f"attributes no caller can read"
        )

    # **Across a line break**: rustfmt wraps a long source path before `for`, and a
    # pattern that required the impl on one line read a generated file as incomplete.
    conversions = re.findall(r"impl From<&[\w:]+>\s+for (Py\w+)\s*\{", text)
    assert len(conversions) == len(re.findall(r"#\[pyclass\(", text)), (
        "every pyclass needs a `From` impl, which is what the transport calls"
    )


def test_the_dataclass_is_frozen_and_carries_no_default() -> None:
    """**No `= None` on an optional field.** The bridge supplies every field explicitly.

    A default would invite a caller to omit one and receive a row the kernel never produced, and
    it would also make the field order load-bearing in a way the spec does not promise: a
    defaulted field ahead of a required one is a `TypeError` when the class is created.
    """
    text = PYTHON_OUT.read_text(encoding="utf-8")
    classes = re.findall(r"@dataclass\((.*?)\)\nclass (\w+):", text)
    assert classes, "no record output dataclass was generated, so there is nothing to hold"

    for decorator, name in classes:
        assert "frozen=True" in decorator, f"`{name}` is not frozen"
        assert "eq=False" in decorator, f"`{name}` compares `pint` quantities with `==`"
    assert "= None" not in text, "an optional field carries a default the bridge always overrides"


def test_the_quantity_import_is_emitted_exactly_where_a_row_needs_one() -> None:
    """`Q` is imported when a row has a dimensioned field, and not otherwise.

    The same rule the transport's `PyQty` follows: a dimensionless quantity crosses as a bare
    number, so an unconditional import would be an `unused_import` ruff refuses, and a missing one
    would be a `NameError` raised while a generated module is imported.
    """
    rows = gen_outputs().record_output_rows()
    dimensioned = any(
        declared.get("unit") not in (None, "dimensionless")
        for declaration in rows.values()
        for declared in declaration["fields"].values()
        if declared.get("type", "quantity") not in ("boolean", "string", "enum")
    )
    text = PYTHON_OUT.read_text(encoding="utf-8")
    assert ("from azoth.core.units import Q" in text) == dimensioned
