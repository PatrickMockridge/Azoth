"""`spec_lint`'s palette rules, and the schema that had never existed.

The palette was the one spec tree with no schema, so its `notes` field grew to five
thousand characters with nothing reading it while `docs/src/architecture/spec-files.md`
stated the format as "a spec file declares ... no prose fields". These tests pin the
three rules that now hold it, and one more that has no runtime effect at all: that the
schema and the Rust struct it mirrors still describe the same document.

That last one is the price of the schema being a second description of a shape whose
first description is `struct UnitOpSpec`. `unit.schema.json` is *generated* from the
vocabulary and cannot drift; this one is hand-written, so a test is what stands in for
the generator.
"""

from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path
from typing import Any

from _toml_fixture import dump

REPO_ROOT = Path(__file__).resolve().parents[2]
SPEC_LINT = REPO_ROOT / "tools" / "spec_lint.py"
SCHEMA = REPO_ROOT / "specs" / "schema" / "unit_ops.schema.json"
UNIT_OP_RS = REPO_ROOT / "crates" / "azoth-process" / "src" / "unit_op.rs"

#: A palette entry that lints clean, and the baseline every mutation starts from.
ENTRY: dict[str, Any] = {
    "id": "unit_ops.widget",
    "name": "Widget",
    "source": {"standard": "NeqSim process/equipment/Widget.java"},
    "parameters": {"duty": {"unit": "W", "required": True, "description": "the duty"}},
    "ports": [
        {
            "name": "inlet",
            "direction": "in",
            "fields": {
                "n": {"dimension": "molar_flow"},
                "z": {"dimension": "dimensionless", "shape": "vector"},
                "P": {"dimension": "pressure"},
                "T": {"dimension": "thermodynamic_temperature"},
                "h": {"dimension": "molar_energy"},
            },
        },
        {
            "name": "outlet",
            "direction": "out",
            "fields": {
                "n": {"dimension": "molar_flow"},
                "P": {"dimension": "pressure"},
            },
        },
    ],
}


def lint(palette_dir: Path) -> subprocess.CompletedProcess[str]:
    """Run spec_lint over a synthetic palette tree."""
    return subprocess.run(
        [sys.executable, str(SPEC_LINT), "--palette-dir", str(palette_dir), "--quiet"],
        capture_output=True,
        text=True,
        check=False,
    )


def palette(tmp_path: Path, **changes: object) -> Path:
    entry = {**ENTRY, **changes}
    directory = tmp_path / "utility"
    directory.mkdir(parents=True, exist_ok=True)
    (directory / "widget.toml").write_text(dump(entry), encoding="utf-8")
    return tmp_path


def test_a_well_formed_entry_lints_clean(tmp_path: Path) -> None:
    """The baseline. Every other case here is a mutation of an entry that passes."""
    assert lint(palette(tmp_path)).returncode == 0


def test_prose_over_the_cap_is_refused(tmp_path: Path) -> None:
    """**The rule the tree had been exempt from by omission.**"""
    result = lint(palette(tmp_path, name="W" * 400))
    assert result.returncode == 1, result.stdout
    assert "over the 300" in result.stderr


def test_two_ports_sharing_a_name_are_refused(tmp_path: Path) -> None:
    """A port's name is how a connection reaches it."""
    ports = [dict(ENTRY["ports"][0]), {**ENTRY["ports"][1], "name": "inlet"}]
    result = lint(palette(tmp_path, ports=ports))
    assert result.returncode == 1, result.stdout
    assert "two ports share a name" in result.stderr


def test_a_dimension_the_vocabulary_lacks_is_refused(tmp_path: Path) -> None:
    """A field's dimension is a vocabulary id, and a name outside it resolves to nothing."""
    ports = [
        {
            **ENTRY["ports"][0],
            "fields": {"n": {"dimension": "enthalpy_of_mixing"}},
        }
    ]
    result = lint(palette(tmp_path, ports=ports))
    assert result.returncode == 1, result.stdout
    assert "vocabulary does not carry" in result.stderr


def test_an_unported_row_beside_a_model_is_refused(tmp_path: Path) -> None:
    """**The exclusivity that kills the twin.**

    `unit_ops.packed_column` has a model spec, so the porting fact belongs to that spec
    and nowhere else - which is how a claim corrected in one file left its twin lying in
    the other.
    """
    entry = {**ENTRY, "id": "unit_ops.packed_column"}
    directory = tmp_path / "column"
    directory.mkdir(parents=True, exist_ok=True)
    (directory / "packed_column.toml").write_text(
        dump({**entry, "unported": {"class": "PackingColumn"}}), encoding="utf-8"
    )
    result = lint(tmp_path)
    assert result.returncode == 1, result.stdout
    assert "the fact belongs there" in result.stderr


def test_an_unported_row_with_no_model_is_allowed(tmp_path: Path) -> None:
    """`unit_ops.simple_absorber` is the entry whose fact has nowhere else to live."""
    entry = {**ENTRY, "id": "unit_ops.simple_absorber"}
    directory = tmp_path / "column"
    directory.mkdir(parents=True, exist_ok=True)
    (directory / "simple_absorber.toml").write_text(
        dump({**entry, "unported": {"class": "AmineKentEisenberg"}}), encoding="utf-8"
    )
    assert lint(tmp_path).returncode == 0


def test_the_schema_and_the_struct_describe_the_same_document() -> None:
    """Every field `UnitOpSpec` declares is a property of the schema, and the reverse.

    `family` is the exception and is stated rather than skipped: it is `#[serde(skip)]`,
    filled from the path the loader was given, so no document carries one and the schema
    must not allow it.
    """
    import json

    properties = set(json.loads(SCHEMA.read_text(encoding="utf-8"))["properties"])
    source = UNIT_OP_RS.read_text(encoding="utf-8")
    body = source.split("pub struct UnitOpSpec {", 1)[1].split("\n}", 1)[0]
    fields = set(re.findall(r"^\s+pub ([a-z_]+):", body, re.M))

    assert fields - properties == {"family"}, (
        f"UnitOpSpec declares {sorted(fields - properties)} and the schema does not"
    )
    assert properties - fields == set(), (
        f"the schema declares {sorted(properties - fields)} and UnitOpSpec does not"
    )
    assert "family" not in properties, "a skipped field is not one a document can carry"
