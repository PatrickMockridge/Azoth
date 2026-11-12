"""The model's declared outputs are generated from the specs, and say what they read.

`gen_model_inputs.py` writes two tables: the inputs a form draws and the outputs a result sheet
explains. The output table is the one this file holds, because a description or an `optional`
dropped there is invisible until a reader meets a value with no sentence beside it.

The generated file itself is checked against the generator's own output by
`python tools/gen_model_inputs.py --check` in CI; what is asserted here is the reading - the
default kind, an enum's values, the dimension derived from the unit, the `optional` flag, and
the refusal of a unit the vocabulary does not carry.
"""

from __future__ import annotations

import importlib
import sys
from pathlib import Path
from typing import Any

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
TOOLS = REPO_ROOT / "tools"


@pytest.fixture(scope="module")
def generator():
    sys.path.insert(0, str(TOOLS))
    return importlib.import_module("gen_model_inputs")


def synthetic() -> tuple[list[dict[str, Any]], list[tuple[dict[str, Any], dict[str, Any]]]]:
    """One palette entry and one model, carrying one output of each kind that matters."""
    palette = [{"id": "unit_ops.demo", "_path": "specs/unit_ops/demo.toml"}]
    model = {
        "id": "process.demo",
        "_path": "specs/models/process/demo.toml",
        "outputs": {
            "outlet_p": {"unit": "Pa", "description": "the outlet pressure.\n"},
            "regime": {
                "type": "enum",
                "values": ["laminar", "turbulent"],
                "description": "the regime.",
            },
            "flag": {"type": "boolean", "description": "a flag."},
            "rows": {
                "type": "record_list",
                "optional": True,
                "description": "one row per fitting.",
            },
        },
    }
    return palette, [(palette[0], model)]


def test_the_table_keeps_every_declaration(generator) -> None:
    palette, rows = synthetic()
    text = generator.emit_outputs_rust(rows, palette)

    assert 'unit_op: "unit_ops.demo"' in text
    assert 'model: "process.demo"' in text
    # Spec order, which is the order a result sheet draws. Read the table only: the banner's
    # own prose uses the same words (`a flag, a name, a record list`).
    table = text.split("static ALL", 1)[1]
    assert (
        table.index("outlet_p") < table.index("regime") < table.index("flag") < table.index("rows")
    )
    # The default kind is a quantity, with its unit and the dimension the vocabulary gives it.
    assert 'kind: "quantity"' in text
    assert 'unit: Some("Pa")' in text
    assert 'dimension: Some("pressure")' in text
    # An enum carries its values; a flag and a record list carry no unit.
    assert 'kind: "enum"' in text and 'values: &["laminar", "turbulent"]' in text
    assert 'kind: "boolean"' in text and "unit: None" in text and "dimension: None" in text
    # The optional flag and the description, with the spec's newline escaped into the literal.
    assert 'kind: "record_list"' in text and "optional: true" in text
    assert "one row per fitting." in text
    assert "the outlet pressure.\\n" in text


def test_an_output_unit_outside_the_vocabulary_is_refused(generator) -> None:
    palette, rows = synthetic()
    rows[0][1]["outputs"]["outlet_p"]["unit"] = "furlong"
    with pytest.raises(SystemExit):
        generator.emit_outputs_rust(rows, palette)


def test_the_committed_table_carries_the_heater_it_was_generated_from(generator) -> None:
    text = (REPO_ROOT / "crates/azoth-process/src/model_outputs_gen.rs").read_text(encoding="utf-8")
    heater = text.split('unit_op: "unit_ops.heater"', 1)[1]
    assert 'model: "process.heater"' in heater
    assert 'name: "outlet_duty"' in heater
    assert 'unit: Some("W")' in heater
    assert 'name: "vfd_efficiency"' in heater
    assert "optional: true" in heater
