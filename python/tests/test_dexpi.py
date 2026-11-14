"""The DEXPI exchange, through the Python binding and an independent parser.

pyDEXPI is a second implementation of the same specification, and it is a dev dependency
rather than an optional extra: a schema that only this crate's own validator understands is a
claim about this crate, while a file a second implementation loads strictly is evidence about
the file. The repo's rule is that a missing check fails rather than skips, so an environment
without pyDEXPI fails this module.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest
from pydexpi.loaders import load_xml_string

from azoth import dexpi

ROOT = Path(__file__).resolve().parents[2]
PALETTE = str(ROOT / "specs" / "unit_ops")
DEMO = (ROOT / "specs" / "flowsheets" / "demo.toml").read_text()

FIXTURE = """
id = "flowsheets.flare"
name = "One flare"

products = ["vent"]

[[inputs]]
name = "feed_1"
components = ["methane"]
n = 1.0
z = [1.0]
P = 5.0e5
T = 300.0

[[instances]]
id = "f1"
unit = "unit_ops.flare"

[[connections]]
from = "feed_1"
to = "f1.inlet"

[[connections]]
from = "f1.product"
to = "vent"
"""


def export_demo() -> str:
    xml, report = dexpi.export_process(
        DEMO,
        PALETTE,
        allow_lossy_mappings=True,
        project_recycles=True,
        with_azoth_extensions=True,
    )
    assert report["steps"] == 7
    assert report["streams"] == 7
    assert report["approval_status"] == "NOT_A_DEXPI_EV_CERTIFICATE"
    return xml


def test_pydexpi_loads_azoths_export_strictly() -> None:
    model = load_xml_string(export_demo(), strict=True)
    steps = model.conceptualModel.processSteps
    assert len(steps) == 7
    kinds = [type(step).__name__ for step in steps]
    assert "MixingSimple" in kinds
    assert "Pumping" in kinds
    assert "SeparatingByGravity" in kinds


def test_the_exchange_round_trips_its_topology() -> None:
    xml = export_demo()
    toml, report = dexpi.import_process(xml, PALETTE)
    assert report["flowsheet"]["id"] == "flowsheets.demo"
    assert [instance["unit"] for instance in report["flowsheet"]["instances"]] == [
        "unit_ops.mixer",
        "unit_ops.pump",
        "unit_ops.heater",
        "unit_ops.separator",
    ]
    assert len(report["flowsheet"]["recycles"]) == 1
    # The scaffold is not runnable, and the report says why: DEXPI carries no composition.
    codes = {diagnostic["code"] for diagnostic in report["diagnostics"]}
    assert "input_composition_missing" in codes
    assert "recycles" in toml or "[[recycles]]" in toml


def test_a_refused_mapping_raises_with_its_report() -> None:
    with pytest.raises(ValueError) as raised:
        dexpi.export_process(FIXTURE, PALETTE)
    report = json.loads(str(raised.value))
    codes = {diagnostic["code"] for diagnostic in report["diagnostics"]}
    assert "refused_mapping" in codes


def test_validate_reports_rather_than_raising() -> None:
    report = dexpi.validate(export_demo())
    assert report["valid"] is True
    broken = dexpi.validate("<Model name='m' uri='u'></Model>")
    assert broken["valid"] is False


def test_the_schema_is_pinned() -> None:
    assert dexpi.schema_sha256() == (
        "0bd6568b7c52fd59ac4ef6d8547ab2b2aa890c009a3a1bfbc5d7cc18c53b36e6"
    )
