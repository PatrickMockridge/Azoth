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
    conceptual = model.conceptualModel
    assert conceptual is not None
    # pyDEXPI's stub types the field as the abstract `ConceptualModel`; the loader returns
    # the concrete `ProcessModel`, which is the one that carries the steps.
    steps = getattr(conceptual, "processSteps", None)
    assert steps is not None
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


PLANT_FIXTURE = (
    ROOT / "crates" / "azoth-dexpi" / "tests" / "fixtures" / "plant_simple.toml"
).read_text()
PROTEUS_FIXTURE = (
    ROOT / "crates" / "azoth-dexpi" / "tests" / "fixtures" / "proteus_simple.xml"
).read_text()


def test_a_plant_export_loads_in_pydexpi() -> None:
    xml, report = dexpi.export_plant(
        PLANT_FIXTURE,
        PALETTE,
        allow_lossy_mappings=True,
        with_azoth_extensions=True,
    )
    assert report["items"] == 3
    assert report["segments"] == 5
    model = load_xml_string(xml, strict=True)
    plant = model.conceptualModel
    assert plant is not None
    items = getattr(plant, "taggedPlantItems", None)
    assert items is not None
    assert len(items) == 3
    functions = getattr(plant, "processInstrumentationFunctions", None)
    assert functions is not None
    assert len(functions) == 2
    switch = next(
        item for item in functions if item.processInstrumentationFunctionModifier == "SHH"
    )
    assert switch.processInstrumentationFunctionCategory == "P"
    assert switch.processInstrumentationFunctionNumber == "101"
    assert switch.safetyRelevanceClass == "SIL 2"
    assert switch.votingSystemRepresentation == "2oo3"
    assert len(switch.signalConveyingFunctions) == 1


def test_the_plant_round_trip_recovers_the_topology() -> None:
    xml, _ = dexpi.export_plant(
        PLANT_FIXTURE,
        PALETTE,
        allow_lossy_mappings=True,
        with_azoth_extensions=True,
    )
    _, report = dexpi.import_plant(xml, PALETTE)
    units = [instance["unit"] for instance in report["flowsheet"]["instances"]]
    assert units == ["unit_ops.separator", "unit_ops.compressor", "unit_ops.cooler"]
    assert len(report["flowsheet"]["connections"]) == 5


def test_a_proteus_plant_model_imports() -> None:
    _, report = dexpi.import_proteus(PROTEUS_FIXTURE, PALETTE)
    units = {instance["id"]: instance["unit"] for instance in report["flowsheet"]["instances"]}
    assert units["P-100"] == "unit_ops.pump"
    assert units["T-100"] == "unit_ops.tank"
    assert units["GlobeValve-1"] == "unit_ops.throttling_valve"
    assert report["flowsheet"]["products"] == ["PRODUCT"]


def test_the_schema_is_pinned() -> None:
    assert dexpi.schema_sha256() == (
        "0bd6568b7c52fd59ac4ef6d8547ab2b2aa890c009a3a1bfbc5d7cc18c53b36e6"
    )


def test_a_dexpi_2_drawing_loads_strictly() -> None:
    """A `Core/Diagram` document is a DEXPI 2.0 document to the independent parser too."""
    fixture = (
        ROOT / "crates" / "azoth-dexpi" / "tests" / "fixtures" / "native_projection_shape.dexpi.xml"
    )
    model = load_xml_string(fixture.read_text(), strict=True)
    assert model is not None


def test_a_pfd_export_carries_a_drawing_and_places_it() -> None:
    xml, report = dexpi.export_process(
        DEMO,
        PALETTE,
        allow_lossy_mappings=True,
        project_recycles=True,
        with_azoth_extensions=True,
        drawing=True,
    )
    assert report["diagram_groups"] > 0
    assert report["diagram_primitives"] > 0
    load_xml_string(xml, strict=True)
    toml, imported = dexpi.import_process(xml, PALETTE)
    assert imported["layout_positions"] > 0
    assert "[layout.instances]" in toml


def test_a_pid_export_carries_a_drawing_and_places_it() -> None:
    xml, report = dexpi.export_plant(
        PLANT_FIXTURE,
        PALETTE,
        allow_lossy_mappings=True,
        with_azoth_extensions=True,
        drawing=True,
    )
    assert report["diagram_groups"] > 0
    load_xml_string(xml, strict=True)
    toml, imported = dexpi.import_plant(xml, PALETTE)
    assert imported["layout_positions"] > 0
    assert "[layout.instances]" in toml


def test_a_drawing_view_reports_coverage_and_layout() -> None:
    report = dexpi.view(PROTEUS_FIXTURE, PALETTE)
    assert report["format"] == "proteus"
    assert report["counts"]["modelled"] > 0
    assert report["scaffold"]["layout_positions"] > 0
    assert report["drawing"]["groups"]


def test_a_valve_exports_inside_a_segment_and_imports_again() -> None:
    valve = (ROOT / "crates/azoth-dexpi/tests/fixtures/plant_valve.toml").read_text()
    xml, _ = dexpi.export_plant(valve, PALETTE, drawing=True)
    assert "Plant/Piping.GlobeValve" in xml
    _, imported = dexpi.import_plant(xml, PALETTE)
    units = {instance["id"]: instance["unit"] for instance in imported["flowsheet"]["instances"]}
    assert units["V-101"] == "unit_ops.throttling_valve"
