"""DEXPI 2.0 exchange: a flowsheet as a DEXPI Process model, and back.

The work is the Rust `azoth-dexpi` crate's; this module is the notebook-shaped door onto it.
A flowsheet goes in as its TOML text, the exchange comes back as XML text beside a JSON report,
and an import comes back as TOML beside its report.

**What an import cannot give you.** DEXPI carries no composition, so an imported file is a
scaffold: the instances, connections and products are there, the feeds are undeclared and no
unit carries parameters. `azoth.process.validate` on the returned text reports exactly that,
and filling in the feeds is the caller's next step.
"""

from __future__ import annotations

import json
from typing import Any

from azoth import _core

__all__ = [
    "export_plant",
    "export_process",
    "import_plant",
    "import_process",
    "import_proteus",
    "schema_sha256",
    "validate",
]


def _report(text: str) -> dict[str, Any]:
    value: dict[str, Any] = json.loads(text)
    return value


def export_process(
    flowsheet: str,
    palette_dir: str = "specs/unit_ops",
    *,
    allow_lossy_mappings: bool = False,
    project_recycles: bool = False,
    with_azoth_extensions: bool = False,
    topology_only: bool = False,
    allow_unconverged: bool = False,
    drawing: bool = False,
    maps: dict[str, str] | None = None,
) -> tuple[str, dict[str, Any]]:
    """Write a flowsheet as a DEXPI 2.0 Process exchange.

    Returns the XML and its report. A refused exchange raises `ValueError` whose message is the
    report JSON - the refusals are the product, not an afterthought.
    """
    xml, report = _core.dexpi_export(
        flowsheet,
        palette_dir,
        allow_lossy_mappings,
        project_recycles,
        with_azoth_extensions,
        topology_only,
        allow_unconverged,
        drawing,
        maps,
    )
    return xml, _report(report)


def import_process(
    xml: str,
    palette_dir: str = "specs/unit_ops",
    *,
    strict: bool = False,
    layout: bool = True,
) -> tuple[str, dict[str, Any]]:
    """Read a DEXPI 2.0 Process exchange as a flowsheet scaffold.

    Returns the flowsheet TOML and its report. A refused import raises `ValueError` whose message
    is the report JSON.
    """
    toml, report = _core.dexpi_import(xml, palette_dir, strict, layout)
    return toml, _report(report)


def export_plant(
    flowsheet: str,
    palette_dir: str = "specs/unit_ops",
    *,
    allow_lossy_mappings: bool = False,
    with_azoth_extensions: bool = False,
    drawing: bool = False,
    maps: dict[str, str] | None = None,
) -> tuple[str, dict[str, Any]]:
    """Write a flowsheet as a DEXPI 2.0 Plant exchange - a P&ID.

    Returns the XML and its report. A refusal raises `ValueError` whose message is the report.
    """
    xml, report = _core.dexpi_export_plant(
        flowsheet, palette_dir, allow_lossy_mappings, with_azoth_extensions, drawing, maps
    )
    return xml, _report(report)


def import_plant(
    xml: str,
    palette_dir: str = "specs/unit_ops",
    *,
    strict: bool = False,
    layout: bool = True,
) -> tuple[str, dict[str, Any]]:
    """Read a DEXPI 2.0 Plant exchange as a flowsheet scaffold."""
    toml, report = _core.dexpi_import_plant(xml, palette_dir, strict, layout)
    return toml, _report(report)


def import_proteus(
    xml: str,
    palette_dir: str = "specs/unit_ops",
    *,
    strict: bool = False,
    layout: bool = True,
) -> tuple[str, dict[str, Any]]:
    """Read a Proteus `PlantModel` - DEXPI 1.x's XML - as a flowsheet scaffold."""
    toml, report = _core.dexpi_import_proteus(xml, palette_dir, strict, layout)
    return toml, _report(report)


def validate(xml: str) -> dict[str, Any]:
    """Check a document against the semantic rules; returns the report."""
    return _report(_core.dexpi_validate(xml))


def schema_sha256() -> str:
    """The digest of the vendored official DEXPI XML schema."""
    return _core.dexpi_schema_sha256()
