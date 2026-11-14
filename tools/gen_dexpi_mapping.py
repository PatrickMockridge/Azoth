#!/usr/bin/env python3
"""Compile the DEXPI Process mapping from its one hand-written table.

`specs/dexpi/process_units.toml` is the only place the azoth palette -> DEXPI 2.0 Process
mapping is stated. This reads it, holds it to the palette and to the vendored schema, and
emits the table the crate compiles: the forward map, the reverse candidate index an importer
needs, and the three counts the documentation claims.

**The completeness rule is the point.** A palette entry with no row, a row naming no palette
entry, a type outside the profile, a refused row that still names a type, or a lossy or refused
row with no reason is a refusal here rather than a gap found at run time. That is the same rule
`gen_palette.py` applies to the palette itself.

    python tools/gen_dexpi_mapping.py            # write
    python tools/gen_dexpi_mapping.py --check    # fail if the tree is not what this emits
"""

from __future__ import annotations

import hashlib
import subprocess
import sys
import tomllib
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SPEC = ROOT / "specs" / "dexpi" / "process_units.toml"
PALETTE_DIR = ROOT / "specs" / "unit_ops"
XSD = ROOT / "crates" / "azoth-dexpi" / "schema" / "DEXPI_XML_Schema.xsd"
OUT = ROOT / "crates" / "azoth-dexpi" / "src" / "process_mapping_gen.rs"

# The official DEXPI XML V2.0.0 schema, unmodified. 2.0.1 publishes a byte-identical file,
# so one digest pins both releases.
XSD_SHA256 = "0bd6568b7c52fd59ac4ef6d8547ab2b2aa890c009a3a1bfbc5d7cc18c53b36e6"

# The reviewed Process types. They come from NeqSim's `Dexpi20ProcessModelWriter.type`
# (the switch at lines 401-449), which is the only reviewed source for this list; a type
# outside it has not been checked against the DEXPI Process model.
SUPPORTED_TYPES = [
    "Process/Process.Compressing",
    "Process/Process.Cooling",
    "Process/Process.Distilling",
    "Process/Process.ExchangingThermalEnergy",
    "Process/Process.HeatingInFurnace",
    "Process/Process.MixingSimple",
    "Process/Process.Pumping",
    "Process/Process.RegulatingFlow",
    "Process/Process.SeparatingByGravity",
    "Process/Process.Sink",
    "Process/Process.Source",
    "Process/Process.SplittingMaterial",
    "Process/Process.StoringFluids",
    "Process/Process.TransportingFluids",
]

DISPOSITIONS = ("faithful", "lossy", "refused")


def format_rust(text: str, generator: str) -> str:
    """Run the emitted Rust through rustfmt, which every generated file is compared against.

    A missing rustfmt is a refusal rather than a fallback: the check compares this output to
    the committed file, so unformatted text would report a drift that is not there.
    """
    try:
        proc = subprocess.run(
            ["rustfmt", "--edition", "2024"],
            input=text,
            capture_output=True,
            text=True,
            check=False,
        )
    except FileNotFoundError:
        sys.exit(
            f"{generator}: rustfmt is not on PATH, and the generated file is compared against "
            f"rustfmt's own output - install the Rust toolchain"
        )
    if proc.returncode != 0:
        sys.exit(f"{generator}: rustfmt refused the output:\n{proc.stderr}")
    return proc.stdout


def palette_ids() -> set[str]:
    ids: set[str] = set()
    for path in sorted(PALETTE_DIR.rglob("*.toml")):
        with path.open("rb") as handle:
            ids.add(tomllib.load(handle)["id"])
    if not ids:
        sys.exit(f"gen_dexpi_mapping: no palette entries under {PALETTE_DIR}")
    return ids


def read_spec() -> list[dict]:
    with SPEC.open("rb") as handle:
        document = tomllib.load(handle)
    if document.get("schema_version") != 1:
        sys.exit("gen_dexpi_mapping: schema_version must be 1")
    rows = document.get("unit", [])
    if not rows:
        sys.exit(f"gen_dexpi_mapping: {SPEC.relative_to(ROOT)} declares no [[unit]] rows")
    return rows


def check_rows(rows: list[dict]) -> None:
    palette = palette_ids()
    seen: dict[str, str] = {}
    problems: list[str] = []
    for row in rows:
        ident = row.get("azoth_id")
        if not ident:
            problems.append("a row has no azoth_id")
            continue
        if ident in seen:
            problems.append(f"{ident}: declared twice")
        seen[ident] = ident
        if ident not in palette:
            problems.append(f"{ident}: no such palette entry")
        disposition = row.get("disposition")
        if disposition not in DISPOSITIONS:
            problems.append(f"{ident}: disposition {disposition!r} is not one of {DISPOSITIONS}")
            continue
        types = row.get("dexpi_types", [])
        note = row.get("note", "").strip()
        if disposition in ("faithful", "lossy") and not types:
            problems.append(f"{ident}: {disposition} row names no DEXPI type")
        if disposition == "refused" and types:
            problems.append(f"{ident}: refused row names a type")
        if disposition in ("lossy", "refused") and not note:
            problems.append(f"{ident}: {disposition} row carries no reason")
        for dexpi_type in types:
            if dexpi_type not in SUPPORTED_TYPES:
                problems.append(f"{ident}: {dexpi_type} is outside the reviewed Process profile")
    for ident in sorted(palette - seen.keys()):
        problems.append(f"{ident}: palette entry has no mapping row")
    if problems:
        for problem in problems:
            print(f"gen_dexpi_mapping: {problem}", file=sys.stderr)
        sys.exit("gen_dexpi_mapping: the mapping table does not cover the palette")


def check_schema() -> None:
    digest = hashlib.sha256(XSD.read_bytes()).hexdigest()
    if digest != XSD_SHA256:
        sys.exit(
            f"gen_dexpi_mapping: {XSD.relative_to(ROOT)} has digest {digest}, not the pinned "
            f"{XSD_SHA256}; the vendored schema is not the DEXPI 2.0 release"
        )


def rust_string(text: str) -> str:
    return '"' + text.replace("\\", "\\\\").replace('"', '\\"') + '"'


def render(rows: list[dict]) -> str:
    rows = sorted(rows, key=lambda row: row["azoth_id"])
    faithful = sum(1 for row in rows if row["disposition"] == "faithful")
    lossy = sum(1 for row in rows if row["disposition"] == "lossy")
    refused = sum(1 for row in rows if row["disposition"] == "refused")

    entries = []
    for row in rows:
        types = ", ".join(rust_string(t) for t in row["dexpi_types"])
        entries.append(
            "    UnitMapping {\n"
            f"        azoth_id: {rust_string(row['azoth_id'])},\n"
            f"        disposition: Disposition::{row['disposition'].capitalize()},\n"
            f"        dexpi_types: &[{types}],\n"
            f"        note: {rust_string(row.get('note', '').strip())},\n"
            "    },\n"
        )

    reverse: dict[str, list[str]] = {}
    for row in rows:
        if row["disposition"] == "refused":
            continue
        for dexpi_type in row["dexpi_types"]:
            reverse.setdefault(dexpi_type, []).append(row["azoth_id"])
    reverse_rows = []
    for dexpi_type in sorted(reverse):
        ids = ", ".join(rust_string(i) for i in sorted(reverse[dexpi_type]))
        reverse_rows.append(f"    ({rust_string(dexpi_type)}, &[{ids}]),\n")

    text = (
        "//! GENERATED FILE - DO NOT EDIT BY HAND.\n"
        "//!\n"
        "//! Generated by `tools/gen_dexpi_mapping.py` from `specs/dexpi/process_units.toml`,\n"
        "//! which is the one hand-written statement of the azoth palette -> DEXPI 2.0 Process\n"
        "//! mapping. Its completeness is checked against `specs/unit_ops/`, and the vendored\n"
        "//! schema's digest is checked against the DEXPI 2.0 release.\n"
        "\n"
        "use crate::mapping::{Disposition, UnitMapping};\n"
        "\n"
        "/// Every palette entry, with what exchanging it as a DEXPI Process step costs.\n"
        "pub static PROCESS_UNITS: &[UnitMapping] = &[\n"
        f"{''.join(entries)}];\n"
        "\n"
        "/// The reverse index an importer reads: a DEXPI Process type and the azoth units it\n"
        "/// may stand for, in palette-id order. More than one is an ambiguity the importer\n"
        "/// reports rather than resolves silently.\n"
        "pub static PROCESS_TYPE_CANDIDATES: &[(&str, &[&str])] = &[\n"
        f"{''.join(reverse_rows)}];\n"
        "\n"
        f"/// Faithful mappings: {faithful}.\n"
        f"pub const FAITHFUL_COUNT: usize = {faithful};\n"
        "\n"
        f"/// Lossy mappings: {lossy}.\n"
        f"pub const LOSSY_COUNT: usize = {lossy};\n"
        "\n"
        f"/// Refused mappings: {refused}.\n"
        f"pub const REFUSED_COUNT: usize = {refused};\n"
    )
    return format_rust(text, "gen_dexpi_mapping")


PLANT_SPEC = ROOT / "specs" / "dexpi" / "plant_units.toml"
PLANT_OUT = ROOT / "crates" / "azoth-dexpi" / "src" / "plant_mapping_gen.rs"

#: The native Plant classes this profile reads, mirroring the Plant half of NeqSim's
#: `Dexpi20SemanticValidator.SUPPORTED_TYPES`. A type outside it is not reviewed here.
PLANT_SUPPORTED_TYPES = [
    "Plant/PlantModel",
    "Plant/ProcessEquipment.ProcessEquipment",
    "Plant/ProcessEquipment.CentrifugalCompressor",
    "Plant/ProcessEquipment.CentrifugalPump",
    "Plant/ProcessEquipment.Separator",
    "Plant/ProcessEquipment.AirCoolingSystem",
    "Plant/ProcessEquipment.TubularHeatExchanger",
    "Plant/ProcessEquipment.FiredHeater",
    "Plant/ProcessEquipment.Tank",
    "Plant/ProcessEquipment.Nozzle",
    "Plant/Piping.PipingNetworkSystem",
    "Plant/Piping.PipingNetworkSegment",
    "Plant/Piping.PipingNode",
    "Plant/Piping.Pipe",
    "Plant/Piping.GlobeValve",
    "Plant/Piping.SpringLoadedGlobeSafetyValve",
    "Plant/Piping.SwingCheckValve",
    "Plant/Piping.FlowInPipeOffPageConnector",
    "Plant/Piping.FlowOutPipeOffPageConnector",
    "Plant/Instrumentation.ProcessInstrumentationFunction",
    "Plant/Instrumentation.ProcessSignalGeneratingFunction",
    "Plant/Instrumentation.MeasuringLineFunction",
    "Plant/Diagram.PlantMetaData",
]


def read_plant_spec() -> list[dict]:
    with PLANT_SPEC.open("rb") as handle:
        document = tomllib.load(handle)
    if document.get("schema_version") != 1:
        sys.exit("gen_dexpi_mapping: plant schema_version must be 1")
    rows = document.get("unit", [])
    if not rows:
        sys.exit(f"gen_dexpi_mapping: {PLANT_SPEC.relative_to(ROOT)} declares no [[unit]] rows")
    return rows


def check_plant_rows(rows: list[dict]) -> None:
    palette = palette_ids()
    seen: dict[str, str] = {}
    reader_owner: dict[str, str] = {}
    problems: list[str] = []
    for row in rows:
        ident = row.get("azoth_id")
        if not ident:
            problems.append("a plant row has no azoth_id")
            continue
        if ident in seen:
            problems.append(f"{ident}: declared twice")
        seen[ident] = ident
        if ident not in palette:
            problems.append(f"{ident}: no such palette entry")
        write = row.get("write")
        if write not in DISPOSITIONS:
            problems.append(f"{ident}: write {write!r} is not one of {DISPOSITIONS}")
            continue
        types = row.get("plant_types", [])
        note = row.get("note", "").strip()
        if write in ("faithful", "lossy") and len(types) != 1:
            problems.append(f"{ident}: a writable row names exactly one Plant type")
        if write == "refused" and types:
            problems.append(f"{ident}: a refused row names a Plant type")
        if write in ("lossy", "refused") and not note:
            problems.append(f"{ident}: {write} row carries no reason")
        placement = row.get("placement", "item")
        if placement not in ("item", "segment"):
            problems.append(f"{ident}: placement {placement!r} is not `item` or `segment`")
        if placement == "segment":
            if write == "refused":
                problems.append(f"{ident}: a segment placement needs a writable row")
            if any(not plant_type.startswith("Plant/Piping.") for plant_type in types):
                problems.append(f"{ident}: a segment placement names a piping component")
        for plant_type in types:
            if plant_type not in PLANT_SUPPORTED_TYPES:
                problems.append(f"{ident}: {plant_type} is outside the reviewed Plant profile")
        for reader_type in row.get("reader_types", []):
            if reader_type in reader_owner:
                problems.append(
                    f"{ident}: reads {reader_type}, already read as {reader_owner[reader_type]}"
                )
            reader_owner[reader_type] = ident
    for ident in sorted(palette - seen.keys()):
        problems.append(f"{ident}: palette entry has no plant mapping row")
    if problems:
        for problem in problems:
            print(f"gen_dexpi_mapping: {problem}", file=sys.stderr)
        sys.exit("gen_dexpi_mapping: the plant table does not cover the palette")


def render_plant(rows: list[dict]) -> str:
    rows = sorted(rows, key=lambda row: row["azoth_id"])
    faithful = sum(1 for row in rows if row["write"] == "faithful")
    lossy = sum(1 for row in rows if row["write"] == "lossy")
    refused = sum(1 for row in rows if row["write"] == "refused")
    readable = sum(1 for row in rows if row.get("reader_types"))

    entries = []
    for row in rows:
        types = ", ".join(rust_string(t) for t in row["plant_types"])
        readers = ", ".join(rust_string(t) for t in row["reader_types"])
        entries.append(
            "    PlantMapping {\n"
            f"        azoth_id: {rust_string(row['azoth_id'])},\n"
            f"        write: Disposition::{row['write'].capitalize()},\n"
            f"        plant_types: &[{types}],\n"
            f"        reader_types: &[{readers}],\n"
            f"        note: {rust_string(row.get('note', '').strip())},\n"
            "    },\n"
        )

    classes = []
    for row in rows:
        for reader_type in row["reader_types"]:
            classes.append((reader_type, row["azoth_id"]))
    class_rows = [
        f"    ({rust_string(name)}, {rust_string(ident)}),\n" for name, ident in sorted(classes)
    ]

    # The native type a projected reader class is spelled as. A valve is a piping component;
    # everything else is equipment. This is what lets the semantic validator accept a
    # projection of a class the writer never emits.
    valve_types = {
        "GlobeValve",
        "ButterflyValve",
        "CheckValve",
        "ControlValve",
        "PressureSafetyValve",
        "PressureReliefValve",
        "PressureReducingValve",
        "BallValve",
        "GateValve",
        "PlugValve",
        "DiaphragmValve",
        "NeedleValve",
        "SwingCheckValve",
        "SpringLoadedGlobeSafetyValve",
        "PilotOperatedSafetyValve",
        "SolenoidValve",
        "PinchValve",
        "RuptureDisk",
        "Valve",
        "Pipe",
    }
    projected = sorted(
        {
            f"Plant/Piping.{name}" if name in valve_types else f"Plant/ProcessEquipment.{name}"
            for name, _ in classes
        }
    )
    projected_block = "".join(f"    {rust_string(name)},\n" for name in projected)
    segment_units = sorted(
        row["azoth_id"] for row in rows if row.get("placement", "item") == "segment"
    )
    segment_block = "".join(f"    {rust_string(unit)},\n" for unit in segment_units)

    text = (
        "//! GENERATED FILE - DO NOT EDIT BY HAND.\n"
        "//!\n"
        "//! Generated by `tools/gen_dexpi_mapping.py` from `specs/dexpi/plant_units.toml`,\n"
        "//! which is the one hand-written statement of the azoth palette -> DEXPI 2.0 Plant\n"
        "//! mapping. Its completeness is checked against `specs/unit_ops/`.\n"
        "\n"
        "use crate::mapping::Disposition;\n"
        "use crate::plant::PlantMapping;\n"
        "\n"
        "/// Every palette entry, with what exchanging it as a tagged plant item costs.\n"
        "pub static PLANT_UNITS: &[PlantMapping] = &[\n"
        f"{''.join(entries)}];\n"
        "\n"
        "/// The units placed inside a segment's `Items`, not among the tagged items.\n"
        "pub static PLANT_SEGMENT_UNITS: &[&str] = &[\n"
        f"{segment_block}];\n"
        "\n"
        "/// The reverse index a reader uses: a class name and the azoth unit it stands for.\n"
        "/// One name, one unit - a class that could mean two is mapped to the general one.\n"
        "pub static PLANT_CLASS_CANDIDATES: &[(&str, &str)] = &[\n"
        f"{''.join(class_rows)}];\n"
        "\n"
        "/// The native types a reader class is spelled as in a projection. A valve is a piping\n"
        "/// component; everything else is equipment. The semantic validator accepts these beside\n"
        "/// the writer's own narrow set, because these are what this reader has reviewed.\n"
        "pub static PLANT_CLASS_TYPES: &[&str] = &[\n"
        f"{projected_block}];\n"
        "\n"
        f"/// Faithful writes: {faithful}.\n"
        f"pub const PLANT_FAITHFUL_COUNT: usize = {faithful};\n"
        "\n"
        f"/// Lossy writes: {lossy}.\n"
        f"pub const PLANT_LOSSY_COUNT: usize = {lossy};\n"
        "\n"
        f"/// Refused writes: {refused}.\n"
        f"pub const PLANT_REFUSED_COUNT: usize = {refused};\n"
        "\n"
        f"/// Rows with at least one readable class name: {readable}.\n"
        f"pub const PLANT_READABLE_COUNT: usize = {readable};\n"
    )
    return format_rust(text, "gen_dexpi_mapping")


def main() -> None:
    check = "--check" in sys.argv
    check_schema()
    rows = read_spec()
    check_rows(rows)
    rendered = render(rows)
    plant_rows = read_plant_spec()
    check_plant_rows(plant_rows)
    plant_rendered = render_plant(plant_rows)

    for out, content, label in ((OUT, rendered, "Process"), (PLANT_OUT, plant_rendered, "Plant")):
        existing = out.read_text(encoding="utf-8") if out.exists() else None
        if existing != content:
            if check:
                print(f"gen_dexpi_mapping: {out.relative_to(ROOT)} is out of date")
                sys.exit(
                    f"gen_dexpi_mapping: generated {label} file is stale; run "
                    "tools/gen_dexpi_mapping.py"
                )
            out.write_text(content, encoding="utf-8")
            print(f"gen_dexpi_mapping: wrote {out.relative_to(ROOT)}")

    faithful = sum(1 for row in rows if row["disposition"] == "faithful")
    lossy = sum(1 for row in rows if row["disposition"] == "lossy")
    refused = sum(1 for row in rows if row["disposition"] == "refused")
    writable = sum(1 for row in plant_rows if row["write"] != "refused")
    print(
        f"gen_dexpi_mapping: {len(rows)} palette entries, {faithful} faithful, "
        f"{lossy} lossy, {refused} refused; Plant {writable} writable, "
        f"{len(plant_rows) - writable} refused" + (", generated files up to date" if check else "")
    )


if __name__ == "__main__":
    main()
