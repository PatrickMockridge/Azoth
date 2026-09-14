"""The two implementations must read the same data tables.

Four docstrings used to claim this test existed - ``fluids.rs``, ``fittings.rs``,
``reference/fittings.py`` and ``_data.py`` - and none of it was true. The test could
not be written either: the compiled extension exposed the calculations and not the
tables they were built from, so a Python test had no Rust-parsed value to compare
against. The claim was corrected in one commit and made true in another, and this is
the second.

# Why there are two kinds of assertion here

**Field by field** catches a parser divergence. Python's `csv` module and the `csv`
crate are different code and do not have to agree about quoting, byte-order marks,
CRLF or exponent notation.

**Byte for byte** catches something the field comparison cannot: two *different files*
that happen to parse to the same values. That is not hypothetical - it is exactly what
a stale copy bundled into a wheel looks like, and it is the one shape of failure that
would leave every numerical test in this repository green while the Rust core and the
Python reference were answering from different tables.

The byte comparison also catches an ordinary working mistake with a clear diagnosis:
editing a CSV without rebuilding the extension leaves the embedded text stale, which
is not a bug in either implementation but is worth failing over rather than leaving
two sides quietly disagreeing in a source tree.
"""

from __future__ import annotations

import hashlib
import importlib
from pathlib import Path
from types import ModuleType
from typing import Any

import pytest

from azoth._data import find
from azoth.eos import components
from azoth.hydraulics.reference.fittings import registry
from azoth.properties import available_fluids, provider_for

pytestmark = pytest.mark.requires_rust

#: Fields compared row by row, per dataset. Written out rather than derived from
#: `dataclasses.fields`, because a field silently dropped from one side is the failure
#: this is here to catch - deriving the list from one of the two sides would let it
#: shrink without complaint.
FITTING_FIELDS = (
    "id",
    "family",
    "name",
    "n_ld",
    "f_t_basis",
    "citation",
    "verify_status",
)

FLUID_FIELDS = (
    "fluid",
    "temperature_c",
    "density_kg_m3",
    "dynamic_viscosity_pa_s",
    "citation",
    "verify_status",
)

#: The columns a calculation reads out of the component databank. Not every column the
#: file carries: a comparison of a field neither side acts on proves nothing about the
#: numbers a flash uses.
COMPONENT_FIELDS = (
    "name",
    "tc_k",
    "pc_pa",
    "acentric_factor",
    "cp_a",
    "cp_b",
    "cp_c",
    "cp_d",
    "cp_e",
)

KIJ_FIELDS = ("component_a", "component_b", "kij_pr")


def extension() -> ModuleType:
    """The compiled extension, imported by name.

    The same accessor `test_cross_impl.py` and `test_registration_completeness.py`
    use, and for the same reason: the module does not exist until the bindings are
    built, and an attribute mypy cannot resolve is a worse trade than a lookup that
    fails clearly.
    """
    try:
        return importlib.import_module("azoth._core")
    except ImportError as exc:  # pragma: no cover - the marker skips these
        raise AssertionError("azoth._core is not built; run `maturin develop`") from exc


def python_fitting_rows() -> list[dict[str, Any]]:
    """The registry as the Python reference parses it."""
    return [
        {
            "id": row.id,
            "family": row.family,
            "name": row.name,
            "n_ld": row.n_ld,
            "f_t_basis": row.f_t_basis,
            "citation": row.citation,
            "verify_status": row.status.value,
        }
        for row in registry()
    ]


def python_fluid_rows(name: str) -> list[dict[str, Any]]:
    """A fluid's table as the Python side parses it."""
    provider = provider_for(name)
    return [
        {
            "fluid": provider.name,
            "temperature_c": point.temperature_c,
            "density_kg_m3": point.density_kg_m3,
            "dynamic_viscosity_pa_s": point.dynamic_viscosity_pa_s,
            "citation": point.citation,
            "verify_status": point.verify_status,
        }
        for point in provider.points
    ]


def python_component_rows() -> list[dict[str, Any]]:
    """The databank as the Python reference parses it, ordered by name.

    `available()` is the public surface and is sorted, which is the one ordering both
    sides can reproduce without agreeing on how the file happens to be laid out.
    """
    out = []
    for name in components.available():
        record = components.entry(name)
        if record.cp is None:  # pragma: no cover - only a keycard-added substance
            continue
        out.append(
            {
                "name": record.name.lower(),
                "tc_k": record.Tc.to("K").magnitude,
                "pc_pa": record.Pc.to("Pa").magnitude,
                "acentric_factor": record.omega,
                "cp_a": record.cp[0],
                "cp_b": record.cp[1],
                "cp_c": record.cp[2],
                "cp_d": record.cp[3],
                "cp_e": record.cp[4],
            }
        )
    return sorted(out, key=lambda row: row["name"])


def python_kij_rows() -> list[dict[str, Any]]:
    """The interaction table as the Python reference parses it, one row per pair.

    `_kij()` rather than a public function: the shipped map is symmetric and keyed both
    ways round, which is the storage a flash wants and not a table. Reaching in here
    keeps a public surface from existing only to be compared.
    """
    return sorted(
        (
            {"component_a": a, "component_b": b, "kij_pr": value}
            for (a, b), value in components._kij().items()
            if a < b
        ),
        key=lambda row: (row["component_a"], row["component_b"]),
    )


def rust_rows(rows: list[Any], fields: tuple[str, ...]) -> list[dict[str, Any]]:
    """Transported rows flattened to plain dicts, in a fixed field order."""
    return [{field: getattr(row, field) for field in fields} for row in rows]


def test_this_build_embeds_something() -> None:
    """Guard against every test below passing vacuously.

    A build that embedded no files, or a table that parsed to no rows, would let the
    comparisons succeed by having nothing to compare.
    """
    files = extension().data_files()
    assert files, "the extension embeds no data files"
    assert extension().fittings_rows(), "the extension parsed no fitting rows"
    assert extension().component_rows(), "the extension parsed no component rows"
    assert extension().kij_rows(), "the extension parsed no interaction rows"
    for name in available_fluids():
        assert extension().fluid_rows(name), f"the extension parsed no rows for {name}"


def test_the_two_sides_agree_on_which_files_exist() -> None:
    """The same set of files, addressed by the same paths.

    A file one side carries and the other does not is the failure that would follow
    from a packaging change made on one side only.
    """
    rust_paths = {data_file.path for data_file in extension().data_files()}
    expected = {
        "data/fittings/crane_k_factors.csv",
        "data/components/components.csv",
        "data/components/kij.csv",
    } | {f"data/fluids/{name}.csv" for name in available_fluids()}
    assert rust_paths == expected, (
        f"the extension embeds {sorted(rust_paths)} but this build of the Python side "
        f"expects {sorted(expected)}"
    )


@pytest.mark.parametrize(
    "path",
    [
        "data/fittings/crane_k_factors.csv",
        "data/fluids/water.csv",
        "data/components/components.csv",
        "data/components/kij.csv",
    ],
)
def test_a_named_file_is_byte_identical_on_both_sides(path: str) -> None:
    """A specific file, so a failure names which table is out of step.

    The parametrized list is short and explicit rather than derived from the
    extension, because this is the test a reader looks at to find out what a passing
    run actually proved.
    """
    by_path = {data_file.path: data_file for data_file in extension().data_files()}
    assert path in by_path, f"the extension does not embed {path}"
    assert_byte_identical(by_path[path])


def test_every_embedded_file_is_byte_identical_on_both_sides() -> None:
    """The same check over whatever the build happens to embed.

    Kept alongside the named cases so a file added to the extension is compared
    without anyone remembering to add it to a list.
    """
    for data_file in extension().data_files():
        assert_byte_identical(data_file)


def assert_byte_identical(data_file: Any) -> None:
    """The bytes this build embedded against the bytes on disk.

    Hashed rather than compared directly so a failure prints two digest prefixes
    rather than a diff of two CSVs, which is unreadable at this size.
    """
    embedded = hashlib.sha256(data_file.text.encode("utf-8")).hexdigest()
    on_disk = hashlib.sha256(find(data_file.path).read_bytes()).hexdigest()
    assert embedded == on_disk, (
        f"{data_file.path}: the Rust core embedded {embedded[:16]}... but this side "
        f"reads {on_disk[:16]}...\n"
        f"  If the file was edited, the extension is stale - rebuild with "
        f"`maturin develop`.\n"
        f"  If it was not, the two sides are reading different files, which is the "
        f"failure this test exists for: parsed values can agree across two different "
        f"files, and a stale copy in a wheel would look exactly like that."
    )


def test_every_fitting_row_agrees_field_by_field() -> None:
    """Every field of every row, compared in order.

    Order matters and is asserted: the registry is parsed in file order on both sides,
    and a divergence in ordering would mean one side is dropping and re-adding rows.
    """
    rust = rust_rows(extension().fittings_rows(), FITTING_FIELDS)
    python = python_fitting_rows()

    assert len(rust) == len(python), (
        f"the two sides parsed {len(rust)} and {len(python)} fitting rows"
    )
    for index, (left, right) in enumerate(zip(rust, python, strict=True)):
        assert left == right, (
            f"fittings row {index} ({left.get('id')}) differs\n  rust:   {left}\n  python: {right}"
        )


@pytest.mark.parametrize("name", sorted(available_fluids()))
def test_every_fluid_row_agrees_field_by_field(name: str) -> None:
    """Every field of every row of a fluid table."""
    rust = rust_rows(extension().fluid_rows(name), FLUID_FIELDS)
    python = python_fluid_rows(name)

    assert len(rust) == len(python), (
        f"{name}: the two sides parsed {len(rust)} and {len(python)} rows"
    )
    for index, (left, right) in enumerate(zip(rust, python, strict=True)):
        assert left == right, f"{name} row {index} differs\n  rust:   {left}\n  python: {right}"


def test_every_component_row_agrees_field_by_field() -> None:
    """Every column a calculation reads, for every substance.

    This is the comparison that makes "both sides read one databank" a checked claim
    rather than an architectural intention. Python's `csv` module and the `csv` crate
    are different code: they do not have to agree about quoting, a byte-order mark or
    exponent notation, and a divergence here is a component whose critical pressure is
    one value in Rust and another in Python.
    """
    rust = rust_rows(extension().component_rows(), COMPONENT_FIELDS)
    python = python_component_rows()

    assert len(rust) == len(python), (
        f"the two sides parsed {len(rust)} and {len(python)} component rows"
    )
    for index, (left, right) in enumerate(zip(rust, python, strict=True)):
        for field in COMPONENT_FIELDS:
            assert left[field] == pytest.approx(right[field], rel=0, abs=0), (
                f"component row {index} ({left['name']}) differs on `{field}`\n"
                f"  rust:   {left[field]!r}\n  python: {right[field]!r}"
            )


def test_every_interaction_row_agrees_field_by_field() -> None:
    """Every pair, in the same order on both sides."""
    rust = rust_rows(extension().kij_rows(), KIJ_FIELDS)
    python = python_kij_rows()

    assert len(rust) == len(python), (
        f"the two sides parsed {len(rust)} and {len(python)} interaction rows"
    )
    for index, (left, right) in enumerate(zip(rust, python, strict=True)):
        assert left["component_a"] == right["component_a"], f"interaction row {index}"
        assert left["component_b"] == right["component_b"], f"interaction row {index}"
        assert left["kij_pr"] == pytest.approx(right["kij_pr"], rel=0, abs=0), (
            f"interaction row {index} ({left['component_a']}/{left['component_b']}): "
            f"rust {left['kij_pr']!r}, python {right['kij_pr']!r}"
        )


def test_an_unknown_fluid_fails_rather_than_returning_nothing() -> None:
    """An empty list for an unknown name would make a typo look like an empty table.

    The Python side raises `PropertyUnavailableError` for the same reason, so the two
    agree about the shape of the mistake as well as about the data.
    """
    with pytest.raises(Exception, match="unknown_fluid"):
        extension().fluid_rows("unknown_fluid")


def test_the_data_files_are_where_the_report_says_they_are() -> None:
    """Each embedded file resolves on this side too.

    The resolution path is the part that broke in an installed wheel, so it is worth
    asserting directly rather than only through the calculations that use it. A
    missing file raises `DataFileNotFoundError`, whose message says the data was not
    packaged.
    """
    for data_file in extension().data_files():
        resolved = find(data_file.path)
        assert isinstance(resolved, Path)
        assert resolved.is_file(), f"{data_file.path} did not resolve to a file"
