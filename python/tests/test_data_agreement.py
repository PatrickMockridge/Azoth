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
        "data/components/UNIFACcomp.csv",
        "data/components/UNIFACGroupParam.csv",
        "data/components/UNIFACInterParam.csv",
        "data/components/mbwr32.csv",
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


def test_the_nrtl_matrices_resolve_from_names() -> None:
    """The NRTL columns, resolved into the matrices `eos.nrtl_activity_coefficients` takes.

    `alpha` is symmetric with a zero diagonal; `dij` is directional (`g_ij != g_ji`),
    so reversing the name order swaps the off-diagonal energy. Row-major, so index 1 is
    `(0, 1)` and index 2 is `(1, 0)`.
    """
    params = components.nrtl_parameters(["methanol", "water"])
    assert params.alpha == (0.0, 0.303, 0.303, 0.0)
    assert params.dij == (0.0, -48.68, 610.6, 0.0)

    reversed_params = components.nrtl_parameters(["water", "methanol"])
    assert reversed_params.dij == (0.0, 610.6, -48.68, 0.0)


def test_the_unifac_parameters_resolve_from_names() -> None:
    """The UNIFAC group tables, resolved into the inputs `unifac_activity_coefficients` takes."""
    groups, group_r, group_q, aij = components.unifac_parameters(["methanol", "water"])
    assert groups == ((1.0, 0.0), (0.0, 1.0))
    assert group_r == (1.4311, 0.92)
    assert group_q == (1.432, 1.4)
    assert aij == ((0.0, -181.0), (289.6, 0.0))


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


# ---------------------------------------------------------------------------
# A keycard, which is the other half of the data
# ---------------------------------------------------------------------------
#
# The comparison above is over the *file*. A keycard is the other thing that decides
# what a calculation runs on, and until an overlay existed in Rust there was nothing to
# compare it against: the merge rule - a card naming only `omega` keeps the shipped `Tc`
# and `Pc` - had one implementation, in `azoth.eos.components`, and a rule with one
# implementation has nothing to be held to.
#
# The cases below are one per *rule* rather than per value, because the rules are what a
# keycard is. Nothing here re-implements the merge: both sides are asked the same
# question and their answers are compared, which is the arrangement the two kernels use.


def _card(**sections: Any) -> Any:
    """A keycard built from a mapping, through the loader's own validation."""
    from azoth import keycard

    return keycard.use({"schema_version": 2, **sections})


def _overlay(card: Any) -> Any:
    """The same card, as the extension reads it."""
    import azoth._rust_bridge as bridge

    return bridge.overlay_from(card)


def python_carded_rows(card: Any) -> list[dict[str, Any]]:
    """Every name the card states, as the Python reference resolves it, by name."""
    rows = []
    for name in sorted(card.components):
        record = components.entry(name, card=card)
        cp = record.cp
        rows.append(
            {
                "name": record.name.lower(),
                "tc_k": record.Tc.to("K").magnitude,
                "pc_pa": record.Pc.to("Pa").magnitude,
                "acentric_factor": record.omega,
                "cp_a": None if cp is None else cp[0],
                "cp_b": None if cp is None else cp[1],
                "cp_c": None if cp is None else cp[2],
                "cp_d": None if cp is None else cp[3],
                "cp_e": None if cp is None else cp[4],
            }
        )
    return rows


def rust_carded_rows(card: Any) -> list[dict[str, Any]]:
    """The same names, as the core resolves them.

    `COMPONENT_FIELDS` and not a set of its own: a card can change `Tc`, `Pc`, `omega`
    and nothing else, and the fields the row carries are the same fields either way. A
    separate tuple would be a second place the row's shape is written down.
    """
    rows = extension().overlay_component_rows(_overlay(card))
    return rust_rows(rows, COMPONENT_FIELDS)


def assert_carded_rows_agree(card: Any) -> None:
    """Both implementations, asked about one card, must answer the same thing."""
    rust = rust_carded_rows(card)
    python = python_carded_rows(card)
    assert rust == python, (
        f"the two merge rules disagree about {sorted(card.components)}\n"
        f"  rust:   {rust}\n  python: {python}"
    )


def test_an_empty_card_resolves_to_nothing() -> None:
    """The guard that makes every case below readable, and the reason it exists.

    A comparison that ignored its argument would pass this case and fail every other
    one. Without it, the others could be green because both sides were shown the same
    empty answer.
    """
    card = _card()
    assert python_carded_rows(card) == []
    assert rust_carded_rows(card) == []


def test_a_parameter_override_keeps_the_parameters_it_does_not_name() -> None:
    """The rule a keycard is built on, compared rather than asserted twice."""
    assert_carded_rows_agree(
        _card(components={"methane": {"omega": {"value": 0.5, "unit": "dimensionless"}}})
    )

    # And both sides kept the shipped `Tc`: a card that replaced the record whole would
    # make a user correcting one value restate the others.
    shipped = components.entry("methane").Tc.to("K").magnitude
    assert (
        python_carded_rows(
            _card(components={"methane": {"omega": {"value": 0.5, "unit": "dimensionless"}}})
        )[0]["tc_k"]
        == shipped
    )


def test_a_full_override_of_a_shipped_substance_agrees() -> None:
    assert_carded_rows_agree(
        _card(
            components={
                "methane": {
                    "Tc": {"value": 190.0, "unit": "K"},
                    "Pc": {"value": 4.599e6, "unit": "Pa"},
                    "omega": {"value": 0.0115, "unit": "dimensionless"},
                }
            }
        )
    )


def test_an_added_substance_has_no_polynomial_on_either_side() -> None:
    """The case that forced `Entry.cp` to be optional in Rust.

    A card supplies the parameters a *cubic* reads, and a heat-capacity polynomial is
    not one of them - so a substance a card adds has none, on both sides, and the
    comparison has to be able to say so.
    """
    card = _card(
        components={
            "unobtainium": {
                "Tc": {"value": 500.0, "unit": "K"},
                "Pc": {"value": 2.0e6, "unit": "Pa"},
                "omega": {"value": 0.3, "unit": "dimensionless"},
            }
        }
    )
    assert_carded_rows_agree(card)

    row = rust_carded_rows(card)[0]
    assert row["cp_a"] is None, "an added substance has no polynomial on the Rust side"
    assert python_carded_rows(card)[0]["cp_a"] is None


def test_a_name_in_neither_source_is_refused_by_the_same_class() -> None:
    """The class matters as much as the refusal.

    Python raises `PropertyUnavailableError` here and Rust has to raise the same thing,
    or one lookup reports itself two ways depending on which backend answered. The
    Python exception is a `LookupError` and the Rust variant is not an `InvalidInput`,
    which is the distinction: "you gave me nonsense" is a different answer from "I do
    not have that substance".
    """
    from azoth.core.errors import PropertyUnavailableError

    card = _card()
    with pytest.raises(PropertyUnavailableError):
        components.entry("unobtainium", card=card)
    with pytest.raises(PropertyUnavailableError):
        extension().overlay_entry_row("unobtainium", _overlay(card))


def test_a_partial_new_substance_is_refused_by_the_same_class() -> None:
    """Completing it from a similar substance would be inventing data, on both sides."""
    from azoth.core.errors import PropertyUnavailableError

    card = _card(components={"unobtainium": {"Tc": {"value": 500.0, "unit": "K"}}})
    with pytest.raises(PropertyUnavailableError):
        components.entry("unobtainium", card=card)
    with pytest.raises(PropertyUnavailableError):
        extension().overlay_entry_row("unobtainium", _overlay(card))


@pytest.mark.parametrize("value", [0.5, 0.0])
def test_a_kij_override_agrees(value: float) -> None:
    """Including a zero, which is the case that reads as "nothing stated".

    Overriding a fitted pair back to ideal mixing is a caller stating something, and
    a lookup that drops the zero silently undoes it. Both sides have to mean the same
    thing by it - Python omits the pair from `kij_for` and Rust resolves it to zero,
    and both mean ideal mixing.
    """
    card = _card(kij=[{"component_a": "methane", "component_b": "n-butane", "value": value}])

    resolved = extension().overlay_kij_rows(_overlay(card))
    assert resolved == [("methane", "n-butane", value)], (
        f"the core resolved the pair to {resolved}, not to the card's {value}"
    )

    # Python's effective value for the pair, with an absent pair meaning ideal mixing.
    python_side = components.kij_for(("methane", "n-butane"), card=card).get((0, 1), 0.0)
    assert python_side == value


def test_the_shipped_pair_is_still_there_without_a_card() -> None:
    """The other direction: a card changes the pair and nothing else does."""
    shipped = components.kij_for(("methane", "n-butane"))[(0, 1)]
    assert shipped != 0.0, "the fixture needs a fitted non-zero pair"
    card = _card(kij=[{"component_a": "methane", "component_b": "n-butane", "value": 0.0}])
    assert components.kij_for(("methane", "n-butane"), card=card).get((0, 1), 0.0) == 0.0
