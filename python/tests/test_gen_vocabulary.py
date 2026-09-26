"""The vocabulary generator, on tables the tree does not have.

`gen_vocabulary` compiles `specs/vocabulary/vocabulary.toml` into four artefacts, and
the drift job fails if any of them is stale. What nothing exercised is the generator's
*refusals* - and a table is exactly the kind of file a person edits by hand, so a
generator whose disagreement between two columns is silently resolved is one that
produces four consistent artefacts describing a unit nobody meant.

Every test here mutates a real copy of the table and asserts the generator refuses it
with a message naming the problem. The `--vocabulary` flag exists for this: it points
the generator at a table that is not the repository's, so a test can compile a wrong
one without the tree ever holding it.

The positive direction is the other way round and lives in `test_units_contract.py`,
which asserts the four generated artefacts are the table that is actually in the tree.
"""

from __future__ import annotations

import copy
import importlib
import sys
import tomllib
from pathlib import Path
from types import ModuleType
from typing import Any

import pytest

from _toml_fixture import dump

REPO_ROOT = Path(__file__).resolve().parents[2]
VOCAB_PATH = REPO_ROOT / "specs" / "vocabulary" / "vocabulary.toml"


def gen_vocabulary() -> ModuleType:
    sys.path.insert(0, str(REPO_ROOT / "tools"))
    try:
        return importlib.import_module("gen_vocabulary")
    finally:
        sys.path.pop(0)


def a_table() -> dict[str, Any]:
    """The repository's table, as a mutable copy."""
    return copy.deepcopy(tomllib.loads(VOCAB_PATH.read_text()))


def compile_table(tmp_path: Path, table: dict[str, Any]) -> None:
    """Write a table and make the generator read it, whatever the outcome.

    The generator is driven through `main` with `--vocabulary`, so the test exercises
    the same path CI does rather than a helper beside it. `sys.argv` is set directly
    because argparse reads it at call time and restoring it is a `finally`; `monkeypatch`
    would work too, but the argument has to be built rather than substituted.

    `--check` and not a bare run: `load_table` runs before either branch, so every
    refusal below is reached identically - and a table that turned out to be *valid*
    would write the generated files over the repository's if this were a bare run.
    """
    path = tmp_path / "vocabulary.toml"
    path.write_text(dump(table))

    module = gen_vocabulary()
    previous = sys.argv
    sys.argv = ["gen_vocabulary.py", "--check", "--vocabulary", str(path)]
    try:
        module.main()
    finally:
        sys.argv = previous


def unit(table: dict[str, Any], name: str) -> dict[str, Any]:
    return next(u for u in table["units"] if u["id"] == name)


def dimension(table: dict[str, Any], name: str) -> dict[str, Any]:
    return next(d for d in table["dimensions"] if d["id"] == name)


# ---------------------------------------------------------------------------
# The refusals
# ---------------------------------------------------------------------------


def test_a_dimension_the_generator_does_not_know_is_refused(tmp_path: Path) -> None:
    """A dimension absent from `UOM_TYPES` stops generation.

    The map is total over the table by construction - the generator exits if it is
    not - which is what makes `None` in it an answer about `uom` rather than an
    oversight. A new dimension therefore forces a statement about whether `uom`
    carries it, instead of defaulting to "carries it" or "does not".
    """
    table = a_table()
    table["dimensions"].append({"id": "made_up", "exponents": [0, 0, 0, 0, 0, 0, 3]})
    with pytest.raises(SystemExit, match="UOM_TYPES"):
        compile_table(tmp_path, table)


def test_a_unit_whose_uom_path_disagrees_with_its_dimension_is_refused(
    tmp_path: Path,
) -> None:
    """The module in `uom:` has to be the dimension's own quantity.

    `uom` derives a module name from a quantity name, so a path whose module is not
    the dimension's is a path into a different quantity - and the generated
    conversion would call a constructor for a dimension the row does not name.
    """
    table = a_table()
    dimension(table, "length")["exponents"] = [2, 0, 0, 0, 0, 0, 0]
    unit(table, "mm")["dimension"] = "length"
    with pytest.raises(SystemExit, match="whose module is"):
        compile_table(tmp_path, table)


def test_a_unit_that_claims_a_uom_path_for_a_dimension_uom_lacks_is_refused(
    tmp_path: Path,
) -> None:
    """`mol/s` may not be given a `uom:` path.

    `uom::si` assembles its `Units` trait from a closed quantity list, so a molar
    flow is not among them. A table row claiming otherwise would generate Rust that
    cannot compile - and the failure would be an obscure type error rather than a
    statement about the vocabulary.
    """
    table = a_table()
    unit(table, "mol/s")["uom"] = "molar_flow::mole_per_second"
    unit(table, "mol/s")["rust_ctor"] = "moles_per_second"
    with pytest.raises(SystemExit, match="maps to None"):
        compile_table(tmp_path, table)


def test_a_unit_with_a_dimension_the_table_does_not_declare_is_refused(
    tmp_path: Path,
) -> None:
    """A dangling `dimension:` reference stops generation."""
    table = a_table()
    unit(table, "Pa")["dimension"] = "pressure_but_misspelled"
    with pytest.raises(SystemExit, match="not declared"):
        compile_table(tmp_path, table)


def test_a_unit_with_no_lean_dimension_is_refused(tmp_path: Path) -> None:
    """A new unit forces the question of what `lean-units` calls it.

    `LEAN_DIMENSIONS` is the hand-written claim about what each unit *is*, and the
    generated theorem is what makes the table agree with it. A unit with no entry
    would get a theorem whose right-hand side came from the table, which proves
    nothing - so the omission is refused rather than defaulted.
    """
    table = a_table()
    table["units"].append(
        {
            "id": "invented",
            "dimension": "length",
            "pint": "furlong",
        }
    )
    with pytest.raises(SystemExit, match="LEAN_DIMENSIONS"):
        compile_table(tmp_path, table)


def test_a_duplicate_unit_id_is_refused(tmp_path: Path) -> None:
    """Two rows for one name is a silent drop, not a merge.

    Every generated map keys on the name, so the second row would win and the first
    would disappear - and the row that disappeared would be the one a reviewer had
    just read.
    """
    table = a_table()
    table["units"].append(dict(unit(table, "mm")))
    with pytest.raises(SystemExit, match="duplicate unit id"):
        compile_table(tmp_path, table)


def test_a_duplicate_dimension_id_is_refused(tmp_path: Path) -> None:
    """The same, for dimensions - where two rows with one name could disagree."""
    table = a_table()
    table["dimensions"].append({"id": "length", "exponents": [9, 0, 0, 0, 0, 0, 0]})
    with pytest.raises(SystemExit, match="duplicate dimension id"):
        compile_table(tmp_path, table)


def test_an_unknown_schema_version_is_refused(tmp_path: Path) -> None:
    """A newer format is refused rather than read as though it were this one.

    The schema owns this rule and the generator does not repeat it: `schema_version`
    is a `const`, so a table declaring anything else fails validation before any
    generator code sees the field. One rule in one place, and the message names the
    field - which is what a duplicated check would have cost, since a second version
    check in the generator could never be reached to say anything.
    """
    table = a_table()
    table["schema_version"] = 99
    with pytest.raises(SystemExit, match="does not validate"):
        compile_table(tmp_path, table)


def test_an_exponent_vector_of_the_wrong_width_is_refused(tmp_path: Path) -> None:
    """One exponent per slot, or the tuple is read against the wrong slots.

    A short vector is worse than a rejection: `[1]` read against seven slots is a
    length as far as anything downstream can tell.
    """
    table = a_table()
    dimension(table, "length")["exponents"] = [1]
    with pytest.raises(SystemExit, match="exponents for"):
        compile_table(tmp_path, table)


# ---------------------------------------------------------------------------
# The shape of the table itself
# ---------------------------------------------------------------------------


def test_a_uom_quantity_with_neither_a_path_nor_a_constructor_is_refused(tmp_path: Path) -> None:
    """A unit uom carries may not fall through to the identity.

    **The identity is a claim, and it is false for every unit uom carries**: it says
    the name is already its own SI base unit, which is true of `m` and `kg` and of
    nothing else in the table. The two ways of being right are uom's own constructor
    (a `uom:` path) and this crate's (a `rust_ctor`, returning the quantity - or, for
    a dimension uom has no quantity for, the base magnitude), so a row names one or
    the generator stops.
    """
    table = a_table()
    row = unit(table, "mm")
    # Both, because either one alone is a legitimate row: `uom` without a
    # constructor is the schema's refusal, and a constructor without `uom` is what
    # `psi` and `hp` are.
    del row["uom"]
    del row["rust_ctor"]
    with pytest.raises(SystemExit, match="names neither a uom path nor a constructor"):
        compile_table(tmp_path, table)


def test_the_repository_table_names_a_constructor_for_every_unit_uom_carries() -> None:
    """The same rule, held against the table that is actually in the tree.

    The schema enforces the half of it that is a `uom` path implying a constructor;
    asserting the whole rule here as well is what makes the other half visible in
    Python, where the generator's own error message is the only other place it
    appears. **`psi` and `hp` are the rows that look like a violation of the old
    rule and are the point of the new one**: uom carries their quantities and defines
    their units at six or seven significant figures, so the constructor is azoth's
    derivation from the definition instead - and the cross-library check is what says
    so, by failing on uom's own.
    """
    module = gen_vocabulary()
    table = a_table()
    exponents = {row["id"]: tuple(row["exponents"]) for row in table["dimensions"]}
    for row in table["units"]:
        if "uom" in row:
            assert row.get("rust_ctor"), f"{row['id']} claims a uom path and names no constructor"
        elif module.UOM_TYPES[exponents[row["dimension"]]] is not None:
            assert row.get("rust_ctor"), (
                f"{row['id']} is a uom quantity and names neither a uom path nor a "
                f"constructor, so the generated conversion would be the identity"
            )


# ---------------------------------------------------------------------------
# The unit sets
# ---------------------------------------------------------------------------


def test_a_unit_set_that_switches_a_dimension_to_another_dimensions_unit_is_refused(
    tmp_path: Path,
) -> None:
    """A set's key and its value's own dimension have to agree.

    The one of the three that would survive a glance: a set naming a *length* for
    pressure is a display showing a number with a unit beside it that does not
    measure it, which is worse than the right number in the wrong unit - a reader
    can convert the second and cannot see the first.
    """
    table = a_table()
    table["unit_sets"][1]["units"]["pressure"] = "ft"
    with pytest.raises(SystemExit, match="which is a length"):
        compile_table(tmp_path, table)


def test_a_unit_set_naming_a_unit_the_table_does_not_have_is_refused(tmp_path: Path) -> None:
    """A set is a choice among declared units and not a second declaration."""
    table = a_table()
    table["unit_sets"][1]["units"]["pressure"] = "torr"
    with pytest.raises(SystemExit, match="which is not in the vocabulary"):
        compile_table(tmp_path, table)


def test_a_unit_set_naming_something_that_is_not_a_dimension_is_refused(tmp_path: Path) -> None:
    """A key no dimension declares is a set with a switch nothing reads."""
    table = a_table()
    table["unit_sets"][1]["units"]["head"] = "m"
    with pytest.raises(SystemExit, match="which are not declared"):
        compile_table(tmp_path, table)


def test_a_display_unit_named_like_a_spec_unit_is_refused(tmp_path: Path) -> None:
    """A name in both tables is one string meaning two conversions.

    `Pa`-style names are what a spec's `unit:` is checked against and a display's names are what a
    switcher offers, so a name in both would convert by a scale for a case file and - where it is
    affine - by a shift for a reader. Refusing the overlap is what lets the front end look a unit
    up by name without asking which table it meant.
    """
    table = a_table()
    table["display_units"].append(
        {
            "id": "K",
            "dimension": "thermodynamic_temperature",
            "pint": "degC",
            "factor": 1.0,
            "offset": 273.15,
        }
    )
    with pytest.raises(SystemExit, match="both a spec-declarable unit and a display unit"):
        compile_table(tmp_path, table)


def test_a_display_unit_with_a_dimension_the_table_lacks_is_refused(tmp_path: Path) -> None:
    """The same dangling reference a unit's `dimension:` is refused for."""
    table = a_table()
    table["display_units"].append(
        {
            "id": "furlong",
            "dimension": "length_of_a_horse",
            "pint": "furlong",
            "factor": 201.168,
            "offset": 0.0,
        }
    )
    with pytest.raises(SystemExit, match="which is not declared"):
        compile_table(tmp_path, table)


def test_a_display_unit_with_a_zero_factor_is_refused(tmp_path: Path) -> None:
    """A factor of zero is not a scale: a display dividing by it has no answer."""
    table = a_table()
    table["display_units"][0]["factor"] = 0.0
    with pytest.raises(SystemExit, match="factor of zero"):
        compile_table(tmp_path, table)


def test_the_two_display_tables_are_one_set_on_both_halves() -> None:
    """The generated Rust and Python carry the table's display units, and the schema does not.

    **The third list is the point of the separation.** `specs/schema/unit.schema.json` is generated
    from `units` alone, so a display unit cannot appear in it - and that is what makes "a spec may
    not declare one" a fact about a generated file rather than a rule to remember. The Lean
    vocabulary is generated from `units` for the same reason: a spec unit makes a dimension claim
    the gate proves, and a display unit makes none.
    """
    module = gen_vocabulary()
    table = a_table()
    declared = sorted(unit["id"] for unit in table["display_units"])

    generated = [line.strip() for line in module.emit_rust(table).splitlines()]
    assert declared, "the table declares no display unit, so this test asserts nothing"
    for name in declared:
        assert f'id: "{name}",' in generated, f"{name} is not in the generated Rust"
        assert f'"{name}"' in module.emit_python(table), f"{name} is not in the generated Python"
        assert f'"{name}"' not in module.emit_unit_schema(table), (
            f"{name} is in the spec-facing enum, so a case file could declare it"
        )
        assert f'"{name}"' not in module.emit_lean(table), (
            f"{name} is in the Lean vocabulary, which makes a dimension claim it cannot prove"
        )


def test_two_unit_sets_covering_different_dimensions_are_refused(tmp_path: Path) -> None:
    """Every set names the same dimensions, so switching back is the same document.

    A dimension one set names and another does not is a display that changes what it
    shows when a reader asks for the *other* set - and there is no unit it could fall
    back to without the reader having chosen it.
    """
    table = a_table()
    del table["unit_sets"][1]["units"]["length"]
    with pytest.raises(SystemExit, match="differs from the first set"):
        compile_table(tmp_path, table)


def test_every_unit_in_the_table_has_a_distinct_pint_name_or_a_distinct_dimension() -> None:
    """No two rows are the same unit twice.

    `mm` and `m` share a dimension and differ in their `pint` name; `m` and `mm` with
    the same name and the same dimension would be one unit under two spellings, and
    the generated enum would offer a spec a choice that means nothing.
    """
    table = a_table()
    seen: dict[tuple[str, str], str] = {}
    for row in table["units"]:
        key = (row["pint"], row["dimension"])
        assert key not in seen, (
            f"{row['id']} and {seen[key]} are the same unit - both {row['pint']!r} of "
            f"dimension {row['dimension']!r}"
        )
        seen[key] = row["id"]
