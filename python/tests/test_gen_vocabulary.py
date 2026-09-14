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


def test_the_repository_table_has_no_uom_path_without_a_constructor() -> None:
    """Every row claiming a `uom` path names a hand-written constructor.

    The schema enforces this too; asserting it here as well is what makes the fact
    visible in Python, where the generator's own error message is the only other
    place it appears.
    """
    table = a_table()
    for row in table["units"]:
        if "uom" not in row:
            assert "rust_ctor" not in row, f"{row['id']} names a constructor with no uom path"
        else:
            assert row.get("rust_ctor"), f"{row['id']} claims a uom path and names no constructor"


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
