"""The Rust result-type index, held to the tree it reads and the specs it must account for.

`tools/rust_index.py` is where the registration generators learn a Rust struct's field names, and
it reads them from the source because no spec carries them: `CalcResult::FIELDS` holds the public
name - a symbol on the `eos` flashes - while the struct spells the word out. A mispaired index
would have every generator downstream emit a field the struct does not have; the compiler would
catch that for the ids that get generated, and these assert it for all of them, before generation.
"""

from __future__ import annotations

import importlib
import sys
from pathlib import Path
from types import ModuleType

from azoth._models_gen import MODELS
from azoth._registry_gen import CALCS

REPO_ROOT = Path(__file__).resolve().parents[2]

#: Every `FIELDS` name that is **not** the struct's own, and the spelling the struct uses.
#: `T` is `temperature` on the heat-specified flashes and the shorter `t` on the three whose
#: struct declares it that way (`pv_reflux_flash`, `pvf_flash`, `vu_flash_single_comp`), so a
#: per-symbol rule would be wrong: this is the measured set of pairs, and the count beside it.
SYMBOL_PAIRS = {("P", "pressure"), ("T", "t"), ("T", "temperature"), ("V", "v")}
SYMBOL_FIELDS = 20

#: The flashes a `macro_rules!` writes rather than a literal `impl CalcResult`.
MACRO_FLASHES = ("eos.th_flash", "eos.ts_flash", "eos.tu_flash")


def rust_index() -> ModuleType:
    sys.path.insert(0, str(REPO_ROOT / "tools"))
    try:
        return importlib.import_module("rust_index")
    finally:
        sys.path.pop(0)


def registry_outputs(calc_id: str) -> set[str]:
    """The `outputs` the spec for `calc_id` declares."""
    for spec in (*CALCS, *MODELS):
        if spec["id"] == calc_id:
            return set(spec["outputs"])
    raise AssertionError(f"{calc_id} is in no registry entry")


def test_the_index_accounts_for_every_registered_id() -> None:
    """Neither a spec without a Rust type nor a Rust type without a spec."""
    types = rust_index().result_types()
    declared = {spec["id"] for spec in (*CALCS, *MODELS)}
    assert {result.calc_id for result in types} == declared


def test_every_type_is_indexed_once() -> None:
    """One type per id, which is what makes the pairing unambiguous."""
    types = rust_index().result_types()
    assert len({result.rust_name for result in types}) == len(types)
    assert len({result.item_path for result in types}) == len(types)


def test_a_symbol_field_is_paired_with_the_word_its_struct_spells() -> None:
    """The `eos` flashes are the whole reason the index reads the source at all, and they are the
    only types in the tree whose `FIELDS` name is not the struct's own."""
    differing = [
        (public, rust)
        for result in rust_index().result_types()
        for public, rust in result.fields
        if public != rust
    ]
    assert len(differing) == SYMBOL_FIELDS, (
        f"{len(differing)} differing field(s), not {SYMBOL_FIELDS}"
    )
    assert set(differing) == SYMBOL_PAIRS


def test_the_macro_generated_flashes_are_indexed() -> None:
    """Three types have no literal `impl` to scan, so they are resolved through the macro."""
    types = {result.calc_id: result for result in rust_index().result_types()}
    for calc_id in MACRO_FLASHES:
        assert calc_id in types, f"{calc_id} is not indexed"
        assert dict(types[calc_id].fields)["P"] == "pressure"


def test_the_public_fields_are_the_specs_own_outputs() -> None:
    """**The pairing is only meaningful if the public half is the spec's.**
    `FIELDS` is the spec's `outputs` plus `warnings`, so anything else in either direction is a
    misread array or a spec the type has not caught up with.
    """
    for result in rust_index().result_types():
        public = set(result.public_fields) - {"warnings"}
        declared = registry_outputs(result.calc_id)
        assert public == declared, f"{result.calc_id}: {sorted(public ^ declared)}"


def test_the_unit_helpers_read_the_vocabulary() -> None:
    """The two accessors the generators convert inputs and outputs with."""
    index = rust_index()
    assert index.rust_ctor("Pa") == "pascals"
    assert index.rust_ctor("dimensionless") is None
    assert index.uom_quantity("Pa") == "Pressure"
    # `uom` has no molar-flow quantity, and no crate outside it can add one.
    assert index.uom_quantity("mol/s") is None
