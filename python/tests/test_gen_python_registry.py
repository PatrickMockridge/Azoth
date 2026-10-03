"""The extension's id tables are generated, and this is the ratchet that says so.

`gen_python_registry` emits `crates/azoth-python/src/registry_tables_gen.rs` - the two
tables `results.rs` used to carry by hand. `result_fields` maps a calc id to its
result's public field names; `calc_ids` is the calculations this extension implements.
Both are derivable, so a hand-maintained pair was a second place to edit for every new
calculation, and forgetting it failed no gate: the extension would register a calc it
could not describe the fields of.

`docs-drift` regenerates and diffs the file, so a spec edited without regenerating is
caught in CI. `test_the_committed_registry_tables_are_current` runs the same `--check`
in the suite, and `test_the_arms_are_the_index_in_order` drives the emitter as a pure
function so an arm that silently reordered or dropped still fails here.
"""

from __future__ import annotations

import importlib
import subprocess
import sys
from pathlib import Path
from types import ModuleType

REPO_ROOT = Path(__file__).resolve().parents[2]


def _tools_module(name: str) -> ModuleType:
    sys.path.insert(0, str(REPO_ROOT / "tools"))
    try:
        return importlib.import_module(name)
    finally:
        sys.path.pop(0)


def test_the_committed_registry_tables_are_current() -> None:
    """The drift gate `docs-drift` runs, run here so it fails in the suite too.

    Through `--check` rather than by comparing the emit: the generator passes the Rust it
    emits through `rustfmt` before writing, so a comparison against the raw emitter output
    fails on formatting rather than on drift.
    """
    result = subprocess.run(
        [sys.executable, str(REPO_ROOT / "tools" / "gen_python_registry.py"), "--check"],
        capture_output=True,
        text=True,
        check=False,
    )

    assert result.returncode == 0, (
        f"the extension's id tables are out of date - run "
        f"`python tools/gen_python_registry.py`\n{result.stdout}{result.stderr}"
    )


def test_the_arms_are_the_index_in_order() -> None:
    """Every indexed result type gets exactly one arm, and the ids are sorted.

    Driving the emitter directly rather than reading the committed file: the point is that
    the *generator* pairs the index with the arms, so a reader change that dropped a type
    or reordered the list shows up here rather than only as a byte diff.
    """
    generator = _tools_module("gen_python_registry")
    rust_index = _tools_module("rust_index")

    source: str = generator.emit()
    results = rust_index.result_types()
    arms = [f'"{result.calc_id}" => {result.item_path}::FIELDS.to_vec(),' for result in results]

    positions = [source.index(arm) for arm in arms]
    assert positions == sorted(positions), "the arms are not in the index's own order"
    assert source.count("::FIELDS.to_vec(),") == len(arms), "the emitter wrote an arm more or less"

    ids = generator.calc_ids()
    assert ids == sorted(ids), "calc_ids is not sorted"
    assert all(f'"{calc_id}".to_string(),' in source for calc_id in ids)
