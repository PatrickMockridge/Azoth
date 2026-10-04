"""The simple wrappers are generated, and the hand-written ones are capped.

`tools/gen_python_wrappers.py` emits `crates/azoth-python/src/wrappers_gen.rs` - one
`#[pyfunction]` per id whose kernel parameter names are the spec's declared input names, which is
the plan's own "the signature *is* the declared input list". `docs-drift` regenerates and diffs
it; the `--check` below runs the same gate in the suite.

The cap is the other half. What this file cannot derive is the mixture-expanded models - 45 `eos`
ids whose kernels take a `Mixture` the spec does not declare - plus the handful whose kernel
signature carries a conversion this file has no rule for. Those stay hand-written in the
namespace modules, and this file holds their number, so a new one is a decision someone made
rather than a line that drifted in.
"""

from __future__ import annotations

import importlib
import re
import subprocess
import sys
from pathlib import Path
from types import ModuleType

REPO_ROOT = Path(__file__).resolve().parents[2]
BINDING = REPO_ROOT / "crates" / "azoth-python" / "src"

#: The hand-written `pub fn`s the namespace modules may carry. Measured when the generator
#: landed; lowering it is the point, and raising it needs a reason in the diff.
KNOWN_HAND_WRITTEN = 64


def _tools_module(name: str) -> ModuleType:
    sys.path.insert(0, str(REPO_ROOT / "tools"))
    try:
        return importlib.import_module(name)
    finally:
        sys.path.pop(0)


def _hand_written() -> list[str]:
    """Every registered id whose wrapper is still written out in a namespace module."""
    rust_index = _tools_module("rust_index")
    registered = {function for _, _, function in rust_index.implementations()}
    out = []
    for module in ("hydraulics", "eos", "thermal", "reactions", "standards", "process"):
        source = (BINDING / f"{module}.rs").read_text(encoding="utf-8")
        text = rust_index._blank_comments(source)
        out.extend(re.findall(r"^pub fn (\w+)\(\s*py: Python<'_>", text, re.M))
    return [name for name in out if name in registered]


def test_the_committed_wrappers_are_current() -> None:
    """The drift gate `docs-drift` runs, run here so it fails in the suite too."""
    result = subprocess.run(
        [sys.executable, str(REPO_ROOT / "tools" / "gen_python_wrappers.py"), "--check"],
        capture_output=True,
        text=True,
        check=False,
    )

    assert result.returncode == 0, (
        f"the wrappers are out of date - run `python tools/gen_python_wrappers.py`\n"
        f"{result.stdout}{result.stderr}"
    )


def test_the_hand_written_wrappers_are_capped() -> None:
    """Every wrapper left in a namespace module is one the generator cannot derive."""
    hand = _hand_written()

    assert len(hand) == KNOWN_HAND_WRITTEN, (
        f"the namespace modules carry {len(hand)} hand-written wrappers against the capped "
        f"{KNOWN_HAND_WRITTEN}: {sorted(set(hand))}. If the new one is genuinely underivable, "
        f"raise the cap here and say why; if its kernel's parameter names are the spec's declared "
        f"inputs, it belongs in the generated file instead."
    )


def test_no_id_is_both_generated_and_hand_written() -> None:
    """The two sets are disjoint. A wrapper emitted twice would not link, and one removed from
    the modules while the generator also refused it would be an id nothing registers."""
    wrappers = _tools_module("gen_python_wrappers")
    generated = {function for _, function in wrappers.covered_ids()}

    assert generated & set(_hand_written()) == set(), "an id is both generated and hand-written"
