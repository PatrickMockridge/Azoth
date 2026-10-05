"""The Rust backend's adapters are generated, and this is the ratchet that says so.

`tools/gen_python_bridge.py` emits `python/src/azoth/_rust_bridge_gen.py` for every id whose
signature is the public wrapper's and whose boundary is the spec's own inputs; `_rust_bridge.py`
re-exports those and keeps the rest. `docs-drift` regenerates and diffs it, and the `--check`
below runs the same gate in the suite.

The hand-written cap is the other half. What the generator cannot derive is the ids whose boundary
is a `Mixture` or a record the spec does not name - and the count of those is a decision someone
made rather than a line that drifted in, so it is held here.
"""

from __future__ import annotations

import importlib
import re
import subprocess
import sys
from pathlib import Path
from types import ModuleType

REPO_ROOT = Path(__file__).resolve().parents[2]
BRIDGE = REPO_ROOT / "python" / "src" / "azoth" / "_rust_bridge.py"

#: The adapters `_rust_bridge.py` may still carry by hand. Every one of them is an id whose
#: boundary is a mixture, a record or a container the spec does not declare; lowering this is the
#: point, and raising it needs a reason in the diff.
KNOWN_HAND_WRITTEN = 22


def _tools_module(name: str) -> ModuleType:
    sys.path.insert(0, str(REPO_ROOT / "tools"))
    try:
        return importlib.import_module(name)
    finally:
        sys.path.pop(0)


def _hand_written() -> list[str]:
    """Every registered id `_rust_bridge.py` still writes out itself."""
    rust_index = _tools_module("rust_index")
    registered = {function for _, _, function in rust_index.implementations()}
    text = rust_index._blank_comments(BRIDGE.read_text(encoding="utf-8"))
    found = re.findall(r"^def (\w+)\(", text, re.M)
    return [name for name in found if name in registered]


def test_the_committed_adapters_are_current() -> None:
    """The drift gate `docs-drift` runs, run here so it fails in the suite too."""
    result = subprocess.run(
        [sys.executable, str(REPO_ROOT / "tools" / "gen_python_bridge.py"), "--check"],
        capture_output=True,
        text=True,
        check=False,
    )

    assert result.returncode == 0, (
        f"the generated adapters are out of date - run `python tools/gen_python_bridge.py`\n"
        f"{result.stdout}{result.stderr}"
    )


def test_the_hand_written_adapters_are_capped() -> None:
    """Every adapter left in `_rust_bridge.py` is one the generator cannot derive."""
    hand = _hand_written()

    assert len(hand) == KNOWN_HAND_WRITTEN, (
        f"`_rust_bridge.py` carries {len(hand)} hand-written adapters against the capped "
        f"{KNOWN_HAND_WRITTEN}: {sorted(set(hand))}. If the new one is genuinely underivable, "
        f"raise the cap here and say why; if its signature is the public wrapper's and its inputs "
        f"are the spec's, it belongs in the generated file instead."
    )


def test_no_id_is_both_generated_and_hand_written() -> None:
    """The two sets are disjoint, and together they are the registry.

    A name defined twice would shadow rather than fail, which is the quiet half of this: the
    generated one would win and the hand-written one would be dead code nobody noticed.
    """
    generator = _tools_module("gen_python_bridge")
    generated = {calc_id.rpartition(".")[2] for calc_id, _, _, _ in generator.covered()[0]}
    registered = {function for _, _, function in _tools_module("rust_index").implementations()}

    both = generated & set(_hand_written())
    assert not both, f"an id is both generated and hand-written: {sorted(both)}"
    assert generated | set(_hand_written()) == registered, (
        "the generated and hand-written sets do not account for the registry: "
        f"{sorted(registered - generated - set(_hand_written()))}"
    )
