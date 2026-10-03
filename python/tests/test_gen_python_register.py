"""The extension's registration is generated, and the hand-written remainder is capped.

`tools/gen_python_register.py` emits `crates/azoth-python/src/register_gen.rs` - one `add_class`
per registered result and one `add_function` per registered id, the namespace and the function
name taken from the spec's own `implementations.rust` path. `docs-drift` regenerates and diffs
it; the `--check` below runs the same gate in the suite.

The cap is the second half. What generated registration cannot reach is genuinely not derivable:
the data-transport classes, the card overlay, the stream and session bindings, and the
introspection helpers. Those stay in `lib.rs`, and this file holds their number - so a new one is
a decision someone made on purpose rather than a line that drifted in.
"""

from __future__ import annotations

import importlib
import re
import subprocess
import sys
from pathlib import Path
from types import ModuleType

REPO_ROOT = Path(__file__).resolve().parents[2]
LIB = REPO_ROOT / "crates" / "azoth-python" / "src" / "lib.rs"

#: The hand-written `add_class`/`add_function` calls `lib.rs` may carry. Measured when the
#: generator landed; lowering it is the point, and raising it needs a reason in the diff.
KNOWN_HAND_WRITTEN = {"add_class": 14, "add_function": 39}


def _tools_module(name: str) -> ModuleType:
    sys.path.insert(0, str(REPO_ROOT / "tools"))
    try:
        return importlib.import_module(name)
    finally:
        sys.path.pop(0)


def test_the_committed_registration_is_current() -> None:
    """The drift gate `docs-drift` runs, run here so it fails in the suite too."""
    result = subprocess.run(
        [sys.executable, str(REPO_ROOT / "tools" / "gen_python_register.py"), "--check"],
        capture_output=True,
        text=True,
        check=False,
    )

    assert result.returncode == 0, (
        f"the registration is out of date - run `python tools/gen_python_register.py`\n"
        f"{result.stdout}{result.stderr}"
    )


def test_the_hand_written_registrations_are_capped() -> None:
    """Every `add_class`/`add_function` left in `lib.rs` is one the generator cannot derive.

    A count rather than a set, because the numbers are the claim: the generated file carries one
    per registered id, and what remains is the list this file's docstring names. A new line here
    is a new hand-maintained registration point, which is what the tranche exists to remove.
    """
    source = LIB.read_text(encoding="utf-8")
    counted = {
        "add_class": len(re.findall(r"add_class::<", source)),
        "add_function": len(re.findall(r"add_function\(", source)),
    }

    assert counted == KNOWN_HAND_WRITTEN, (
        f"lib.rs carries {counted} hand-written registrations against the capped "
        f"{KNOWN_HAND_WRITTEN}. If the new one is genuinely underivable, raise the cap here and "
        f"say why; if it is a calculation or a result class, it belongs in the spec instead."
    )


def test_every_namespace_carries_a_registered_id() -> None:
    """The six-crate map is exercised, so a namespace that lost its ids is visible.

    `implementations()` refuses a seventh crate; this is the other direction, and it is what
    stops the map from quietly becoming a list of five with one entry no spec names.
    """
    generator = _tools_module("gen_python_register")
    used = {module for _, module, _ in generator.implementations()}

    assert used == set(generator.NAMESPACES.values()), (
        f"namespaces with no registered id: {sorted(set(generator.NAMESPACES.values()) - used)}"
    )
