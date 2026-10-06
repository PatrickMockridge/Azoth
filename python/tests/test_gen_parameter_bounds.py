"""The domain a card parameter has to lie in is declared once, and the table is total.

`tools/gen_parameter_bounds.py` reads `specs/schema/component.schema.json` and emits the bounds both
readers check - `crates/azoth-core/src/parameter_bounds_gen.rs` and
`python/src/azoth/core/_bounds_gen.py`. `docs/src/calculus/values.md` states the claim it serves:
*a datum a bank holds is inside the domain the layer declares for it, or the lookup is refused*.

Three assertions, and the third is the one this repository always needs: a table that decided
nothing would pass every check by having nothing to check.
"""

from __future__ import annotations

import importlib
import json
import subprocess
import sys
from pathlib import Path
from types import ModuleType
from typing import Any

REPO_ROOT = Path(__file__).resolve().parents[2]
SCHEMA = REPO_ROOT / "specs" / "schema" / "component.schema.json"
RUST_TABLE = REPO_ROOT / "crates" / "azoth-core" / "src" / "parameter_bounds_gen.rs"


def generator_decisions() -> list[dict[str, Any]]:
    """The schema's parameters as the generator itself reads them.

    Bound to a name rather than returned directly: the tools tree is untyped, so the call
    is `Any` and `--strict` refuses to return one from a typed function.
    """
    decisions: list[dict[str, Any]] = _tools_module("gen_parameter_bounds").decisions()
    return decisions


def _tools_module(name: str) -> ModuleType:
    sys.path.insert(0, str(REPO_ROOT / "tools"))
    try:
        return importlib.import_module(name)
    finally:
        sys.path.pop(0)


def test_the_committed_tables_are_current() -> None:
    """The gate `docs-drift` runs, run here too so a push is not the first thing to see it."""
    result = subprocess.run(
        [sys.executable, str(REPO_ROOT / "tools" / "gen_parameter_bounds.py"), "--check"],
        capture_output=True,
        text=True,
        check=False,
    )

    assert result.returncode == 0, (
        f"the bounds tables are out of date - run `python tools/gen_parameter_bounds.py`\n"
        f"{result.stdout}{result.stderr}"
    )


def test_every_parameter_states_its_domain_or_why_it_has_none() -> None:
    """**A parameter nothing checks is the thing this exists to stop.**

    Read from the schema rather than from the generated tables, so a parameter added to the format
    and forgotten here fails even if both tables were regenerated from it - the generator refuses
    it, and this says the refusal is still about the rule rather than about a table that happens to
    be empty.
    """
    schema = json.loads(SCHEMA.read_text(encoding="utf-8"))
    numeric = {
        name: spec
        for name, spec in schema["properties"].items()
        if isinstance(spec, dict) and "$ref" in spec
    }
    assert numeric, "no numeric parameter was read - this is looking in the wrong place"

    undecided = sorted(
        name
        for name, spec in numeric.items()
        if not any(key.startswith("x-azoth-min") or key.startswith("x-azoth-max") for key in spec)
        and "x-azoth-unbounded" not in spec
    )
    assert not undecided, (
        f"{undecided} state neither a bound nor `x-azoth-unbounded`, so nothing says whether the "
        f"silence is a decision. A bound is `x-azoth-min`/`x-azoth-max` with a rationale; no bound "
        f"is `x-azoth-unbounded` with the reason there is none."
    )

    bounded = sorted(name for name, spec in numeric.items() if "x-azoth-min" in spec)
    assert bounded, (
        "no parameter carries a bound at all, so the layer `values.md` specifies has no claim in "
        "it - which is the state this tranche exists to leave"
    )


def test_the_two_readers_carry_the_same_bounds() -> None:
    """Both readers are generated from one schema, and this says the emission did not diverge.

    The Python table is imported and the Rust one parsed, because the two are what the readers
    actually hold: comparing the *generator's* two outputs would pass whatever it emitted.
    """
    bounds = importlib.import_module("azoth.core._bounds_gen").PARAMETER_BOUNDS
    unbounded = importlib.import_module("azoth.core._bounds_gen").UNBOUNDED
    assert generator_decisions(), "the generator read no parameter"

    rust = RUST_TABLE.read_text(encoding="utf-8")
    for name, bound in bounds.items():
        # **Read the alias and the bound, not the layout.** The Rust table is emitted through
        # `rustfmt`, which decides for itself whether a tuple fits on one line; an assertion on
        # the pre-rustfmt shape passed until the generator started formatting its output.
        assert f'"{name}"' in rust, f"{name} is in the Python table and not the Rust one"
        # `0` and `0.0` are the same bound and different spellings; the value is what is compared.
        assert f"min: Some({float(bound.min)!r})" in rust, (
            f"{name}: the minimum differs between the two readers"
        )

    declared = len(bounds) + len(unbounded)
    schema = json.loads(SCHEMA.read_text(encoding="utf-8"))
    numeric = sum(
        1 for spec in schema["properties"].values() if isinstance(spec, dict) and "$ref" in spec
    )
    assert declared == numeric, (
        f"the tables account for {declared} parameter(s) and the schema declares {numeric} - one "
        f"is in neither table, which is the omission the generator refuses"
    )
