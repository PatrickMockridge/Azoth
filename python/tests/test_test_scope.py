"""The test-selection tool's own rules.

`tools/test_scope.py` decides what a push runs, so a rule that drifts here is a rule that
lets something through - and in the direction that costs an hour, because the fallback when
the tool cannot be trusted is to run everything. Nothing else in the suite checks it.

What is pinned is the *shape* of the selection, not the wording of any one command: which
crates a change reaches, that the crate graph is read rather than listed, that a case file is
attributed to its model rather than to its own name, and that a path no rule claims is an
error rather than a shrug.
"""

from __future__ import annotations

import importlib
import sys
from pathlib import Path
from types import ModuleType

TOOLS = Path(__file__).resolve().parents[2] / "tools"


def tool() -> ModuleType:
    """`tools/test_scope.py`, imported by name - the tools are not a package."""
    sys.path.insert(0, str(TOOLS))
    try:
        return importlib.import_module("test_scope")
    finally:
        sys.path.pop(0)


def commands_for(paths: list[str]) -> list[str]:
    """Every command the paths reach."""
    module = tool()
    commands, unclaimed = module.scope(paths, module.crate_edges())
    assert not unclaimed, f"unclaimed: {unclaimed}"
    return [command for _, _, command in commands]


def test_a_hydraulics_change_does_not_reach_the_equation_of_state():
    """**The property the tool exists for.**

    A push that touches the hydraulics of produced water must not run the equation of state's
    736 tests, and it must still run the crates that use hydraulics - `azoth-process` among
    them. This is the pair that decides whether the tool is worth having: too narrow and it
    misses a break, too wide and it is the sweep again.
    """
    joined = " ".join(commands_for(["crates/azoth-hydraulics/src/darcy_weisbach.rs"]))
    assert "-p azoth-hydraulics" in joined
    assert "-p azoth-process" in joined, "a crate that depends on it has to be tested too"
    assert "-p azoth-eos" not in joined, "and a crate that does not is not"


def test_a_shared_crate_reaches_everything_that_reads_it():
    """`azoth-core` is the root of the graph, so a change there is as wide as it gets.

    Asserted rather than assumed because the graph is *read* from the manifests: a new member
    that depends on core has to appear here without anyone editing a list.
    """
    joined = " ".join(commands_for(["crates/azoth-core/src/units.rs"]))
    module = tool()
    for crate in module.members():
        assert f"-p {crate}" in joined, crate


def test_a_case_file_is_attributed_to_its_model():
    """A case has an `id` of its own, and it is not a registered id.

    Reading that `id` first put case names into the `-k` expression, where they select nothing
    and the run looks like a pass.
    """
    joined = " ".join(commands_for(["specs/cases/eos/co2_brine_phase_equilibrium.toml"]))
    assert "eos.co2_brine_phase_equilibrium" in joined
    assert "the_40_bara" not in joined, "a case's own name is not an id"


def test_a_path_nothing_claims_is_an_error():
    """The rule that keeps the scoped set trustworthy.

    A tool that shrugged at an unknown directory would be one whose silence means either
    "nothing to run" or "I do not understand this" - and the two are not the same answer.
    """
    module = tool()
    _, unclaimed = module.scope(["somewhere/new_thing.txt"], module.crate_edges())
    assert unclaimed == ["somewhere/new_thing.txt"]


def test_a_docs_change_runs_no_tests_but_the_doc_checks():
    joined = " ".join(commands_for(["docs/src/safety.md"]))
    assert "check_doc_claims" in joined
    assert "pytest" not in joined, "a page cannot break a model"


def test_the_nightly_lane_is_the_whole_suite_and_is_not_gated():
    """`tools/gated.sh` does not exist on a runner, so the lane that CI runs must not name it."""
    lines = tool().nightly(gate=False)
    assert any("--workspace --all-targets" in line for line in lines)
    assert any("pytest python/tests" in line for line in lines)
    assert not any("tools/gated.sh" in line for line in lines)
