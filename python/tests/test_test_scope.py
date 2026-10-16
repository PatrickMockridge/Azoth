"""The test-selection tool's own rules.

`tools/test_scope.py` decides what a push runs, so a rule that drifts here is a rule that
lets something through - and in the direction that costs an hour, because the fallback when
the tool cannot be trusted is to run everything. Nothing else in the suite checks it.

What is pinned is the *shape* of the selection, not the wording of any one command: which
crates a change reaches, that the crate graph is read rather than listed, that a case file is
attributed to its model rather than to its own name, and that a path no rule claims is an
error rather than a shrug. The last section pins the pair of failures that are worse than a
refusal: a path a rule *claims and then places nowhere*, and an id a `-k` can never select.
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


def test_a_hydraulics_change_does_not_reach_the_equation_of_state() -> None:
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


def test_a_shared_crate_reaches_everything_that_reads_it() -> None:
    """`azoth-core` is the root of the graph, so a change there is as wide as it gets.

    Asserted rather than assumed because the graph is *read* from the manifests: a new member
    that depends on core has to appear here without anyone editing a list.
    """
    joined = " ".join(commands_for(["crates/azoth-core/src/units.rs"]))
    module = tool()
    for crate in module.members():
        assert f"-p {crate}" in joined, crate


def test_a_case_file_is_attributed_to_its_model() -> None:
    """A case has an `id` of its own, and it is not a registered id.

    Reading that `id` first put case names into the `-k` expression. A lone wrong id is loud -
    pytest exits 5 when `-k` deselects everything, which was measured rather than assumed - but
    a wrong id *inside an `or`* is silent, because the expression still selects the real ones.
    That is why the ids are checked against the registry, and why this pins the case's model
    rather than the case's own name.
    """
    joined = " ".join(commands_for(["specs/cases/eos/co2_brine_phase_equilibrium.toml"]))
    assert "eos.co2_brine_phase_equilibrium" in joined
    assert "the_40_bara" not in joined, "a case's own name is not an id"


def test_a_path_nothing_claims_is_an_error() -> None:
    """The rule that keeps the scoped set trustworthy.

    A tool that shrugged at an unknown directory would be one whose silence means either
    "nothing to run" or "I do not understand this" - and the two are not the same answer.
    """
    module = tool()
    _, unclaimed = module.scope(["somewhere/new_thing.txt"], module.crate_edges())
    assert unclaimed == ["somewhere/new_thing.txt"]


def test_a_docs_change_runs_no_tests_but_the_doc_checks() -> None:
    joined = " ".join(commands_for(["docs/src/safety.md"]))
    assert "check_doc_claims" in joined
    assert "pytest" not in joined, "a page cannot break a model"


def test_the_nightly_lane_is_the_whole_suite_and_is_not_gated() -> None:
    """`tools/gated.sh` does not exist on a runner, so the lane that CI runs must not name it."""
    lines = tool().nightly(gate=False)
    assert any("--workspace --all-targets" in line for line in lines)
    assert any("pytest python/tests" in line for line in lines)
    assert not any("tools/gated.sh" in line for line in lines)


# ------------------------------------------------------------ the map places, or it refuses


def test_a_namespace_without_a_test_directory_still_reaches_its_tests() -> None:
    """**The silent zero this tool shipped with.**

    `process` has no `python/tests/process/`, so a change under `python/src/azoth/process/`
    emitted *nothing* and the run printed "nothing reaches a test" and exited 0 - a green push
    that tested nothing. Its tests live at the root of the test tree and are named here.
    """
    joined = " ".join(commands_for(["python/src/azoth/process/kernels.py"]))
    assert "test_process.py" in joined
    assert "-k " in joined, "and it reaches the namespace's ids, not just the files"


def test_a_change_inside_a_namespace_reaches_that_namespaces_ids() -> None:
    """Only a `specs/` path used to contribute ids, so a standards change reached nothing."""
    joined = " ".join(commands_for(["python/src/azoth/standards/kernels.py"]))
    assert "standards." in joined, "the namespace's own ids are what verifies it"


def test_a_namespace_the_table_does_not_know_is_refused() -> None:
    """The guard that makes the table complete rather than merely written down.

    A namespace with no entry is a namespace whose tests nothing runs, and the failure has to
    be loud: the fallback when the tool cannot be trusted is to run everything.
    """
    module = tool()
    try:
        module.namespace_targets("definitely_not_a_namespace")
    except SystemExit as exc:
        assert "NAMESPACE_TESTS" in str(exc)
    else:
        raise AssertionError("a namespace with no mapping has to be refused")


def test_the_namespace_table_covers_the_tree() -> None:
    """`discovered_namespaces` is the union of three sources; the table must match it."""
    module = tool()
    assert module.table_problems() == []
    assert set(module.NAMESPACE_TESTS) == set(module.discovered_namespaces())


def test_an_id_that_is_not_registered_is_refused() -> None:
    """A `-k` naming a non-id selects nothing, and that is the one thing that looks like a pass.

    The real hazard is not a lone wrong id - pytest exits 5 for that - but a wrong id inside an
    `or` expression, where the expression still selects the valid ones.
    """
    module = tool()
    assert "eos.pt_flash" in module.registered_ids()
    assert "dimensionless" not in module.registered_ids()


def test_the_whole_tracked_tree_is_placed() -> None:
    """`--audit`, as a test: every path a rule claims, and every claim reaches something.

    This is what makes the map *complete* rather than merely plausible - a `skills/` or `ui/`
    or `lean/` change used to be refused outright, which failed the CI job that ran the tool.
    """
    assert tool().audit() == 0
