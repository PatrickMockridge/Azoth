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
import subprocess
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
    return [entry.command for entry in commands]


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


def test_a_target_that_selects_nothing_is_caught() -> None:
    """`--verify-selection`, measured rather than reasoned about.

    `-k` deselecting everything exits 5, so this is the check that turns "the gate ran" into
    "the gate ran something" - the distinction the map's silent zeros lived in.
    """
    module = tool()
    target = f"{module.TESTS}/test_test_scope.py"
    assert module.collects_nothing(f"{target} -k definitely_not_a_test") == "collects nothing"
    assert module.collects_nothing(target) is None


# ------------------------------------------------------ what the iterate lane may narrow


def test_the_iterate_lane_deselects_by_marker_rather_than_by_name() -> None:
    """**The list went stale; the marker cannot.** A node id in this tool goes stale the moment
    a sweep is renamed, and the 462-second layer diff was missing from the list entirely."""
    module = tool()
    command = module.Command("x", module.PYTEST, f"{module.TESTS}/eos")
    line = module.render(command, "iterate", gate=False)
    assert f'-m "not {module.MARKER}"' in line
    assert "--deselect" not in line
    push = module.render(command, "push", gate=False)
    assert module.MARKER not in push, "the push lane runs the sweeps: that is the cost it buys"


def test_a_wholly_marked_file_is_skipped_rather_than_filtered_to_nothing() -> None:
    """A filter over a file where every test is marked leaves an empty collection, and pytest
    calls an empty collection exit 5 - which the tool would read as a failure. So the iterate
    lane skips such a file whole and says why, rather than turning a skip into a red."""
    module = tool()
    line = module.render(
        module.Command("x", module.PYTEST, f"{module.TESTS}/test_dispatch.py"),
        "iterate",
        gate=False,
    )
    assert line.startswith("#"), line
    assert "skips" in line
    assert module.module_is_a_sweep(f"{module.TESTS}/test_dispatch.py")
    assert not module.module_is_a_sweep(f"{module.TESTS}/eos")


def test_the_tool_and_pytest_agree_on_which_files_sweep() -> None:
    """The tool names them and pytest deselects them, so the two must describe one set.

    Checked against pytest rather than against a second list, which is the drift this replaced.
    """
    module = tool()
    named = module.marked_sweeps([module.Command("the suite", module.PYTEST, module.TESTS)])
    files = sorted({entry.split("::")[0] for entry in named})
    assert files, "the suite has slow sweeps; if this is empty the markers were dropped"
    result = subprocess.run(
        [
            module.python_runner(),
            "-m",
            "pytest",
            "--collect-only",
            "-q",
            "-m",
            module.MARKER,
            "-o",
            "addopts=",
            "-p",
            "no:cacheprovider",
            *files,
        ],
        cwd=module.ROOT,
        capture_output=True,
        text=True,
        env=module.pytest_env(),
    )
    collected = {line.split("::", 1)[0] for line in result.stdout.splitlines() if "::" in line}
    assert collected == set(files), "the marker and the tool disagree about which files sweep"


# ------------------------------------------------------------ which jobs a diff has to run


def test_the_plan_reads_what_the_commands_actually_are() -> None:
    """**`--plan` is what every job reads**, so it has to agree with the commands.

    A Python-only change runs the Python suite and the drift job (a registered calc's reference
    source is in its provenance block) and does **not** recompile the workspace. That pair is
    the whole point of scoping twelve jobs by the diff instead of running ten of them flat.
    """
    module = tool()
    commands, _ = module.scope(["python/src/azoth/eos/flash.py"], module.crate_edges())
    jobs = dict(module.plan(commands, full=False))
    assert jobs["python-test"] and jobs["docs-drift"]
    assert not jobs["rust-test"], "a Python change does not recompile the workspace"
    assert jobs["spec-validate"], "and the cheap flat ones still run"


def test_a_command_no_job_runs_is_refused() -> None:
    """`jobs_of` refuses rather than returning nothing, because nothing is the silent answer."""
    module = tool()
    commands, _ = module.scope(["crates/azoth-eos/src/pr_kappa.rs"], module.crate_edges())
    declared = set(module.known_jobs())
    for entry in commands:
        named = set(module.jobs_of(entry))
        assert named, entry
        assert named <= declared, f"{entry.command} names a job ci.yml does not declare"
    assert "rust-test" in {job for entry in commands for job in module.jobs_of(entry)}, (
        "a crate change is what the Rust job is for"
    )


# ------------------------------------------------- a kernel change reaches its Python mirror


def test_a_rust_kernel_change_reaches_its_python_mirror() -> None:
    """**The hole that made the two-implementation promise not hold.**

    Python targets were produced only from `python/` and `specs/` paths, so a change to a Rust
    kernel ran cargo, clippy and `maturin develop` and **not one pytest**. The
    cross-implementation agreement test - the whole reason the library is written twice - never
    ran for the change most likely to break it.
    """
    joined = " ".join(commands_for(["crates/azoth-eos/src/pr_kappa.rs"]))
    assert "test_cross_impl.py" in joined, "the agreement test is the point"
    assert "test_model_cases.py" in joined
    assert '-k "eos."' in joined, "and it is the namespace, not a hundred ids nobody can check"


def test_a_crate_with_no_namespace_says_what_it_reaches() -> None:
    """Five members name no namespace; each has to say, or a change reaches nothing in silence."""
    module = tool()
    assert module.crate_claim("azoth-wasm").jobs == ("ui",)
    assert f"{module.TESTS}" in module.crate_claim("azoth-python").files
    assert len(module.crate_claim("azoth-core").namespaces) == len(module.crate_namespaces())
    assert module.crate_claim("azoth-test-support").empty(), "its users' tests carry it"
    assert module.table_problems() == []
