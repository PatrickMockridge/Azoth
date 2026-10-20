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
import re
import subprocess
import sys
import time
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


def test_the_working_lanes_deselect_by_marker_rather_than_by_name() -> None:
    """**The list went stale; the marker cannot.** A node id in this tool goes stale the moment
    a sweep is renamed, and the 462-second layer diff was missing from the list entirely. The
    iterate and push lanes deselect; the nightly runs everything."""
    module = tool()
    command = module.Command("x", module.PYTEST, f"{module.TESTS}/eos")
    for lane in ("iterate", "push", "tranche"):
        line = module.render(command, lane, gate=False)
        assert f'-m "not {module.MARKER}"' in line, lane
        assert "--deselect" not in line, lane
    nightly = module.render(command, "nightly", gate=False)
    assert module.MARKER not in nightly, "the nightly runs the sweeps whole"


def test_a_wholly_marked_file_is_skipped_rather_than_filtered_to_nothing() -> None:
    """A filter over a file where every test is marked leaves an empty collection, and pytest
    calls an empty collection exit 5 - which the tool would read as a failure. So the working
    lanes skip such a file whole and say why, rather than turning a skip into a red."""
    module = tool()
    sweep = f"{module.TESTS}/test_degenerate_sweep.py"
    line = module.render(module.Command("x", module.PYTEST, sweep), "iterate", gate=False)
    assert line.startswith("#"), line
    assert "skips" in line
    assert module.module_is_a_sweep(sweep)
    assert not module.module_is_a_sweep(f"{module.TESTS}/test_dispatch.py"), (
        "the dispatch module's slow loop is one function, not the whole file"
    )


def test_the_tranche_narrows_a_reached_sweep_to_the_selected_ids() -> None:
    """The close-out runs the reached ids inside a sweep, not the sweep.

    A diff that reaches both a spec and core has real terms (`eos.`, `eos.pt_flash`) and a
    marked file whose tests are parametrized over ids (`test_registration_completeness`); a
    core-only diff has neither, so its sweeps are the nightly's.
    """
    module = tool()
    edges = module.crate_edges()
    mixed, _ = module.scope(
        ["specs/models/eos/ge_nrtl_phase.toml", "python/src/azoth/core/result.py"], edges
    )
    selection = module.selection_terms(mixed)
    assert "eos." in selection and "eos.ge_nrtl_phase" in selection
    entry = next(
        entry
        for entry in mixed
        if entry.kind == module.PYTEST
        and entry.command.startswith(f"{module.TESTS}/test_registration_completeness.py")
    )
    line = module.sweep_line(entry, selection, gate=False, budget_seconds=None)
    assert line is not None and f"-m {module.MARKER}" in line and '-k "' in line

    core, _ = module.scope(["python/src/azoth/core/result.py"], edges)
    core_entry = next(
        entry
        for entry in core
        if entry.kind == module.PYTEST
        and entry.command.startswith(f"{module.TESTS}/test_registration_completeness.py")
    )
    assert (
        module.sweep_line(core_entry, module.selection_terms(core), gate=False, budget_seconds=None)
        is None
    )


def test_a_narrowed_sweep_is_emitted_only_where_it_collects(tmp_path: Path) -> None:
    """A marked file whose tests cannot match the selection must not emit a `-k` line.

    `test_batch.py` carries slow-sweep tests that are unparametrized and name no spec id,
    so the tranche lane's `-m slow_sweep -k "eos."` collected nothing and pytest exited 5
    - the lane went red over an optimisation with nothing to optimise. A parametrized test
    is still emitted, because its collected ids carry the parametrize ids the selection
    terms are drawn from.
    """
    module = tool()
    unrelated = tmp_path / "test_unrelated.py"
    unrelated.write_text(
        "import pytest\n\n"
        "pytestmark = pytest.mark.slow_sweep\n\n\n"
        "def test_nothing_to_do():\n"
        "    pass\n",
        encoding="utf-8",
    )
    assert not module.file_has_selected_marked_tests(str(unrelated), ("eos.",))
    parametrized = tmp_path / "test_parametrized.py"
    parametrized.write_text(
        "import pytest\n\n"
        "@pytest.mark.slow_sweep\n"
        "@pytest.mark.parametrize('x', [1])\n"
        "def test_case(x):\n"
        "    pass\n",
        encoding="utf-8",
    )
    assert module.file_has_selected_marked_tests(str(parametrized), ("eos.",))
    named = tmp_path / "test_named.py"
    named.write_text(
        "import pytest\n\n@pytest.mark.slow_sweep\ndef test_the_case():\n    pass\n",
        encoding="utf-8",
    )
    # A term is matched against the *identifier*, and a spec id's dot cannot appear in
    # one: the parametrized path is how a dotted term reaches a test, which is the case
    # above. A plain word matches the name.
    assert module.file_has_selected_marked_tests(str(named), ("case",))
    assert not module.file_has_selected_marked_tests(str(named), ("eos.",))


def test_the_smoke_ids_do_not_open_a_registry_sweep() -> None:
    """One id per namespace is not a reached set; narrowing a sweep to it is the sweep in a
    disguise. The smoke lives in the cross-implementation command, not in `-k` over sweeps."""
    module = tool()
    commands, _ = module.scope(["python/src/azoth/core/result.py"], module.crate_edges())
    assert module.selection_terms(commands) == ""
    joined = " ".join(entry.command for entry in commands if entry.kind == module.PYTEST)
    assert "eos.pt_flash" in joined


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
    extension = module.crate_claim("azoth-python")
    assert extension.files == module.EXTENSION_TESTS
    assert extension.ids == module.SMOKE_IDS
    core = module.crate_claim("azoth-core")
    assert core.files == module.CORE_TESTS
    assert core.ids == module.SMOKE_IDS
    assert not core.namespaces, "the core suites and the smoke, not every namespace"
    assert module.crate_claim("azoth-test-support").empty(), "its users' tests carry it"
    assert module.table_problems() == []


def test_a_generated_file_at_any_depth_is_its_generators_output() -> None:
    """**The hour this tranche exists for.**

    `core/result_gen.py` is generated, but the old short-circuit fired only for a module
    directly under `python/src/azoth/`; a nested one fell through to `NAMESPACE_TESTS["core"]`,
    which named the whole suite. The scope now asks for the generator's `--check` alone.
    """
    module = tool()
    claim = module.claims("python/src/azoth/core/result_gen.py")
    assert claim.generated
    assert not claim.namespaces and not claim.files
    commands = commands_for(["python/src/azoth/core/result_gen.py"])
    assert module.TESTS not in commands, "the whole suite is not a generator check"
    assert any("gen_python_result.py --check" in command for command in commands)


def test_a_core_source_change_runs_the_core_suites_and_the_smoke() -> None:
    """The narrowed core claim: its own suites plus one id per namespace, never the suite."""
    module = tool()
    commands = commands_for(["python/src/azoth/core/result.py"])
    joined = " ".join(commands)
    for target in module.CORE_TESTS:
        assert target in joined, target
    assert "test_cross_impl.py" in joined
    for smoke in module.SMOKE_IDS:
        assert smoke in joined, smoke
    assert module.TESTS not in commands, "the whole suite is what this replaced"
    assert "test_degenerate_sweep" not in joined


def test_the_public_surface_change_runs_the_surface_set() -> None:
    module = tool()
    commands = commands_for(["python/src/azoth/__init__.py"])
    joined = " ".join(commands)
    for target in module.SURFACE_TESTS:
        assert target in joined, target
    assert module.TESTS not in commands
    assert "test_degenerate_sweep" not in joined


def test_a_binding_change_runs_the_binding_suites() -> None:
    module = tool()
    commands = commands_for(["crates/azoth-python/src/eos.rs"])
    joined = " ".join(commands)
    for target in module.EXTENSION_TESTS:
        assert target in joined, target
    assert "eos.pt_flash" in joined, (
        "the smoke, because the binding is every namespace's front door"
    )
    assert module.TESTS not in commands


# ------------------------------------------------- the tool is not defeated by its own caller


def test_a_workflow_writes_the_command_list_outside_the_tree() -> None:
    """**The redirect is a path in the diff, and the refusal lands inside it.**

    An untracked file is treated as a change, deliberately: a new kernel or model is exactly the
    file that must not escape the map. A shell redirect creates its target *before* the command
    runs, so `--job python-test > scope.sh` handed the tool an unclaimed path named after its own
    output - and the message explaining the refusal was written into the file being refused. The
    job failed with exit 2 and nothing in the log, on every push. Nothing about the tool's rules
    was wrong; the caller was writing into the tree it was asking about.
    """
    for workflow in sorted((TOOLS.parent / ".github" / "workflows").glob("*.yml")):
        for number, line in enumerate(workflow.read_text(encoding="utf-8").splitlines(), 1):
            if "scope.sh" in line and not line.lstrip().startswith("#"):
                assert "$RUNNER_TEMP/scope.sh" in line, f"{workflow.name}:{number}: {line}"


# ------------------------------------------------- a full run is a property of the commands


def test_full_emits_the_whole_surface_and_not_a_diff() -> None:
    """**`mode: full` used to change only the plan.**

    It marked every job `true`, so every job started - and each one then ran `--job <name>`
    against the diff. On a schedule `github.event.before` is empty, so the base became `HEAD~1`:
    the nightly ran four and a half minutes, executed no Python tests, and failed at its coverage
    step with `No data to report`. A job starting is not a test running, which is the same
    distinction the whole tool exists to hold.
    """
    module = tool()
    commands = module.full_surface()

    def of(job: str) -> str:
        return " ".join(entry.command for entry in commands if job in module.jobs_of(entry))

    # A pytest `Command` carries its *target*, and `render` adds the word `pytest`, so this
    # compares targets rather than rendered lines.
    python = [entry.command for entry in commands if "python-test" in module.jobs_of(entry)]
    assert module.TESTS in python, "the whole suite, not the diff's corner"
    assert "skills" in python
    assert any(command.endswith("/maturin develop") for command in python), (
        "the extension both backends need"
    )
    rust = of("rust-test")
    assert "cargo test --workspace --all-targets" in rust
    assert "cargo test --workspace --doc" in rust
    assert "cargo clippy --workspace" in rust
    drift = [entry for entry in commands if "docs-drift" in module.jobs_of(entry)]
    assert drift, "the generator sweep is what the drift job runs"
    assert all("gen_" in entry.command and "--check" in entry.command for entry in drift)


def test_the_nightly_lane_renders_the_same_list() -> None:
    """One list, so `--full --job x` and the nightly lane cannot disagree."""
    module = tool()
    rendered = module.nightly(gate=False)
    assert len(rendered) == len(module.full_surface())
    assert any("pytest python/tests" in line for line in rendered)


# ------------------------------------------------- the local gate is the push


def test_every_unconditional_ci_step_is_in_the_close_out_table() -> None:
    """**The gap three mypy errors shipped through.**

    `--run` emitted only what the diff reaches, while `spec-validate`, `rust-lint`, `python-lint`
    and `provenance` run whatever the diff says. The table is hand-kept because the commands are,
    and this is the guard that stops it drifting - in both directions, so a step added to `ci.yml`
    fails here and a command in the table that CI does not run fails too.
    """
    module = tool()
    workflow = (module.ROOT / ".github" / "workflows" / "ci.yml").read_text(encoding="utf-8")
    declared = {entry.command for entry in module.CLOSE_OUT}
    expected: set[str] = set()
    for job in module.ALWAYS_JOBS:
        block = re.search(rf"^  {re.escape(job)}:\n(.*?)(?=^  [\w-]+:|\Z)", workflow, re.M | re.S)
        assert block, f"{job} is in ALWAYS_JOBS and not in ci.yml"
        # `- run:` inline, or a `run:` line under a named step; a `run: |` block is not a
        # one-line check and is skipped by the `|` test below.
        for command in re.findall(r"^\s+(?:- )?run: (.+)$", block.group(1), re.M):
            command = command.strip()
            if command == "|":
                continue
            if not ("tools/" in command or ".venv/bin/" in command or command.startswith("cargo ")):
                continue
            if re.search(r"tools/gen_\w+\.py --check", command):
                # `generator_sweep()` reads these out of the same workflow; a second copy here
                # would be the list this table exists to avoid.
                continue
            expected.add(command)
    assert declared == expected


def test_the_cheap_checks_run_before_the_expensive_targets() -> None:
    """A stale generated file, or a line to reflow, should cost seconds and not the suite.

    Measured three times in one session: the provenance hash carried past a reflow, the palette
    fixtures, and a `cargo fmt` failure that arrived after a ninety-minute scoped run. Each was
    found at the end of the expensive half; the check that would have said so in seconds had not
    run yet.
    """
    result = subprocess.run(
        [
            sys.executable,
            str(TOOLS / "test_scope.py"),
            "--paths",
            "tools/test_scope.py",
            "--no-gate",
        ],
        capture_output=True,
        text=True,
        check=True,
    )
    lines = result.stdout.splitlines()
    assert lines[0].startswith("python tools/gen_") and lines[0].endswith("--check"), lines[0]
    # The close-out block is cheap too, and it comes before any target.
    audit = lines.index("python tools/test_scope.py --audit")
    targets = [index for index, line in enumerate(lines) if " -m pytest " in line]
    assert targets, "the diff reaches a test"
    assert audit < targets[0], "the close-out checks run before the targets"

    module = tool()
    ordered = module.cheap_first(module.full_surface())
    assert ordered[0].command.startswith("python tools/gen_")
    cheap = [entry for entry in ordered if module.is_cheap(entry)]
    assert len(cheap) > len([e for e in ordered if module.is_generator_check(e)])
    assert ordered[: len(cheap)] == cheap


# ------------------------------------------------------------ the lanes and their budgets


def _cli(*args: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, str(TOOLS / "test_scope.py"), *args],
        capture_output=True,
        text=True,
        check=False,
    )


def test_the_timeout_sits_inside_the_memory_cap() -> None:
    """`timeout` signals its own process group; putting it outside the cgroup would signal
    the shell nobody is waiting for."""
    module = tool()
    entry = module.Command("x", module.PYTEST, f"{module.TESTS}/eos")
    gated = module.render(entry, "push", gate=True, budget_seconds=900)
    assert gated.startswith("tools/gated.sh timeout --signal=TERM"), gated
    bare = module.render(entry, "push", gate=False, budget_seconds=900)
    assert bare.startswith("timeout --signal=TERM"), bare
    assert module.render(entry, "push", gate=False, budget_seconds=0).startswith("env "), (
        "zero disables the cap"
    )


def test_the_watchdog_kills_the_process_group() -> None:
    """A shell-only kill leaves xdist workers and servers running; this spawns a process and
    expects the whole group to go."""
    module = tool()
    started = time.monotonic()
    code = module.run_with_watchdog(
        [("a process that would outlive the lane", "sleep 30")], 1, "iterate"
    )
    elapsed = time.monotonic() - started
    assert code == 1
    assert elapsed < 15, elapsed


def test_the_estimate_names_the_lane_budget() -> None:
    result = _cli(
        "--paths",
        f"{tool().TESTS}/test_guard_bash.py",
        "--lane",
        "iterate",
        "--estimate",
        "--no-gate",
    )
    assert result.returncode == 0, result.stderr
    assert "estimated:" in result.stdout
    assert "iterate" in result.stdout


def test_full_forces_the_nightly_render() -> None:
    """`mode: full` must not inherit the push lane's deselection, or the nightly silently
    stops running every sweep."""
    result = _cli("--paths", f"{tool().TESTS}/test_guard_bash.py", "--full", "--no-gate")
    assert result.returncode == 0, result.stderr
    assert f'-m "not {tool().MARKER}"' not in result.stdout


def test_a_working_lane_refuses_the_nightly_run() -> None:
    result = _cli("--lane", "nightly", "--run")
    assert result.returncode == 2
    assert "print-only" in result.stderr


# ------------------------------------------------------------ the generator subset


def test_a_test_change_asks_no_generator() -> None:
    joined = " ".join(commands_for([f"{tool().TESTS}/test_guard_bash.py"]))
    assert "gen_" not in joined


def test_a_docs_change_asks_no_generator() -> None:
    joined = " ".join(commands_for(["docs/src/safety.md"]))
    assert "gen_" not in joined


def test_a_generated_file_asks_its_own_generator_and_no_other() -> None:
    joined = " ".join(commands_for(["python/src/azoth/core/result_gen.py"]))
    assert "gen_python_result.py --check" in joined
    assert "gen_python_bridge.py --check" not in joined


def test_a_reference_module_asks_only_the_provenance_generator() -> None:
    joined = " ".join(commands_for(["python/src/azoth/eos/reference/ge_nrtl_phase.py"]))
    assert "gen_registry.py --check" in joined
    assert "gen_python_bridge.py --check" not in joined


def test_the_root_conftest_is_not_the_whole_suite() -> None:
    """Every test imports it, so the nearest-directory rule named the whole suite - the
    exact 81-minute answer the lanes exist to remove. The fixtures it defines are exercised
    by the extension, backend and selector tests, and the smoke covers the namespaces."""
    module = tool()
    commands = commands_for([f"{module.TESTS}/conftest.py"])
    joined = " ".join(commands)
    assert module.TESTS not in commands
    assert f"{module.TESTS}/test_extension.py" in joined
    assert "eos.pt_flash" in joined
