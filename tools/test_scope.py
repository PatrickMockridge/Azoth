#!/usr/bin/env python3
"""What a change reaches, as commands.

    python tools/test_scope.py                     # against origin/main...HEAD
    python tools/test_scope.py --staged            # against the index
    python tools/test_scope.py --explain           # and say why each command is there
    python tools/test_scope.py --run               # run them, in order
    python tools/test_scope.py --lane nightly      # the whole thing, either language
    python tools/test_scope.py --audit             # every tracked path, through the table

**Why this exists.** The rule "test what the change touches" has been written down three
times - 2026-09-21, 09-22, 10-13 - and enforced once by a PreToolUse hook, and it still
degraded, because nothing in the tree answered the question it asks. There was no command
that said *this diff reaches these tests*; there was only a hook saying no to the obvious
wrong answer. So every session re-derived the boundary by hand, got it slightly wrong once,
and fell back to running everything - an hour, killed by the system about half the time.

The rule is one function, [`claims`]. **A path that matches no rule is an error**, not a
shrug - and so is a path a rule *claims and then declines to place*, which is worse, because
it comes back green. Seven namespaces did exactly that until `--audit` walked the tree.

# The lanes

`iterate` is the inner loop, and is what `--run` uses when no lane is named: the scoped set
with the `slow_sweep`-marked sweeps deselected **by marker**. A marked test is a registry-wide
sweep, so `-k` cannot reach inside it; the marker lives on the test rather than in a list here
because a list went stale, and could not see the 462-second layer diff at all until someone
timed the suite.

`push` is what CI runs for a scoped diff. It also deselects the marked sweeps - the nightly is
where they run - and adds the five unconditional jobs' commands.

`tranche` is the close-out: the push set with each reached sweep narrowed to the ids the scope
selected (`-m slow_sweep -k "<ids>"`). A core or binding change selects one smoke id per
namespace, not a reached set, so it skips the sweeps and leaves them to the nightly.

`nightly` is everything in both languages, run by `.github/workflows/nightly.yml` on a
24-hour cycle. The full set is the nightly's job and is not something anyone types locally;
`--full` forces this render.

# What this costs, measured 2026-10-16

    gen_python_wrappers --check    7m27s -> 22.9s   (`record_fields` memoised)
    gen_python_register --check    7m26s -> 24.1s   (it calls the same)
    gen_python_bridge  --check     0m59s
    gen_python_result  --check     0m07s
    gen_models / gen_registry / gen_docs / gen_stub / gen_python_transport /
    gen_python_registry            0m01s - 0m04s each
    the whole generator sweep      ~2m40s

So the sweep runs whenever a generated file could have moved, which is exact and is minutes.
"""

from __future__ import annotations

import argparse
import ast
import contextlib
import os
import re
import signal
import subprocess
import sys
import time
import tomllib
from dataclasses import dataclass
from functools import cache, lru_cache
from pathlib import Path
from typing import NamedTuple

ROOT = Path(__file__).resolve().parent.parent
CRATES = ROOT / "crates"
SPECS = ROOT / "specs"
TOOLS = ROOT / "tools"

# The crate/namespace map and the id registry are *read*, never re-typed: `rust_index` is what
# the binding generator itself uses, so a namespace added there is covered here with no second
# edit, and an id list copied into this file would be a second list to drift.
sys.path.insert(0, str(TOOLS))
import rust_index  # noqa: E402

#: Where the Python tests live, as a path prefix and as a pytest argument.
TESTS = "python/tests"

#: The marker a test carries when it loops over every registered id and so cannot be narrowed
#: by `-k`. **Registered in `pyproject.toml`, and the iterate lane deselects it by marker** -
#: not by a list of node ids here, which goes stale in silence and missed the 462-second layer
#: diff entirely until someone timed the suite.
MARKER = "slow_sweep"

#: A command's kind, which decides how it is wrapped and whether the lane can narrow it.
CARGO, CLIPPY, PYTEST, SCRIPT = "cargo", "clippy", "pytest", "script"

#: Where a tool can live: the generators and checkers, and the hooks.
TOOL_DIRS = (ROOT / "tools", ROOT / ".claude" / "hooks")


class Command(NamedTuple):
    """One command, **and the CI job that runs it**.

    The job is a field rather than a second table because the skip decision and the run decision
    have to be one computation: `--plan` says which jobs run, `--job` says what one of them runs,
    and both read this. A second list would be a second thing to keep in step, and the two
    drifting apart is how a job ends up skipped while a command waits for it.
    """

    reason: str
    kind: str
    command: str
    jobs: tuple[str, ...] = ()


#: The jobs that run whatever the diff says, and why each is allowed to. Everything else is
#: reached by a tag on a command, so a job the diff does not reach is *skipped* - visibly, and
#: with the reason in the `scope` job's own summary.
ALWAYS_JOBS = (
    # Reads the tree and writes nothing; stdlib only, no build.
    "spec-validate",
    # `cargo fmt --all --check` is seconds and has no inputs narrower than the tree. The
    # workspace clippy lives in the nightly; the diff's own clippy comes from `rust-test`.
    "rust-lint",
    # ruff and mypy: no build, and mypy is whole-package by nature.
    "python-lint",
    # Runs `python tools/provenance.py`, which is git and stdlib.
    "provenance",
    # The job that reads this file's `--plan`; it is what decides all the others.
    "scope",
)


#: **What CI runs on every push, as commands a local `--run` can execute.**
#:
#: `--run` emits only what the diff reaches, but the five `ALWAYS_JOBS` run whatever it says -
#: and until this table existed they contributed no local command at all. A change could pass
#: `--run` and fail `python-lint` on the push; three mypy errors shipped through exactly that gap.
#:
#: The `gen_*.py --check` lines in `spec-validate` are **not** repeated here - `generator_sweep()`
#: already reads them out of the workflow. The table is hand-kept because the commands are, and the
#: test beside it reads `ci.yml` and refuses a step in either direction, the same hand-list-plus-
#: guard shape `generator_sweep()` and `discovered_namespaces()` use.
CLOSE_OUT: tuple[Command, ...] = (
    Command("the schemas and the semantic rules", SCRIPT, "python tools/spec_lint.py"),
    Command("the prose standard", SCRIPT, "python tools/prose_lint.py"),
    Command("duplicate JSON keys", SCRIPT, "python tools/check_json_keys.py"),
    Command("the numerical-safety rule", SCRIPT, "python tools/check_numerics.py"),
    Command(
        "every executable claim names its falsifier",
        SCRIPT,
        "python tools/check_claims.py",
    ),
    Command("the front end's wire shapes", SCRIPT, "python tools/check_ui.py"),
    Command("the vendoring record", SCRIPT, "python tools/check_manifest.py"),
    Command(
        "the manifest's unread rows, against the tree",
        SCRIPT,
        "python tools/audit_manifest.py --check",
    ),
    Command("the two version numbers", SCRIPT, "python tools/check_versions.py"),
    Command("the vendored MCP schemas", SCRIPT, "python tools/check_mcp_schema.py"),
    Command("the numbers the normative pages state", SCRIPT, "python tools/check_doc_claims.py"),
    Command("the porting refusals, in all three places", SCRIPT, "python tools/check_unported.py"),
    Command("the skill catalog", SCRIPT, "python tools/validate_skills.py"),
    Command("the guard theorems' ownership", SCRIPT, "python tools/check_guards.py"),
    Command("the pages against the Lean tree", SCRIPT, "python tools/audit_lean_claims.py --check"),
    # `.venv/bin/` is rewritten to this box's venv by `render`, so the command reads as CI
    # writes it.
    Command("the workspace's formatting", SCRIPT, "cargo fmt --all --check", jobs=("rust-lint",)),
    Command("ruff", SCRIPT, ".venv/bin/ruff check", jobs=("python-lint",)),
    Command("ruff's formatting", SCRIPT, ".venv/bin/ruff format --check", jobs=("python-lint",)),
    Command("the strict type check", SCRIPT, ".venv/bin/mypy", jobs=("python-lint",)),
    Command("the provenance record", SCRIPT, "python tools/provenance.py --out provenance.json"),
    Command("every tracked path claims something", SCRIPT, "python tools/test_scope.py --audit"),
)


# ------------------------------------------------------------------------ what a path claims
#
# **One table, one dispatch, one answer per path.** The first version classified a path in
# `scope()` and classified it a *second* time in `python_scope()`, and the second chain was
# weaker than the first: it emitted a namespace's tests only `if python/tests/<namespace>/`
# existed. Seven namespaces have no such directory - `process`, `reactions`, `standards`,
# `characterization`, `batch`, `properties`, `core` - so a change to any of them was claimed
# by the first chain, dropped by the second, and came back **green having tested nothing**.
# A path no rule claims is at least an error; a path a rule claims and then declines to place
# is worse, because it looks like a pass. That asymmetry is why `claims()` exists and why
# `--audit` walks the whole tree through it.


@dataclass(frozen=True)
class Claim:
    """What one changed path reaches. Every field defaults to "nothing"."""

    crates: tuple[str, ...] = ()
    namespaces: tuple[str, ...] = ()
    ids: tuple[str, ...] = ()
    files: tuple[str, ...] = ()
    tools: tuple[str, ...] = ()
    jobs: tuple[str, ...] = ()
    commands: tuple[Command, ...] = ()
    generated: bool = False
    data: bool = False
    docs: bool = False
    config: bool = False
    #: A path that reaches no test **and that no test should reach** - a licence, a gitignore.
    #: The distinction from an empty claim is the whole point: this one has to say why.
    inert: str = ""

    def empty(self) -> bool:
        """True when the claim places the path nowhere and names no reason."""
        return not (
            self.crates
            or self.namespaces
            or self.ids
            or self.files
            or self.tools
            or self.jobs
            or self.commands
            or self.generated
            or self.data
            or self.docs
            or self.config
            or self.inert
        )


#: **The wide claims, narrowed.** `core`, `__init__`, `azoth-core` and `azoth-python` used to
#: name `python/tests` or every namespace, which made an ordinary diff run the whole suite
#: (81-116 minutes). These are the tests that actually name the changed surface; the smoke
#: command below keeps one cross-implementation case per namespace, and the whole tree still
#: runs nightly.
CORE_TESTS: tuple[str, ...] = (
    f"{TESTS}/test_units_contract.py",
    f"{TESTS}/test_units_conversion.py",
    f"{TESTS}/test_units_cross_library.py",
    f"{TESTS}/test_dispatch.py",
    f"{TESTS}/test_solver_contract.py",
    f"{TESTS}/test_solver.py",
    f"{TESTS}/test_range.py",
    f"{TESTS}/test_core_numerics.py",
    f"{TESTS}/test_result_json.py",
    f"{TESTS}/test_provenance.py",
    f"{TESTS}/test_registration_completeness.py",
    f"{TESTS}/test_gen_python_result.py",
    f"{TESTS}/test_gen_python_record_outputs.py",
    f"{TESTS}/test_extension.py",
)

#: `python/src/azoth/__init__.py` is the public surface: registration, dispatch, the units
#: gate and the extension.
SURFACE_TESTS: tuple[str, ...] = (
    f"{TESTS}/test_extension.py",
    f"{TESTS}/test_dispatch.py",
    f"{TESTS}/test_registration_completeness.py",
    f"{TESTS}/test_registry_contract.py",
    f"{TESTS}/test_units_contract.py",
)

#: The root `conftest.py` is imported by every test, so a change to it is genuinely as wide as
#: the suite - and "as wide as the suite" is the 81-minute answer the lanes exist to remove. It
#: gets the tests that exercise what it defines (the extension guard, the backend fixtures, the
#: selector's own rules) plus the smoke; the whole tree is the nightly's.
CONFTEST_TESTS: tuple[str, ...] = (
    f"{TESTS}/test_extension.py",
    f"{TESTS}/test_dispatch.py",
    f"{TESTS}/test_test_scope.py",
)

#: The binding crate's own tests: the extension, the generated wrappers and the registries.
EXTENSION_TESTS: tuple[str, ...] = (
    f"{TESTS}/test_extension.py",
    f"{TESTS}/test_dispatch.py",
    f"{TESTS}/test_registration_completeness.py",
    f"{TESTS}/test_batch.py",
    f"{TESTS}/test_units_contract.py",
    f"{TESTS}/test_solver_contract.py",
    f"{TESTS}/test_gen_python_bridge.py",
    f"{TESTS}/test_gen_python_register.py",
    f"{TESTS}/test_gen_python_registry.py",
    f"{TESTS}/test_gen_python_transport.py",
    f"{TESTS}/test_gen_python_wrappers.py",
)

#: **One registered id per namespace**, the smoke that a core or binding change has not broken
#: a namespace end to end. `-k` over these names the cases and the cross-implementation
#: comparison, in one command each, without opening the namespace's directory.
SMOKE_IDS: tuple[str, ...] = (
    "eos.pt_flash",
    "eos.bubble_pressure",
    "hydraulics.darcy_weisbach",
    "hydraulics.control_valve_cv",
    "process.compressor",
    "thermal.conduction_plane_wall",
    "reactions.chemical_equilibrium",
    "characterization.delumping",
    "standards.iso6976",
)

#: The wall-clock budget each lane is held to, in minutes. The tool enforces the local run
#: with a process-group watchdog and renders a matching `timeout` for CI; `--max-minutes`
#: overrides. The numbers are the point of the lanes: an edit is minutes, not hours.
LANE_BUDGET_MINUTES: dict[str, int] = {
    "iterate": 5,
    "push": 15,
    "tranche": 20,
    "nightly": 45,
}

#: **Which generators a diff can make stale.** Conservative by design: a miss here is caught
#: by CI's unconditional `spec-validate` and the drift job, so the cost is one stale file
#: found on the push rather than one merged. The alternative - the full 23-generator sweep on
#: every edit - spends 2m40s of a 5-minute budget before a test runs.
SPEC_GENERATORS = frozenset(
    {
        "gen_dispatch",
        "gen_docs",
        "gen_guard_theorems",
        "gen_model_inputs",
        "gen_models",
        "gen_palette",
        "gen_parameter_bounds",
        "gen_python_bridge",
        "gen_python_record_outputs",
        "gen_python_records",
        "gen_python_register",
        "gen_python_registry",
        "gen_python_result",
        "gen_python_transport",
        "gen_python_wrappers",
        "gen_reaction_data",
        "gen_registry",
        "gen_rust_model_records",
        "gen_stub",
        "gen_vocabulary",
    }
)
RUST_GENERATORS = SPEC_GENERATORS
#: A hand-written Python reference module is hashed into the provenance block `gen_registry`
#: emits; nothing else generated reads it.
PYTHON_GENERATORS = frozenset({"gen_registry"})
DATA_GENERATORS = frozenset({"gen_databank", "gen_keycard"})
CAPTURE_GENERATORS = frozenset({"gen_neqsim_cases"})

#: A tree that cannot make any generated file stale.
QUIET_TREES = (
    "docs/",
    f"{TESTS}/",
    "CLAUDE.md",
    "README.md",
    "ROADMAP.md",
    ".claude/",
    "skills/",
    "agents/",
    "ui/",
    "lean/",
    ".github/",
    # The Python toolchain's own configuration: it changes how the tests run, not what the
    # generators emit. `Cargo.toml`/`Cargo.lock` stay out of this list - a member or a
    # dependency change can move the registry the generators read.
    "pyproject.toml",
    "uv.lock",
)

#: namespace -> the tests a change **in that namespace** runs. Every entry names files that
#: exist, checked by [`namespace_targets`]; an entry that is empty (`standards`) is one whose
#: coverage is its registered ids, which is stated here rather than discovered later.
NAMESPACE_TESTS: dict[str, tuple[str, ...]] = {
    "eos": (f"{TESTS}/eos",),
    "hydraulics": (f"{TESTS}/hydraulics",),
    "thermal": (f"{TESTS}/thermal",),
    "process": (
        f"{TESTS}/test_process.py",
        f"{TESTS}/test_process_layer_diff.py",
        f"{TESTS}/test_designer_report.py",
        f"{TESTS}/test_mechanical_report.py",
    ),
    "reactions": (f"{TESTS}/test_reaction_tables.py",),
    "characterization": (
        f"{TESTS}/test_characterization_front_door.py",
        f"{TESTS}/test_characterize_to_reference.py",
        f"{TESTS}/test_delumping.py",
    ),
    # Nothing imports `azoth.standards` by name: its ids reach their tests through the
    # registration sweeps, which is what a standards change is verified by.
    "standards": (),
    "batch": (f"{TESTS}/test_batch.py",),
    "properties": (f"{TESTS}/test_properties.py", f"{TESTS}/test_data_agreement.py"),
    # `azoth.core` is the Python mirror of `azoth-core`: units, dispatch and the record
    # machinery that every other namespace imports. The tests that name it are the unit,
    # dispatch, solver and record suites; the namespaces are covered by the smoke ids.
    "core": CORE_TESTS,
}

#: module directly under `python/src/azoth/` (not a namespace) -> the tests that name it.
MODULE_TESTS: dict[str, tuple[str, ...]] = {
    "__init__": SURFACE_TESTS,
    "keycard": (
        f"{TESTS}/test_keycard_loader.py",
        f"{TESTS}/test_keycard_writer.py",
        f"{TESTS}/test_keycard_schema.py",
    ),
    "_data": (f"{TESTS}/test_data.py", f"{TESTS}/test_data_agreement.py"),
    "_dispatch": (f"{TESTS}/test_dispatch.py", f"{TESTS}/test_cross_impl.py"),
    "_rust_bridge": (f"{TESTS}/test_registration_completeness.py",),
}

#: The `specs/` trees that carry no id, and what a change in each reaches.
SPEC_TREES: tuple[tuple[str, Claim], ...] = (
    (
        "specs/schema/",
        Claim(
            commands=(
                Command(
                    "a schema is what the specs are checked against",
                    SCRIPT,
                    "python tools/spec_lint.py",
                ),
                Command(
                    "a key the schema forbids is still a key",
                    SCRIPT,
                    "python tools/check_json_keys.py",
                ),
            ),
            files=(f"{TESTS}/test_keycard_schema.py", f"{TESTS}/test_spec_lint.py"),
            generated=True,
        ),
    ),
    (
        "specs/vocabulary/",
        Claim(
            files=(
                f"{TESTS}/test_units_contract.py",
                f"{TESTS}/test_units_conversion.py",
                f"{TESTS}/test_units_cross_library.py",
                f"{TESTS}/test_gen_vocabulary.py",
            ),
            generated=True,
        ),
    ),
    (
        "specs/unit_ops/",
        Claim(
            commands=(Command("the palette the UI reads", SCRIPT, "python tools/check_ui.py"),),
            files=(f"{TESTS}/test_check_ui.py",),
            generated=True,
        ),
    ),
    (
        "specs/flowsheets/",
        Claim(
            commands=(
                Command(
                    "the fixture a flowsheet produces, which only the CLI can make",
                    SCRIPT,
                    "python tools/gen_ui_fixtures.py --check",
                ),
            ),
            files=(f"{TESTS}/test_gen_ui_fixtures.py",),
            generated=True,
        ),
    ),
)

#: Top-level paths, each with what it reaches. The audit refuses any tracked path that is not
#: here and not matched by an earlier rule - which is the list a diff used to be refused by.
ROOT_CLAIMS: dict[str, Claim] = {
    "skills.toml": Claim(
        commands=(
            Command(
                "the skills are validated, not just present",
                SCRIPT,
                "python tools/validate_skills.py",
            ),
        ),
        files=(f"{TESTS}/test_validate_skills.py",),
    ),
    "keycard.example.toml": Claim(
        files=(
            f"{TESTS}/test_keycard_loader.py",
            f"{TESTS}/test_keycard_writer.py",
            f"{TESTS}/test_user_data_template.py",
        ),
    ),
    "SPEC.md": Claim(
        commands=(
            Command(
                "a normative page states counts and paths the tree must match",
                SCRIPT,
                "python tools/check_doc_claims.py",
            ),
        ),
    ),
    "skills/": Claim(
        commands=(
            Command(
                "the skills are validated, not just present",
                SCRIPT,
                "python tools/validate_skills.py",
            ),
            Command("the skills' own smoke checks", PYTEST, "skills"),
        ),
    ),
    "tools/gated.sh": Claim(inert="the memory wrapper the hook requires; no test reads it"),
    # **These run somewhere this box is not.** They are placed in a job rather than silently
    # dropped, so a diff that reaches them is a job that runs rather than a shrug.
    "ui/": Claim(jobs=("ui",)),
    "lean/": Claim(jobs=("lean-verify",)),
    "AGENTS.md": Claim(inert="prompt-layer prose; no test reads it"),
    "agents/": Claim(inert="prompt-layer prose; no test reads it"),
    ".opencode/": Claim(inert="an export target for the agent definitions; no test reads it"),
    ".editorconfig": Claim(inert="editor settings; no test reads it"),
    ".gitignore": Claim(inert="which files git ignores; no test reads it"),
    "LICENSE": Claim(inert="a licence; no test reads it"),
    "LICENSE-CC-BY-4.0": Claim(inert="a licence; no test reads it"),
    "CODE_OF_CONDUCT.md": Claim(inert="community prose; no test reads it"),
    "SECURITY.md": Claim(inert="community prose; no test reads it"),
    "NOTICE": Claim(inert="attribution prose; no test reads it"),
}

#: Files directly under `python/src/azoth/` that are not modules and reach no test.
INERT_MODULES: dict[str, str] = {
    "py.typed": "a PEP 561 marker; no test reads it",
}

#: Generators that are not run as a bare `--check`, and why. Every other `gen_*.py` has to
#: appear in CI's list or [`generator_sweep`] refuses - an unchecked generator that nothing
#: names by hand is how a stale generated file ships, which is not hypothetical: two did.
UNCHECKABLE = {
    "gen_user_data": "it takes a target path, so it is not a `--check` generator",
}

#: Generators that are checked, but **not by the sweep**: each needs something a plain
#: `--check` runner does not have, so a claim that reaches one emits it explicitly.
CARGO_GENERATORS = ("gen_ui_fixtures",)

#: Checkers that cannot be run **here**, and why. Each needs something this box does not have -
#: a token, or an artefact an instrumented run writes - so a local invocation could only fail.
#: Their own test is what a change to one reaches.
UNRUNNABLE: dict[str, str] = {
    "check_nightly": "it asks the GitHub API; .github/workflows/liveness.yml runs it",
    "check_coverage": "it reads a coverage document the instrumented run writes; CI runs it",
}


@lru_cache(maxsize=1)
def registered_ids() -> frozenset[str]:
    """Every id the specs register - **the same reader the binding generator uses**.

    Not a second list: `rust_index.spec_ids()` is what the extension is built from, so the two
    cannot disagree. An id that is not one of these cannot select a test, and the one thing a
    `-k` that selects nothing looks like is a pass.
    """
    return frozenset(rust_index.spec_ids())


@lru_cache(maxsize=1)
def crate_namespaces() -> dict[str, str]:
    """Cargo package -> namespace: `azoth-eos` -> `eos`, from the map the binding keeps."""
    return {
        crate.replace("_", "-"): namespace for crate, namespace in rust_index.NAMESPACES.items()
    }


#: The members that name no namespace, kept as a tuple so [`table_problems`] can refuse a new
#: crate that is in neither this list nor the namespace map. A crate that reaches nothing has to
#: say so, rather than reach nothing in silence.
NAMELESS_CRATES = ("azoth-cli", "azoth-core", "azoth-python", "azoth-test-support", "azoth-wasm")


def crate_claim(crate: str) -> Claim:
    """The Python side a change to **this crate** reaches.

    **The mirror, not the sweep.** A change to `crates/azoth-eos` runs the `eos` namespace's
    Python tests, because those are the other half of the same calculation - and until this
    existed, a kernel change ran no pytest at all, so the cross-implementation agreement test,
    which is the whole reason the library is written twice, never ran for the change most likely
    to break it.

    A crate that merely *depends* on the changed one is covered on the Rust side by [`reaches`]
    and does not drag its Python mirror in. That is the difference between targeted and
    everything.
    """
    namespace = crate_namespaces().get(crate)
    if namespace is not None:
        return Claim(namespaces=(namespace,))
    if crate == "azoth-core":
        # Units, the shared types and the databank: the core suites, plus one registered id per
        # namespace so a break in a namespace path is not silent. It used to name every
        # namespace, which ran the whole Python suite for a change to a shared helper.
        return Claim(files=CORE_TESTS, ids=SMOKE_IDS)
    if crate == "azoth-python":
        # The extension itself: its own suites and the smoke ids, not the whole surface.
        return Claim(files=EXTENSION_TESTS, ids=SMOKE_IDS)
    if crate == "azoth-cli":
        return Claim(files=(f"{TESTS}/validation/test_mcp_conformance.py",))
    if crate == "azoth-wasm":
        return Claim(jobs=("ui",))
    if crate == "azoth-test-support":
        # A dev-dependency: `reaches` puts every crate that uses it in the Rust set, and their
        # namespaces are theirs to claim.
        return Claim()
    raise SystemExit(
        f"test_scope: the member {crate!r} names no namespace and is in no crate entry, so a "
        "change to it would reach no Python test. Add one to NAMELESS_CRATES and to "
        "`crate_claim`; seven namespaces once reached nothing this way."
    )


@lru_cache(maxsize=1)
def discovered_namespaces() -> frozenset[str]:
    """Every namespace three independent sources know about, unioned.

    The package directories, the crates the specs' Rust paths name, and the trees under
    `specs/`. Any one of them alone would miss a namespace the other two carry, and the union
    is what [`table_problems`] checks the mapping against.
    """
    package = {
        path.name
        for path in (ROOT / "python" / "src" / "azoth").iterdir()
        if path.is_dir() and not path.name.startswith("__")
    }
    specs = {
        path.name
        for kind in ("models", "calcs")
        for path in (SPECS / kind).iterdir()
        if path.is_dir()
    }
    return frozenset(package | specs | set(crate_namespaces().values()))


def namespace_targets(namespace: str) -> tuple[str, ...]:
    """The tests a namespace change runs, **refused rather than empty**.

    A namespace the table does not know is an error. That is the same rule [`generator_sweep`]
    applies to a generator nobody checks, and for the same reason: the thing that is not there
    is invisible, and the fallback when the tool cannot be trusted is to run everything.
    """
    if namespace not in NAMESPACE_TESTS:
        raise SystemExit(
            f"test_scope: the namespace {namespace!r} has no NAMESPACE_TESTS entry, so a change "
            f"in it reaches nothing and would go green. Known: {sorted(NAMESPACE_TESTS)}. Add a "
            "mapping; a namespace with no tests is a namespace whose tests nothing runs."
        )
    targets = NAMESPACE_TESTS[namespace]
    for target in targets:
        if not (ROOT / target).exists():
            raise SystemExit(
                f"test_scope: the namespace {namespace!r} names {target}, which is not in the "
                "tree. Fix the mapping or delete the entry."
            )
    if not targets and not ids_of(namespace):
        raise SystemExit(
            f"test_scope: the namespace {namespace!r} maps to no test file and declares no id"
        )
    return targets


def table_problems() -> list[str]:
    """Everything wrong with the tables, as messages. Empty is the only good answer."""
    problems: list[str] = []
    known = discovered_namespaces()
    for namespace in sorted(known - set(NAMESPACE_TESTS)):
        problems.append(
            f"the namespace {namespace!r} is in the tree and in no NAMESPACE_TESTS entry"
        )
    for namespace in sorted(set(NAMESPACE_TESTS) - known):
        problems.append(
            f"NAMESPACE_TESTS names {namespace!r}, which is not a namespace in the tree"
        )
    for namespace, targets in sorted(NAMESPACE_TESTS.items()):
        for target in targets:
            if not (ROOT / target).exists():
                problems.append(f"NAMESPACE_TESTS[{namespace!r}] names {target}, which is absent")
    for module, targets in sorted(MODULE_TESTS.items()):
        for target in targets:
            if not (ROOT / target).exists():
                problems.append(f"MODULE_TESTS[{module!r}] names {target}, which is absent")
    for label, targets in (
        ("CORE_TESTS", CORE_TESTS),
        ("SURFACE_TESTS", SURFACE_TESTS),
        ("CONFTEST_TESTS", CONFTEST_TESTS),
        ("EXTENSION_TESTS", EXTENSION_TESTS),
    ):
        for target in targets:
            if not (ROOT / target).exists():
                problems.append(f"{label} names {target}, which is absent")
    unknown_smoke = sorted(set(SMOKE_IDS) - registered_ids())
    if unknown_smoke:
        problems.append(
            f"SMOKE_IDS names {unknown_smoke}, which are not registered ids - a `-k` on one "
            "selects nothing, and that looks like a pass"
        )
    triggered = SPEC_GENERATORS | DATA_GENERATORS | CAPTURE_GENERATORS
    for entry in generator_sweep():
        stem = Path(entry.command.split()[1]).stem
        if stem not in triggered:
            problems.append(
                f"the generator {stem!r} is in the sweep and in no trigger set, so no diff "
                "asks for it - add it to SPEC_GENERATORS, DATA_GENERATORS or CAPTURE_GENERATORS"
            )
    for crate in members():
        if crate not in crate_namespaces() and crate not in NAMELESS_CRATES:
            problems.append(
                f"the member {crate!r} names no namespace and is in no crate entry, so a change "
                "to it would reach no Python test"
            )
    return problems


def claims(path: str) -> Claim:
    """What one path reaches. **The only place a path is classified.**

    Raises `SystemExit` for a path no rule places, and for a path the tables are not complete
    enough to place - the two are the same failure and both have to stop the run.
    """
    parts = path.split("/")

    if parts[0] == "crates" and len(parts) > 1:
        # **A generated file is its generator's output, not a source edit.** The file carries
        # "DO NOT EDIT BY HAND", its content is a function of the specs, and the `--check` that
        # reproduces it is the thing that proves it current - so this asks for that check rather
        # than for the crate's tests and every dependent's. `provenance_gen.rs` changing would
        # otherwise name all twelve members, because every crate reads it.
        if parts[-1].endswith("_gen.rs"):
            return Claim(generated=True)
        return Claim(crates=(parts[1],))

    if parts[0] == "python":
        return python_claim(path, parts)

    if parts[0] == "specs":
        return specs_claim(path, parts)

    if path.endswith(".py") and (parts[0] == "tools" or parts[:2] == [".claude", "hooks"]):
        return Claim(tools=(Path(path).stem,))

    if parts[0] in ("databank", "data"):
        # The databank is packaged into the wheel, which is what `wheel-data` builds and checks.
        return Claim(data=True, jobs=("wheel-data",))

    if parts[0] == "validation":
        # The NeqSim captures: read by `python/tests/validation`, which its own job runs.
        return Claim(data=True, jobs=("validation",))

    if parts[0] == "docs" or path in ("README.md", "ROADMAP.md"):
        return Claim(docs=True)

    if path in ("Cargo.lock", "rust-toolchain.toml"):
        # **A dependency bump or a toolchain change recompiles every member.** Naming one crate
        # would be a lie, and naming none is how a lockfile edit slips through green.
        return Claim(crates=tuple(members()))

    if path in ("Cargo.toml", "pyproject.toml", "uv.lock"):
        # What a wheel is built from, so `wheel-data` is what proves the change packaged.
        return Claim(config=True, jobs=("wheel-data",))

    if parts[0] == ".github":
        # **`generator_sweep` reads `ci.yml`** to learn which generators to check, so a change
        # to a workflow changes what the sweep runs - which is the drift job's business.
        return Claim(config=True, jobs=("docs-drift",))

    if parts[0] == ".claude" or path == "CLAUDE.md":
        # Nothing in the suite reads a prompt or a policy page. Claimed, placed nowhere, and it
        # says so rather than emitting a `true` no job runs.
        return Claim(config=True)

    for prefix, claim in ROOT_CLAIMS.items():
        if path == prefix or path.startswith(prefix):
            return claim

    raise SystemExit(
        f"test_scope: no rule claims {path}, so this tool cannot say what it reaches. Add a "
        "rule to tools/test_scope.py - or, if it reaches nothing and nothing should reach it, "
        "an inert entry with the reason. A path nothing places is how the scoped set stops "
        "being trustworthy, and an untrusted tool is a sweep waiting to happen."
    )


def nearest_test_dir(path: str) -> str | None:
    """The closest directory holding tests, for a fixture under `python/tests/`.

    A vendored schema is not a test and must not be handed to pytest as one, but the tests
    that read it are the thing a change to it reaches - so walk up to them rather than either
    guessing a name or running nothing.
    """
    parent = (ROOT / path).parent
    root = ROOT / TESTS
    while True:
        if any(parent.glob("test_*.py")):
            return str(parent.relative_to(ROOT))
        if parent == root:
            return None
        parent = parent.parent


#: The header a generated file carries; the Python generators write one of these two.
GENERATED_MARK = re.compile(r"GENERATED FILE|Generated by [`\"']?tools/")


def generated_file(path: str) -> bool:
    """Whether a path is generated - by its name or by its own header, **at any depth**.

    The old check fired only for a module directly under `python/src/azoth/`
    (`len(parts) == 4`), so `core/result_gen.py` - a generated file if there ever was one -
    fell through to the `core` namespace and selected the whole suite. A generated file's
    proof is its generator's `--check`: it is not asked to have run a line under coverage.
    """
    name = Path(path).name
    if name.endswith(".pyi") or name.removesuffix(".py").endswith("_gen"):
        return True
    if not name.endswith(".py"):
        return False
    try:
        with (ROOT / path).open(encoding="utf-8", errors="replace") as handle:
            head = "".join(handle.readline() for _ in range(5))
    except OSError:
        return False
    return bool(GENERATED_MARK.search(head))


GENERATOR_BY_HEADER = re.compile(r"Generated by [`\"']?tools/(gen_\w+)\.py")


def generator_of(path: str) -> str | None:
    """The generator whose output this file is, read from its own header."""
    file = ROOT / path
    if not file.is_file():
        return None
    try:
        with file.open(encoding="utf-8", errors="replace") as handle:
            head = "".join(handle.readline() for _ in range(5))
    except OSError:
        return None
    found = GENERATOR_BY_HEADER.search(head)
    return found.group(1) if found else None


def triggered_generators(paths: list[str]) -> frozenset[str] | None:
    """The generators this diff can make stale, or None for "ask all of them".

    **None is the conservative answer, not the empty one.** A path no rule places takes the
    full sweep: this map is an optimisation, and an optimisation that misses is a stale file
    CI finds rather than a file nothing checked.
    """
    names: set[str] = set()
    for path in paths:
        if generated_file(path):
            generator = generator_of(path)
            if generator is None:
                return None
            names.add(generator)
            continue
        if path.startswith("tools/gen_") and path.endswith(".py"):
            names.add(Path(path).stem)
            continue
        matched = False
        for prefix, generators in (
            ("specs/", SPEC_GENERATORS),
            ("crates/", RUST_GENERATORS),
            ("python/src/azoth/", PYTHON_GENERATORS),
            ("data/", DATA_GENERATORS),
            ("databank/", DATA_GENERATORS),
            ("validation/neqsim/captures/", CAPTURE_GENERATORS),
        ):
            if path.startswith(prefix):
                names |= generators
                matched = True
                break
        if matched:
            continue
        if path.startswith(QUIET_TREES) or path in QUIET_TREES:
            continue
        return None
    return frozenset(names)


def python_claim(path: str, parts: list[str]) -> Claim:
    """A path under `python/`: a test file, a fixture, a module, or a namespace package."""
    if path.startswith(f"{TESTS}/"):
        name = parts[-1]
        if name.startswith("test_") and name.endswith(".py"):
            return Claim(files=(path,))
        if path == f"{TESTS}/conftest.py":
            # **The root conftest is not a directory's conftest.** It is imported by every
            # test, so the nearest-directory rule would name the whole suite; the fixtures it
            # defines are exercised by the extension, backend and selector tests, and the
            # smoke keeps a namespace path honest.
            return Claim(files=CONFTEST_TESTS, ids=SMOKE_IDS)
        # **`_helpers.py` and a vendored fixture are not test files.** Handing one to pytest
        # collects nothing, and the tests it reaches are the ones beside it - which for a
        # nested `conftest` is a whole directory, because that is what a conftest is for.
        directory = nearest_test_dir(path)
        if directory is None:
            return Claim(inert="a fixture under the test tree with no test beside it")
        return Claim(files=(directory,))
    if not path.startswith("python/src/azoth/"):
        if path.startswith("python/azoth/"):
            # **Outside the package and on no import path.** `pythonpath` is `python/src`, so
            # this file is imported by nothing; its own docstring says it belongs under
            # `python/src/azoth`. Recorded rather than refused, because a refusal here is a gate
            # that is red for ever - and a gate nobody can green is scenery. Whoever moves it
            # should give `agents` a NAMESPACE_TESTS entry on the way.
            return Claim(inert="outside the package and on no import path; nothing reads it")
        raise SystemExit(
            f"test_scope: {path} is under python/ but outside {TESTS} and python/src/azoth/"
        )
    if len(parts) < 4:
        raise SystemExit(f"test_scope: {path} names no module under python/src/azoth/")
    if generated_file(path):
        return Claim(generated=True)
    if len(parts) == 4:
        if parts[3] in INERT_MODULES:
            return Claim(inert=INERT_MODULES[parts[3]])
        stem = parts[3].removesuffix(".py").removesuffix(".pyi")
        targets = MODULE_TESTS.get(stem)
        if targets is None:
            raise SystemExit(
                f"test_scope: python/src/azoth/{parts[3]} is a module with no MODULE_TESTS entry,"
                f" so a change to it reaches nothing. Known: {sorted(MODULE_TESTS)}."
            )
        return Claim(files=targets)
    return Claim(namespaces=(parts[3],))


def specs_claim(path: str, parts: list[str]) -> Claim:
    """A path under `specs/`: a model, a calc, a case, or one of the id-less trees."""
    if len(parts) >= 3 and parts[1] in ("models", "calcs"):
        found = declared_id(ROOT / path)
        return Claim(ids=tuple(found), namespaces=(parts[2],))
    if len(parts) >= 3 and parts[1] == "cases":
        found = declared_id(ROOT / path)
        return Claim(
            ids=tuple(found),
            namespaces=tuple(sorted({i.split(".")[0] for i in found})),
            jobs=("validation",),
        )
    for prefix, claim in SPEC_TREES:
        if path.startswith(prefix):
            return claim
    raise SystemExit(f"test_scope: {path} is under specs/ and no rule places it")


def venv_bin() -> str:
    """The venv's `bin`, for the tools that live there and are on nobody's PATH.

    `maturin` is the one: `uv sync --extra dev` installs it into the venv and this box has no
    `maturin` on `PATH`, so a bare `maturin develop` fails with "command not found" - on the
    one step a Rust change cannot skip, which is how it would go unnoticed until a crate
    change ran it.
    """
    return python_runner().rsplit("/", 1)[0]


def python_runner() -> str:
    """The interpreter the suite runs under, as a path relative to the checkout.

    **A worktree has no `.venv` of its own**, and the one it should use is the checkout's, a
    directory or two up - a worktree that grew its own would be testing a different install
    while looking right. So this is searched rather than assumed, and the relative form is
    kept wherever it is correct: CI runs from the workspace root, where `.venv` is, and an
    absolute path baked in here would not exist on a runner.

    `AZOTH_VENV` overrides, for a checkout whose venv is somewhere neither of these.
    """
    override = os.environ.get("AZOTH_VENV")
    if override:
        return f"{override}/bin/python"
    if (ROOT / ".venv" / "bin" / "python").exists():
        return ".venv/bin/python"
    # Walk up, bounded: a worktree sits under `<checkout>/.claude/worktrees/<name>`, so the
    # venv is four levels above it and a fixed number of parents would be a guess.
    candidate = ROOT.parent
    for _ in range(6):
        if (candidate / ".venv" / "bin" / "python").exists():
            return f"{candidate}/.venv/bin/python"
        if candidate.parent == candidate:
            break
        candidate = candidate.parent
    return ".venv/bin/python"


def members() -> list[str]:
    """The workspace members, from the root manifest's own list."""
    manifest = tomllib.loads((ROOT / "Cargo.toml").read_text(encoding="utf-8"))
    return [Path(entry).name for entry in manifest["workspace"]["members"]]


def crate_edges() -> dict[str, set[str]]:
    """Each crate's workspace dependencies, **lib and dev alike**.

    A test-only edge is an edge: `azoth-eos`'s tests use `azoth-test-support`, so a change
    there breaks `azoth-eos`'s tests without breaking anything that ships.
    """
    known = set(members())
    edges: dict[str, set[str]] = {}
    for crate in known:
        manifest = tomllib.loads((CRATES / crate / "Cargo.toml").read_text(encoding="utf-8"))
        depends: set[str] = set()
        for section in ("dependencies", "dev-dependencies"):
            depends |= {name for name in manifest.get(section, {}) if name in known}
        edges[crate] = depends
    return edges


def reaches(crate: str, edges: dict[str, set[str]]) -> set[str]:
    """`crate` and every crate that depends on it, transitively."""
    out = {crate}
    while True:
        grown = {name for name, depends in edges.items() if depends & out and name not in out}
        if not grown:
            return out
        out |= grown


def ids_of(namespace: str) -> list[str]:
    """Every registered id in a namespace, read from the specs rather than the registries."""
    out: list[str] = []
    for kind in ("models", "calcs"):
        for path in sorted((SPECS / kind / namespace).glob("*.toml")):
            text = path.read_text(encoding="utf-8")
            match = re.search(r'^id = "([^"]+)"', text, re.M)
            if match:
                out.append(match.group(1))
    return sorted(out)


def declared_id(path: Path) -> list[str]:
    """The **model or calc** a spec file is about.

    A case file under `specs/cases/` has an `id` of its own - the *case's* name, which is not
    a registered id at all - and names its model on a `model = ` line. So that line is asked
    for first: reading `id` first put case names into the `-k` expression, where they select
    nothing.
    """
    text = path.read_text(encoding="utf-8")
    match = re.search(r'^model = "([^"]+)"', text, re.M)
    if match:
        return [match.group(1)]
    match = re.search(r'^id = "([^"]+)"', text, re.M)
    return [match.group(1)] if match else []


def generator_sweep() -> list[Command]:
    """Every `--check` generator CI runs, **read from the workflow rather than listed again**.

    `.github/workflows/ci.yml`'s `spec-validate` job is where the hand-written decision lives
    about which generators can be checked by running them - most take `--check` and one takes
    a target path and cannot. Reading it means a generator added there is swept here with no
    second edit, and the completeness check below means one that is *missing* everywhere is
    refused rather than skipped in silence.
    """
    workflow = (ROOT / ".github" / "workflows" / "ci.yml").read_text(encoding="utf-8")
    # **CI checks generators two ways and both are read.** `spec-validate` runs
    # `gen_X.py --check`, which reproduces the file and diffs it; the drift job runs
    # `gen_X.py` and then asks git whether anything moved, which checks the same thing from
    # the other end. A generator in either list is checkable, and `--check` is how to ask
    # without writing.
    named = set(re.findall(r"tools/(gen_[\w]+)\.py", workflow))
    # **`gen_ui_fixtures` is not in the sweep**: it is the one generator that runs the CLI
    # (`cargo run -p azoth-cli`) rather than reading the specs, so it cannot run in a job
    # without a cargo build - and tagging it for the Rust job would make every sweep run that
    # job. It is emitted by the claims that reach it instead, and reaches `rust-test` there.
    commands = sorted(
        f"python tools/{stem}.py --check" for stem in named if stem not in CARGO_GENERATORS
    )
    unaccounted = sorted(
        path.stem
        for path in (ROOT / "tools").glob("gen_*.py")
        if path.stem not in named and path.stem not in UNCHECKABLE
    )
    if unaccounted:
        raise SystemExit(
            "test_scope: these generators are checked by nothing - neither CI's `--check` list "
            f"nor UNCHECKABLE: {unaccounted}. Add one or the other; a generator nobody checks "
            "is a stale generated file waiting to be pushed."
        )
    return [
        Command("a generated file embeds the specs or the sources", SCRIPT, command)
        for command in commands
    ]


def is_generator_check(entry: Command) -> bool:
    """A `gen_*.py --check`: the cheap half of the sweep, and the half that fails most often."""
    return "tools/gen_" in entry.command and entry.command.endswith("--check")


def is_cheap(entry: Command) -> bool:
    """A check that costs seconds: the generator sweep, and every unconditional CI step.

    `CLOSE_OUT` is flat-cost by construction - that is why those five jobs run on every push -
    so a formatting or lint failure belongs at the front of the run with the generators.
    """
    return is_generator_check(entry) or entry.command in {step.command for step in CLOSE_OUT}


def cheap_first(commands: list[Command]) -> list[Command]:
    """**The cheap checks before everything else, because they are seconds and the tests are
    an hour.**

    They catch the most common failures - a generated file gone stale, a line to reflow, a type
    the strict check refuses - and running them last means each is found after the expensive
    half has run. Measured three times in one session: a provenance hash carried past a reflow,
    the palette fixtures, and a `cargo fmt` failure that cost a ninety-minute scoped run. The
    order is otherwise unchanged, so the plan still reads in the order it did.
    """
    return [e for e in commands if is_cheap(e)] + [e for e in commands if not is_cheap(e)]


def rust_scope(crates: set[str], edges: dict[str, set[str]]) -> list[Command]:
    """The cargo commands for a set of changed crates, plus the extension rebuild."""
    touched = sorted({crate for crate in crates} | {n for c in crates for n in reaches(c, edges)})
    if not touched:
        return []
    packages = " ".join(f"-p {crate}" for crate in touched)
    out: list[Command] = [
        Command(
            "changed crates and everything that depends on them",
            CARGO,
            f"cargo test {packages} --all-targets -j 10 -- --test-threads 10",
        ),
        Command(
            "their doc tests, which `--all-targets` cannot carry",
            CARGO,
            f"cargo test {packages} --doc -j 10",
        ),
        Command(
            "the same crates' lints",
            CLIPPY,
            f"cargo clippy {packages} --all-targets -j 10 -- -D warnings",
        ),
    ]
    if set(touched) & reaches("azoth-python", edges):
        out.append(
            Command(
                "a crate the extension links, so it must be rebuilt before any Rust-backed "
                "Python test means anything",
                SCRIPT,
                f"{venv_bin()}/maturin develop",
            )
        )
    return out


def id_named_files(ids: set[str]) -> tuple[str, ...]:
    """The test files named after an id, which `-k` over the two sweeps does not reach.

    `python/tests/models/test_pt_flash.py` is `eos.pt_flash`'s own test and lives in
    `python/tests/models/`, which is a slice of no namespace - so without this it ran nightly
    and not on the push that changed it.
    """
    found: set[str] = set()
    for identifier in sorted(ids):
        namespace, _, leaf = identifier.partition(".")
        for candidate in (
            f"{TESTS}/models/test_{leaf}.py",
            f"{TESTS}/{namespace}/test_{leaf}.py",
        ):
            if (ROOT / candidate).exists():
                found.add(candidate)
                break
    return tuple(sorted(found))


def python_commands(namespaces: set[str], ids: set[str], files: set[str]) -> list[Command]:
    """The pytest invocations for those namespaces, ids and files."""
    out: list[Command] = []
    for path in sorted(files):
        out.append(Command(f"the test file itself: {path}", PYTEST, path))

    #: A `-k` term per namespace, plus the ids a spec named outright.
    terms: set[str] = set(ids)
    covered: list[str] = []
    named: set[str] = set()
    for namespace in sorted(namespaces):
        # **The core namespace reaches the smoke, not every directory.** Its own table names
        # the unit/dispatch/solver/record suites; one cross-implementation id per namespace is
        # what keeps a shared-code break from being invisible in a namespace nobody opened.
        if namespace == "core":
            terms |= set(SMOKE_IDS)
        for target in namespace_targets(namespace):
            out.append(Command(f"the {namespace} namespace's own tests", PYTEST, target))
            if (ROOT / target).is_dir():
                covered.append(str(ROOT / target) + "/")
        # **A change inside a namespace reaches that whole namespace**, whether or not a spec
        # moved with it: only a `specs/` path contributed ids before, so `azoth/standards/...`
        # reached nothing - and the file it changed was real code. The term is the namespace and
        # a dot, which pytest accepts (`-k "process."` selects 106 of 431), rather than a hundred
        # ids that nobody can check by reading.
        namespace_ids = set(ids_of(namespace))
        if namespace_ids:
            terms.add(f"{namespace}.")
            named |= set(id_named_files(namespace_ids))
    named = {path for path in named if not any(str(ROOT / path).startswith(d) for d in covered)}

    if terms:
        unknown = sorted(ids - registered_ids())
        if unknown:
            raise SystemExit(
                f"test_scope: {unknown} are not registered ids, so a `-k` naming them selects "
                "nothing - and the one thing that looks like is a pass. `declared_id` reads an "
                "`id =` line out of any spec file, including the ones whose id is not a model."
            )
        expression = " or ".join(sorted(terms))
        out.append(
            Command(
                f"the cases of {sorted(terms)}",
                PYTEST,
                f'{TESTS}/models/test_model_cases.py -k "{expression}"',
            )
        )
        out.append(
            Command(
                f"{sorted(terms)} on both kernels",
                PYTEST,
                f'{TESTS}/test_cross_impl.py -k "{expression}"',
            )
        )
        if named:
            out.append(
                Command(
                    "the files that test these ids by name",
                    PYTEST,
                    " ".join(sorted(named)),
                )
            )
    return out


def scope(paths: list[str], edges: dict[str, set[str]]) -> tuple[list[Command], list[str]]:
    """Every command the paths reach, and the paths nothing claimed."""
    crates: set[str] = set()
    namespaces: set[str] = set()
    ids: set[str] = set()
    files: set[str] = set()
    tools: set[str] = set()
    jobs: set[str] = set()
    commands: list[Command] = []
    data = docs = config = generated = False
    unclaimed: list[str] = []

    for path in paths:
        try:
            claim = claims(path)
        except SystemExit:
            # **Collected, not raised.** One path the map cannot place must not hide the rest of
            # the diff: the caller reports every one of them together and exits 2.
            unclaimed.append(path)
            continue
        crates |= set(claim.crates)
        namespaces |= set(claim.namespaces)
        ids |= set(claim.ids)
        files |= set(claim.files)
        tools |= set(claim.tools)
        jobs |= set(claim.jobs)
        commands.extend(claim.commands)
        data |= claim.data
        docs |= claim.docs
        config |= claim.config
        generated |= claim.generated

    # **A crate change reaches its Python mirror.** Python targets came only from `python/` and
    # `specs/` paths, so a kernel edit ran cargo and not one pytest - and `test_cross_impl.py`,
    # the Rust-against-Python agreement test the library is built around, never ran for the
    # change most likely to break it.
    for crate in sorted(crates):
        crate_side = crate_claim(crate)
        namespaces |= set(crate_side.namespaces)
        files |= set(crate_side.files)
        ids |= set(crate_side.ids)
        jobs |= set(crate_side.jobs)

    commands = rust_scope(crates, edges) + commands
    commands += python_commands(namespaces, ids, files)

    for tool in sorted(tools):
        script = next(
            (
                directory / f"{tool}.py"
                for directory in TOOL_DIRS
                if (directory / f"{tool}.py").exists()
            ),
            None,
        )
        if script is not None:
            # **Each kind of tool is run the way it is run.** A generator has a `--check` that
            # reproduces its output and diffs; a checker *is* the check and is run plainly; and
            # anything else - a hook, a helper - is covered by its own test and by nothing
            # here. Asking a tool for a flag it does not take is an argparse failure dressed as
            # a test, and a rule that cannot be satisfied teaches a reader to ignore the tool.
            if tool.startswith("gen_"):
                commands.append(
                    Command(
                        f"the generator that changed: {tool}",
                        SCRIPT,
                        f"python tools/{tool}.py --check",
                    )
                )
            elif tool in UNRUNNABLE:
                # **Writing the command is not running it.** These need a token or an artefact
                # this box does not have, so a change to one is covered by its own test and by
                # the workflow that runs it - not by an invocation that could only fail.
                pass
            elif tool.startswith(("check_", "audit_")) or tool.endswith("_lint"):
                commands.append(
                    Command(f"the checker that changed: {tool}", SCRIPT, f"python tools/{tool}.py")
                )
        test = ROOT / TESTS / f"test_{tool}.py"
        if test.exists():
            commands.append(
                Command(f"the test that names {tool}", PYTEST, f"{TESTS}/test_{tool}.py")
            )
        if tool == "guard_bash":
            commands.append(
                Command(
                    "the hook's own rule table",
                    SCRIPT,
                    "python .claude/hooks/guard_bash.py --selftest",
                )
            )

    # **The sweep runs for a Python change too.** A registered calc's reference source is in its
    # provenance block, so editing `python/src/azoth/eos/reference/...` makes `gen_registry`'s
    # output stale - which is not a guess: it happened in this very tree, twice. The *subset*
    # is the diff's: a docs or test change asks no generator, a generated file asks its own,
    # and a path the map cannot place falls back to the whole sweep.
    if crates or ids or tools or data or generated or namespaces or files:
        selected = triggered_generators(paths)
        sweep = generator_sweep()
        if selected is not None:
            sweep = [entry for entry in sweep if Path(entry.command.split()[1]).stem in selected]
        commands += sweep
    if data:
        for path in ("test_data.py", "test_data_agreement.py", "test_manifest_tally.py"):
            commands.append(
                Command("the databank is every model's input", PYTEST, f"{TESTS}/{path}")
            )
    if docs:
        commands += [
            Command(
                "a page states counts and paths the tree must match",
                SCRIPT,
                "python tools/check_doc_claims.py",
            ),
            Command("a link can rot", SCRIPT, "python tools/check_links.py"),
        ]
    for job in sorted(jobs):
        # **A job this box cannot run is still placed.** `ui/`, `lean/` and the wheel are jobs
        # rather than commands - the work happens on a runner, so all this has to say is *run
        # that job*. It renders as a comment and `--plan` reads it as a yes.
        commands.append(
            Command(f"this reaches the {job} job, which runs on a runner", SCRIPT, "true", (job,))
        )

    seen = set()
    unique = []
    for entry in commands:
        # Keyed on the job too, so two `true` placeholders for two different jobs both survive.
        if (entry.command, entry.jobs) not in seen:
            seen.add((entry.command, entry.jobs))
            unique.append(entry)
    return unique, unclaimed


def script_jobs(command: str) -> set[str]:
    """Which job runs a non-cargo, non-pytest command, read off the command itself."""
    if command.endswith("/maturin develop"):
        return {"python-test"}
    if "tools/gen_ui_fixtures" in command:
        # The one generator that *runs* the CLI rather than reading the specs, so it belongs in
        # the job that has a cargo build - which is why `ci.yml` keeps it beside `rust-test`.
        return {"rust-test"}
    if "tools/gen_" in command:
        # The `--check` sweeps and the drift job check the same thing from opposite ends.
        return {"docs-drift"}
    if "check_links" in command:
        return {"docs-build"}
    if "tools/provenance" in command:
        return {"provenance"}
    if "tools/test_scope" in command:
        # The tool's own audit, which the `scope` job runs on every push. `--plan` is not here:
        # it is the job's decision, not a check a local run repeats.
        return {"scope"}
    if command in ("true",):
        return set()
    if command.startswith(("python tools/", "python .claude/")):
        # `spec_lint`, `check_json_keys`, `check_ui`, `validate_skills`, `check_doc_claims` and
        # the hook's own selftest are the stdlib-only checkers, and `spec-validate` is where the
        # tree keeps them.
        return {"spec-validate"}
    return set()


def jobs_of(entry: Command) -> tuple[str, ...]:
    """Which CI jobs run a command. **The command is what runs, so the job is read from it.**

    A tag stored beside the command would be a second answer to the same question, and the way
    the two drift apart is a job that is skipped while a command waits for it. A command that
    names no job is a refusal: a command no job runs is a command nobody runs.
    """
    jobs = set(entry.jobs)
    if entry.kind in (CARGO, CLIPPY):
        jobs.add("rust-test")
    elif entry.kind == PYTEST:
        jobs.add("python-test")
    else:
        jobs |= script_jobs(entry.command)
    if not jobs:
        raise SystemExit(
            f"test_scope: no CI job runs {entry.command!r}, so it would never run in CI. Add "
            "one to `script_jobs`, or tag the claim that emitted it."
        )
    return tuple(sorted(jobs))


def known_jobs() -> tuple[str, ...]:
    """Every job `.github/workflows/ci.yml` declares, read from the workflow."""
    workflow = (ROOT / ".github" / "workflows" / "ci.yml").read_text(encoding="utf-8")
    jobs = re.search(r"^jobs:\n(.*?)(?=^\S|\Z)", workflow, re.M | re.S)
    if jobs is None:
        raise SystemExit("test_scope: ci.yml has no `jobs:` block to read")
    return tuple(sorted(re.findall(r"^  ([\w-]+):", jobs.group(1), re.M)))


def needs_extension(commands: list[Command]) -> bool:
    """Whether `python-test` must build the extension for this set.

    A skills-only diff never touches `python/tests` and needs no `maturin develop`; a Rust
    change emits the rebuild command itself. The job's separate build step is guarded on this,
    so a docs or skills push does not pay a Rust link.
    """
    return any(
        entry.command.endswith("/maturin develop")
        or (entry.kind == PYTEST and entry.command.split()[0].startswith(TESTS))
        for entry in commands
    )


def plan(commands: list[Command], full: bool) -> list[tuple[str, bool]]:
    """Which jobs run, as `(job, runs)`.

    **One computation decides both.** `--plan` writes these for `$GITHUB_OUTPUT` and the jobs
    read their own line, so the decision to skip and the decision to run cannot disagree - which
    is the only way a skipped job is a statement rather than a silence.
    """
    tagged = {job for entry in commands for job in jobs_of(entry)}
    known = set(known_jobs())
    unknown = sorted(tagged - known)
    if unknown:
        raise SystemExit(
            f"test_scope: {unknown} are named as jobs and are not in ci.yml, so naming them "
            f"decides nothing. ci.yml has {sorted(known)}."
        )
    return [(job, full or job in ALWAYS_JOBS or job in tagged) for job in sorted(known)]


def with_timeout(line: str, seconds: int | None) -> str:
    """Prefix a command with GNU `timeout`, **inside** the memory cap.

    `subprocess.run(timeout=)` kills only the shell; xdist workers, `cargo` and a server
    survive it. `timeout` signals its own process group, and placing it inside
    `tools/gated.sh` means the cgroup receives the command it is timing.
    """
    if not seconds or line.startswith("#"):
        return line
    cap = f"timeout --signal=TERM --kill-after=30s {seconds}s"
    if line.startswith("tools/gated.sh "):
        return f"tools/gated.sh {cap} {line[len('tools/gated.sh ') :]}"
    return f"{cap} {line}"


@lru_cache(maxsize=1)
def pytest_workers() -> int:
    """How many xdist workers a pytest target gets.

    `AZOTH_PYTEST_WORKERS` overrides; a runner gets four (ubuntu-latest is 4 vCPU), a
    developer box eight, and `0` disables xdist. The suite has no session-scoped fixtures, and
    the three module-scoped ones carry `xdist_group`, so `--dist loadgroup` keeps each of them
    on one worker while the 6,800-parametrized sweeps distribute individually - which is the
    difference between a sweep landing on one worker and an 81-minute suite becoming minutes.
    """
    override = os.environ.get("AZOTH_PYTEST_WORKERS")
    if override is not None:
        try:
            return max(0, int(override))
        except ValueError:
            return 0
    if os.environ.get("CI"):
        return 4
    return min(8, os.cpu_count() or 1)


def pytest_line(target: str, gate: bool, budget_seconds: int | None) -> str:
    """A pytest target as the shell line the gate runs."""
    workers = pytest_workers()
    if workers and target.split()[0].startswith(TESTS):
        target = f"{target} -n {workers} --dist loadgroup --durations=25"
    if os.environ.get("AZOTH_COVERAGE") == "1":
        # **One instrumented run, not two.** `--cov-append` accumulates across the several
        # pytest invocations a scoped set is, and a job that ran the suite twice - once to
        # test and once to measure - would pay the whole cost again for the measurement.
        # The document is written once afterwards by `coverage json`.
        target = f"{target} --cov=azoth --cov-append --cov-report="
    # **`PYTHONPATH` is set, not inherited.** This box has ROS on the path, which puts
    # `/opt/ros/...` in front of the interpreter and makes pytest's plugin autoload import
    # a `launch` module that needs a `yaml` the venv does not have - so a bare `pytest`
    # fails here and passes on a runner, which is the worst kind of difference. Stating it
    # removes the ambient dependency in both places.
    # **`env`, not a bare assignment.** `tools/gated.sh` ends in `exec systemd-run ... "$@"`,
    # which is not a shell - so `PYTHONPATH=x python` reaches it as a *program named
    # `PYTHONPATH=x`* and fails with "No such file or directory". `env` is a real program
    # and passes the assignment through.
    line = f"env PYTHONPATH={ROOT / 'python' / 'src'} {python_runner()} -m pytest {target} -q"
    line = f"tools/gated.sh {line}" if gate else line
    return with_timeout(line, budget_seconds)


def render(
    entry: Command,
    lane: str,
    gate: bool = True,
    budget_seconds: int | None = None,
) -> str:
    """One command as a shell line.

    `gate` wraps the heavy ones in `tools/gated.sh`, which caps memory - **and does not exist
    on a CI runner**, so the nightly calls this with `--no-gate` and runs the bare commands.

    `iterate`, `push` and `tranche` deselect the marked global sweeps; the tranche lane's
    narrowed half is [`sweep_line`]'s. `nightly` runs everything.
    """
    reason, kind, command = entry.reason, entry.kind, entry.command
    if kind == SCRIPT and command == "true":
        # **A claim that reaches no command says so on the line, not only under `--explain`.**
        # `ui/` and `lean/` reach a job this box cannot run, and config reaches nothing at all;
        # both used to print a bare `true`, which is legible to nobody.
        return f"# {reason}"
    if kind == PYTEST:
        target = command
        if lane != "nightly":
            wanted = bare_targets(command)
            if wanted and all(module_is_a_sweep(name) for name in wanted):
                # Filtering a wholly-marked file leaves an empty collection, which pytest calls
                # exit 5 and this tool calls a failure. Skip it whole, and say so.
                return f"# the {lane} lane skips {target}: every test in it is a marked sweep"
            # **Deselected by marker, and unconditionally.** `-m` composes with `-k` by AND, so
            # a narrowed target is narrowed further rather than left carrying the sweeps.
            target = f'{target} -m "not {MARKER}"'
        return pytest_line(target, gate, budget_seconds)
    if kind in (CARGO, CLIPPY) or command.endswith("/maturin develop"):
        line = f"tools/gated.sh {command}" if gate else command
        return with_timeout(line, budget_seconds)
    # **The venv's own tools.** `CLOSE_OUT` stores `ruff` and `mypy` as `ci.yml` writes them, so
    # the table can be read straight against the workflow; here the checkout's venv is substituted,
    # which is what makes them runnable on a box with nothing of theirs on `PATH`.
    return command.replace(".venv/bin/", f"{venv_bin()}/")


@cache
def file_has_marked_tests(path: str) -> bool:
    """True when a *file* carries the marker anywhere - module-level or on a function."""
    file = ROOT / path
    if file.suffix != ".py" or not file.is_file():
        return False
    try:
        tree = ast.parse(file.read_text(encoding="utf-8"))
    except SyntaxError:
        return False
    return MARKER in ast.dump(tree)


def sweep_line(
    entry: Command, selection: str, gate: bool, budget_seconds: int | None
) -> str | None:
    """The tranche lane's marked half: the reached ids inside a sweep, and nothing else.

    `None` when the entry is not pytest, when it names no marked test, or when the scope
    selected no id (a core or binding change - its sweeps are the nightly's).
    """
    if entry.kind != PYTEST or not selection:
        return None
    targets = bare_targets(entry.command)
    files: list[str] = []
    for target in targets:
        path = ROOT / target
        if path.is_dir():
            files.extend(str(p.relative_to(ROOT)) for p in path.rglob("test_*.py"))
        else:
            files.append(target)
    if not any(file_has_marked_tests(path) for path in files):
        return None
    return pytest_line(f'{entry.command} -m {MARKER} -k "{selection}"', gate, budget_seconds)


def selection_terms(commands: list[Command]) -> str:
    """The ids and namespace terms the scope selected, minus the smoke.

    The smoke ids are one per namespace by construction, not a reached set: narrowing a
    registry-wide sweep to them would be the full sweep wearing nine names.
    """
    terms: set[str] = set()
    for entry in commands:
        if entry.kind != PYTEST:
            continue
        found = re.search(r'-k "([^"]+)"', entry.command)
        if found:
            terms |= {term.strip() for term in found.group(1).split(" or ")}
    terms -= set(SMOKE_IDS)
    return " or ".join(sorted(terms))


#: Advisory seconds per target, anchored on the measured numbers the module docstring states.
#: The watchdog and the CI job caps are the enforcement; this only says what is coming.
ESTIMATE_SECONDS: dict[str, int] = {
    f"{TESTS}/test_degenerate_sweep.py": 836,
    f"{TESTS}/test_layer_diff.py": 462,
    f"{TESTS}/test_units_conversion.py": 120,
    f"{TESTS}/test_units_contract.py": 60,
    f"{TESTS}/test_result_json.py": 90,
    f"{TESTS}/test_batch.py": 120,
    f"{TESTS}/test_registry_contract.py": 90,
    f"{TESTS}/models/test_model_contract.py": 120,
    f"{TESTS}/models/test_model_cases.py": 60,
    f"{TESTS}/test_cross_impl.py": 180,
}
ESTIMATE_GENERATORS: dict[str, int] = {
    "gen_python_wrappers": 23,
    "gen_python_register": 24,
    "gen_python_bridge": 59,
    "gen_python_result": 7,
}
ESTIMATE_PYTEST_DEFAULT = 8
ESTIMATE_CARGO = 45


def estimate_seconds(entry: Command) -> int:
    """A rough cost for one command, in seconds."""
    if entry.kind == PYTEST:
        targets = bare_targets(entry.command)
        return max(
            (ESTIMATE_SECONDS.get(target, ESTIMATE_PYTEST_DEFAULT) for target in targets),
            default=ESTIMATE_PYTEST_DEFAULT,
        )
    if is_generator_check(entry):
        stem = Path(entry.command.split()[1]).stem
        return ESTIMATE_GENERATORS.get(stem, 2)
    if entry.kind in (CARGO, CLIPPY) or entry.command.endswith("/maturin develop"):
        return ESTIMATE_CARGO
    return 3


def run_with_watchdog(
    labelled: list[tuple[str, str]], budget_seconds: int | None, lane: str
) -> int:
    """Run the lines, killing the whole process group when the lane's budget is spent."""
    start = time.monotonic()
    for reason, line in labelled:
        remaining = None if not budget_seconds else budget_seconds - (time.monotonic() - start)
        if remaining is not None and remaining <= 0:
            print(f"\ntest_scope: the {lane} budget of {budget_seconds // 60} min is spent")
            return 1
        print(f"\n== {reason}\n== {line}", flush=True)
        process = subprocess.Popen(line, shell=True, cwd=ROOT, start_new_session=True)
        try:
            code = process.wait(timeout=remaining)
        except subprocess.TimeoutExpired:
            _terminate(process)
            print(f"\ntest_scope: the {lane} budget is spent at: {line}")
            return 1
        if code != 0:
            print(f"\ntest_scope: failed: {line}")
            return 1
    return 0


def _terminate(process: subprocess.Popen) -> None:
    """Kill the process group, so xdist workers and servers die with the shell."""
    try:
        os.killpg(os.getpgid(process.pid), signal.SIGTERM)
    except ProcessLookupError:
        return
    try:
        process.wait(timeout=10)
    except subprocess.TimeoutExpired:
        with contextlib.suppress(ProcessLookupError):
            os.killpg(os.getpgid(process.pid), signal.SIGKILL)


def bare_targets(command: str) -> list[str]:
    """The path arguments of a pytest command, past any `-k`/`-m` that follows them."""
    return command.split(" -")[0].split()


@cache
def module_is_a_sweep(path: str) -> bool:
    """True when a *file* carries the marker for its whole module.

    Needed because `-m "not slow_sweep"` over a file like this collects **nothing**, and pytest
    calls that exit 5 - so the iterate lane has to skip the file outright rather than filter it
    to zero and read the empty collection as a failure.
    """
    file = ROOT / path
    if file.suffix != ".py" or not file.is_file():
        return False
    try:
        tree = ast.parse(file.read_text(encoding="utf-8"))
    except SyntaxError:
        return False
    return any(
        isinstance(entry, ast.Assign)
        and any(getattr(name, "id", "") == "pytestmark" for name in entry.targets)
        and MARKER in ast.dump(entry.value)
        for entry in tree.body
    )


def marked_sweeps(commands: list[Command]) -> list[str]:
    """The marked sweeps a command set reaches, by node id.

    **For naming the cost, not for selecting the tests.** Selection is pytest's own
    (`-m "not slow_sweep"` or its absence); this only says *which* minutes a reader is about to
    spend, which is the difference between a gate that is trusted and one that is tolerated. It
    reads the tree with `ast` rather than asking pytest, because it runs on every push.
    """
    found: set[str] = set()
    for entry in commands:
        if entry.kind != PYTEST:
            continue
        for target in bare_targets(entry.command):
            path = ROOT / target
            files = sorted(path.rglob("test_*.py")) if path.is_dir() else [path]
            for file in files:
                if file.suffix != ".py" or not file.exists():
                    continue
                try:
                    tree = ast.parse(file.read_text(encoding="utf-8"))
                except SyntaxError:
                    continue
                relative = file.relative_to(ROOT)
                if module_is_a_sweep(str(relative)):
                    # One line for a whole file, not one per function: a module-level mark means
                    # the file, and listing its forty tests buries the ones that are marked
                    # individually.
                    found.add(str(relative))
                    continue
                for entry in tree.body:
                    if not isinstance(entry, (ast.FunctionDef, ast.AsyncFunctionDef)):
                        continue
                    if any(MARKER in ast.dump(decorator) for decorator in entry.decorator_list):
                        found.add(f"{relative}::{entry.name}")
    return sorted(found)


def changed_paths(args: argparse.Namespace) -> list[str]:
    """The changed paths, from the source the caller named.

    **The default is the *working tree* against the base, not the committed range.** The
    first version diffed the three-dot range, which sees only commits - so a developer
    asking what their in-progress change reaches was told "nothing reaches a test". That is
    the worst answer this tool can give: a silent pass that means "I compared the wrong
    thing". Two dots include the index and the working tree, so the default is what the
    person typing it is actually holding.
    """
    if args.paths:
        return args.paths
    command = ["git", "diff", "--name-only"]
    if args.staged:
        command.append("--cached")
    command.append(args.base)
    result = subprocess.run(command, cwd=ROOT, capture_output=True, text=True, check=True)
    paths = [line for line in result.stdout.split("\n") if line.strip()]
    if not args.staged:
        # **A new file is a change.** `diff` reports tracked files only, so a file that is not
        # committed is invisible to it - and a *new* kernel, model or directory is exactly the
        # change most likely to reach something nothing else reached. That is the hole the
        # "no rule claims this path" error exists to close, so it cannot be left open on the
        # one kind of path that is always new.
        untracked = subprocess.run(
            ["git", "ls-files", "--others", "--exclude-standard"],
            cwd=ROOT,
            capture_output=True,
            text=True,
            check=True,
        )
        paths += [line for line in untracked.stdout.split("\n") if line.strip()]
    return sorted(set(paths))


def tracked_paths() -> list[str]:
    """Every path git tracks, plus the untracked ones a diff would carry."""
    listed = subprocess.run(
        ["git", "ls-files"], cwd=ROOT, capture_output=True, text=True, check=True
    ).stdout
    untracked = subprocess.run(
        ["git", "ls-files", "--others", "--exclude-standard"],
        cwd=ROOT,
        capture_output=True,
        text=True,
        check=True,
    ).stdout
    return sorted({line for line in f"{listed}\n{untracked}".split("\n") if line.strip()})


def pytest_env() -> dict[str, str]:
    """The environment the suite runs under: `PYTHONPATH` stated rather than inherited."""
    return {**os.environ, "PYTHONPATH": str(ROOT / "python" / "src")}


def collects_nothing(command: str) -> str | None:
    """Why a pytest target is empty, or `None` when it collects something.

    `-k` deselecting everything exits 5 rather than 0, so a lone wrong id is already loud. The
    silent case is a wrong id *inside an `or`*, where the expression still selects the valid
    ones - and a target whose file moved. Both are what this catches, without running a test.
    """
    target, _, expression = command.partition(" -k ")
    argv = [
        python_runner(),
        "-m",
        "pytest",
        "--collect-only",
        "-q",
        "-o",
        "addopts=",
        "-p",
        "no:cacheprovider",
    ]
    argv += target.split()
    if expression:
        argv += ["-k", expression.strip().strip('"')]
    result = subprocess.run(argv, cwd=ROOT, capture_output=True, text=True, env=pytest_env())
    if result.returncode == 5:
        return "collects nothing"
    if result.returncode != 0:
        tail = [line for line in result.stdout.splitlines() if line.strip()][-1:]
        return f"collection failed ({' '.join(tail)})"
    if not result.stdout.strip():
        return "collects nothing"
    return None


def verify_selection(commands: list[Command]) -> int:
    """**Every emitted pytest target has to collect at least one test.**

    A target that collects nothing is a green run over nothing, which is the failure this whole
    tool exists to make impossible. Checking it costs a collection and no execution, and it is
    the only thing that can catch an id that is registered but not parametrised.
    """
    bad: list[str] = []
    for entry in commands:
        if entry.kind != PYTEST:
            continue
        problem = collects_nothing(entry.command)
        if problem:
            bad.append(f"{entry.command}: {problem}")
    if bad:
        print(f"test_scope: {len(bad)} pytest target(s) select nothing:")
        for entry in bad:
            print(f"    {entry}")
        return 2
    print("# every pytest target collects something")
    return 0


def audit() -> int:
    """**Every path in the tree, through the table.**

    A diff can only show the holes a diff happens to touch. This walks everything, so the day
    a directory is added and no rule knows it - or a namespace grows a package and no mapping
    does - the gate says so, instead of the map quietly acquiring another path that comes back
    green. It is the guard for every table above, which is why it runs in CI's `scope` job.
    """
    problems = table_problems()
    for path in tracked_paths():
        try:
            claim = claims(path)
        except SystemExit as exc:
            problems.append(f"{path}: {exc}")
            continue
        if claim.inert:
            if not claim.inert.strip():
                problems.append(f"{path}: marked inert with no reason")
        elif claim.empty():
            problems.append(f"{path}: claimed, and reaches nothing")
    if problems:
        print(f"test_scope: {len(problems)} path(s) the map does not place:")
        for problem in problems:
            print(f"    {problem}")
        return 2
    print("# every tracked path claims something")
    return 0


def full_surface() -> list[Command]:
    """**Every command a full run executes, whatever the diff says.**

    `mode: full` used to change only the *plan*: every job was marked `true`, and each job then
    ran `--job <name>` against the diff. On a schedule `github.event.before` is empty, so the
    base is `HEAD~1` - one commit - and the nightly ran four and a half minutes, executed no
    Python tests at all, and failed at its coverage step with `No data to report`. **A full run
    is a property of the commands, not of a flag that decides whether a job starts.**

    The job each command belongs to is read off the command by `jobs_of`, so this list carries no
    second copy of that mapping - the nightly lane and `--full --job <name>` are the same list.
    """
    return [
        Command("the whole Rust workspace", CARGO, "cargo test --workspace --all-targets"),
        Command(
            "its doc tests, which `--all-targets` cannot carry",
            CARGO,
            "cargo test --workspace --doc",
        ),
        Command(
            "the lints over the same workspace",
            CLIPPY,
            "cargo clippy --workspace --all-targets -- -D warnings",
        ),
        Command("the extension both backends need", SCRIPT, f"{venv_bin()}/maturin develop"),
        Command("the whole Python suite", PYTEST, TESTS),
        Command("the skills' smoke checks", PYTEST, "skills"),
        *generator_sweep(),
        *CLOSE_OUT,
    ]


def nightly(gate: bool) -> list[str]:
    """The full surface, rendered: what a developer types to reproduce the nightly by hand."""
    return [render(entry, "nightly", gate) for entry in full_surface()]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--base", default="origin/main", help="what to diff against")
    parser.add_argument("--staged", action="store_true", help="diff the index instead")
    parser.add_argument("--paths", nargs="*", help="explicit paths, rather than a diff")
    parser.add_argument(
        "--lane",
        choices=("iterate", "push", "tranche", "nightly"),
        default=None,
        help="the wall-time budget and sweep policy; --run defaults to iterate, printing to push",
    )
    parser.add_argument(
        "--max-minutes",
        type=int,
        default=None,
        help="override the lane budget; 0 disables the watchdog",
    )
    parser.add_argument(
        "--estimate",
        action="store_true",
        help="print the advisory minutes and whether they fit the lane's budget",
    )
    parser.add_argument("--explain", action="store_true")
    parser.add_argument("--run", action="store_true")
    parser.add_argument(
        "--audit",
        action="store_true",
        help="walk every tracked path through the table, rather than the diff",
    )
    parser.add_argument(
        "--verify-selection",
        action="store_true",
        help="collect every pytest target, and refuse one that selects nothing",
    )
    parser.add_argument(
        "--no-gate",
        action="store_true",
        help="print the bare commands, for a runner that has no tools/gated.sh",
    )
    parser.add_argument(
        "--only",
        choices=("rust", "python"),
        help="one language's half of the set, for a CI job that has only that toolchain",
    )
    parser.add_argument("--job", help="one CI job's commands, as `--plan` names them")
    parser.add_argument(
        "--plan",
        action="store_true",
        help="print `job=<true|false>` for every job in ci.yml, for $GITHUB_OUTPUT",
    )
    parser.add_argument(
        "--full",
        action="store_true",
        help="the whole command set rather than the diff's; with --plan, every job runs",
    )
    args = parser.parse_args()

    if args.audit:
        return audit()

    # **`--run` is the edit lane unless told otherwise.** A local run that silently reaches the
    # whole suite is the hour this tool exists to prevent; the close-out is explicit.
    lane = args.lane or ("iterate" if args.run else "push")
    if args.full:
        # The full surface is the nightly's render: no deselection, whichever lane was typed.
        lane = "nightly"

    if args.lane == "nightly":
        if args.run:
            print(
                "test_scope: --lane nightly is print-only. CI runs it; a local whole-tree run "
                "is `--full --run` with AZOTH_FULL_GATE=1 (and no lane).",
                file=sys.stderr,
            )
            return 2
        for line in nightly(gate=not args.no_gate):
            print(line)
        return 0

    if args.full:
        # **`--full` is the command set, not a permission for a job to start.** `--job` then
        # filters this list, so `--full --job python-test` is the whole Python surface and not
        # the diff's corner of it.
        paths: list[str] = []
        commands, unclaimed = full_surface(), []
    else:
        paths = changed_paths(args)
        commands, unclaimed = scope(paths, crate_edges())
        if lane in ("push", "tranche"):
            # **The push and the close-out also run what the diff does not reach**, because CI
            # does: the five `ALWAYS_JOBS` are flat-cost and unconditional. Deduplicated by
            # command, since a docs change already claims `check_doc_claims` itself.
            seen = {entry.command for entry in commands}
            commands += [entry for entry in CLOSE_OUT if entry.command not in seen]
    commands = cheap_first(commands)
    if args.explain:
        if args.full:
            print("# the full surface, not the diff")
        else:
            print(f"# {len(paths)} path(s) changed against {args.base}")
            for path in sorted(paths):
                print(f"#   {path}")
    if unclaimed:
        # **To stderr, because the caller may be redirecting stdout into a file.** CI writes the
        # command list to one, and the shell creates that file before this runs - so it is an
        # untracked path this tool then refuses, and the reason went into the file it refused.
        # The message has to survive the redirect that caused it.
        print(
            "test_scope: no rule claims these paths, so this tool cannot say what they reach:",
            file=sys.stderr,
        )
        for path in unclaimed:
            print(f"    {path}", file=sys.stderr)
        print(
            "Add a rule to tools/test_scope.py. A path nothing claims is how the scoped set "
            "stops being trustworthy, and an untrusted tool is a sweep waiting to happen.",
            file=sys.stderr,
        )
        return 2

    if args.plan:
        # **The skip decision and the run decision are one computation.** These lines go to
        # `$GITHUB_OUTPUT` and each job reads its own, so a job that does not run is a job this
        # said `false` about - visible in the summary rather than absent without explanation.
        print(f"extension={'true' if needs_extension(commands) else 'false'}")
        for job, runs in plan(commands, full=args.full):
            print(f"{job}={'true' if runs else 'false'}")
        if args.explain:
            for job, runs in plan(commands, full=args.full):
                tags = sorted({entry.reason for entry in commands if job in jobs_of(entry)})
                print(f"# {job}: {'; '.join(tags) if runs and tags else 'nothing reaches it'}")
        return 0

    if args.only == "rust":
        commands = [entry for entry in commands if entry.kind in (CARGO, CLIPPY)]
    elif args.only == "python":
        # Everything that is not a cargo invocation: the pytest targets, the generator
        # checks, the doc checkers and the extension rebuild.
        commands = [entry for entry in commands if entry.kind not in (CARGO, CLIPPY)]
    if args.job:
        commands = [entry for entry in commands if args.job in jobs_of(entry)]
    if args.verify_selection:
        return verify_selection(commands)
    if not commands:
        print(f"# nothing in {len(paths)} changed path(s) reaches a test")
        return 0

    budget = args.max_minutes if args.max_minutes is not None else LANE_BUDGET_MINUTES[lane]
    budget_seconds = budget * 60 if budget else None
    selection = selection_terms(commands) if lane == "tranche" else ""

    sweeps = marked_sweeps(commands)
    if sweeps:
        if args.full:
            print(f"# the full surface runs {len(sweeps)} marked slow sweep(s):")
        else:
            verb = (
                "narrows to the reached ids"
                if lane == "tranche"
                else "has deselected"
                if lane in ("iterate", "push")
                else "reaches"
            )
            print(f"# this diff {verb} {len(sweeps)} marked slow sweep(s):")
        for sweep in sweeps:
            print(f"#   {sweep}")

    labelled: list[tuple[str, str]] = []
    for entry in commands:
        line = render(entry, lane, gate=not args.no_gate, budget_seconds=budget_seconds)
        labelled.append((entry.reason, line))
        marked = sweep_line(entry, selection, not args.no_gate, budget_seconds)
        if marked:
            labelled.append((f"the reached ids inside {entry.command}", marked))

    if args.estimate:
        minutes = sum(estimate_seconds(entry) for entry in commands) / 60
        print(f"# estimated: {minutes:.1f} min against a {budget}-minute {lane} lane")
        if budget and minutes > budget:
            print(f"# WARNING: over the {lane} budget; the tranche lane is the close-out")

    for reason, line in labelled:
        print(f"# {reason}\n{line}" if args.explain else line)

    if args.run:
        return run_with_watchdog(labelled, budget_seconds, lane)
    return 0


if __name__ == "__main__":
    sys.exit(main())
