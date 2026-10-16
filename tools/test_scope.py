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

`iterate` is the inner loop - the push set with the known-slow global sweeps deselected by
node id. They need deselecting by node id because a single test function loops over every
registered id internally, so `-k` cannot reach inside it.

`push` is the default and is what the pre-push gate runs.

`nightly` is everything in both languages, run by `.github/workflows/nightly.yml` on a
24-hour cycle. The full set is the nightly's job and is not something anyone types.

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
import os
import re
import subprocess
import sys
import tomllib
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

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

#: The global sweeps whose expensive functions cannot be deselected by id, because one
#: function loops over every registered model internally. `iterate` deselects them.
#:
#: Named rather than made selectable: parametrising them is a second, larger piece of work,
#: and until it is done this tuple is the honest statement of what `-k` cannot reach.
SLOW_SWEEPS = (
    f"{TESTS}/models/test_model_contract.py::test_every_declared_range_check_runs",
    f"{TESTS}/test_registry_contract.py::test_every_declared_range_check_runs",
    f"{TESTS}/test_registration_completeness.py",
    f"{TESTS}/test_result_json.py",
    f"{TESTS}/test_dispatch.py",
    f"{TESTS}/test_units_contract.py",
    f"{TESTS}/test_units_conversion.py",
    f"{TESTS}/test_solver_contract.py",
    f"{TESTS}/test_batch.py::test_batch_of_one_equals_scalar_including_warnings",
)

#: A command's kind, which decides how it is wrapped and whether the lane can narrow it.
CARGO, CLIPPY, PYTEST, SCRIPT = "cargo", "clippy", "pytest", "script"

#: Where a tool can live: the generators and checkers, and the hooks.
TOOL_DIRS = (ROOT / "tools", ROOT / ".claude" / "hooks")


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
    commands: tuple[tuple[str, str, str], ...] = ()
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
    # machinery that every other namespace imports. A change there is as wide as it gets, so
    # it says so instead of picking a directory that would look narrower than it is.
    "core": (TESTS,),
}

#: module directly under `python/src/azoth/` (not a namespace) -> the tests that name it.
MODULE_TESTS: dict[str, tuple[str, ...]] = {
    "__init__": (TESTS,),
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
                (
                    "a schema is what the specs are checked against",
                    SCRIPT,
                    "python tools/spec_lint.py",
                ),
                (
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
            commands=(("the palette the UI reads", SCRIPT, "python tools/check_ui.py"),),
            files=(f"{TESTS}/test_check_ui.py",),
            generated=True,
        ),
    ),
    (
        "specs/flowsheets/",
        Claim(files=(f"{TESTS}/test_gen_ui_fixtures.py",), generated=True),
    ),
)

#: Top-level paths, each with what it reaches. The audit refuses any tracked path that is not
#: here and not matched by an earlier rule - which is the list a diff used to be refused by.
ROOT_CLAIMS: dict[str, Claim] = {
    "skills.toml": Claim(
        commands=(
            (
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
            (
                "a normative page states counts and paths the tree must match",
                SCRIPT,
                "python tools/check_doc_claims.py",
            ),
        ),
    ),
    "skills/": Claim(
        commands=(
            (
                "the skills are validated, not just present",
                SCRIPT,
                "python tools/validate_skills.py",
            ),
            ("the skills' own smoke checks", PYTEST, "skills"),
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

    if parts[0] in ("databank", "data", "validation"):
        return Claim(data=True)

    if parts[0] == "docs" or path in ("README.md", "ROADMAP.md"):
        return Claim(docs=True)

    if path in ("Cargo.lock", "rust-toolchain.toml"):
        # **A dependency bump or a toolchain change recompiles every member.** Naming one crate
        # would be a lie, and naming none is how a lockfile edit slips through green.
        return Claim(crates=tuple(members()))

    if parts[0] in (".github", ".claude") or path in (
        "Cargo.toml",
        "pyproject.toml",
        "uv.lock",
        "CLAUDE.md",
    ):
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


def python_claim(path: str, parts: list[str]) -> Claim:
    """A path under `python/`: a test file, a fixture, a module, or a namespace package."""
    if path.startswith(f"{TESTS}/"):
        if path.endswith(".py"):
            return Claim(files=(path,))
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
    if len(parts) == 4:
        if parts[3] in INERT_MODULES:
            return Claim(inert=INERT_MODULES[parts[3]])
        stem = parts[3].removesuffix(".py").removesuffix(".pyi")
        if parts[3].endswith(".pyi") or stem.endswith("_gen"):
            return Claim(generated=True)
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
        return Claim(ids=tuple(found), namespaces=tuple(sorted({i.split(".")[0] for i in found})))
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


def generator_sweep() -> list[tuple[str, str, str]]:
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
    commands = sorted(f"python tools/{stem}.py --check" for stem in named)
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
        ("a generated file embeds the specs or the sources", SCRIPT, command)
        for command in commands
    ]


def rust_scope(crates: set[str], edges: dict[str, set[str]]) -> list[tuple[str, str, str]]:
    """The cargo commands for a set of changed crates, plus the extension rebuild."""
    touched = sorted({crate for crate in crates} | {n for c in crates for n in reaches(c, edges)})
    if not touched:
        return []
    packages = " ".join(f"-p {crate}" for crate in touched)
    out: list[tuple[str, str, str]] = [
        (
            "changed crates and everything that depends on them",
            CARGO,
            f"cargo test {packages} --all-targets -j 10 -- --test-threads 10",
        ),
        (
            "their doc tests, which `--all-targets` cannot carry",
            CARGO,
            f"cargo test {packages} --doc -j 10",
        ),
        (
            "the same crates' lints",
            CLIPPY,
            f"cargo clippy {packages} --all-targets -j 10 -- -D warnings",
        ),
    ]
    if set(touched) & reaches("azoth-python", edges):
        out.append(
            (
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


def python_commands(
    namespaces: set[str], ids: set[str], files: set[str]
) -> list[tuple[str, str, str]]:
    """The pytest invocations for those namespaces, ids and files."""
    out: list[tuple[str, str, str]] = []
    for path in sorted(files):
        out.append((f"the test file itself: {path}", PYTEST, path))
    for namespace in sorted(namespaces):
        for target in namespace_targets(namespace):
            out.append((f"the {namespace} namespace's own tests", PYTEST, target))
        # **A change inside a namespace reaches that namespace's ids**, whether or not a spec
        # moved with it. Only a `specs/` path contributed ids before, so `azoth/standards/...`
        # with no spec edit reached nothing - and the file it changed was real code.
        ids = ids | set(ids_of(namespace))
    if ids:
        unknown = sorted(ids - registered_ids())
        if unknown:
            raise SystemExit(
                f"test_scope: {unknown} are not registered ids, so a `-k` naming them selects "
                "nothing - and the one thing that looks like is a pass. `declared_id` reads an "
                "`id =` line out of any spec file, including the ones whose id is not a model."
            )
        expression = " or ".join(sorted(ids))
        out.append(
            (
                f"the cases of {sorted(ids)}",
                PYTEST,
                f'{TESTS}/models/test_model_cases.py -k "{expression}"',
            )
        )
        out.append(
            (
                f"{sorted(ids)} on both kernels",
                PYTEST,
                f'{TESTS}/test_cross_impl.py -k "{expression}"',
            )
        )
        named = id_named_files(ids)
        if named:
            out.append(
                (
                    "the files that test these ids by name",
                    PYTEST,
                    " ".join(named),
                )
            )
    return out


def scope(
    paths: list[str], edges: dict[str, set[str]]
) -> tuple[list[tuple[str, str, str]], list[str]]:
    """Every command the paths reach, and the paths nothing claimed."""
    crates: set[str] = set()
    namespaces: set[str] = set()
    ids: set[str] = set()
    files: set[str] = set()
    tools: set[str] = set()
    jobs: set[str] = set()
    commands: list[tuple[str, str, str]] = []
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
                    (
                        f"the generator that changed: {tool}",
                        SCRIPT,
                        f"python tools/{tool}.py --check",
                    )
                )
            elif tool.startswith(("check_", "audit_")) or tool.endswith("_lint"):
                commands.append(
                    (f"the checker that changed: {tool}", SCRIPT, f"python tools/{tool}.py")
                )
        test = ROOT / TESTS / f"test_{tool}.py"
        if test.exists():
            commands.append((f"the test that names {tool}", PYTEST, f"{TESTS}/test_{tool}.py"))
        if tool == "guard_bash":
            commands.append(
                (
                    "the hook's own rule table",
                    SCRIPT,
                    "python .claude/hooks/guard_bash.py --selftest",
                )
            )

    if crates or ids or tools or data or generated:
        commands += generator_sweep()
    if data:
        for path in ("test_data.py", "test_data_agreement.py", "test_manifest_tally.py"):
            commands.append(("the databank is every model's input", PYTEST, f"{TESTS}/{path}"))
    if docs:
        commands += [
            (
                "a page states counts and paths the tree must match",
                SCRIPT,
                "python tools/check_doc_claims.py",
            ),
            ("a link can rot", SCRIPT, "python tools/check_links.py"),
        ]
    if config:
        commands.append(("configuration changed; nothing in the suite reads it", SCRIPT, "true"))

    if jobs:
        # **A job this box cannot run is still placed.** `ui/` and `lean/` reach a job rather
        # than a command, and saying so is what keeps their diff from being a shrug - the same
        # honesty the `INERT` entries carry, in the other direction.
        commands.append(
            (
                f"these reach the {', '.join(sorted(jobs))} job(s), which run on a runner",
                SCRIPT,
                "true",
            )
        )

    seen = set()
    unique = []
    for entry in commands:
        if entry[2] not in seen:
            seen.add(entry[2])
            unique.append(entry)
    return unique, unclaimed


def render(entry: tuple[str, str, str], lane: str, gate: bool = True) -> str:
    """One command as a shell line.

    `gate` wraps the heavy ones in `tools/gated.sh`, which caps memory - **and does not exist
    on a CI runner**, so the nightly calls this with `--no-gate` and runs the bare commands.
    """
    reason, kind, command = entry
    if kind == SCRIPT and command == "true":
        # **A claim that reaches no command says so on the line, not only under `--explain`.**
        # `ui/` and `lean/` reach a job this box cannot run, and config reaches nothing at all;
        # both used to print a bare `true`, which is legible to nobody.
        return f"# {reason}"
    if kind == PYTEST:
        target = command
        # A target already narrowed by `-k` is the id itself and needs no further narrowing;
        # anything wider - a file, a namespace directory - carries the slow sweeps with it.
        if lane == "iterate" and " -k " not in target:
            deselects = " ".join(f"--deselect {node}" for node in SLOW_SWEEPS)
            target = f"{target} {deselects}"
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
        return (
            f"tools/gated.sh {line}"
            if gate
            else f"PYTHONPATH={ROOT / 'python' / 'src'} {python_runner()} -m pytest {target} -q"
        )
    if kind in (CARGO, CLIPPY) or command.endswith("/maturin develop"):
        return f"tools/gated.sh {command}" if gate else command
    return command


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


def nightly(gate: bool) -> list[str]:
    """The whole set: every test in both languages, and every generator reproduced.

    **`.github/workflows/nightly.yml` runs this**, which is what makes the two agree - the
    scheduled job is the only thing that runs what no diff claimed, so a second copy of the
    list in YAML would be a second thing to keep in step.
    """
    entries: list[tuple[str, str, str]] = [
        ("the whole Rust workspace", CARGO, "cargo test --workspace --all-targets"),
        (
            "its doc tests, which `--all-targets` cannot carry",
            CARGO,
            "cargo test --workspace --doc",
        ),
        (
            "the lints over the same workspace",
            CLIPPY,
            "cargo clippy --workspace --all-targets -- -D warnings",
        ),
        ("the extension both backends need", SCRIPT, f"{venv_bin()}/maturin develop"),
        ("the whole Python suite", PYTEST, TESTS),
        ("the skills' smoke checks", PYTEST, "skills"),
    ]
    return [render(entry, "nightly", gate) for entry in entries + generator_sweep()]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--base", default="origin/main", help="what to diff against")
    parser.add_argument("--staged", action="store_true", help="diff the index instead")
    parser.add_argument("--paths", nargs="*", help="explicit paths, rather than a diff")
    parser.add_argument("--lane", choices=("iterate", "push", "nightly"), default="push")
    parser.add_argument("--explain", action="store_true")
    parser.add_argument("--run", action="store_true")
    parser.add_argument(
        "--audit",
        action="store_true",
        help="walk every tracked path through the table, rather than the diff",
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
    args = parser.parse_args()

    if args.audit:
        return audit()

    if args.lane == "nightly":
        for line in nightly(gate=not args.no_gate):
            print(line)
        return 0

    paths = changed_paths(args)
    commands, unclaimed = scope(paths, crate_edges())
    if args.explain:
        print(f"# {len(paths)} path(s) changed against {args.base}")
        for path in sorted(paths):
            print(f"#   {path}")
    if args.only == "rust":
        commands = [entry for entry in commands if entry[1] in (CARGO, CLIPPY)]
    elif args.only == "python":
        # Everything that is not a cargo invocation: the pytest targets, the generator
        # checks, the doc checkers and the extension rebuild.
        commands = [entry for entry in commands if entry[1] not in (CARGO, CLIPPY)]
    if unclaimed:
        print("test_scope: no rule claims these paths, so this tool cannot say what they reach:")
        for path in unclaimed:
            print(f"    {path}")
        print(
            "Add a rule to tools/test_scope.py. A path nothing claims is how the scoped set "
            "stops being trustworthy, and an untrusted tool is a sweep waiting to happen."
        )
        return 2
    if not commands:
        print(f"# nothing in {len(paths)} changed path(s) reaches a test")
        return 0

    for reason, kind, command in commands:
        line = render((reason, kind, command), args.lane, gate=not args.no_gate)
        print(f"# {reason}\n{line}" if args.explain else line)

    if args.run:
        for reason, kind, command in commands:
            line = render((reason, kind, command), args.lane, gate=not args.no_gate)
            print(f"\n== {reason}\n== {line}", flush=True)
            if subprocess.run(line, shell=True, cwd=ROOT).returncode != 0:
                print(f"\ntest_scope: failed: {line}")
                return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
