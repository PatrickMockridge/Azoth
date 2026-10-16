#!/usr/bin/env python3
"""What a change reaches, as commands.

    python tools/test_scope.py                     # against origin/main...HEAD
    python tools/test_scope.py --staged            # against the index
    python tools/test_scope.py --explain           # and say why each command is there
    python tools/test_scope.py --run               # run them, in order
    python tools/test_scope.py --lane nightly      # the whole thing, either language

**Why this exists.** The rule "test what the change touches" has been written down three
times - 2026-09-21, 09-22, 10-13 - and enforced once by a PreToolUse hook, and it still
degraded, because nothing in the tree answered the question it asks. There was no command
that said *this diff reaches these tests*; there was only a hook saying no to the obvious
wrong answer. So every session re-derived the boundary by hand, got it slightly wrong once,
and fell back to running everything - an hour, killed by the system about half the time.

The rules are the two functions below, [`rust_scope`] and [`python_scope`]. **A path that
matches no rule is an error**, not a shrug: a new directory that silently escapes the map is
how "run everything" becomes the safe answer again.

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
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
CRATES = ROOT / "crates"
SPECS = ROOT / "specs"

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


#: Generators that are not run as a bare `--check`, and why. Every other `gen_*.py` has to
#: appear in CI's list or [`generator_sweep`] refuses - an unchecked generator that nothing
#: names by hand is how a stale generated file ships, which is not hypothetical: two did.
UNCHECKABLE = {
    "gen_user_data": "it takes a target path, so it is not a `--check` generator",
}


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
                "maturin develop",
            )
        )
    return out


def python_scope(paths: list[str]) -> tuple[set[str], set[str], set[str]]:
    """The namespaces, ids and test files a set of Python-side paths reaches."""
    namespaces: set[str] = set()
    ids: set[str] = set()
    files: set[str] = set()
    for path in paths:
        parts = path.split("/")
        if path.startswith("python/src/azoth/") and len(parts) > 3:
            namespaces.add(parts[3].removesuffix(".py"))
        elif path.startswith(f"{TESTS}/"):
            files.add(path)
        elif path.startswith("specs/"):
            for found in declared_id(ROOT / path):
                ids.add(found)
                namespaces.add(found.split(".")[0])
    return namespaces, ids, files


def python_commands(
    namespaces: set[str], ids: set[str], files: set[str]
) -> list[tuple[str, str, str]]:
    """The pytest invocations for those namespaces, ids and files."""
    out: list[tuple[str, str, str]] = []
    for path in sorted(files):
        out.append((f"the test file itself: {path}", PYTEST, path))
    for namespace in sorted(namespaces):
        directory = ROOT / TESTS / namespace
        if directory.is_dir():
            out.append((f"the {namespace} namespace's own tests", PYTEST, f"{TESTS}/{namespace}"))
    if ids:
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
    return out


def scope(
    paths: list[str], edges: dict[str, set[str]]
) -> tuple[list[tuple[str, str, str]], list[str]]:
    """Every command the paths reach, and the paths nothing claimed."""
    crates: set[str] = set()
    tools: set[str] = set()
    data = docs = config = generated = False
    python_paths: list[str] = []
    unclaimed: list[str] = []

    for path in paths:
        parts = path.split("/")
        if parts[0] == "crates" and len(parts) > 1 and parts[-1].endswith("_gen.rs"):
            # **A generated file is its generator's output, not a source edit.** The file
            # carries "DO NOT EDIT BY HAND", its content is a function of the specs, and the
            # `--check` that reproduces it is the thing that proves it current - so this asks
            # for that check rather than for the crate's tests and every dependent's.
            # `crates/azoth-core/src/provenance_gen.rs` changing would otherwise name all
            # twelve members, because every crate reads it; it changed here because the specs
            # did, and the specs are what the sweep reproduces from.
            generated = True
        elif parts[0] == "crates" and len(parts) > 1:
            crates.add(parts[1])
        elif path.startswith("python/") or path.startswith("specs/"):
            python_paths.append(path)
        elif path.endswith(".py") and (parts[0] == "tools" or parts[:2] == [".claude", "hooks"]):
            tools.add(Path(path).stem)
        elif parts[0] == "databank" or path.startswith("data/") or path.startswith("validation/"):
            data = True
        elif parts[0] == "docs" or path in ("README.md", "ROADMAP.md"):
            docs = True
        elif parts[0] in (".github", ".claude") or path in (
            "Cargo.toml",
            "pyproject.toml",
            "uv.lock",
            "CLAUDE.md",
        ):
            config = True
        else:
            unclaimed.append(path)

    commands = rust_scope(crates, edges)
    namespaces, ids, files = python_scope(python_paths)
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
    _, kind, command = entry
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
    if kind in (CARGO, CLIPPY) or command.startswith((".venv/bin/maturin", "maturin")):
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
        ("the extension both backends need", SCRIPT, "maturin develop"),
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
