#!/usr/bin/env python3
"""PreToolUse guard for Bash: block uncapped and repo-wide build/test commands.

Wired from .claude/settings.json. Reads the hook JSON on stdin and prints a
PreToolUse decision on stdout.

The box has 47 GiB RAM, 2 GiB swap and no OOM guard, so a runaway thrashes the
desktop instead of failing one process. Two rules keep that from happening:

  - every heavy command runs through tools/gated.sh, which puts it in a
    memory-capped cgroup scope (and caps the JVM heap); **but the wrapper is
    not a license to run the whole suite** - its inner command is checked too;
  - the whole-tree sweep is nightly-only and needs the literal token
    AZOTH_FULL_GATE=1, which is visible in the transcript.

**A shape matcher has to resolve, not just match.** The first version read the program
token and compared it: `cargo`, `pytest`. That is defeated by anything that puts a word in
front - `uv run pytest`, `bash -c '…'`, `xargs cargo test`, a subshell, a command
substitution - and every one of those runs the same sweep the guard exists to refuse, in a
spelling a developer is likely to type by accident. So a segment is *resolved* past its
launchers, and what a subshell or an `xargs` passes on is checked as its own command.

Fail-open on any internal error, so a broken hook cannot brick every Bash call. The one
place that is **not** fail-open is a command `shlex` cannot read at all: unbalanced quotes
used to mean "nothing to check", which made an unterminated quote the cheapest way past the
guard. It now falls back to a regex over the raw string.

Run `guard_bash.py --selftest` to check the rule table.
"""

from __future__ import annotations

import json
import re
import shlex
import sys
from pathlib import Path

WRAPPER = "tools/gated.sh"
FULL_GATE_TOKEN = "AZOTH_FULL_GATE=1"

# Directory roots that mean "the whole suite" when handed to pytest.
WIDE_PYTEST_ROOTS = {
    "python/tests",
    "python/tests/",
    "python/tests/models",
    "python/tests/models/",
    "skills",
    "skills/",
}

# pytest flags that consume the next token, so it is not a positional path.
PYTEST_VALUE_FLAGS = {
    "-k",
    "-m",
    "-p",
    "-n",
    "-c",
    "-o",
    "-W",
    "-x",
    "--maxfail",
    "--ignore",
    "--ignore-glob",
    "--rootdir",
    "--deselect",
    "--junitxml",
    "--tb",
    "--capture",
    "--only",
    "--job",
    "--paths",
    "--base",
}

HEAVY_CARGO_SUBCOMMANDS = {"build", "test", "clippy", "run", "bench", "install"}
HEAVY_PROGRAMS = {"maturin", "lake", "wasm-pack", "mdbook"}
JVM_PROGRAMS = {"java", "javac"}

#: Programs that run whatever they are told to. Stripped before the real program is read:
#: `uv run pytest python/tests` is pytest, and it was the spelling that walked straight
#: through the first version of this guard.
LAUNCHERS = {"env", "nice", "nohup", "time", "timeout", "stdbuf", "command", "exec", "sudo"}
RUNNERS = {"uv", "uvx", "poetry", "pipenv", "pdm", "hatch", "rye"}
#: Package runners whose next word is the program they run (`npx cargo test`).
PACKAGE_RUNNERS = {"npx"}
#: Flags of a package runner that consume the following word.
PACKAGE_RUNNER_VALUE_FLAGS = {"-p", "--package", "--prefix", "--npm", "--node-options", "--cache"}
#: The longest a single heavy command may run, in seconds. The tranche lane is 20 minutes.
CEILING_SECONDS = 1800
#: A sentinel a forwarder's segment carries when the script it runs names a heavy program.
HEAVY_SCRIPT = "the-script-runs-heavy"
#: Runners whose next word is the subcommand to drop (`uv run`, `poetry run`).
RUNNER_VERBS = {"run", "tool", "exec"}

#: Flags that consume the word after them, so stripping the flag does not leave its value
#: standing where the program should be. `nice -n 5 cargo test` is the case that found this.
LAUNCHER_VALUE_FLAGS = {
    "-n",
    "-w",
    "-u",
    "-g",
    "-C",
    "-S",
    "-s",
    "-k",
    "--signal",
    "--kill-after",
    "--adjustment",
    "--user",
    "--group",
}
RUNNER_VALUE_FLAGS = {"--with", "-w", "--python", "-p", "--index", "--extra"}

#: Programs that pass the rest of their arguments on as a command, or run a recipe nobody here
#: can read. `xargs` and `bash -c` are *resolved*; a build tool is refused outright, because
#: its default target is the sweep and this repository has no Makefile to read.
FORWARDERS = {"bash", "sh", "zsh", "dash", "ksh", "xargs"}
BUILD_TOOLS = {"make", "just", "ninja", "task", "cargo-make"}

#: What the guard looks for when it cannot parse the command at all.
HEAVY_WORDS = re.compile(
    r"\b(cargo\s+(?:test|clippy|build|run|bench|install)|pytest|maturin|lake|wasm-pack|"
    r"mdbook|java|javac)\b"
)


def segments(cmd: str) -> list[list[str]]:
    """Split a shell command on its operators, **respecting quotes**.

    The first version split the raw string on `|`, `;`, `&&` and `||` and then tried to
    tokenise the pieces. That reads a quoted pattern as an operator: `grep -n -E "cargo
    test|pytest" ci.yml` was cut down the middle of its own pattern, and the `pytest`
    fragment was then refused as a repo-wide run - for a command that runs no tests at all.

    `shlex` with `punctuation_chars` emits an operator as its own token **only where it is
    not quoted**, which is the distinction that was missing. Subshell parentheses, a command
    substitution's `$(`, and backticks are operators too, so `(cargo test --workspace)` and
    `$(cargo test --workspace)` are the same command as the unadorned one - both of which
    the first version read as starting with `(cargo` and passed.

    Unbalanced quotes mean the command cannot be read: an empty result, which the caller
    treats as a signal to fall back rather than as nothing to check.
    """
    lexer = shlex.shlex(cmd, posix=True, punctuation_chars="|&;$()`")
    lexer.whitespace_split = True
    lexer.commenters = ""
    operators = "|&;$()`"
    out: list[list[str]] = []
    current: list[str] = []
    try:
        for token in lexer:
            if token and all(character in operators for character in token):
                if current:
                    out.append(current)
                current = []
            else:
                current.append(token)
    except ValueError:
        return []
    if current:
        out.append(current)
    return out


def strip_launchers(tokens: list[str]) -> list[str]:
    """Past the assignments, the launchers and the runners, to the program that runs.

    **A launcher's own arguments are not the program.** `nice -n 5 cargo test` left `-n 5` in
    front of the command, so the program token was `-n` - a word no rule has an opinion about,
    which is the same hole in a different costume.
    """
    while tokens:
        first = tokens[0]
        if re.fullmatch(r"[A-Za-z_][A-Za-z0-9_]*=.*", first):
            tokens = tokens[1:]
            continue
        name = Path(first).name
        if name in LAUNCHERS:
            tokens = tokens[1:]
            while tokens and (tokens[0].startswith("-") or tokens[0] in RUNNER_VERBS):
                tokens = tokens[2:] if tokens[0] in LAUNCHER_VALUE_FLAGS else tokens[1:]
            # **`timeout` takes a DURATION before its command.** Without this the program
            # token would be `300s`, and `timeout 300 pytest python/tests` would be a word
            # no rule has an opinion about - the same hole as `nice -n 5`.
            if name == "timeout" and tokens:
                tokens = tokens[1:]
            continue
        if name in PACKAGE_RUNNERS:
            tokens = tokens[1:]
            while tokens and tokens[0].startswith("-"):
                tokens = tokens[2:] if tokens[0] in PACKAGE_RUNNER_VALUE_FLAGS else tokens[1:]
            continue
        if name in RUNNERS:
            tokens = tokens[1:]
            while tokens and (tokens[0].startswith("-") or tokens[0] in RUNNER_VERBS):
                tokens = tokens[2:] if tokens[0] in RUNNER_VALUE_FLAGS else tokens[1:]
            continue
        break
    return tokens


#: systemd-run options that consume the following word, for unwrapping a scope.
SYSTEMD_VALUE_FLAGS = {
    "-p",
    "--property",
    "-u",
    "--unit",
    "--slice",
    "--description",
    "-E",
    "--setenv",
    "--working-directory",
    "--uid",
    "--gid",
    "--nice",
    "-M",
    "--machine",
}


def wrapper_inner(tokens: list[str]) -> list[str] | None:
    """The command a wrapper segment runs, or None when the segment is not a wrapper.

    **This is what the wide check was skipping.** A wrapper was treated as a segment to
    ignore, so `tools/gated.sh pytest python/tests` - the spelling the guard's own message
    recommends - ran the whole suite unchecked. A wrapper satisfies the *wrap* rule; it does
    not satisfy the *wide* rule.
    """
    name = Path(tokens[0]).name
    if name in ("gated.sh", WRAPPER):
        return strip_launchers(tokens[1:])
    if "systemd-run" in tokens and "--scope" in tokens:
        rest = tokens[tokens.index("systemd-run") + 1 :]
        i = 0
        while i < len(rest):
            if rest[i].startswith("-"):
                i += 2 if rest[i] in SYSTEMD_VALUE_FLAGS else 1
                continue
            break
        return strip_launchers(rest[i:])
    return None


def script_runs_heavy(arguments: list[str]) -> bool:
    """Whether a shell script the command names contains a heavy program.

    A file that cannot be read is not refused here: the shell will fail on it, and guessing
    at an unreadable path would block a command for a reason the caller cannot act on.
    """
    for argument in arguments:
        if argument.startswith("-"):
            continue
        path = Path(argument)
        try:
            if not path.is_file():
                continue
            text = path.read_text(encoding="utf-8", errors="ignore")[:65536]
        except OSError:
            continue
        if HEAVY_WORDS.search(text):
            return True
    return False


def expand(
    tokens: list[str], wrapped_so_far: bool, gate: bool
) -> list[tuple[list[str], bool, bool]]:
    """One resolved command and everything it passes on, with its wrap and gate flags."""
    is_wrapped = wrapped_so_far or wrapped(tokens)
    out: list[tuple[list[str], bool, bool]] = [(tokens, is_wrapped, gate)]
    name = Path(tokens[0]).name
    rest = tokens[1:]
    if name in FORWARDERS:
        if name in ("bash", "sh", "zsh", "dash", "ksh") and "-c" in rest:
            out.extend(candidates(" ".join(rest[rest.index("-c") + 1 :]), is_wrapped, gate))
        elif name in ("bash", "sh", "zsh", "dash", "ksh"):
            # **A script file is a command the guard cannot otherwise read.** If it names a
            # heavy program, the segment is refused; a script with nothing heavy in it passes.
            if script_runs_heavy(rest):
                out.append(([HEAVY_SCRIPT], is_wrapped, gate))
        elif name == "xargs":
            out.extend(candidates(" ".join(rest), is_wrapped, gate))
    elif name in BUILD_TOOLS:
        # **A recipe this guard cannot read is the sweep until proven otherwise.** There is
        # no Makefile in this tree, so an unreadable target is a command nobody can check.
        out.append(([name, "the sweep"], is_wrapped, gate))
    if name == "find" and "-exec" in tokens and HEAVY_WORDS.search(" ".join(tokens)):
        out.append(([HEAVY_SCRIPT], is_wrapped, gate))
    inner = wrapper_inner(tokens)
    if inner:
        out.extend(expand(inner, wrapped(tokens), gate))
    return out


def candidates(
    cmd: str, inherited_wrapped: bool = False, inherited_gate: bool = False
) -> list[tuple[list[str], bool, bool]]:
    """Every token list worth checking: each segment, resolved, and what its shell passes on.

    Each entry is `(tokens, wrapped, gate)`: whether a memory-capped wrapper encloses it, and
    whether the full-gate token was in its own segment. The gate is read **before** the
    leading assignments are stripped, or `AZOTH_FULL_GATE=1 tools/gated.sh ...` would be
    checked without the token it carries.
    """
    out: list[tuple[list[str], bool, bool]] = []
    for raw in segments(cmd):
        gate = inherited_gate or FULL_GATE_TOKEN in raw
        tokens = strip_launchers(raw)
        if not tokens:
            continue
        out.extend(expand(tokens, inherited_wrapped, gate))
    return out


#: A `timeout` duration: seconds by default, with an optional s/m/h suffix.
DURATION = re.compile(r"^(\d+(?:\.\d+)?)([smh]?)$")


def _seconds(text: str) -> float | None:
    match = DURATION.fullmatch(text)
    if match is None:
        return None
    value = float(match.group(1))
    return value * {"": 1.0, "s": 1.0, "m": 60.0, "h": 3600.0}[match.group(2)]


def timeout_refusal(cmd: str) -> str | None:
    """Refuse a `timeout` longer than the per-command ceiling in front of a heavy program."""
    for raw in segments(cmd):
        for index, token in enumerate(raw):
            if Path(token).name != "timeout" or index + 1 >= len(raw):
                continue
            seconds = _seconds(raw[index + 1])
            if seconds is None or seconds <= CEILING_SECONDS:
                continue
            rest = strip_launchers(raw[index + 2 :])
            if rest and is_heavy(rest):
                return (
                    f"No single heavy command runs longer than {CEILING_SECONDS // 60} "
                    f"minutes, and this one asks for {seconds:.0f} seconds. Lanes are bounded "
                    "by their own budgets; the longer surface is CI's nightly."
                )
    return None


def test_scope_refusal(cmd: str) -> str | None:
    """Refuse a local run of the whole tree through the scoped runner."""
    for tokens, _, gate in candidates(cmd):
        if not any(t.endswith("test_scope.py") for t in tokens):
            continue
        if gate or "--run" not in tokens:
            continue
        lane = ""
        if "--lane" in tokens:
            i = tokens.index("--lane")
            if i + 1 < len(tokens):
                lane = tokens[i + 1]
        lane = lane or next((t.split("=", 1)[1] for t in tokens if t.startswith("--lane=")), "")
        if "--full" in tokens or lane == "nightly":
            return (
                "`test_scope.py --run` with `--full` or `--lane nightly` is the whole tree. "
                "It runs nightly in .github/workflows/nightly.yml (budget 45 min). Run the "
                "edit lane (budget 5 min) or `--lane tranche --run` at close-out; for the "
                f"whole tree, prefix `{FULL_GATE_TOKEN}`."
            )
        if "--run" in tokens:
            budget = None
            for index, token in enumerate(tokens):
                if token == "--max-minutes" and index + 1 < len(tokens):
                    try:
                        budget = int(tokens[index + 1])
                    except ValueError:
                        budget = None
                elif token.startswith("--max-minutes="):
                    try:
                        budget = int(token.split("=", 1)[1])
                    except ValueError:
                        budget = None
            if budget is not None:
                lane_budget = {"iterate": 5, "push": 15, "tranche": 20, "nightly": 45}.get(
                    lane or "iterate", 5
                )
                if budget <= 0:
                    return (
                        "`--max-minutes 0` disables the watchdog and every per-command "
                        "timeout. A lane runs to its budget or not at all."
                    )
                if budget > lane_budget:
                    return (
                        f"The {lane or 'iterate'} lane's budget is {lane_budget} minutes and "
                        f"`--max-minutes` may only lower it. A diff too wide for it is "
                        "`--lane tranche` (20 min) or CI, never a longer run."
                    )
    return None


def program(tokens: list[str]) -> str:
    """The program token, past any leading env assignments."""
    for tok in tokens:
        if re.fullmatch(r"[A-Za-z_][A-Za-z0-9_]*=.*", tok):
            continue
        return Path(tok).name
    return ""


def is_pytest(tokens: list[str]) -> bool:
    prog = program(tokens)
    if prog == "pytest":
        return True
    if prog.startswith("python") and "-m" in tokens:
        i = tokens.index("-m")
        if i + 1 < len(tokens) and tokens[i + 1] == "pytest":
            return True
    return False


def pytest_positionals(tokens: list[str]) -> list[str]:
    idx = tokens.index("pytest") if "pytest" in tokens else -1
    args = tokens[idx + 1 :] if idx >= 0 else []
    positional, skip = [], False
    for tok in args:
        if skip:
            skip = False
            continue
        if tok in PYTEST_VALUE_FLAGS:
            skip = True
            continue
        if tok.startswith("-"):
            continue
        positional.append(tok)
    return positional


def cargo_has_package(tokens: list[str]) -> bool:
    return any(t in ("-p", "--package") or t.startswith("--package=") for t in tokens)


def is_wide(tokens: list[str]) -> bool:
    prog = program(tokens)
    if prog == "cargo":
        if "test" in tokens:
            if any(t in ("--workspace", "--all-targets") for t in tokens):
                return True
            return not cargo_has_package(tokens)
        if "clippy" in tokens and "--workspace" in tokens:
            return True
    if is_pytest(tokens):
        pos = pytest_positionals(tokens)
        if not pos:
            return True
        return all(p in WIDE_PYTEST_ROOTS for p in pos)
    if prog == "mypy":
        return not any(not t.startswith("-") for t in tokens[1:])
    return prog in BUILD_TOOLS


def is_heavy(tokens: list[str]) -> bool:
    if HEAVY_SCRIPT in tokens:
        return True
    prog = program(tokens)
    if prog == "cargo":
        return bool(HEAVY_CARGO_SUBCOMMANDS & set(tokens))
    if is_pytest(tokens):
        return True
    if prog in HEAVY_PROGRAMS or prog in JVM_PROGRAMS or prog in BUILD_TOOLS:
        return True
    # **The sweep is heavy because of what it spawns, not because of the word.** Denying the
    # mention was denying `grep -n oracle_sweep tools/`, which runs no JVM at all.
    return prog.startswith("python") and any(t.endswith("oracle_sweep.py") for t in tokens[1:])


def wrapped(tokens: list[str]) -> bool:
    """A segment that *is* the wrapper, or a `systemd-run --scope` scope that caps memory.

    **A scope without a `MemoryMax` is not the wrapper.** `systemd-run --user --scope java`
    runs uncapped, which is the failure the wrapper exists to prevent; only a scope carrying
    the property the wrapper passes counts.
    """
    if program(tokens) in ("gated.sh", WRAPPER):
        return True
    if "systemd-run" not in tokens or "--scope" not in tokens:
        return False
    for index, token in enumerate(tokens):
        if "MemoryMax=" in token:
            return True
        if (
            token in ("-p", "--property")
            and index + 1 < len(tokens)
            and tokens[index + 1].startswith("MemoryMax=")
        ):
            return True
    return False


def decide(cmd: str) -> tuple[bool, str]:
    """Return (allow, reason)."""
    # **An unreadable command is checked by its words.** `shlex` cannot parse it, so there are
    # no segments; treating that as "nothing to check" made an unterminated quote the cheapest
    # bypass in the file. The fail-open property is for a broken *hook*, not for a command the
    # guard cannot read.
    if not segments(cmd):
        if HEAVY_WORDS.search(cmd):
            return False, (
                "That command cannot be parsed (unbalanced quotes), and it names a heavy "
                "program - so it is refused rather than assumed harmless. Close the quote, "
                f"and wrap what it runs: `{WRAPPER} <cmd>`."
            )
        return True, ""

    if re.search(r"\bcargo\s+fmt\b", cmd) and re.search(r"(^|\s)--all(\s|$)", cmd):
        return False, (
            "`cargo fmt --all` reformats other sessions' files in this shared "
            "checkout. Format only your own files: "
            "`rustfmt --edition 2024 <files>`."
        )
    refusal = test_scope_refusal(cmd) or timeout_refusal(cmd)
    if refusal:
        return False, refusal
    for tokens, is_wrapped, full_gate in candidates(cmd):
        # **The token has to be in the same segment as the command it attests to.** Read off
        # the whole line, `echo AZOTH_FULL_GATE=1; cargo test --workspace` attested for a
        # command it had nothing to do with, and so did `echo tools/gated.sh && cargo test`.
        if is_wide(tokens) and not full_gate:
            return False, (
                "That is the whole suite. The lanes are iterate (an edit, budget 5 min) and "
                "tranche (close-out, 20 min); push CI is scoped (15 min) and the whole tree "
                "runs nightly (45 min). Run `python tools/test_scope.py --run` for the set "
                "this change reaches, `--lane tranche --run` at close-out, or for the whole "
                f"tree prefix: {FULL_GATE_TOKEN} {WRAPPER} <cmd>"
            )
        if is_heavy(tokens) and not is_wrapped:
            return False, (
                "Heavy command: it must run inside the memory-capped scope. "
                f"Wrap it: `{WRAPPER} <cmd>` "
                "(add -j 10 to cargo, --test-threads 10 after `--` to tests)."
            )
    return True, ""


def emit_deny(reason: str) -> None:
    print(
        json.dumps(
            {
                "hookSpecificOutput": {
                    "hookEventName": "PreToolUse",
                    "permissionDecision": "deny",
                    "permissionDecisionReason": reason,
                },
                "systemMessage": "Resource guard blocked a heavy command.",
            }
        )
    )


SELFTEST: list[tuple[str, bool]] = [
    # **the wrapper satisfies the wrap rule, not the wide rule.** These were the bypass.
    ("tools/gated.sh cargo test --workspace -j 10", False),
    # **a bare scope is not the wrapper**: without a MemoryMax it runs uncapped
    ("systemd-run --user --scope -q cargo build -j 10", False),
    ("systemd-run --user --scope -q -p MemoryMax=12G cargo build -j 10", True),
    ("systemd-run --user --scope -q pytest python/tests", False),
    # wide -> denied
    ("cargo test --workspace --all-targets", False),
    ("cargo test", False),
    ("cargo clippy --workspace --all-targets -- -D warnings", False),
    ("python -m pytest", False),
    ("pytest skills", False),
    (".venv/bin/python -m pytest python/tests", False),
    ("pytest -k flash", False),
    ("mypy", False),
    # heavy but targeted, unwrapped -> denied (every heavy command is wrapped)
    ("cargo test -p azoth-eos --test flash", False),
    ("cargo clippy -p azoth-eos", False),
    ("python -m pytest python/tests/models/test_model_cases.py -k flash", False),
    ("cargo build -j 10", False),
    ("maturin develop", False),
    ("python tools/oracle_sweep.py", False),
    ("javac Probe.java", False),
    ("lake build", False),
    # the same, wrapped -> allowed
    ("tools/gated.sh cargo test -p azoth-eos --test flash -j 10 -- --test-threads 10", True),
    ("tools/gated.sh python -m pytest python/tests/models/test_model_cases.py -k flash", True),
    ("tools/gated.sh python tools/oracle_sweep.py", True),
    # not heavy -> allowed
    (".venv/bin/ruff check", True),
    # fmt --all -> denied
    ("cargo fmt --all", False),
    # full-gate token (unwrapped) still must be wrapped
    ("AZOTH_FULL_GATE=1 cargo test --workspace", False),
    # unrelated -> allowed
    ("git status", True),
    ("ls -la tools", True),
    ("ruff check python/src", True),
    # the scoped runner is what the deny message tells you to use, so it must be allowed
    ("python tools/test_scope.py --run", True),
    ("python tools/test_scope.py --lane nightly --no-gate", True),
    ("python tools/test_scope.py --plan --explain", True),
    # **a quoted pattern is not an operator.** This was refused as a repo-wide run before the
    # split learned about quotes, and it runs no tests at all.
    ('grep -n -E "cargo test|pytest" .github/workflows/ci.yml', True),
    # and a real pipe still separates
    ("cargo test -p azoth-eos | tail -5", False),
    # **the launcher spellings.** Each runs the same sweep the guard exists to refuse, and
    # each read as a program the rules had never heard of.
    ("uv run pytest python/tests", False),
    ("uv run cargo test --workspace", False),
    ("uv run --with numpy pytest python/tests", False),
    ("poetry run pytest", False),
    ("env RUSTFLAGS=x cargo test --workspace", False),
    ("nice -n 5 cargo test --workspace", False),
    # subshells and command substitution
    ("(cargo test --workspace)", False),
    ("$(cargo test --workspace)", False),
    ("`cargo test --workspace`", False),
    ("echo $(cargo test --workspace)", False),
    # the rest of what a forwarder passes on
    ("bash -c 'cargo test --workspace'", False),
    ("echo x | xargs cargo test -p azoth-eos", False),
    ("make test", False),
    # ...and the wrapped forms of the same are still wide
    ("tools/gated.sh uv run pytest python/tests", False),
    ("tools/gated.sh bash -c 'cargo test --workspace'", False),
    ("tools/gated.sh env PYTHONPATH=/x python -m pytest python/tests -q", False),
    ("tools/gated.sh timeout 300 pytest python/tests", False),
    ("uv run ruff check", True),
    # **reading a file named oracle_sweep is not running it**
    ("grep -n oracle_sweep tools/", True),
    ("sed -n '1,20p' tools/oracle_sweep.py", True),
    # an unreadable command with no heavy word in it is still allowed
    ("echo 'unterminated", True),
    ("cargo test --workspace -- -k 'foo", False),
    # a token that attests for a different segment does not count
    ("echo AZOTH_FULL_GATE=1; cargo test --workspace", False),
    ("echo tools/gated.sh; cargo test --workspace", False),
    # **the gate token is read before its assignment is stripped**, so the documented
    # whole-tree escape still works on the wrapped forms
    ("AZOTH_FULL_GATE=1 tools/gated.sh cargo test --workspace -j 10", True),
    ("AZOTH_FULL_GATE=1 tools/gated.sh pytest python/tests", True),
    ("AZOTH_FULL_GATE=1 systemd-run --user --scope -q cargo test --workspace", False),
    (
        "AZOTH_FULL_GATE=1 systemd-run --user --scope -q -p MemoryMax=12G cargo test --workspace",
        True,
    ),
    # **the spellings that ran on 2026-10-22.** `timeout` must not hide the JVM, and the
    # per-command ceiling refuses an hour even when it is written down.
    ("timeout 900 java -cp .:neqsim-f0c7436.jar Probe", False),
    ("timeout 3600 java -cp .:neqsim-f0c7436.jar Probe", False),
    ("tools/gated.sh timeout 3600 java -cp .:neqsim-f0c7436.jar Probe", False),
    # a bare scope is uncapped, so it is not the wrapper
    ("systemd-run --user --scope java -cp .:neqsim-f0c7436.jar Probe", False),
    (
        "systemd-run --user --scope -p MemoryMax=12G -p MemorySwapMax=0 "
        "java -cp .:neqsim-f0c7436.jar Probe",
        True,
    ),
    # package runners resolve to the program they run
    ("npx cargo test", False),
    ("npx -y pytest python/tests", False),
    # the lane budget may only be lowered, and zero is not a budget
    ("python tools/test_scope.py --lane push --run --max-minutes 35", False),
    ("python tools/test_scope.py --lane push --run --max-minutes 0", False),
    ("python tools/test_scope.py --lane push --run --max-minutes 15", True),
    ("python tools/test_scope.py --lane tranche --run --max-minutes 21", False),
    ("python tools/test_scope.py --lane tranche --run --max-minutes 20", True),
    # `timeout` is resolved past its duration, wrapped or not
    ("timeout 300 pytest python/tests", False),
    (
        "timeout 300 tools/gated.sh python -m pytest "
        "python/tests/models/test_model_cases.py -k flash",
        True,
    ),
    # the scoped runner's whole-tree modes
    ("python tools/test_scope.py --lane nightly --run", False),
    ("python tools/test_scope.py --full --run", False),
    ("python tools/test_scope.py --lane tranche --run", True),
    ("AZOTH_FULL_GATE=1 python tools/test_scope.py --full --run", True),
]


def selftest() -> int:
    bad = 0
    for cmd, expected in SELFTEST:
        allow, _ = decide(cmd)
        if allow != expected:
            bad += 1
            print(f"FAIL: {cmd!r} -> allow={allow}, expected {expected}", file=sys.stderr)
    print(f"{len(SELFTEST) - bad}/{len(SELFTEST)} passed")
    return 1 if bad else 0


def main() -> int:
    if "--selftest" in sys.argv:
        return selftest()
    try:
        payload = json.load(sys.stdin)
    except Exception as exc:
        print(f"guard_bash: unreadable hook input ({exc}); allowing", file=sys.stderr)
        return 0
    if payload.get("tool_name") != "Bash":
        return 0
    cmd = (payload.get("tool_input") or {}).get("command") or ""
    try:
        allow, reason = decide(cmd)
    except Exception as exc:
        print(f"guard_bash: internal error ({exc}); allowing", file=sys.stderr)
        return 0
    if not allow:
        emit_deny(reason)
    return 0


if __name__ == "__main__":
    sys.exit(main())
