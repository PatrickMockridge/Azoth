#!/usr/bin/env python3
"""PreToolUse guard for Bash: block uncapped and repo-wide build/test commands.

Wired from .claude/settings.json. Reads the hook JSON on stdin and prints a
PreToolUse decision on stdout.

The box has 47 GiB RAM, 2 GiB swap and no OOM guard, so a runaway thrashes the
desktop instead of failing one process. Two rules keep that from happening:

  - every heavy command runs through tools/gated.sh, which puts it in a
    memory-capped cgroup scope (and caps the JVM heap);
  - the repo-wide sweep is close-out only, and needs the literal token
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
LAUNCHERS = {"env", "nice", "nohup", "time", "stdbuf", "command", "exec", "sudo"}
RUNNERS = {"uv", "uvx", "poetry", "pipenv", "pdm", "hatch", "rye"}
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
            continue
        if name in RUNNERS:
            tokens = tokens[1:]
            while tokens and (tokens[0].startswith("-") or tokens[0] in RUNNER_VERBS):
                tokens = tokens[2:] if tokens[0] in RUNNER_VALUE_FLAGS else tokens[1:]
            continue
        break
    return tokens


def candidates(cmd: str) -> list[list[str]]:
    """Every token list worth checking: each segment, resolved, and what its shell passes on."""
    out: list[list[str]] = []
    for raw in segments(cmd):
        tokens = strip_launchers(raw)
        if not tokens:
            continue
        out.append(tokens)
        name = Path(tokens[0]).name
        rest = tokens[1:]
        if name in FORWARDERS:
            if name in ("bash", "sh", "zsh", "dash", "ksh") and "-c" in rest:
                script = rest[rest.index("-c") + 1 :]
                out.extend(candidates(" ".join(script)))
            elif name == "xargs":
                out.extend(candidates(" ".join(rest)))
        elif name in BUILD_TOOLS:
            # **A recipe this guard cannot read is the sweep until proven otherwise.** There is
            # no Makefile in this tree, so an unreadable target is a command nobody can check.
            out.append([name, "the sweep"])
    return out


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
    if prog in BUILD_TOOLS:
        return True
    return False


def is_heavy(tokens: list[str]) -> bool:
    prog = program(tokens)
    if prog == "cargo":
        return bool(HEAVY_CARGO_SUBCOMMANDS & set(tokens))
    if is_pytest(tokens):
        return True
    if prog in HEAVY_PROGRAMS or prog in JVM_PROGRAMS or prog in BUILD_TOOLS:
        return True
    # **The sweep is heavy because of what it spawns, not because of the word.** Denying the
    # mention was denying `grep -n oracle_sweep tools/`, which runs no JVM at all.
    if prog.startswith("python") and any(t.endswith("oracle_sweep.py") for t in tokens[1:]):
        return True
    return False


def wrapped(tokens: list[str]) -> bool:
    """A segment that *is* the wrapper, or a `systemd-run --scope` scope of its own."""
    if program(tokens) in ("gated.sh", WRAPPER):
        return True
    return "systemd-run" in tokens and "--scope" in tokens


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
    for tokens in candidates(cmd):
        if wrapped(tokens):
            continue
        # **The token has to be in the same segment as the command it attests to.** Read off
        # the whole line, `echo AZOTH_FULL_GATE=1; cargo test --workspace` attested for a
        # command it had nothing to do with, and so did `echo tools/gated.sh && cargo test`.
        full_gate = FULL_GATE_TOKEN in tokens
        if is_wide(tokens) and not full_gate:
            return False, (
                "That is the whole suite (1,549 Rust tests / 1,165 Python tests). It runs "
                "nightly, in .github/workflows/nightly.yml, and is not what a push needs: run "
                "`python tools/test_scope.py --run` for the set this change reaches, or if you "
                f"do want everything, prefix it: {FULL_GATE_TOKEN} {WRAPPER} <cmd>"
            )
        if is_heavy(tokens):
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
    # wrapped -> allowed
    ("tools/gated.sh cargo test --workspace -j 10", True),
    ("systemd-run --user --scope -q cargo build -j 10", True),
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
    # ...and the wrapped forms of the same
    ("tools/gated.sh uv run pytest python/tests", True),
    ("tools/gated.sh bash -c 'cargo test --workspace'", True),
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
