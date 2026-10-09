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

Fail-open: on any internal error the command is allowed, so a broken hook
cannot brick every Bash call.

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
}

HEAVY_CARGO_SUBCOMMANDS = {"build", "test", "clippy", "run", "bench", "install"}
HEAVY_PROGRAMS = {"maturin", "lake", "wasm-pack", "mdbook"}
JVM_PROGRAMS = {"java", "javac"}


def segments(cmd: str) -> list[list[str]]:
    """Split a shell command on operators and tokenise each part."""
    out = []
    for part in re.split(r"&&|\|\||\||;|\n", cmd):
        part = part.strip()
        if not part:
            continue
        try:
            out.append(shlex.split(part))
        except ValueError:
            continue
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
    return False


def is_heavy(tokens: list[str]) -> bool:
    prog = program(tokens)
    if prog == "cargo":
        return bool(HEAVY_CARGO_SUBCOMMANDS & set(tokens))
    if is_pytest(tokens):
        return True
    if prog in HEAVY_PROGRAMS or prog in JVM_PROGRAMS:
        return True
    return any(t.endswith("oracle_sweep.py") for t in tokens)


def decide(cmd: str) -> tuple[bool, str]:
    """Return (allow, reason)."""
    if WRAPPER in cmd or ("systemd-run" in cmd and "--scope" in cmd):
        return True, ""
    if re.search(r"\bcargo\s+fmt\b", cmd) and re.search(r"(^|\s)--all(\s|$)", cmd):
        return False, (
            "`cargo fmt --all` reformats other sessions' files in this shared "
            "checkout. Format only your own files: "
            "`rustfmt --edition 2024 <files>`."
        )
    full_gate = FULL_GATE_TOKEN in cmd
    for tokens in segments(cmd):
        if not tokens:
            continue
        if is_wide(tokens) and not full_gate:
            return False, (
                "That is the full gate set (197 Rust test files / ~111 Python "
                "test files). It runs once per tranche, at close-out, not per "
                "edit. Test what the change touches instead, or if this really "
                f"is close-out, prefix it: {FULL_GATE_TOKEN} {WRAPPER} <cmd>"
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
