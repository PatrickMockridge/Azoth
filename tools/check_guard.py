#!/usr/bin/env python3
"""The guard in the working tree is the tracked one, and the settings point at it.

**This check exists because a stale copy silently disabled the guard on 2026-10-22.** The
hook that ran was an untracked file from the checkout's earliest day, 8 KB against the
tracked 22 KB, and it let `timeout 900 java ...` through where the tracked hook refuses it.
The guard is the only thing enforcing the memory cap and the lane budgets, so a checkout
whose copy differs from `HEAD` must fail a gate rather than fail open.

    python tools/check_guard.py            # write nothing; exit 2 on a stale guard
"""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

HOOK = Path(".claude/hooks/guard_bash.py")
SETTINGS = Path(".claude/settings.json")


def _git(root: Path, *arguments: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        ["git", *arguments], cwd=root, capture_output=True, text=True, check=False
    )


def main() -> int:
    top = _git(Path.cwd(), "rev-parse", "--show-toplevel")
    if top.returncode != 0:
        print("check_guard: not inside a git checkout", file=sys.stderr)
        return 2
    root = Path(top.stdout.strip())

    tracked = _git(root, "ls-files", "--error-unmatch", str(HOOK))
    if tracked.returncode != 0:
        print(
            f"check_guard: {HOOK} is untracked, so the guard that runs here is not the one "
            "the tree carries. Remove it and restore the tracked file: "
            f"`rm {HOOK} && git checkout -- {HOOK}`.",
            file=sys.stderr,
        )
        return 2

    head = _git(root, "show", f"HEAD:{HOOK}")
    if head.returncode != 0:
        print(f"check_guard: HEAD carries no {HOOK}", file=sys.stderr)
        return 2
    working = (root / HOOK).read_text(encoding="utf-8")
    if head.stdout != working:
        print(
            f"check_guard: {HOOK} differs from HEAD. The guard is the enforcement; an "
            "edited copy is not the rule the tree states. Restore it: "
            f"`git checkout -- {HOOK}`.",
            file=sys.stderr,
        )
        return 2

    settings = json.loads((root / SETTINGS).read_text(encoding="utf-8"))
    commands = [
        hook.get("command", "")
        for entry in settings.get("hooks", {}).get("PreToolUse", [])
        for hook in entry.get("hooks", [])
    ]
    if not any("guard_bash.py" in command for command in commands):
        print(
            f"check_guard: {SETTINGS} does not wire {HOOK} into a PreToolUse hook",
            file=sys.stderr,
        )
        return 2

    print(f"check_guard: OK ({HOOK} matches HEAD and is wired)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
