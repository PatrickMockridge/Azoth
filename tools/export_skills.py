#!/usr/bin/env python3
"""Export azoth's skill catalog into the native skill format of agent tools.

Reads `skills.toml` - the single source of truth - and each `SKILL.md` body, and
writes one skill per azoth skill, as YAML frontmatter (`name`, `description`) over
the body, into each target's discovery directory. The frontmatter and body are
identical for every target; only the directory and file layout differ.

Targets:

    dsh      <root>/.dsh/skills/<name>.md           DeepSeek Harness (flat)
    claude   <root>/.claude/skills/<name>/SKILL.md  Claude Code
    opencode <root>/.opencode/skills/<name>/SKILL.md  OpenCode

Usage:
    python tools/export_skills.py --target all
    python tools/export_skills.py --target claude
"""

from __future__ import annotations

import argparse
import sys
import tomllib
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
CATALOG = ROOT / "skills.toml"

#: Target -> (output subdirectory under the root, flat-file layout). A flat target
#: writes `<name>.md`; a directory target writes `<name>/SKILL.md`.
TARGETS: dict[str, tuple[str, bool]] = {
    "dsh": (".dsh/skills", True),
    "claude": (".claude/skills", False),
    "opencode": (".opencode/skills", False),
}


def _quote(value: str) -> str:
    """A YAML double-quoted scalar, escaping what YAML escapes."""
    return '"' + value.replace("\\", "\\\\").replace('"', '\\"') + '"'


def render(name: str, description: str, body: str) -> str:
    """One skill: `name` and `description` frontmatter over the azoth body."""
    return f"---\nname: {name}\ndescription: {_quote(description)}\n---\n{body}"


def export_one(out: Path, flat: bool, name: str, description: str, body: str) -> None:
    """Write one skill, as a flat file or a `SKILL.md` under its own directory."""
    target = out / f"{name}.md" if flat else out / name / "SKILL.md"
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(render(name, description, body), encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--target",
        choices=["all", *TARGETS],
        default="all",
        help="which tool to write for (default: all)",
    )
    parser.add_argument(
        "--root", type=Path, default=ROOT, help="output root (default: the repository root)"
    )
    args = parser.parse_args()

    if not CATALOG.is_file():
        sys.exit(f"export_skills: {CATALOG} does not exist")

    document = tomllib.loads(CATALOG.read_text(encoding="utf-8"))
    entries = document.get("skill", [])

    targets = list(TARGETS) if args.target == "all" else [args.target]

    for target in targets:
        subdir, flat = TARGETS[target]
        out = args.root / subdir
        written = 0
        for entry in entries:
            name = str(entry["name"])
            body = ROOT / str(entry["path"])
            if not body.is_file():
                print(f"  skip  {name}: {body} does not exist", file=sys.stderr)
                continue
            export_one(out, flat, name, str(entry["description"]), body.read_text(encoding="utf-8"))
            written += 1
        print(f"export_skills: wrote {written} skill(s) to {out}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
