#!/usr/bin/env python3
"""Export azoth's skill catalog into DeepSeek Harness's skill format.

Reads `skills.toml` - the single source of truth - and each `SKILL.md` body, and
writes one DeepSeek Harness skill per azoth skill: YAML frontmatter (`name`,
`description`) prepended to the azoth body, into `.dsh/skills/`, which dsh's local
provider scans at project rank.

# Why this exists

azoth's skills carry no frontmatter: their metadata lives only in `skills.toml`.
DeepSeek Harness requires `name` and `description` frontmatter on every skill. This
script is the adapter between the two, and the only place in the tree that knows
dsh's format - so a change to dsh's dev-preview format is a change here and nowhere
else.

Usage:
    python tools/export_dsh_skills.py [--out .dsh/skills]
"""

from __future__ import annotations

import argparse
import sys
import tomllib
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
CATALOG = ROOT / "skills.toml"


def _quote(value: str) -> str:
    """A YAML double-quoted scalar, escaping what YAML escapes."""
    return '"' + value.replace("\\", "\\\\").replace('"', '\\"') + '"'


def render(name: str, description: str, body: str) -> str:
    """One dsh skill: `name` and `description` frontmatter over the azoth body."""
    return f"---\nname: {name}\ndescription: {_quote(description)}\n---\n{body}"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--out", default=str(ROOT / ".dsh" / "skills"), help="output directory")
    args = parser.parse_args()

    if not CATALOG.is_file():
        sys.exit(f"export_dsh_skills: {CATALOG} does not exist")

    document = tomllib.loads(CATALOG.read_text(encoding="utf-8"))
    entries = document.get("skill", [])

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    written = 0
    for entry in entries:
        name = str(entry["name"])
        body = ROOT / str(entry["path"])
        if not body.is_file():
            print(f"  skip  {name}: {body} does not exist", file=sys.stderr)
            continue
        (out / f"{name}.md").write_text(
            render(name, str(entry["description"]), body.read_text(encoding="utf-8")),
            encoding="utf-8",
        )
        written += 1

    print(f"export_dsh_skills: wrote {written} skill(s) to {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
