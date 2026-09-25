#!/usr/bin/env python3
"""Check `skills.toml` against the skill folders it catalogues.

# Why this exists

A skill is a `SKILL.md` plus examples and tests, and `skills.toml` is the one place
its name, trigger and basis are written down. That record is worth having only if it
cannot rot, and rot here is silent in both directions: a skill added without a
catalog entry is undiscoverable, and an entry whose `path` points at a `SKILL.md`
that does not exist is a trigger that loads nothing. Neither breaks the build; both
leave an agent unable to find the skill.

# What it checks

In one line each: the catalog parses; every entry carries the required fields, and a
name, version, description (with `USE WHEN:`), basis and path that satisfy the
standard; every `SKILL.md` under `skills/` is listed exactly once; and every listed
`SKILL.md` carries the required sections.

Usage:
    python tools/validate_skills.py
"""

from __future__ import annotations

import re
import sys
import tomllib
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
CATALOG = ROOT / "skills.toml"
SKILLS_DIR = ROOT / "skills"

#: The calculation basis a skill may declare, and what each means for how its
#: numbers are produced. Mirrors NeqSim's enum, with `azoth` where NeqSim has
#: `neqsim-java`: the skill drives the validated library rather than a placeholder.
BASIS = ("azoth", "screening", "advisory", "data-retrieval", "hybrid")

#: What a `screening` skill is backed by. `azoth`-basis skills omit it.
#:
#: The P-number tranches are `docs/src/architecture/specification.md`'s order. The
#: named families are the ones that page puts **beyond** P12 - "mechanical design,
#: safety, cost, electrical, automation, `standards/`, `statistics/`" - because a
#: tranche numbering cannot name something that has no tranche.
#:
#: **The field is what would back the skill, and it is required to be honest.** A
#: skill used to have to name a P-number, and the ones whose physics no tranche
#: covers named the nearest - which is how a noise screening and a PSV orifice
#: calculation came to read `P11` after P11 had closed. A named family is the
#: honest answer where the physics belongs to a tree the port does not reach.
TRANCHES = frozenset(f"P{i}" for i in range(13)) | frozenset(
    (
        "mechanical-design",
        "safety",
        "statistics",
        "cost",
        "electrical",
        "automation",
        "standards",
        "fluidmechanics",
        "pvtsimulation",
    )
)

NAME = re.compile(r"^azoth-[a-z0-9]+(-[a-z0-9]+)*$")
VERSION = re.compile(r"^\d+\.\d+\.\d+$")

#: Required keys on every catalog entry, and the required `## ` sections on every
#: `SKILL.md`. Kept in one place so the standard and the checker cannot disagree.
REQUIRED_FIELDS = ("name", "version", "description", "calculation_basis", "path", "tags")
REQUIRED_SECTIONS = (
    "When to Use",
    "Inputs",
    "Outputs",
    "How a calculation runs",
    "Python usage pattern",
    "The keycard",
    "Validation checklist",
    "Common mistakes",
    "Limitations",
    "Related Azoth functionality",
    "References",
)


def problems(entries: list[dict[str, object]]) -> list[str]:
    """The catalog's errors, in the order a reader meets them."""
    found: list[str] = []
    names: dict[str, str] = {}
    for index, entry in enumerate(entries):
        where = f"skills.toml skill[{index}]"
        missing = [field for field in REQUIRED_FIELDS if field not in entry]
        if missing:
            found.append(f"{where}: missing {missing}")
            continue

        name = str(entry["name"])
        if not NAME.match(name):
            found.append(f"{where}: name {name!r} does not match azoth-<kebab>")
        if name in names:
            found.append(f"{where}: name {name!r} duplicates {names[name]}")
        names[name] = where

        version = str(entry["version"])
        if not VERSION.match(version):
            found.append(f"{where}: version {version!r} is not semver x.y.z")

        description = str(entry["description"])
        if "USE WHEN:" not in description:
            found.append(f"{where}: description has no USE WHEN: trigger")

        basis = entry.get("calculation_basis")
        if basis not in BASIS:
            found.append(f"{where}: calculation_basis {basis!r} is not one of {BASIS}")

        tranche = entry.get("tranche")
        if tranche is not None and tranche not in TRANCHES:
            found.append(
                f"{where}: tranche {tranche!r} is not one of P0..P12 or the beyond-P12 "
                f"families {sorted(TRANCHES - {f'P{i}' for i in range(13)})}"
            )
        if basis == "screening" and tranche is None:
            found.append(f"{where}: a screening skill names what will back it")

        path = str(entry.get("path", ""))
        resolved = ROOT / path
        if not resolved.is_file() or not path.startswith("skills/"):
            found.append(f"{where}: path {path!r} is not a file under skills/")
        else:
            found.extend(_section_problems(resolved, where))

        tags = entry.get("tags")
        if not isinstance(tags, list) or not tags:
            found.append(f"{where}: tags must be a non-empty list")

    # Every SKILL.md on disk must be listed, and listed exactly once.
    listed = {ROOT / str(e["path"]) for e in entries if "path" in e}
    on_disk = {p for p in SKILLS_DIR.rglob("SKILL.md")}
    for orphan in sorted(on_disk - listed):
        found.append(f"{orphan.relative_to(ROOT)} is a SKILL.md not listed in skills.toml")
    for stale in sorted(listed - on_disk):
        found.append(f"{stale.relative_to(ROOT)} is listed but does not exist")
    return found


def _section_problems(path: Path, where: str) -> list[str]:
    text = path.read_text(encoding="utf-8")
    headings = {line[3:].strip() for line in text.splitlines() if line.startswith("## ")}
    missing = [section for section in REQUIRED_SECTIONS if section not in headings]
    if missing:
        return [f"{where}: {path.relative_to(ROOT)} is missing section(s) {missing}"]
    return []


def main() -> int:
    if not CATALOG.is_file():
        sys.exit(f"validate_skills: {CATALOG} does not exist")

    try:
        document = tomllib.loads(CATALOG.read_text(encoding="utf-8"))
    except tomllib.TOMLDecodeError as exc:
        sys.exit(f"validate_skills: {CATALOG} does not parse: {exc}")

    if "catalog_version" not in document:
        sys.exit(f"validate_skills: {CATALOG} has no catalog_version")

    entries = document.get("skill", [])
    found = problems(entries)
    if found:
        for problem in found:
            print(f"  ERROR  {problem}", file=sys.stderr)
        print(f"\nvalidate_skills: FAILED with {len(found)} problem(s)", file=sys.stderr)
        return 1

    on_disk = sum(1 for _ in SKILLS_DIR.rglob("SKILL.md"))
    print(f"validate_skills: OK ({len(entries)} skill(s), {on_disk} SKILL.md file(s))")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
