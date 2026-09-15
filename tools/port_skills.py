#!/usr/bin/env python3
"""Port NeqSim's community skills into azoth as screening placeholders.

One-time migration. Reads each `SKILL.md` under the NeqSim community-skills checkout,
extracts its name and `USE WHEN:` description from the YAML frontmatter (by regex,
not a YAML parser), and writes a minimal placeholder under
`skills/<domain>/<skill>/SKILL.md` whose `Related Azoth functionality` section names
the tranche that will back it. The NeqSim content is Apache-2.0 and is credited in
`NOTICE`.

The placeholders have no examples, tests or `src/`: the method they point at is not
ported yet, so there is no code to run. They satisfy the same eleven-section standard
as an `azoth`-basis skill, so an agent always knows the shape.

Usage:
    python tools/port_skills.py /tmp/neqsim-check/neqsim-community-skills
"""

from __future__ import annotations

import re
import sys
import tomllib
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
CATALOG = ROOT / "skills.toml"

#: Domain -> (calculation_basis, port tranche). A `None` tranche means no port is
#: needed — the skill is advisory or data-retrieval, not a calculation azoth will
#: back. `screening` skills must name one, which the validator enforces.
DOMAINS: dict[str, tuple[str, str | None]] = {
    "process": ("screening", "P11"),
    "safety": ("screening", "P11"),
    "flow-assurance": ("screening", "P9"),
    "pvt": ("screening", "P4"),
    "subsurface": ("screening", "P6"),
    "environment": ("advisory", None),
    "field-development": ("advisory", None),
    "subsea": ("data-retrieval", None),
    "reporting": ("advisory", None),
    "engineering-data": ("data-retrieval", None),
}

#: `name: neqsim-hydrate-screening` and `description: "..."` in the frontmatter.
_NAME = re.compile(r"^name:\s*(\S+)\s*$")
_DESCRIPTION = re.compile(r'^description:\s*"(.*)"\s*$')

TEMPLATE = """# {title}

{description}

This is an `{basis}` skill and, until azoth ports the calculation, a placeholder: it
produces no azoth number. {backing}

## When to Use

{description}

## Inputs

To be declared by the calculation that backs this skill.

## Outputs

To be declared by the calculation that backs this skill.

## How a calculation runs

No azoth calculation runs yet. The skill exists so that the task, and its place in
the roadmap, is named rather than lost.

## Python usage pattern

None — the calculation is not ported.

## The keycard

Not applicable until the calculation is ported.

## Validation checklist

- [ ] Confirm the task actually needs this calculation.
- [ ] Use NeqSim's validated engine for real work until the tranche lands.

## Common mistakes

| Symptom | Cause | Fix |
|---|---|---|
| A placeholder used as a result | The calculation is not ported | Use NeqSim's engine |

## Limitations

Azoth has not ported this calculation; nothing here produces a number.

## Related Azoth functionality

None yet. {backing} The validated engine today is NeqSim's, credited in `NOTICE`.

## References

- NeqSim: https://github.com/equinor/neqsim
- NeqSim community skills: https://github.com/equinor/neqsim-community-skills
"""


def title_of(skill: str) -> str:
    """`hydrate-screening` -> `Hydrate screening`."""
    return skill.replace("-", " ").capitalize()


def backing_of(tranche: str | None) -> str:
    """The one-line `backing` clause for a placeholder."""
    if tranche is None:
        return "No azoth backing is planned — this is advisory or data, not a calculation."
    return f"Azoth backs this at tranche {tranche}."


def frontmatter(text: str) -> tuple[str, str]:
    """The `name` and `description` from a NeqSim SKILL.md, with safe defaults."""
    name, description = "", ""
    for line in text.splitlines():
        if (match := _NAME.match(line)) and not name:
            name = match.group(1)
        elif (match := _DESCRIPTION.match(line)) and not description:
            description = match.group(1)
        if name and description:
            break
    return name, description


def write_placeholder(
    domain: str, skill: str, description: str, basis: str, tranche: str | None
) -> Path:
    target = ROOT / "skills" / domain / skill / "SKILL.md"
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(
        TEMPLATE.format(
            title=title_of(skill),
            description=description or f"A {title_of(skill)} placeholder.",
            basis=basis,
            backing=backing_of(tranche),
        ),
        encoding="utf-8",
    )
    return target


def catalog_entry(
    name: str, description: str, basis: str, path: str, domain: str, tranche: str | None
) -> str:
    lines = [
        "[[skill]]",
        f'name = "{name}"',
        'version = "0.1.0"',
        f'description = "{description}"',
        f'calculation_basis = "{basis}"',
        f'path = "{path}"',
    ]
    if tranche is not None:
        lines.append(f'tranche = "{tranche}"')
    lines.append(f'tags = ["{domain}"]')
    return "\n".join(lines)


def main() -> int:
    if len(sys.argv) != 2:
        sys.exit("usage: python tools/port_skills.py <neqsim-community-skills checkout>")

    source = Path(sys.argv[1])
    skills_dir = source / "skills"
    if not skills_dir.is_dir():
        sys.exit(f"port_skills: {skills_dir} is not a skills directory")

    existing = {
        str(entry["name"])
        for entry in tomllib.loads(CATALOG.read_text(encoding="utf-8")).get("skill", [])
    }

    entries: list[str] = []
    written = 0
    for path in sorted(skills_dir.rglob("SKILL.md")):
        relative = path.relative_to(skills_dir)
        domain, skill = relative.parts[0], relative.parts[1]
        if domain not in DOMAINS:
            print(f"  skip  {relative}: unknown domain", file=sys.stderr)
            continue

        name, description = frontmatter(path.read_text(encoding="utf-8"))
        azoth_name = name.replace("neqsim-", "azoth-", 1) if name.startswith("neqsim-") else name
        if azoth_name in existing:
            print(f"  skip  {relative}: {azoth_name} already hand-written", file=sys.stderr)
            continue

        basis, tranche = DOMAINS[domain]
        write_placeholder(domain, skill, description, basis, tranche)
        entries.append(
            catalog_entry(
                azoth_name,
                description,
                basis,
                f"skills/{domain}/{skill}/SKILL.md",
                domain,
                tranche,
            )
        )
        written += 1

    out = ROOT / "skills" / "placeholders.toml"
    out.write_text("\n\n".join(entries) + "\n", encoding="utf-8")
    print(f"port_skills: wrote {written} placeholder SKILL.md file(s)")
    print(f"port_skills: catalog entries in {out.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
