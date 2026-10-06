#!/usr/bin/env python3
"""Does the tree read the columns `databank/manifest.toml` says it does not?

# Why this exists

A manifest row is a claim about the tree: this column was carried across, this one was carried
and nothing reads it, this one was left because the model does not exist yet. `tools/manifest.py`
decides everything about a row that the manifest alone can decide - that the file exists, that the
compiled header agrees, that a `not-yet` column's `consumer` names a tranche or an id that
resolves.

**The one thing it cannot decide is whether *nothing reads it* is true**, and that is the claim
the disposition makes. `consumer_problems` records the first half of this: the `consumer` pointer
was unvalidated, eleven columns pointed at `roadmap:2b`, and the prose beside them was free to go
stale. The pointer resolves now. **Whether the column is read is still nobody's question**, and it
has been answered wrongly twice by hand:

* `PARACHOR` says so in its own reason - "It was dispositioned `not-ported`, and the rate-based
  tranche is what made that false".
* `LJEPS` and `LIQUIDCONDUCTIVITY1`/`2`/`3` are the same, and were found by reading rather than by
  any check: `eos.phase_transport` names its components, and its kernel reads the Lennard-Jones
  pair and the conductivity polynomial straight off the databank entry.

# What it decides, and how

**The alias-to-field mapping is read, not tabled.** A column is carried under an alias and lands in
a Rust field or a Python one, and both readers say so on the line the value is assigned:

    lennard_jones_energy: number(&record, index["ljeps"], "ljeps", row)?,
    liquid_conductivity: [ number(&record, index["liquidconductivity1"], ... ]   # the array form
    lennard_jones_energy=float(row["ljeps"]),

so this file walks up from each `index["<alias>"]` to the field that encloses it. A table here
would be a second declaration of what the readers already declare, and the second one is the one
that goes stale.

A field is **consumed** when something outside the readers *reads a member of it* - `.<field>` on a
value, in Rust or Python. A bare mention is not enough: a field's own declaration in a struct is a
name, and a name is not a reader. That is the same distinction `check_manifest.py` draws between
*the data is here* and *a model consumes it*.

# What it cannot decide, and says so

* **A field read only inside a reader** - a default, a fallback in `databank.rs` itself - reports as
  unconsumed. The tool errs toward *not* reporting, which is the direction `check_manifest.py`
  already chose when it refuses to guess whether a vendored slice is current.
* **Whether a reason is *right* in any richer sense.** A reason is prose, and prose is a reader's
  job. The mechanical question is exactly one: does the tree read this.

    python tools/audit_manifest.py --report    # every unread-shaped row, with its verdict
    python tools/audit_manifest.py --check     # fail on a row the tree reads
"""

from __future__ import annotations

import re
import sys
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parent))
import manifest as manifest_module

ROOT = Path(__file__).resolve().parent.parent

#: The readers whose assignments are *carrying* a value rather than consuming it, and the
#: generated tables, which restate a spec rather than read a field.
NOT_CONSUMERS = (
    "crates/azoth-eos/src/databank.rs",
    "python/src/azoth/eos/components.py",
)

#: The reason prefixes that say **nothing reads this column**. `not-ported` is work azoth has
#: not done, the three `*-upstream` words are readings of NeqSim's own reachability, and each is
#: a claim this tool is the reader of. `used` is *not* in this list: a `used` column is one
#: something reads, whatever its reason's prose says - except that `LJEPS` is `used` with a
#: `not-ported` reason, which is why the prefix is read from the reason and not the disposition.
UNREAD_PREFIXES = (
    "not-ported",
    "unreachable-upstream",
    "uncalled-upstream",
    "unread-upstream",
)

#: `index["ljeps"]` in the Rust reader and `row["ljeps"]` in the Python one.
_ALIAS = re.compile(r'(?:index|row)\["(\w+)"\]')

#: A struct field or a keyword argument opening a block: `lennard_jones_energy: [` or
#: `        liquid_conductivity: [`.
_FIELD = re.compile(r"^\s*(\w+):\s")

#: The readers' own vocabulary, which the walk up from an alias meets before the field does:
#: a function parameter (`row: usize,`) sits at the same indentation as a struct field and
#: `associationscheme` is read inside a function whose arguments are three of these.
_NOT_A_FIELD = ("index", "number", "record", "row", "fn", "let", "mut")

#: A field whose **name is a word the tree uses for its own objects**. `.name` and
#: `.reference` appear in seventy-odd files that have nothing to do with the column being
#: asked about, so a match on one is not evidence and reporting it would be the noisy check
#: `CONTRIBUTING.md` refuses.
#:
#: **This is a language fact and not a list of this repository's fields** - the same kind of
#: datum as `_NOT_A_FIELD` above. It costs a possible finding: `UNIFACGroupParam.csv.Name` is
#: `not-ported` and its claim may be false, and this cannot tell. The alternative is a check
#: that reports seventy files for one column, which teaches a reader to skim.
_COMMON_WORD = (
    "name",
    "reference",
    "value",
    "class",
    "type",
    "id",
    "main",
    "secondary",
    "index",
    "row",
    "record",
    "source",
)


def _walk_up_to_field(lines: list[str], index: int) -> str | None:
    """The field an alias's assignment belongs to, looking upward from it.

    The single-line form names the field on the same line; the array form opens the field with
    `field: [` and the assignments follow inside it, so the search has to walk up rather than
    pattern-match one line. `index` and `number` are excluded because they are the reader's own
    vocabulary and a search that stopped on them would attribute every column to it.
    """
    for back in range(index, max(index - 30, -1), -1):
        match = _FIELD.match(lines[back])
        if match and match.group(1) not in _NOT_A_FIELD:
            return match.group(1)
    return None


def _python_field(line: str, alias: str) -> str | None:
    """The field on a Python assignment line: `lennard_jones_energy=float(row["ljeps"])`."""
    match = re.search(rf"^\s*(\w+)\s*=\s*[^=]*row\[\"{re.escape(alias)}\"\]", line)
    return match.group(1) if match else None


def carried_as() -> dict[str, str]:
    """Every alias, and the field the readers put it in."""
    out: dict[str, str] = {}

    rust = (ROOT / "crates" / "azoth-eos" / "src" / "databank.rs").read_text(encoding="utf-8")
    lines = rust.split("\n")
    for number, line in enumerate(lines):
        for alias in _ALIAS.findall(line):
            field = _walk_up_to_field(lines, number)
            if field:
                out.setdefault(alias, field)

    python = (ROOT / "python" / "src" / "azoth" / "eos" / "components.py").read_text(
        encoding="utf-8"
    )
    for line in python.split("\n"):
        for alias in _ALIAS.findall(line):
            field = _python_field(line, alias)
            if field:
                out.setdefault(alias, field)
    return out


def _searchable() -> list[Path]:
    """Every source in the tree except the readers and the generated tables."""
    files: list[Path] = []
    for pattern in ("crates/**/*.rs", "python/src/**/*.py"):
        for path in ROOT.glob(pattern):
            relative = path.relative_to(ROOT).as_posix()
            if relative in NOT_CONSUMERS:
                continue
            # A generated table restates a spec; it is not a reader of anything.
            if path.stem.endswith("_gen"):
                continue
            files.append(path)
    return files


def consumers(field: str) -> list[str]:
    """The files that *read a member* of `field`, rather than naming it.

    `\\.<field>` is a member access, which is what a reader is. A bare mention is excluded
    because a struct's own declaration of a field is a name and not a reader - and an assignment
    to it is a write, so `self.field = ` and `.field =` are excluded too.
    """
    pattern = re.compile(rf"\.{re.escape(field)}\b(?!\s*=[^=])")
    if field in _COMMON_WORD:
        return []
    out: list[str] = []
    for path in _searchable():
        text = path.read_text(encoding="utf-8")
        if pattern.search(text):
            out.append(path.relative_to(ROOT).as_posix())
    return out


def verdicts() -> tuple[list[dict[str, Any]], int]:
    """`(rows with a verdict, how many were read)`."""
    parsed = manifest_module.read()
    document, problems = parsed[0], parsed[1]
    if problems:
        raise SystemExit(
            "audit_manifest: the manifest does not read\n  " + "\n  ".join(problems[:5])
        )

    mapping = carried_as()
    if not mapping:
        raise SystemExit("audit_manifest: no alias was read out of the readers")

    out: list[dict[str, Any]] = []
    read = 0
    for entry in document.files():
        for column in entry.columns:
            if column.prefix not in UNREAD_PREFIXES:
                continue
            read += 1
            # **A column of a file that is not compiled reaches no field.** `COMP_EXT.csv` is
            # vendored whole, so its `NAME` is not the compiled table's `name` - keying by the
            # alias alone borrowed the compiled file's field and reported a `used` column as
            # unread.
            alias = (column.as_field or column.name).lower()
            field = mapping.get(alias) if entry.compiled_to else None
            row = {
                "where": f"{entry.id}.{column.name}",
                "alias": alias,
                "field": field,
                "disposition": column.disposition,
                "prefix": column.prefix,
                "consumer": column.consumer,
            }
            if field is None:
                row["verdict"] = "unparsed"
                row["readers"] = []
            else:
                row["readers"] = consumers(field)
                row["verdict"] = "READ" if row["readers"] else "unread"
            out.append(row)
    return out, read


def main() -> None:
    rows, read = verdicts()
    if not read:
        raise SystemExit(
            "audit_manifest: no column claims to be unread, so this checked nothing - the "
            "manifest is not where this looks, or every column is read"
        )

    false_rows = [row for row in rows if row["verdict"] == "READ"]
    if "--check" in sys.argv:
        if false_rows:
            raise SystemExit(
                "audit_manifest: the manifest says nothing reads these, and the tree does\n  "
                + "\n  ".join(
                    f"{row['where']}: `{row['alias']}` is read by {', '.join(row['readers'])}"
                    for row in false_rows
                )
            )
        print(f"audit_manifest: {read} unread-shaped row(s), every one unread")
        return

    counts: dict[str, int] = {}
    for row in rows:
        counts[row["verdict"]] = counts.get(row["verdict"], 0) + 1
    print(f"{read} unread-shaped row(s): {counts}")
    for row in rows:
        detail = f"  read by {', '.join(row['readers'])}" if row["readers"] else ""
        print(
            f"  {row['verdict']:9} {row['where']:34} [{row['disposition']}, {row['prefix']}]"
            f"{detail}"
        )


if __name__ == "__main__":
    main()
