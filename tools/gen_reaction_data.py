#!/usr/bin/env python3
"""Generate `data/reactions/` from NeqSim's element, stoichiometry and reaction tables.

    python tools/gen_reaction_data.py --check                  # against the vendored copies
    python tools/gen_reaction_data.py /tmp/neqsim-check/neqsim # against a fresh checkout

Sibling of `gen_databank.py`, and for the same reason it exists rather than adding these
to that tool's `VENDORED_FILES`: those tables are the component and interaction data, and
these are the reaction data, which a different crate reads. The transformation is the
same shape, though - re-render with the manifest's `as` names as the header, convert no
units and filter no rows - so the two tools share their conventions even though they
write to different directories.

# Why these are compiled at all

`databank/sources/` is the vendored upstream copy, kept so the derivation is reproducible
without NeqSim installed, and it is not shipped: a wheel carries `data/`. So a table a
kernel reads at runtime has to exist under `data/`, which is what this writes. Five
files, and each maps to one upstream file rather than to a merged table, so the manifest
can hold one `compiled_to` per vendor entry and the Pitzer source's extra
`ValidationStatus` column does not have to be a hole in the other two.

# What is not here

No unit conversion, because NeqSim converts none of these either: the columns are the
coefficients of fitted correlations and the identifiers that select a row, and the only
one with a stated unit is the reference temperature, which the correlation does not use.
No row filtering, so `compiled_rows` equals `source_rows` for all five.
"""

from __future__ import annotations

import argparse
import csv
import io
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import manifest as manifest_module

ROOT = Path(__file__).resolve().parent.parent
SOURCES = ROOT / "databank" / "sources" / "neqsim"
OUT_DIR = ROOT / "data" / "reactions"

#: `(manifest id, source file name, compiled file name)`.
COMPILED = (
    ("neqsim/element.csv", "element.csv", "elements.csv"),
    ("neqsim/STOCCOEFDATA.csv", "STOCCOEFDATA.csv", "stoichiometry.csv"),
    ("neqsim/REACTIONDATA.csv", "REACTIONDATA.csv", "REACTIONDATA.csv"),
    (
        "neqsim/REACTIONDATAPITZER.csv",
        "REACTIONDATAPITZER.csv",
        "REACTIONDATAPITZER.csv",
    ),
    (
        "neqsim/REACTIONDATAKENTEISENBERG.csv",
        "REACTIONDATAKENTEISENBERG.csv",
        "REACTIONDATAKENTEISENBERG.csv",
    ),
)


def columns_for(file_id: str) -> tuple[tuple[str, str], ...]:
    """The upstream-to-compiled name pairs for one file, from the manifest.

    Read rather than written out here, so the compiled header and the manifest cannot
    disagree: the manifest is what says a column is read under the name the kernel looks
    for, and a second copy of that mapping is a second thing to keep in step.
    """
    found, problems = manifest_module.read()
    if problems:
        raise SystemExit("gen_reaction_data: " + "; ".join(problems))
    for entry in found.files():
        if entry.id == file_id:
            pairs: list[tuple[str, str]] = []
            for column in entry.columns:
                # A carried column names the field it becomes, and every column of these
                # five files is carried. Refused rather than defaulted, because a
                # column that reaches here without one has no compiled name to be read
                # under and would silently leave a hole in the header.
                if column.as_field is None:
                    raise SystemExit(
                        f"gen_reaction_data: {file_id}.{column.name} is "
                        f"`{column.disposition}` and names no compiled field"
                    )
                pairs.append((column.name, column.as_field))
            return tuple(pairs)
    raise SystemExit(f"gen_reaction_data: {file_id} is not declared in the manifest")


def read_rows(path: Path) -> list[dict[str, str]]:
    """One upstream CSV, as rows keyed by its own header.

    `utf-8-sig`, because NeqSim's files are written on Windows and one of them carries a
    byte-order mark. Left unread it becomes part of the first column's name, and the
    column then looks absent rather than misspelled.
    """
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def build(source: Path, file_id: str, name: str) -> tuple[tuple[str, ...], list[dict[str, str]]]:
    """One file's compiled header and rows."""
    pairs = columns_for(file_id)
    rows = read_rows(source / name)
    if not rows:
        raise SystemExit(f"gen_reaction_data: {name} has no rows")

    missing = [upstream for upstream, _ in pairs if upstream not in rows[0]]
    if missing:
        raise SystemExit(f"gen_reaction_data: {name} has no column {missing}")

    compiled = tuple(as_field for _, as_field in pairs)
    renamed = [{as_field: row[upstream].strip() for upstream, as_field in pairs} for row in rows]
    return compiled, renamed


def render(header: tuple[str, ...], rows: list[dict[str, str]]) -> str:
    """A file's whole contents, as it will be written."""
    buffer = io.StringIO()
    writer = csv.DictWriter(buffer, fieldnames=list(header), lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return buffer.getvalue()


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "neqsim",
        type=Path,
        nargs="?",
        default=SOURCES,
        help=(
            "path to a NeqSim checkout, or to a directory holding the reaction tables "
            f"directly. Defaults to {SOURCES.relative_to(ROOT)}, the vendored copies, "
            f"which is what makes --check runnable without a checkout."
        ),
    )
    parser.add_argument(
        "--check",
        action="store_true",
        help="compare against what is committed instead of writing it",
    )
    args = parser.parse_args(argv)

    source = args.neqsim
    if (source / "src" / "main" / "resources" / "data").is_dir():
        source = source / "src" / "main" / "resources" / "data"

    found, problems = manifest_module.read()
    if problems:
        print("gen_reaction_data: " + "; ".join(problems), file=sys.stderr)
        return 1
    declared = {entry.id: entry for entry in found.files()}

    stale: list[str] = []
    for file_id, name, compiled_name in COMPILED:
        header, rows = build(source, file_id, name)

        entry = declared[file_id]
        if len(rows) != entry.source_rows:
            print(
                f"gen_reaction_data: {name} has {len(rows)} rows and the manifest says "
                f"{entry.source_rows}",
                file=sys.stderr,
            )
            return 1

        text = render(header, rows)
        target = OUT_DIR / compiled_name
        if args.check:
            if not target.exists() or target.read_text(encoding="utf-8") != text:
                stale.append(str(target.relative_to(ROOT)))
        else:
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_text(text, encoding="utf-8")

    if args.check:
        if stale:
            print(
                "gen_reaction_data: stale, re-run without --check: " + ", ".join(stale),
                file=sys.stderr,
            )
            return 1
        print(f"gen_reaction_data: OK ({len(COMPILED)} file(s) match their specs)")
    else:
        print(f"gen_reaction_data: wrote {len(COMPILED)} file(s) to {OUT_DIR.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
