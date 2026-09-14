#!/usr/bin/env python3
"""Compile a checked user data file into the repository's data files.

`keycard.toml` is a *source*, not a runtime input. This turns its `fittings` and
`fluids` sections - the ones `specs/schema/keycard.schema.json` marks `compiled` - into
the canonical CSVs both languages read, exactly as `specs/calcs/*.toml` is compiled into
two registries. Run it, rebuild, and your values are what the calculations use.

    python tools/check_user_data.py keycard.toml     # check first
    python tools/gen_user_data.py   keycard.toml     # then write
    python tools/gen_user_data.py   keycard.toml --check

It writes **in place**, over the committed placeholders, because the Rust core embeds
the file with `include_str!` and needs it at the path the code names. If the values came
from a standard you licensed, committing them redistributes them - read the warning this
prints, and do not run `git add`.

A fluid name that no implementation can name is refused rather than written: both sides
hardcode the fluid list, so `data/fluids/<name>.csv` for an unregistered name would be a
file nothing reads and nothing complains about.
"""

from __future__ import annotations

import argparse
import csv
import io
import sys
import tomllib
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parent))

# The validator, not a second copy of it. A file this tool writes must be a file
# that tool passes, and restating the rules here would make that two claims.
from check_user_data import SCHEMA_VERSION, Report, check_document

REPO_ROOT = Path(__file__).resolve().parent.parent
FITTINGS_PATH = Path("data/fittings/crane_k_factors.csv")
FLUIDS_DIR = Path("data/fluids")

#: Column order, which is the loader's contract on both sides. `id` becomes
#: `fitting_id` because the CSV's header names the thing, not the field.
FITTING_COLUMNS = (
    "fitting_id",
    "family",
    "name",
    "n_ld",
    "f_t_basis",
    "citation",
    "verify_status",
)
FLUID_COLUMNS = (
    "temperature_c",
    "density_kg_m3",
    "dynamic_viscosity_pa_s",
    "citation",
    "verify_status",
)

#: Where a new fluid has to be registered before its table can be read. Listed
#: rather than described, because "somewhere in the Rust core" is not a thing a
#: person can act on.
FLUID_REGISTRATION = (
    "crates/azoth-hydraulics/src/fluids.rs: a `const <NAME>_CSV` beside WATER_CSV",
    "crates/azoth-hydraulics/src/fluids.rs: `embedded_csv` and `embedded_path`",
    "crates/azoth-hydraulics/src/fluids.rs: a `pub fn <name>()` and `provider_for`",
    "crates/azoth-hydraulics/src/fluids.rs: `available_fluids`",
    "python/src/azoth/properties/__init__.py: a TableProvider and `_BUILTINS`",
)


BOX = "=" * 76

#: The loud banner: every row is a placeholder. Its wording is the claim the file makes
#: about itself, and the line a reader must not skim past.
DUMMY_WARNING = """\
#  WARNING: EVERY VALUE IN THIS FILE IS AN ESTIMATED DUMMY VALUE.
#
#  These numbers were NOT taken from Crane TP-410, or from any other standard. They are
#  placeholders of a plausible magnitude. THEY ARE NOT ENGINEERING DATA. DO NOT USE THEM
#  TO SIZE ANYTHING: a pressure drop computed from this file can be wrong by a factor of
#  two or more and will still look reasonable, and no test here can detect that.
#
#  Replace every row with values read from the primary standard before this library is
#  used for design work, and set its verify_status to `verified` with the source in
#  `citation`."""

#: The banner for a file that is not all placeholders, which must not make the claim
#: above. The water and air tables are real published values that were never near this
#: tool, so it claims neither that every row is a placeholder nor that the file was
#: generated: read the column, and do not commit licensed values.
NOT_A_PLACEHOLDER_WARNING = """\
#  WARNING: NOT EVERY VALUE IN THIS FILE IS A PLACEHOLDER.
#
#  Rows marked `unverified` or `verified` are real values; rows marked
#  `estimated_dummy` are not. Which is which is in the verify_status column, and
#  it is the column to read before trusting anything here.
#
#  IF YOU RAN tools/gen_user_data.py TO MAKE THIS FILE: do not commit it. Values
#  that came from a standard you licensed are redistributed by committing them.
#  `git checkout -- data/` brings the shipped placeholders back."""

#: How a number reaches the file, which is the one thing about it a reader cannot recover
#: afterwards: `0.0000200` and `2e-05` are the same float and not the same statement about
#: precision.
NUMBER_NOTE = """\
# Numbers in these rows are rendered from their values, not from the text written in
# the source: `0.0000200` in a keycard becomes `2e-05` here. The value is the same and
# the significant figures are not, so a table whose trailing zeros are meant has to be
# edited in this file rather than generated into it."""

STATUS_BODY = """\
# verify_status: estimated_dummy (placeholder) | unverified (cited, unconfirmed) |
# verified (confirmed against an authoritative copy)."""

FITTINGS_BODY = (
    """\
# K      = f_t * sum(n_ld)   resistance coefficient of the fittings
# n_ld   = L_eq / D          equivalent length ratio, fully turbulent flow
# f_t                        Darcy friction factor for fully turbulent flow
# See specs/calcs/hydraulics/crane_k_factors.toml for the equation and its status.
#
"""
    + STATUS_BODY
)

FLUIDS_BODY = (
    """\
# Read by `azoth pipe` and `azoth.properties`: temperature into density and viscosity,
# by linear interpolation between the points, with no extrapolation past either end.
# Ascending temperature order, because the checker refuses a table that is not.
#
"""
    + STATUS_BODY
)

LICENCE_BLOCK = """\
# Licence: CC-BY-4.0. This data is deliberately not AGPL: the reference
# values are meant to be reusable and citable without a copyleft
# obligation. See LICENSE-CC-BY-4.0."""


def render_header(title: str, *, statuses: dict[str, int], body: str) -> str:
    """The comment block above the data rows.

    Derived from the rows rather than fixed, because the loudest line in it is a
    claim about the file: "every value here is a placeholder" is true of the
    repository's shipped data and false the moment a user supplies real numbers.
    A banner that said it anyway would be worse than none - it would teach a
    reader to skip the line that matters most.
    """
    dummy = statuses.get("estimated_dummy", 0)
    total = sum(statuses.values())
    warning = DUMMY_WARNING if dummy == total else NOT_A_PLACEHOLDER_WARNING
    counts = ", ".join(f"{n} {status}" for status, n in sorted(statuses.items()))

    lines = [
        f"# azoth {title}",
        "#",
        f"# {BOX}",
        *warning.splitlines(),
        f"# {BOX}",
        "#",
        *NUMBER_NOTE.splitlines(),
        "#",
        *body.splitlines(),
        "#",
        f"# Rows in this file: {counts}.",
        "#",
        *LICENCE_BLOCK.splitlines(),
    ]
    return "\n".join(lines) + "\n"


def _status(row: dict[str, Any]) -> str:
    """The `verify_status` column value for a row, derived rather than asked for.

    The shipped CSVs carry this column and the loaders read it, because a result
    computed from placeholder data has to say so at the point of use - the file's
    banner reaches a reader who opens the file, and the column reaches one who never
    does.

    **It is derived, not supplied.** A keycard row has no `verify_status`: the field
    was a form nobody could check, required of every user, and it taught people to
    fill it in. The one thing worth carrying is whether a value is a placeholder, and
    a citation that says DUMMY already says that.
    """
    citation = str(row.get("citation") or "")
    return "estimated_dummy" if "DUMMY" in citation.upper() else "unverified"


def _statuses(rows: list[dict[str, Any]]) -> dict[str, int]:
    counts: dict[str, int] = {}
    for row in rows:
        status = _status(row)
        counts[status] = counts.get(status, 0) + 1
    return counts


def _write_rows(columns: tuple[str, ...], rows: list[dict[str, Any]]) -> str:
    """The data rows as CSV text.

    `csv.writer` rather than string joining, so a citation containing a comma is
    quoted rather than silently becoming an extra column - which would shift every
    value after it into the wrong field and still read as a well-formed file.
    `lineterminator` is set because the default is ``\\r\\n``, which would make a
    file whose line endings depend on the platform that wrote it.
    """
    buffer = io.StringIO()
    writer = csv.writer(buffer, lineterminator="\n")
    writer.writerow(columns)
    for row in rows:
        writer.writerow(
            [
                _cell(_status(row) if column == "verify_status" else row.get(column))
                for column in columns
            ]
        )
    return buffer.getvalue()


def _cell(value: Any) -> str:
    """One value as it will appear in the file.

    A number is rendered from its *value*, so `0.0000200` reaches the file as `2e-05`
    - the same number and not the same claim about precision, which is why every file
    this writes says so in its banner (`NUMBER_NOTE`). A text value is passed through
    untouched, which is what a citation is.
    """
    if value is None:
        return ""
    return value if isinstance(value, str) else str(value)


def render_fittings(rows: list[dict[str, Any]]) -> str:
    """The fittings registry, as bytes."""
    header = render_header("fittings registry", statuses=_statuses(rows), body=FITTINGS_BODY)
    rows_out = [{**row, "fitting_id": row["id"]} for row in rows]
    return header + _write_rows(FITTING_COLUMNS, rows_out)


def render_fluid(fluid: str, rows: list[dict[str, Any]]) -> str:
    """One fluid property table, as bytes."""
    header = render_header(
        f"fluid property table: {fluid}", statuses=_statuses(rows), body=FLUIDS_BODY
    )
    return header + _write_rows(FLUID_COLUMNS, rows)


def built_in_fluids() -> tuple[str, ...]:
    """The fluid names the implementations can actually resolve.

    Read from the Python side rather than listed here: `available_fluids` is the
    one place that already exists, and a second list would be a second answer to
    "which fluids are built in" - the drift this project removes wherever it finds
    it.

    The two languages' lists are asserted equal by a test, so reading one is
    enough.
    """
    sys.path.insert(0, str(REPO_ROOT / "python" / "src"))
    from azoth.properties import available_fluids

    return available_fluids()


def plan(document: dict[str, Any]) -> list[tuple[Path, str]]:
    """What would be written, as (path, text) pairs, without writing anything.

    Separate from writing so that `--check` and the real run cannot disagree about
    what the output is: there is one function that produces it and one that puts
    it on disk.
    """
    outputs: list[tuple[Path, str]] = []

    if "fittings" in document:
        outputs.append((FITTINGS_PATH, render_fittings(document["fittings"])))

    known = built_in_fluids()
    for fluid, rows in document.get("fluids", {}).items():
        name = str(fluid).strip().lower()
        if name not in known:
            raise UnregisteredFluid(name, known)
        outputs.append((FLUIDS_DIR / f"{name}.csv", render_fluid(name, rows)))

    return outputs


class UnregisteredFluid(Exception):
    """A fluid table that no implementation could read.

    Raised rather than written out, and the message carries the fix. A generated
    file nothing loads is worse than a refusal: it looks like data in use, and the
    next person to debug "why is my fluid missing" has to find that out by
    reading Rust.
    """

    def __init__(self, name: str, known: tuple[str, ...]) -> None:
        self.name = name
        self.known = known
        super().__init__(name)


def validate(document: dict[str, Any], path: Path) -> Report:
    """Run the repository's own rules over the document.

    The same functions `check_user_data.py` uses, on the same parsed document, so
    a file that passes there generates here and a file that fails there is refused
    here with the same words.
    """
    report = Report()
    check_document(report, document, str(path))
    return report


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("path", type=Path, help="the user data file to compile")
    parser.add_argument(
        "--check",
        action="store_true",
        help="report what would change and exit non-zero if anything would; write nothing",
    )
    args = parser.parse_args(argv)

    if not args.path.is_file():
        print(f"gen_user_data: {args.path} does not exist", file=sys.stderr)
        return 1

    try:
        document = tomllib.loads(args.path.read_text(encoding="utf-8"))
    except tomllib.TOMLDecodeError as exc:
        print(f"gen_user_data: {args.path} does not parse as TOML: {exc}", file=sys.stderr)
        return 1

    if not isinstance(document, dict):
        print(
            f"gen_user_data: {args.path} should be a mapping with a schema_version key",
            file=sys.stderr,
        )
        return 1

    # Compared as text rather than coerced. The version is an integer in a well-formed
    # document, so the comparison is exact where it matters, and a version written as
    # `1.5` or `"two"` is refused rather than rounded into acceptance.
    if str(document.get("schema_version")) != str(SCHEMA_VERSION):
        print(
            f"gen_user_data: {args.path} declares schema_version "
            f"{document.get('schema_version')!r}; this tool understands {SCHEMA_VERSION}. "
            f"Refusing rather than guessing at what a newer shape meant.",
            file=sys.stderr,
        )
        return 1

    report = validate(document, args.path)
    for error in report.errors:
        print(f"  ERROR    {error}", file=sys.stderr)
    if report.errors:
        print(
            f"\ngen_user_data: {args.path.name} FAILED its checks with "
            f"{len(report.errors)} error(s). Run tools/check_user_data.py for the "
            f"full report. Nothing was written.",
            file=sys.stderr,
        )
        return 1

    try:
        outputs = plan(document)
    except UnregisteredFluid as exc:
        print(
            f"\ngen_user_data: `{exc.name}` is not a built-in fluid, so a table for it "
            f"would be a file nothing reads.\n  Built in: {', '.join(exc.known)}\n\n"
            f"  To add one, register it in:\n"
            + "\n".join(f"    - {where}" for where in FLUID_REGISTRATION)
            + "\n\n  Nothing was written.",
            file=sys.stderr,
        )
        return 1

    changed = []
    for relative, text in outputs:
        target = REPO_ROOT / relative
        current = target.read_text(encoding="utf-8") if target.is_file() else None
        if current != text:
            changed.append((relative, target, text, current is not None))

    if args.check:
        for relative, _, _, existed in changed:
            print(f"gen_user_data: {relative} would {'change' if existed else 'be created'}")
        if changed:
            print(
                f"\ngen_user_data: {len(changed)} file(s) out of date. Run without "
                f"--check to write them.",
                file=sys.stderr,
            )
            return 1
        print(f"gen_user_data: {len(outputs)} file(s) up to date")
        return 0

    for relative, target, text, _ in changed:
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(text, encoding="utf-8")
        print(f"gen_user_data: wrote {relative}")

    total = sum(report.by_status.values())
    print(
        f"gen_user_data: {total} row(s) from {args.path.name} "
        f"({', '.join(f'{n} {s}' for s, n in sorted(report.by_status.items()))})"
    )
    print("\nRebuild so Rust picks the new bytes up: `maturin develop`", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
