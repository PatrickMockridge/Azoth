#!/usr/bin/env python3
"""Check a user data file before it is used to generate the repository's data.

`azoth-data.example.yaml` at the repository root describes the format. A user
fills in the values they are entitled to use, and this checks the result.

It checks *provenance*, not values. Nothing here can tell whether a coefficient
is right - only whether the row says where it came from in a form a tool can
fetch, and whether its status is one the rest of the pipeline understands. That
is the same division of labour the specs use: the machine enforces that a claim
is checkable, and a person decides whether it is true.

# Why the rules are imported rather than restated

`VALID_VERIFY_STATUS` and `SOURCE_REF_PATTERNS` come from `tools/spec_lint.py`,
which applies them to the repository's own data files. A user's file and the
repository's files are the same kind of thing, so they are held to the same
rules - and a second copy of those rules here would be a second definition of
what a valid citation is, which is exactly the drift this project organises
against.

# What this does not do

It does not write anything. Generating the repository's data files from a
checked file is a separate tool, so that checking and writing cannot become one
step that does both when one of them fails.

Usage:
    python tools/check_user_data.py azoth-data.yaml

Exit status is non-zero if any error is found.

There is deliberately no `--quiet`. This tool's job is to say what it found, and
the one thing it must never do quietly is report that a file is full of
placeholders - that note is the reason it prints at all. A flag to suppress it
would be a flag to hide the only thing here that a user needs to be told.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Any

try:
    import yaml
except ImportError:  # pragma: no cover
    sys.exit("check_user_data requires PyYAML: pip install pyyaml")

sys.path.insert(0, str(Path(__file__).resolve().parent))

# The repository's own rules for a data row, shared so there is one definition of
# what a fetchable citation and a legal status are. `check_source` is the rule
# itself, not a restatement of it: this file used to carry its own copy, and the
# two had already diverged - the repository's rejected a `verified` row whose
# citation still said DUMMY, and this one accepted it.
from spec_lint import check_source

SCHEMA_VERSION = 1

#: Required keys per row, by section. Taken from the columns of the data files
#: the generator will produce, so a row that validates here has everything the
#: loader will later ask for.
FITTING_FIELDS = (
    "id",
    "family",
    "name",
    "n_ld",
    "f_t_basis",
    "citation",
    "verify_status",
    "source_ref",
    "source_locator",
)

FLUID_FIELDS = (
    "temperature_c",
    "density_kg_m3",
    "dynamic_viscosity_pa_s",
    "citation",
    "verify_status",
)


class Report:
    """Collected problems, plus the counts a user wants to see."""

    def __init__(self) -> None:
        self.errors: list[str] = []
        self.by_status: dict[str, int] = {}

    def error(self, where: str, message: str) -> None:
        self.errors.append(f"{where}: {message}")

    def count(self, status: str) -> None:
        self.by_status[status] = self.by_status.get(status, 0) + 1


def check_numeric(report: Report, where: str, row: dict[str, Any], field: str) -> float | None:
    """One numeric field, present and parseable."""
    if field not in row:
        report.error(where, f"is missing '{field}'")
        return None
    try:
        return float(row[field])
    except (TypeError, ValueError):
        report.error(where, f"'{field}' is {row[field]!r}, which is not a number")
        return None


def check_fittings(report: Report, raw: Any) -> None:
    """The fitting registry: equivalent-length ratios, looked up by id."""
    if not isinstance(raw, list) or not raw:
        report.error("fittings", "must be a non-empty list of rows")
        return

    seen: set[str] = set()
    for index, row in enumerate(raw):
        if not isinstance(row, dict):
            report.error(f"fittings[{index}]", "is not a mapping")
            continue
        fitting_id = str(row.get("id") or "")
        where = f"fittings[{fitting_id or index}]"

        missing = [field for field in FITTING_FIELDS if field not in row]
        if missing:
            report.error(where, f"is missing {missing}")
            continue

        if not fitting_id:
            report.error(where, "has an empty id")
        elif fitting_id in seen:
            report.error(where, "is a duplicate id; ids are looked up by name")
        seen.add(fitting_id)

        for field in ("family", "name", "f_t_basis"):
            if not str(row.get(field) or "").strip():
                report.error(where, f"has an empty '{field}'")

        n_ld = check_numeric(report, where, row, "n_ld")
        if n_ld is not None and n_ld <= 0:
            report.error(
                where,
                f"has n_ld={n_ld}; a non-positive equivalent length would give a "
                f"negative or zero fitting loss",
            )

        status = check_source(report, where, row)
        if status is not None:
            report.count(status)


def check_fluids(report: Report, raw: Any) -> None:
    """Fluid property tables, one per fluid, interpolated over temperature."""
    if not isinstance(raw, dict) or not raw:
        report.error("fluids", "must be a non-empty mapping of fluid name to table")
        return

    for fluid, rows in raw.items():
        if not isinstance(rows, list) or len(rows) < 2:
            report.error(
                f"fluids.{fluid}",
                "must be a list of at least two rows; one point cannot be interpolated "
                "between, and this provider does not extrapolate",
            )
            continue

        previous_temperature: float | None = None
        for index, row in enumerate(rows):
            if not isinstance(row, dict):
                report.error(f"fluids.{fluid}[{index}]", "is not a mapping")
                continue
            where = f"fluids.{fluid}[{index}]"

            missing = [field for field in FLUID_FIELDS if field not in row]
            if missing:
                report.error(where, f"is missing {missing}")
                continue

            temperature = check_numeric(report, where, row, "temperature_c")
            density = check_numeric(report, where, row, "density_kg_m3")
            viscosity = check_numeric(report, where, row, "dynamic_viscosity_pa_s")

            if density is not None and density <= 0:
                report.error(where, f"has density {density}; it must be positive")
            if viscosity is not None and viscosity <= 0:
                report.error(where, f"has viscosity {viscosity}; it must be positive")
            if temperature is not None:
                if previous_temperature is not None and temperature <= previous_temperature:
                    report.error(
                        where,
                        f"has temperature {temperature}, which does not increase. "
                        f"Rows are interpolated in order, so an unsorted table would "
                        f"give answers that depend on how it happened to be written.",
                    )
                previous_temperature = temperature

            status = check_source(report, where, row)
            if status is not None:
                report.count(status)


def check(path: Path) -> int:
    """Load and check a user data file. Returns a process exit status."""
    try:
        document = yaml.safe_load(path.read_text(encoding="utf-8"))
    except yaml.YAMLError as exc:
        print(f"check_user_data: {path} does not parse: {exc}", file=sys.stderr)
        return 1

    if not isinstance(document, dict):
        print(
            f"check_user_data: {path} should be a mapping with a schema_version key",
            file=sys.stderr,
        )
        return 1

    version = document.get("schema_version")
    if version != SCHEMA_VERSION:
        print(
            f"check_user_data: {path} declares schema_version {version!r}; this tool "
            f"understands {SCHEMA_VERSION}. Refusing rather than guessing at what a "
            f"newer shape meant.",
            file=sys.stderr,
        )
        return 1

    report = Report()
    if "fittings" in document:
        check_fittings(report, document["fittings"])
    if "fluids" in document:
        check_fluids(report, document["fluids"])

    if not report.by_status:
        report.error("document", "has neither a 'fittings' nor a 'fluids' section")

    return report_result(report, path)


def report_result(report: Report, path: Path) -> int:
    """Print what was found, and fail if anything was wrong."""
    for error in report.errors:
        print(f"  ERROR    {error}", file=sys.stderr)

    counts = ", ".join(f"{count} {status}" for status, count in sorted(report.by_status.items()))
    if report.errors:
        print(
            f"\ncheck_user_data: FAILED with {len(report.errors)} error(s) in {path.name}",
            file=sys.stderr,
        )
        return 1

    print(f"check_user_data: OK ({path.name}: {counts})")

    # Placeholders are not an error - a file full of them is a valid starting
    # point - but they are the thing a user most needs to be told about, so they
    # are stated rather than left in a count.
    dummy = report.by_status.get("estimated_dummy", 0)
    if dummy:
        print(
            f"\n  NOTE  {dummy} row(s) are placeholders. They are not engineering "
            f"data; every result computed from them carries an ESTIMATED_DATA "
            f"warning. Replace them before using this for design work.",
            file=sys.stderr,
        )
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("path", type=Path, help="the user data file to check")
    args = parser.parse_args()

    if not args.path.is_file():
        print(f"check_user_data: {args.path} does not exist", file=sys.stderr)
        return 1
    return check(args.path)


if __name__ == "__main__":
    raise SystemExit(main())
