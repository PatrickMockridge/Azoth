#!/usr/bin/env python3
"""Check a user data file before it is used to generate the repository's data.

`keycard.example.yaml` at the repository root describes the format. A user
fills in the values they are entitled to use, and this checks the result.

It checks *provenance*, not values. Nothing here can tell whether a coefficient
is right - only whether the row says where it came from in a form a tool can
fetch, and whether its status is one the rest of the pipeline understands. That
is the same division of labour the specs use: the machine enforces that a claim
is checkable, and a person decides whether it is true.

# Why the rules are imported rather than restated

`Report`, `check_fittings` and `check_fluids` are the row rules, and they live in
`tools/spec_lint.py`, which applies them to the repository's own data files. A
user's file and the repository's files are the same kind of thing, so they are
held to the same rules - and a second copy of those rules here was a second
definition of what a valid citation is, which is exactly the drift this project
organises against.

That is not hypothetical: the two copies had already diverged. The repository's
rejected a row marked `verified` whose citation still said DUMMY, and this one
accepted it, so a file could pass the checker a user is told to run and fail the
one CI runs.

# What this does not do

It does not write anything. Generating the repository's data files from a
checked file is a separate tool, so that checking and writing cannot become one
step that does both when one of them fails.

Usage:
    python tools/check_user_data.py keycard.yaml

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

# The repository's own rules for a data row, imported rather than restated. This
# file used to carry its own copies and they had already diverged: the
# repository's `spec_lint` rejected a `verified` row whose citation still said
# DUMMY, and this one accepted it. One definition now serves both, and the
# constants and functions live beside the rules they belong to.
import schema_registry
from spec_lint import Report, check_fittings, check_fluids

#: The version of the *keycard format*. It lives here rather than in `spec_lint`
#: because it is this format's version, not a calc spec's - and it is imported by
#: `gen_user_data.py` from here for the same reason. Version 2 added `components`,
#: `kij`, `coefficients` and `models`.
SCHEMA_VERSION = 2

#: The sections a keycard may carry, and the row rule each one is checked by here.
#: Paired rather than listed because a section this file does not know about is
#: rejected below, and a list that could drift from the things actually checked
#: would make that rejection wrong in the direction that matters.
KNOWN_SECTIONS = (
    ("fittings", check_fittings),
    ("fluids", check_fluids),
)

#: Sections that carry no row rule *here*, because something else owns them: the
#: JSON Schema checks their shape and `azoth.keycard` checks whether this build can
#: actually use them. Listed rather than skipped silently, so that a section this
#: checker does not know about is still an error.
OTHER_SECTIONS = ("keyholder", "components", "kij", "coefficients", "models")


def check_document(report: Report, document: dict[str, Any], where: str) -> None:
    """Every section of a parsed data file, and nothing that is not one.

    **Unknown top-level sections are an error, not something to skip past.** The
    checker used to read `schema_version`, `fittings` and `fluids` and ignore every
    other key, which meant a misspelled `fitting:` was accepted and its rows were
    dropped - a file that looks like data in use and is read by nothing, which is
    the failure this whole mechanism exists to prevent. It is also what makes the
    format's version mean anything: without it, a file written for a newer shape
    passes the older checker and quietly loses whatever it added.

    Both the checker and the generator call this, on the same parsed document, so
    a file that passes one generates in the other with the same words.
    """
    known = {"schema_version"} | {name for name, _ in KNOWN_SECTIONS} | set(OTHER_SECTIONS)
    unknown = sorted(set(document) - known)
    if unknown:
        report.error(
            where,
            f"has section(s) {unknown} that this format does not define. Known "
            f"sections: {sorted(known)}. Refused rather than ignored: a section "
            f"nothing reads is data that looks in use and is not.",
        )

    present = [name for name in known if name != "schema_version" and name in document]
    for name, check in KNOWN_SECTIONS:
        if name in document:
            check(report, document[name])

    if not present and not unknown:
        report.error(
            where,
            "has no data sections at all. A keycard with nothing in it is most often "
            "a file that failed to save rather than a deliberate empty one.",
        )


def check_against_schema(report: Report, document: dict[str, Any], where: str) -> None:
    """Validate the document against `keycard.schema.json`, the written contract.

    This is where the format is *defined*. The row rules above and the loader in
    `azoth.keycard` are two implementations of parts of it, and this is the
    document itself - so a shape the schema forbids is an error here even when both
    of those would have coped.

    The `$ref`s into `calc.schema.json` are resolved by putting every schema under
    `specs/schema/` in one registry keyed by `$id`; they are separate files and the
    citation and unit definitions are deliberately shared rather than copied.
    """
    try:
        import jsonschema
    except ImportError:  # pragma: no cover
        report.warn(where, "jsonschema is not installed, so the schema was not checked")
        return

    # Every schema in specs/schema/, because the keycard schema `$ref`s the calc
    # schema for what a citation and a unit are, and the calc schema `$ref`s the
    # unit schema for the enum. See `schema_registry`.
    registry = schema_registry.registry()
    validator = jsonschema.Draft202012Validator(
        schema_registry.load("keycard.schema.json"), registry=registry
    )

    for error in sorted(validator.iter_errors(document), key=lambda e: list(e.path)):
        location = ".".join(str(part) for part in error.path) or "<document>"
        report.error(where, f"`{location}`: {error.message}")


def check_with_loader(report: Report, document: dict[str, Any], where: str) -> None:
    """Load the document through `azoth.keycard`, the way a caller would.

    The schema says what shape a keycard may take; the loader says whether *this
    build* can use it - whether an enum admits a model the code implements, whether
    a component parameter is one something reads. A document can satisfy the schema
    and still be one the library would refuse, so both are run and the loader's
    refusal is reported as an error rather than discovered at first use.
    """
    try:
        from azoth import keycard
    except ImportError:  # pragma: no cover
        report.warn(where, "azoth is not importable, so the loader was not run")
        return
    try:
        keycard.use(document, path=Path(where))
    except Exception as exc:  # the loader raises KeycardError; anything else is a bug
        report.error(where, f"the loader refuses it: {exc}")
    finally:
        keycard.clear()


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
    check_against_schema(report, document, path.name)
    check_document(report, document, path.name)
    check_with_loader(report, document, path.name)
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
