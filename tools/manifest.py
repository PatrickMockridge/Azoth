#!/usr/bin/env python3
"""The databank manifest: what is vendored from each upstream, and why the rest is not.

# What this answers

`databank/sources/` holds upstream files this project does not own - NeqSim's
component and interaction tables, and in time Crane TP-410's fitting coefficients
and the fluid property tables. What gets carried across from them is a *decision*,
and the decision is invisible in the output: a reader looking at a compiled CSV
cannot tell whether a column is absent because nobody needed it, because the
physics is out of scope, or because nobody looked.

`databank/manifest.yaml` records the decision per column. This module is what
stops it becoming prose. It is checked, not read.

# Two steps, and the manifest keeps them apart

    upstream file  ->  databank/sources/  ->  the compiled file both languages read

The first step vendors a file whole, revision and all, so the derivation is
reproducible without the upstream project installed. The second step is where
columns are chosen, converted and renamed. Only the second step is a decision, and
only the second step is what `columns:` describes - but the first is what makes the
description checkable in CI, because a manifest listing 170 columns can be compared
against a 170-column file that is actually present.

# The reason vocabulary, and why it is closed

A free-text reason cannot be counted, and counting is the point: the interesting
number is how many columns are absent because a model does not exist yet, because
that number is a roadmap. A closed prefix makes it `grep -c 'reason: not-yet'`.

  no-implementation   a physical property, and the physics it belongs to is out of
                      scope for a cubic-EOS library
  not-yet             a physical property, the physics is in scope, and a model or
                      its wiring is missing. Names the `consumer` that would read it.
  not-a-cubic-input   not a property of a substance at all - an identifier, an
                      index, or a selector naming a model
  empty-upstream      the column carries nothing over the rows this project keeps,
                      so there was no decision to make
  licence             a source that may not be redistributed
  superseded-by       another column in the same file says the same thing
"""

from __future__ import annotations

import csv
import io
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import yaml

ROOT = Path(__file__).resolve().parent.parent
MANIFEST = ROOT / "databank" / "manifest.yaml"
SOURCES = Path("databank") / "sources"

SCHEMA_VERSION = 1

#: The closed reason vocabulary. A reason is `<prefix>: <free text>`, and the prefix
#: is what makes the list countable. See the module docstring.
REASON_PREFIXES = (
    "no-implementation",
    "not-yet",
    "not-a-cubic-input",
    "empty-upstream",
    "licence",
    "superseded-by",
)

#: What may be done with a column. `used` carries a value into the compiled file;
#: `filter` selects rows rather than carrying a value; `dropped` does neither.
DISPOSITIONS = ("used", "filter", "dropped")

#: How a `not-yet` column names what would read it: a spec id, or a roadmap tranche.
#: The two are told apart by shape, so the check knows whether to look the name up.
SPEC_ID = re.compile(r"^[a-z]+\.[a-z_]+$")


@dataclass(frozen=True, slots=True)
class Column:
    """One upstream column, and what was decided about it."""

    name: str
    disposition: str
    reason: str
    #: The compiled column name, for `used`. Checked against the compiled header, so
    #: a rename that skips the manifest is a failure rather than a silent shift.
    as_field: str | None = None
    #: For `not-yet`, what would read it once it exists.
    consumer: str | None = None

    @property
    def prefix(self) -> str:
        return self.reason.split(":", 1)[0].strip()


@dataclass(frozen=True, slots=True)
class RowGroup:
    """Rows of an upstream file that are excluded, and why."""

    selector: str
    count: int
    reason: str


@dataclass(frozen=True, slots=True)
class UpstreamFile:
    id: str
    upstream_path: str
    #: The upstream file, vendored whole. Empty for a source that is quoted rather
    #: than copied - the Crane coefficients, which may not be redistributed.
    vendored_source: str
    #: The file the two languages read, derived from the source.
    compiled_to: str
    source_rows: int
    compiled_rows: int
    columns: tuple[Column, ...]
    rows: tuple[RowGroup, ...] = ()
    #: Columns the compiler writes that no upstream column produces - the per-row
    #: `citation` on components.csv is the one today. Listed rather than left
    #: implicit, because the compiled header is checked against the manifest in both
    #: directions and an unlisted column would read as one nobody decided about.
    generated: tuple[str, ...] = ()


@dataclass(frozen=True, slots=True)
class Upstream:
    id: str
    project: str
    url: str
    version: str
    commit: str
    licence: str
    retrieved: str
    notice: str
    files: tuple[UpstreamFile, ...] = ()


@dataclass(frozen=True, slots=True)
class Manifest:
    schema_version: int
    upstreams: tuple[Upstream, ...]
    not_vendored: tuple[dict[str, Any], ...] = field(default_factory=tuple)

    def files(self) -> tuple[UpstreamFile, ...]:
        return tuple(f for upstream in self.upstreams for f in upstream.files)

    def used(self, file_id: str) -> dict[str, str]:
        """The compiled column each `used` upstream column becomes."""
        return {
            column.name: column.as_field or ""
            for column in self.file(file_id).columns
            if column.disposition == "used"
        }

    def file(self, file_id: str) -> UpstreamFile:
        for candidate in self.files():
            if candidate.id == file_id:
                return candidate
        raise KeyError(file_id)


def _column(raw: dict[str, Any], where: str, problems: list[str]) -> Column:
    name = str(raw.get("name") or "")
    disposition = str(raw.get("disposition") or "")
    reason = str(raw.get("reason") or "")
    as_field = raw.get("as")
    consumer = raw.get("consumer")

    if not name:
        problems.append(f"{where}: a column has no `name`")
    if disposition not in DISPOSITIONS:
        problems.append(
            f"{where}.{name}: disposition {disposition!r} is not one of {list(DISPOSITIONS)}"
        )

    prefix = reason.split(":", 1)[0].strip()
    if disposition != "used" and prefix not in REASON_PREFIXES:
        problems.append(
            f"{where}.{name}: reason {reason!r} starts with {prefix!r}, which is not a "
            f"known prefix. Known: {list(REASON_PREFIXES)}. The prefix is what makes the "
            f"list countable, so a reason without one is a comment."
        )

    if prefix == "not-yet" and not consumer:
        problems.append(
            f"{where}.{name}: a `not-yet` column names the `consumer` that would read "
            f"it, or it is a parking space rather than a plan."
        )
    if disposition == "used" and not as_field:
        problems.append(
            f"{where}.{name}: a `used` column names the compiled column it becomes, in `as:`"
        )
    if disposition != "used" and as_field:
        problems.append(
            f"{where}.{name}: only a `used` column carries `as:`, and this one is {disposition!r}"
        )

    return Column(
        name=name,
        disposition=disposition,
        reason=reason,
        as_field=str(as_field) if as_field else None,
        consumer=str(consumer) if consumer else None,
    )


def _upstream_file(raw: dict[str, Any], problems: list[str]) -> UpstreamFile:
    file_id = str(raw.get("id") or "?")
    return UpstreamFile(
        id=file_id,
        upstream_path=str(raw.get("upstream_path") or ""),
        vendored_source=str(raw.get("vendored_source") or ""),
        compiled_to=str(raw.get("compiled_to") or ""),
        source_rows=int(raw.get("source_rows") or 0),
        compiled_rows=int(raw.get("compiled_rows") or 0),
        columns=tuple(_column(c, file_id, problems) for c in (raw.get("columns") or ())),
        rows=tuple(
            RowGroup(selector=str(r["selector"]), count=int(r["count"]), reason=str(r["reason"]))
            for r in (raw.get("rows") or ())
        ),
        generated=tuple(str(g) for g in (raw.get("generated") or ())),
    )


def read(path: Path = MANIFEST, root: Path = ROOT) -> tuple[Manifest, list[str]]:
    """Parse the manifest and return it with every problem found, parse or rule.

    One function rather than two, because a shape error and a rule error are the same
    kind of news to the person running the check, and splitting them would mean a
    malformed manifest reported half its faults and stopped.
    """
    problems: list[str] = []
    try:
        document = yaml.safe_load(path.read_text(encoding="utf-8"))
    except OSError as error:
        return Manifest(SCHEMA_VERSION, ()), [f"{path}: {error}"]

    if not isinstance(document, dict):
        return Manifest(SCHEMA_VERSION, ()), [f"{path}: the manifest must be a mapping"]

    version = document.get("schema_version")
    if version != SCHEMA_VERSION:
        return Manifest(SCHEMA_VERSION, ()), [
            f"{path}: schema_version is {version!r}, and this reader implements "
            f"{SCHEMA_VERSION}. Bumped rather than extended in place, because a reader "
            f"that ignores what it does not know would accept the newer file and drop "
            f"whatever it added."
        ]

    manifest = Manifest(
        schema_version=version,
        upstreams=tuple(
            Upstream(
                id=str(raw.get("id") or ""),
                project=str(raw.get("project") or ""),
                url=str(raw.get("url") or ""),
                version=str(raw.get("version") or ""),
                commit=str(raw.get("commit") or ""),
                licence=str(raw.get("licence") or ""),
                retrieved=str(raw.get("retrieved") or ""),
                notice=str(raw.get("notice") or ""),
                files=tuple(_upstream_file(f, problems) for f in (raw.get("files") or ())),
            )
            for raw in (document.get("upstreams") or ())
        ),
        not_vendored=tuple(document.get("not_vendored") or ()),
    )

    problems.extend(vendoring_problems(manifest, root))
    return manifest, problems


def _body(path: Path) -> io.StringIO:
    """The file's text from its first record onward, with any leading banner removed.

    Only *leading* comment and blank lines are stripped, not every line beginning
    with `#`. The compiled files carry a banner and need it removed; the upstream
    files carry none, and stripping throughout would corrupt a quoted field whose
    text happens to start with a hash.

    Read as text and split with `keepends`, so a record with a newline inside a
    quoted field survives unchanged - `INTER.csv` has five of them.
    """
    lines = path.read_text(encoding="utf-8").splitlines(keepends=True)
    start = 0
    for index, line in enumerate(lines):
        if line.strip() and not line.lstrip().startswith("#"):
            start = index
            break
    return io.StringIO("".join(lines[start:]))


def reasons(manifest: Manifest) -> dict[str, int]:
    """How many columns carry each reason prefix, in the vocabulary's order.

    Reported by the check rather than left to `grep`, because a column's reason is a
    quoted flow mapping - `reason: "not-yet: ..."` - so `grep 'reason: not-yet'`
    matches nothing. The count is the point of the vocabulary: it is how many columns
    are absent because a model does not exist yet, which is a roadmap.
    """
    counts = dict.fromkeys(REASON_PREFIXES, 0)
    for entry in manifest.files():
        for column in entry.columns:
            if column.disposition == "used":
                continue
            if column.prefix in counts:
                counts[column.prefix] += 1
    return {prefix: count for prefix, count in counts.items() if count}


def header(path: Path) -> tuple[str, ...]:
    """The fields of the first record."""
    first = next(iter(csv.reader(_body(path))), None)
    if not first:
        raise ValueError(f"{path}: no header line found")
    return tuple(field.strip() for field in first)


def data_rows(path: Path) -> int:
    """CSV records after the header, blank records excluded.

    Parsed rather than counted by line, because `INTER.csv` carries five records
    with an embedded newline inside a quoted field - so `wc -l` reports 1309 lines
    for 1304 records, and a line-based count is wrong by five. A compiled file here
    has no embedded newlines, but one helper for both is the only way the two counts
    can be compared.
    """
    records = [r for r in csv.reader(_body(path)) if any(field.strip() for field in r)]
    return max(len(records) - 1, 0)


def spec_ids(root: Path = ROOT) -> set[str]:
    """Every registered id, read from the specs rather than a generated file.

    Read here rather than imported so this check keeps working while the generators
    are being changed, and so it needs nothing beyond PyYAML.
    """
    found: set[str] = set()
    for kind in ("calcs", "models"):
        for path in (root / "specs" / kind).rglob("*.yaml"):
            for line in path.read_text(encoding="utf-8").splitlines():
                if line.startswith("id:"):
                    found.add(line.split(":", 1)[1].strip().strip("'\""))
                    break
    return found


def vendoring_problems(manifest: Manifest, root: Path = ROOT) -> list[str]:
    """Every rule this manifest breaks, as messages. Empty means it holds.

    1. A column list and the file it describes agree, both ways. The manifest is the
       only place a column's fate is written down, so a column in one and not the
       other is a decision nobody made - and a one-way check passes on a rename that
       skipped one side.
    2. The `used` columns and the compiled header agree, both ways, for the same
       reason one step further down.
    3. The row counts are what the manifest says they are, so a refreshed slice
       fails here rather than silently.
    4. A `not-yet` column naming a spec id names one that exists. This is the rule
       that keeps that list a plan rather than a parking space.
    5. Every file under the sources directory is declared, and every declared source
       is present. A file dropped in with no entry is a file nobody decided about.
    """
    messages: list[str] = []

    for entry in manifest.files():
        if entry.vendored_source:
            source = root / entry.vendored_source
            if not source.is_file():
                messages.append(
                    f"{entry.id}: vendored_source {entry.vendored_source!r} does not exist"
                )
            else:
                _agree(
                    messages,
                    f"{entry.id}: {entry.vendored_source}",
                    declared=[c.name for c in entry.columns],
                    actual=header(source),
                    subject="upstream column",
                )
                rows = data_rows(source)
                if rows != entry.source_rows:
                    messages.append(
                        f"{entry.id}: {entry.vendored_source} has {rows} rows and the "
                        f"manifest says {entry.source_rows}"
                    )

        if entry.compiled_to:
            compiled = root / entry.compiled_to
            if not compiled.is_file():
                messages.append(f"{entry.id}: compiled_to {entry.compiled_to!r} does not exist")
            else:
                _agree(
                    messages,
                    f"{entry.id}: {entry.compiled_to}",
                    declared=[
                        *(c.as_field or "" for c in entry.columns if c.disposition == "used"),
                        *entry.generated,
                    ],
                    actual=header(compiled),
                    subject="compiled column",
                )
                rows = data_rows(compiled)
                if rows != entry.compiled_rows:
                    messages.append(
                        f"{entry.id}: {entry.compiled_to} has {rows} rows and the manifest "
                        f"says {entry.compiled_rows}"
                    )

    declared_sources = {
        entry.vendored_source for entry in manifest.files() if entry.vendored_source
    }
    present = {
        str(path.relative_to(root))
        for path in (root / SOURCES).rglob("*")
        if path.is_file() and path.name != "README.md"
    }
    for path in sorted(present - declared_sources):
        messages.append(f"{path}: present under {SOURCES} but not declared in the manifest")
    for path in sorted(declared_sources - present):
        messages.append(f"{path}: declared in the manifest but not present under {SOURCES}")

    known = spec_ids(root)
    for entry in manifest.files():
        for column in entry.columns:
            if column.prefix != "not-yet" or not column.consumer:
                continue
            if not SPEC_ID.match(column.consumer) or column.consumer in known:
                continue
            messages.append(
                f"{entry.id}.{column.name}: `not-yet` names consumer {column.consumer!r}, "
                f"which is not a registered id"
            )

    return messages


def _agree(
    messages: list[str], where: str, *, declared: list[str], actual: tuple[str, ...], subject: str
) -> None:
    """Report both directions of a set comparison, naming the direction that failed."""
    missing = sorted(set(declared) - set(actual))
    if missing:
        messages.append(
            f"{where}: {missing} declared but absent. The manifest and the file disagree "
            f"about what is there."
        )
    undeclared = sorted(set(actual) - set(declared))
    if undeclared:
        messages.append(
            f"{where}: {undeclared} present but declared by no {subject}. A {subject} with "
            f"no entry is one nobody decided about."
        )
