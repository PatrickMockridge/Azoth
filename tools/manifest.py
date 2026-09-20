#!/usr/bin/env python3
"""The databank manifest: what is vendored from each upstream, and why the rest is not.

# What this answers

`databank/sources/` holds upstream files this project does not own - NeqSim's
component and interaction tables, and in time Crane TP-410's fitting coefficients
and the fluid property tables. What gets carried across from them is a *decision*,
and the decision is invisible in the output: a reader looking at a compiled CSV
cannot tell whether a column is absent because nobody needed it, because the
physics is out of scope, or because nobody looked.

`databank/manifest.toml` records the decision per column. This module is what
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
number is how many columns are absent because a model is not ported, because that
number is the porting backlog. `tools/check_manifest.py` prints the tally.

  not-ported          a physical property whose model exists in NeqSim and is not
                      ported here. Names the class or package that would close it,
                      so the count is a work list rather than a boundary.
  not-yet             the model exists here and the data does not reach it yet.
                      Names the `consumer` - a registered id - that would read it.
  not-a-value         not a property of a substance at all: an identifier, an index,
                      or a selector naming a model
  empty-upstream      the column carries nothing over the rows this project keeps,
                      so there was no decision to make
  licence             a source that may not be redistributed
  superseded-by       another column in the same file says the same thing
"""

from __future__ import annotations

import csv
import io
import re
import tomllib
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent.parent
MANIFEST = ROOT / "databank" / "manifest.toml"
SOURCES = Path("databank") / "sources"

SCHEMA_VERSION = 1

#: The closed reason vocabulary. A reason is `<prefix>: <free text>`, and the prefix
#: is what makes the list countable. See the module docstring.
REASON_PREFIXES = (
    "not-ported",
    "not-yet",
    "not-a-value",
    "empty-upstream",
    "licence",
    "superseded-by",
    # Carried data whose upstream model **cannot be run at all**: the class that would
    # read it is never constructed anywhere in the checkout. Distinct from `not-ported`,
    # which is work azoth has not done: this is work NeqSim has not done either, so
    # counting it as backlog overstates the port.
    "unreachable-upstream",
)

#: What may be done with a column.
#:
#: `used` carries a value into the compiled file and something reads it. `vendored`
#: carries it and nothing reads it yet - the distinction that stops this file
#: conflating *the data is here* with *a model consumes it*, which is the question
#: "can every calculation's inputs be supplied?" and a different question from
#: "which calculations are ported?". `filter` selects rows rather than carrying a
#: value; `dropped` does neither.
DISPOSITIONS = ("used", "vendored", "filter", "dropped")

#: The dispositions that put a value in the compiled file, and so must agree with its
#: header, carry a field name, and state a unit.
CARRIED = ("used", "vendored")

#: The unit a value is stored in, or `unknown`.
#:
#: **A unit is recorded, never guessed.** Where NeqSim's reader states one - a setter
#: taking a unit, a division by 1000, an addition of 273.15 - it is written here as the
#: unit the compiled file holds. Where the reader states nothing, the value is carried
#: exactly as NeqSim stores it and this says `neqsim-internal`: not a guess, and
#: countable, so the set with no stated unit is visible rather than assumed away.
NO_STATED_UNIT = "neqsim-internal"

#: How a `not-yet` column names what would read it: a spec id, or a roadmap tranche.
#: The two are told apart by shape, so the check knows whether to look the name up.
SPEC_ID = re.compile(r"^[a-z]+\.[a-z_]+$")


@dataclass(frozen=True, slots=True)
class Column:
    """One upstream column, and what was decided about it."""

    name: str
    disposition: str
    reason: str
    #: The compiled column name, for `used` and `vendored`. Checked against the compiled
    #: header, so a rename that skips the manifest is a failure rather than a silent
    #: shift.
    as_field: str | None = None
    #: The unit the compiled value is in. Required of every carried column, because a
    #: number whose unit is not written down is a number nobody can use.
    unit: str | None = None
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
    #: Whether `columns:` was written at all.
    #:
    #: A file that is not a table - `fiscal_parameters.json` - has no columns to
    #: disposition, and one whose format defeats the parser has no columns this file can
    #: honestly state. Both are vendored whole, and the difference between "declared and
    #: present" and "declared and present, columns not yet read" has to be visible or the
    #: second reads as an oversight. Absent `columns:` means the second.
    columns_declared: bool = True
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
            if column.disposition in CARRIED
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
    unit = raw.get("unit")
    consumer = raw.get("consumer")

    if not name:
        problems.append(f"{where}: a column has no `name`")
    if disposition not in DISPOSITIONS:
        problems.append(
            f"{where}.{name}: disposition {disposition!r} is not one of {list(DISPOSITIONS)}"
        )

    prefix = reason.split(":", 1)[0].strip()
    if disposition == "vendored":
        # A vendored column is carried and unread. Its reason says what would read it,
        # or that nothing does - it is not a prefix from the dropped vocabulary, because
        # a column that is present is not one this file decided against.
        pass
    elif disposition != "used" and prefix not in REASON_PREFIXES:
        problems.append(
            f"{where}.{name}: reason {reason!r} starts with {prefix!r}, which is not a "
            f"known prefix. Known: {list(REASON_PREFIXES)}. The prefix is what makes the "
            f"list countable, so a reason without one is a comment."
        )

    if disposition != "vendored" and prefix == "not-yet" and not consumer:
        problems.append(
            f"{where}.{name}: a `not-yet` column names the `consumer` that would read "
            f"it, or it is a parking space rather than a plan."
        )
    if disposition != "vendored" and prefix == "not-ported" and "`" not in reason:
        problems.append(
            f"{where}.{name}: a `not-ported` reason names the NeqSim class or package "
            f"that would close it, in backticks. Without one this list stops being a "
            f"porting backlog and becomes somewhere to put a column."
        )
    if disposition in ("used", "vendored") and not as_field:
        problems.append(
            f"{where}.{name}: a `{disposition}` column names the compiled column it "
            f"becomes, in `as:`"
        )
    if disposition in ("used", "vendored") and not unit:
        problems.append(
            f"{where}.{name}: a carried column states the unit its compiled value is in, "
            f"or `{NO_STATED_UNIT}` where NeqSim's reader states none. A number whose "
            f"unit is not written down is one nobody can use."
        )
    if disposition not in ("used", "vendored") and as_field:
        problems.append(
            f"{where}.{name}: only a carried column carries `as:`, and this one is {disposition!r}"
        )

    return Column(
        name=name,
        disposition=disposition,
        reason=reason,
        as_field=str(as_field) if as_field else None,
        unit=str(unit) if unit else None,
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
        columns_declared="columns" in raw,
        rows=tuple(
            RowGroup(selector=str(r["selector"]), count=int(r["count"]), reason=str(r["reason"]))
            for r in (raw.get("rows") or ())
        ),
        generated=tuple(str(g) for g in (raw.get("generated") or ())),
    )


def read(path: Path = MANIFEST) -> tuple[Manifest, list[str]]:
    """Parse the manifest and report what is wrong with the file itself.

    Parsing and cross-checking are separate steps, and this is the first. The second,
    :func:`validate`, compares the manifest against the files on disk - which a
    generator cannot do before it has written them, and would otherwise make it
    impossible to add a column: the file would be stale until regenerated, and
    regeneration would refuse until the file was current.
    """
    problems: list[str] = []
    try:
        document = tomllib.loads(path.read_text(encoding="utf-8"))
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

    return manifest, problems


def validate(manifest: Manifest, root: Path = ROOT) -> list[str]:
    """Where the manifest and the files on disk disagree.

    The second step of :func:`read`'s work, and the one `tools/check_manifest.py` runs:
    every vendored file present and declared, every column of it accounted for in both
    directions, and the row counts what the manifest says.
    """
    return vendoring_problems(manifest, root)


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


#: Which of NeqSim's `COMPTYPE` values the compiled table keeps.
#:
#: Read from NeqSim's own `COMPTYPE` column rather than guessed. Of `COMP.csv`'s 389
#: rows, 62 are `ion`, 19 `GEN`, 18 carry no type, and 4 are `ice`, `salt`, `seawater`
#: and `asphaltene`; the nine types below are the other 286.
#:
#: **The ion rows are kept for what an electrolyte model reads, not for their critical
#: constants.** Those columns are filler: `Pc = 290.89 bar` on 27 of the 62 rows,
#: `omega = 0.344` on 29, `Vc = 99.0 cm3/mol` on 38 - one shared default per column,
#: over different row sets, because a cubic has no notion of an ion. The remainder
#: inherit a neutral parent's numbers (`MDEA+` carries MDEA's `Tc`, `Pc` and `omega`
#: unchanged) or carry values with no stated source. 13 distinct `TC`, 26 `PC` and 25
#: `ACSFACT` values across 62 rows is what the mixture looks like.
#:
#: What `PhasePitzer`, `PhaseDesmukhMather` and the Furst variants actually read is
#: `IONICCHARGE`, `DeshMatIonicDiameter`, `MOLARMASS` and `DIELECTRICPARAMETER1..5`, so
#: the rows are carried and `mixture_of` refuses a cubic over one rather than computing
#: with the filler. That refusal is what makes keeping them safe; without it, vendoring
#: these rows would be shipping plausible-looking wrong numbers, which is the failure
#: this project is organised against.
#:
#: Here rather than in `gen_databank`, which used to own it, because the `empty-upstream`
#: claims are about *the kept rows* and this is the rule that decides which those are.
#: One definition, imported by the generator.
KEEP_TYPES = frozenset(
    {
        "HC",
        "inert",
        "other",
        "glycol",
        "acid",
        "alcohol",
        "amine",
        "chlorine",
        "water",
        "ion",
    }
)


def _records(path: Path) -> list[dict[str, str]]:
    """The file as dictionaries, blanks stripped, parsed the way `header` does."""
    reader = csv.DictReader(_body(path))
    return [
        {(k or "").strip(): (v or "").strip() for k, v in record.items()}
        for record in reader
        if any((v or "").strip() for v in record.values())
    ]


def kept_component_names(root: Path = ROOT) -> set[str]:
    """The substances of NeqSim's `COMP.csv` the vendored slice carries, lower case."""
    rows = _records(root / SOURCES / "neqsim" / "COMP.csv")
    return {r["NAME"].lower() for r in rows if r.get("COMPTYPE", "") in KEEP_TYPES}


#: The pages that define the port's tranches: the specification's table and the two
#: roadmap pages. Read rather than listed, so a tranche renamed in a page is a tranche
#: this stops accepting - the same "decide the claim" move `empty_upstream_problems`
#: makes for the other prefix.
ROADMAP_PAGES = (
    "ROADMAP.md",
    "docs/src/agentic/roadmap.md",
    "docs/src/architecture/specification.md",
)

#: How a tranche is spelled in those pages: `P7`, or `Tier 2`.
TRANCHES = (
    re.compile(r"\bP(\d{1,2})\b"),
    re.compile(r"\bTier (\d)\b"),
)


def roadmap_tranches(root: Path = ROOT) -> set[str]:
    """Every tranche name the roadmap pages carry, as `P7` and `Tier 2` are spelled."""
    found: set[str] = set()
    for page in ROADMAP_PAGES:
        path = root / page
        if not path.is_file():
            continue
        text = path.read_text(encoding="utf-8")
        found.update(f"P{m.group(1)}" for m in TRANCHES[0].finditer(text))
        found.update(f"Tier {m.group(1)}" for m in TRANCHES[1].finditer(text))
    return found


def consumer_problems(manifest: Manifest, root: Path = ROOT) -> list[str]:
    """Check every `consumer` a column names, whatever its reason prefix.

    A `consumer` says *what would read this column*, and it is the field that makes an
    unwritten column a plan rather than an omission. Two forms are allowed - a registered
    spec id, or `roadmap:<tranche>` - and both are decidable, so both are decided here.

    The `roadmap:` form used to be skipped entirely: `SPEC_ID` cannot match a string
    containing `:`, so the check `continue`d and eleven columns pointed at
    `roadmap:2b`, a tranche neither roadmap page defines and nothing compared them
    against. The prose beside them was free to be stale as a result - it claimed
    `eos.viscosity` was "planned and unimplemented" while `specs/models/eos/viscosity.toml`
    was registered - which is what made this worth closing rather than widening the skip.

    The check is on *every* prefix rather than on `not-yet` alone, because the field means
    the same thing wherever it appears and a rule that only fires on one prefix is a rule
    that stops firing the moment the prefix changes.
    """
    messages: list[str] = []
    known = spec_ids(root)
    tranches = roadmap_tranches(root)
    for entry in manifest.files():
        for column in entry.columns:
            consumer = column.consumer
            if not consumer:
                continue
            where = f"{entry.id}.{column.name}"
            if consumer.startswith("roadmap:"):
                tranche = consumer.removeprefix("roadmap:")
                if tranche not in tranches:
                    messages.append(
                        f"{where}: consumer {consumer!r} names tranche {tranche!r}, which "
                        f"none of {', '.join(ROADMAP_PAGES)} defines. A pointer nobody can "
                        f"resolve is a pointer nothing checks."
                    )
                continue
            if not SPEC_ID.match(consumer):
                messages.append(
                    f"{where}: consumer {consumer!r} is neither a registered id nor "
                    f"`roadmap:<tranche>`"
                )
            elif consumer not in known:
                messages.append(f"{where}: consumer {consumer!r} is not a registered id")
    return messages


#: How a document cites an upstream revision: the full object name, wherever it sits in the
#: line. An abbreviation is not matched, because the ambiguous prefix is exactly what a
#: reader cannot fetch, and a bare 40-hex token is a git object name and nothing else -
#: `databank/README.md` writes one in parentheses, two lines after the word `commit`.
COMMIT_CITATION = re.compile(r"\b([0-9a-f]{40})\b")

#: A line that is nothing but one word, which is how `NOTICE` heads each upstream - and
#: how it does *not* head a TOML key or a Markdown heading. **Three characters or more**,
#: because `databank/README.md` draws its stages with a lone `v` on a line of its own and
#: treating that as a heading would silently stop that page being read.
_SECTION_HEADING = re.compile(r"^[A-Za-z][\w.-]{2,}$")


def neqsim_citations(text: str) -> list[str]:
    """Every commit named by a passage about NeqSim, in order, with duplicates.

    **A document may cite several upstreams.** `NOTICE` gives NeqSim and `lean-units` a
    section each, and only NeqSim's revision is this manifest's business, so a passage is
    NeqSim's when it follows a bare-word heading naming NeqSim or when the line itself names
    it. A document with no such heading - a spec, or `databank/README.md` - is one passage
    throughout, because everything it cites an upstream for is NeqSim.
    """
    found: list[str] = []
    section = "neqsim"
    for line in text.splitlines():
        stripped = line.strip()
        if _SECTION_HEADING.match(stripped):
            section = stripped.lower()
        if section == "neqsim" or "neqsim" in line.lower():
            found.extend(COMMIT_CITATION.findall(line))
    return found


def citation_problems(manifest: Manifest, root: Path = ROOT) -> list[str]:
    """Check the NeqSim revision the specs cite against the one the manifest records.

    A spec's `references` names the upstream revision it was ported from, and the manifest
    names the revision vendored under `databank/sources/`. Those are the same revision, so a
    document naming a different one cites a file this project does not hold - and nothing
    compared them until now.

    **They had come apart.** The refresh that moved the manifest to `805cf0f` left seventeen
    specs and `NOTICE` citing `dedba873`, and neither hash is wrong on its own terms, so no
    gate could see it: the difference was between two documents, not between a document and a
    file. Measured afterwards, the two revisions differ in one of the twenty-eight cited
    sources, so the stale citation was not harmless - it named a revision whose source says
    something else.

    The check is on the *shape* rather than on a known list, so a citation to a third revision
    is caught the same way, and a document citing none is left alone.
    """
    recorded = next((u.commit for u in manifest.upstreams if u.id == "neqsim"), "")
    if not recorded:
        return []
    documents = [*(root / "specs").rglob("*.toml"), root / "NOTICE", MANIFEST.parent / "README.md"]
    messages: list[str] = []
    for path in sorted(documents):
        if not path.is_file():
            continue
        for found in sorted(set(neqsim_citations(path.read_text(encoding="utf-8")))):
            if found != recorded:
                messages.append(
                    f"{path.relative_to(root)}: cites NeqSim commit {found}, and the manifest "
                    f"records {recorded} for the revision vendored under {SOURCES}. A citation "
                    f"to a revision the project does not hold describes a file nobody can check "
                    f"this port against."
                )
    return messages


def empty_upstream_problems(manifest: Manifest, root: Path = ROOT) -> list[str]:
    """Check every `empty-upstream` claim against the vendored file.

    `empty-upstream` is a *claim about data* - this column carries no values over the
    rows we keep - and a claim like that is decidable, so it is decided here rather than
    believed. It is not a stylistic preference: a populated column recorded as empty is
    a column dropped out of the roadmap `manifest.reasons` counts, and the one that was
    found this way was `REFERENCESTATETYPE`, which NeqSim's `ComponentGE.fugcoef` reads
    to choose between the Raoult and Henry reference states.

    The row scope is the row scope the compiled file has: for `COMP.csv`, the substances
    `KEEP_TYPES` selects; for `INTER.csv`, the pairs with both ends among them.
    """
    messages: list[str] = []
    for entry in manifest.files():
        if not entry.vendored_source:
            continue
        claims = [
            c for c in entry.columns if c.disposition == "dropped" and c.prefix == "empty-upstream"
        ]
        if not claims:
            continue
        source = root / entry.vendored_source
        if not source.is_file():
            continue
        rows = _records(source)
        if entry.id.endswith("COMP.csv") and not entry.id.endswith("COMP_EXT.csv"):
            names = kept_component_names(root)
            in_scope = [r for r in rows if r.get("NAME", "").lower() in names]
        elif entry.id.endswith("INTER.csv"):
            names = kept_component_names(root)
            in_scope = [
                r
                for r in rows
                if r.get("COMP1", "").lower() in names and r.get("COMP2", "").lower() in names
            ]
        else:  # pragma: no cover - no other vendored file makes this claim
            in_scope = rows
        for column in claims:
            populated = [
                r for r in in_scope if r.get(column.name, "") not in ("", "0", "0.0", "-1")
            ]
            if populated:
                sample = sorted({r[column.name] for r in populated})[:4]
                messages.append(
                    f"{entry.id}.{column.name}: recorded `empty-upstream`, but "
                    f"{len(populated)} of the {len(in_scope)} in-scope rows carry a "
                    f"value - {sample}. The column is not empty, so either its "
                    f"disposition or its reason is wrong."
                )
    return messages


def reasons(manifest: Manifest) -> dict[str, int]:
    """How many columns carry each reason prefix, in the vocabulary's order.

    Reported by the check rather than left to `grep`, because a column's reason is a
    quoted flow mapping - `reason: "not-yet: ..."` - so `grep 'reason: not-yet'`
    matches nothing. The count is the point of the vocabulary: it is how many columns
    are held back because a model does not exist yet, which is a roadmap.

    **Carried columns are counted, and that is the fix to an inversion.** This skipped
    them, so the one number a reader wanted - how many vendored columns are waiting on a
    model - was the one it did not print, while the printed `not-ported` counted only
    *dropped* columns, which are not waiting on anything.
    """
    counts = dict.fromkeys(REASON_PREFIXES, 0)
    for entry in manifest.files():
        for column in entry.columns:
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
    """Every registered id, read from the specs rather than a generated file."""
    found: set[str] = set()
    for kind in ("calcs", "models"):
        for path in (root / "specs" / kind).rglob("*.toml"):
            for line in path.read_text(encoding="utf-8").splitlines():
                if line.startswith("id = "):
                    found.add(line.split("=", 1)[1].strip().strip("'\""))
                    break
    return found


#: The crate file that parses the compiled component table, and the anchors that find the
#: list of column names it indexes. Both are names in the source, so a rename is something
#: this reports rather than a check that quietly stops firing.
COMPONENT_PARSER = Path("crates") / "azoth-eos" / "src" / "databank.rs"
COMPONENTS_CONSTANT = "COMPONENTS_CSV"
PARSER_FUNCTION = "fn parse_components()"


def indexed_column_problems(manifest: Manifest, root: Path = ROOT) -> list[str]:
    """A column the component parser indexes is `used`, whatever the manifest says.

    `used` and `vendored` differ by one claim - something reads it - and that claim is the
    one a manifest cannot check for itself. The parser can: its allow-list is the set of
    names it looks up in every record, so being in it *is* the claim. `LJDIAMETER` and
    `SCHWARTZENTRUBER1`-`3` were `vendored` while all four were indexed into every `Entry`,
    and the tally printed four columns a registered model reads as carried-and-unread. It
    was wrong for three tranches, and nothing but reading the Rust could have told anyone.

    **One file, one direction.** The allow-list belongs to `parse_components`, whose
    compiled output is one table; `COMP_EXT.csv` and the reaction tables carry columns of
    the same names with no reader at all, so applying this to every vendored file would
    fail on sixty-five rows that are right. The reverse direction is false for the same
    reason in miniature: `cas`, `liquid_density_kg_per_m3` and `viscosity_correction_factor`
    are `used` and read by something other than this parser.
    """
    base = root.resolve()
    try:
        source = (base / COMPONENT_PARSER).read_text(encoding="utf-8")
    except OSError as error:
        return [f"{COMPONENT_PARSER}: {error}"]

    include = re.search(rf'const {COMPONENTS_CONSTANT}: &str = include_str!\("([^"]+)"\)', source)
    if include is None:
        return [
            f"{COMPONENT_PARSER}: no `{COMPONENTS_CONSTANT}` include, so the table this rule "
            f"is about cannot be named"
        ]
    try:
        compiled = str(
            (base / COMPONENT_PARSER.parent / include.group(1)).resolve().relative_to(base)
        )
    except ValueError:
        return [
            f"{COMPONENT_PARSER}: `{COMPONENTS_CONSTANT}` points outside the tree, at "
            f"{include.group(1)!r}"
        ]

    entry = next((item for item in manifest.files() if item.compiled_to == compiled), None)
    if entry is None:
        return [f"{compiled}: parsed by {COMPONENT_PARSER}, and no manifest file compiles to it"]

    start = source.find(PARSER_FUNCTION)
    opened = source.find("for name in [", start) if start >= 0 else -1
    if start < 0 or opened < 0:
        return [
            f"{COMPONENT_PARSER}: no `{PARSER_FUNCTION}` with a `for name in [`, so the "
            f"columns it indexes cannot be read and this rule decides nothing. It reports "
            f"rather than passes, because a rule that cannot find its subject is the one "
            f"case where silence looks like agreement."
        ]
    indexed = set(re.findall(r'"([A-Za-z0-9_]+)"', source[opened : source.find("]", opened)]))

    return [
        f"{entry.id}.{column.name}: disposition is {column.disposition!r} and "
        f"{COMPONENT_PARSER} indexes it into every record. A column the parser reads is "
        f"`used`; `vendored` says nothing reads it, and the tally prints it as carried."
        for column in entry.columns
        if column.disposition == "vendored" and column.as_field in indexed
    ]


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
            elif entry.columns_declared:
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
                        *(c.as_field or "" for c in entry.columns if c.disposition in CARRIED),
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

    messages.extend(empty_upstream_problems(manifest, root))

    messages.extend(consumer_problems(manifest, root))

    messages.extend(indexed_column_problems(manifest, root))

    messages.extend(citation_problems(manifest, root))

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
