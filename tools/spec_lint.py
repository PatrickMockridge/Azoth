#!/usr/bin/env python3
"""Validate every calc spec against the JSON Schema, then against the project's
semantic rules.

The schema catches malformed specs. It cannot catch a spec that is well formed but
incoherent - a range check on a quantity the calc never computes, a worked example
that does not match the function signature, a fitting id that does not exist, a
value attributed to a source nobody verified. Those are the failure modes that
matter here, because they look like validation while doing nothing. This tool
exists to catch them.

Usage:
    python tools/spec_lint.py            # check everything
    python tools/spec_lint.py --quiet    # only report problems
    python tools/spec_lint.py --spec-dir <dir>   # check a different tree of specs

Exit status is non-zero if any error-level problem is found. Warnings do not fail
the run, but they are printed, because a warning nobody reads is the same as no
check at all.
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
import tomllib
from collections.abc import Iterator
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

try:
    from jsonschema import Draft202012Validator
except ImportError:  # pragma: no cover
    sys.exit("spec_lint requires jsonschema: pip install jsonschema")


ROOT = Path(__file__).resolve().parent.parent
SPEC_DIR = ROOT / "specs" / "calcs"
MODEL_DIR = ROOT / "specs" / "models"
CASE_DIR = ROOT / "specs" / "cases"
SCHEMA_PATH = ROOT / "specs" / "schema" / "calc.schema.json"
MODEL_SCHEMA_PATH = ROOT / "specs" / "schema" / "model.schema.json"
CASE_SCHEMA_PATH = ROOT / "specs" / "schema" / "case.schema.json"

# Quantities that name a state rather than a number. These have no units and are
# not function parameters, so they are allowed in range checks without appearing
# in `inputs`.
STATE_KEYWORDS = {"phase"}

# The tolerance a worked example may declare before it stops being a meaningful
# check. Relative tolerance, so 1e-3 means "agrees to 0.1%".
MAX_SENSIBLE_TOLERANCE = 0.05

# Permitted provenance states for a data row. **Optional**: a row may carry one
# to say what the library knows about its own shipped data, and a user's row need
# not carry one at all.
#
# This vocabulary is for the CSVs in `data/` and for nothing else. **There is no
# equivalent for a spec file, and there must not be one.** A spec declaring whether
# it is verified would be the library grading the engineer's judgement, which is not
# its job: the choice of which equation of state and which data apply to a situation
# belongs to the engineer and to their keycard. See
# `docs/src/architecture/specification.md`, "What the library owes instead of a
# status field". A schema edit is the only way a new
# spec field can appear at all - every spec schema sets `additionalProperties: false`
# - so if you are here to add a `status` or a `provenance` key to one, that edit is
# the thing to not make.
VALID_VERIFY_STATUS = {"verified", "unverified", "estimated_dummy"}

# A source reference used to have to be machine-fetchable - `arweave:<txid>`,
# `doi:`, `https://` - with a separate locator giving the place inside the
# document. That was built to make a *single transcribed value* auditable, which is
# the case it fits. Against a vendored databank of 258 rows whose provenance is
# institutional it stops fitting: there is no per-row URL, and inventing one would
# be a citation-shaped thing that is not a citation. A `citation` naming the source
# in prose is the whole record now.


@dataclass
class Report:
    """Collected problems, in the order they were found."""

    errors: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    # Every provenance label seen, counted - including `unstated`, which is now the
    # common case. Used by the user-data checker to report what it found rather than
    # only what was wrong, and to tell an empty file from an unchecked one.
    by_status: dict[str, int] = field(default_factory=dict)

    def error(self, spec: str, msg: str) -> None:
        self.errors.append(f"{spec}: {msg}")

    def warn(self, spec: str, msg: str) -> None:
        self.warnings.append(f"{spec}: {msg}")

    def count(self, status: str) -> None:
        self.by_status[status] = self.by_status.get(status, 0) + 1


def check_citation(report: Report, where: str, row: dict[str, Any]) -> str | None:
    """A data row's optional citation, and the one thing that *is* checked.

    **Nothing here is required.** A row may say where its values came from and may
    say they are placeholders; it does not have to do either. The engineer supplying
    data is responsible for its provenance and for their right to use it, which is a
    professional responsibility this tool cannot discharge and does not pretend to.

    The single rule kept is the one that is not about trust but about contradiction:
    a row whose `verify_status` says `estimated_dummy` while its citation says nothing
    about it is a row whose two halves disagree, and a reader skimming the citation
    will not realise the number is a placeholder. That check is about the library
    understanding *its own* shipped data, which really does contain placeholders -
    not about auditing a user's.

    Returns the status when one is given, `None` when none is, so a caller can count
    what it found.
    """
    status = row.get("verify_status")
    citation = str(row.get("citation") or "")

    if status is None:
        # Counted under its own label rather than skipped, so a caller can tell "no
        # rows were checked" from "rows were checked and none declared a status".
        report.count("unstated")
        return None
    if status not in VALID_VERIFY_STATUS:
        report.error(
            where,
            f"verify_status is {status!r}; expected one of {sorted(VALID_VERIFY_STATUS)}",
        )
        return None
    if status == "estimated_dummy" and "DUMMY" not in citation.upper():
        report.error(
            where,
            "marked estimated_dummy but its citation does not say so. The marker and "
            "the citation are the same claim in two places, and a reader skimming the "
            "citation will not realise the number is a placeholder.",
        )
    return status


#: Required keys per row, by section. Taken from the columns of the data files
#: the generator produces, so a row that validates has everything the loader
#: will later ask for.
FITTING_FIELDS = (
    "id",
    "family",
    "name",
    "n_ld",
    "f_t_basis",
)

FLUID_FIELDS = (
    "temperature_c",
    "density_kg_m3",
    "dynamic_viscosity_pa_s",
)


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
    """The fitting registry: equivalent-length ratios, looked up by id.

    Shared, not restated: `check_user_data.py` applies these to a user's file and
    this file applies them to the repository's, and there is one definition of
    what a valid row is rather than two that can drift.
    """
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

        for column in ("family", "name", "f_t_basis"):
            if not str(row.get(column) or "").strip():
                report.error(where, f"has an empty '{column}'")

        n_ld = check_numeric(report, where, row, "n_ld")
        if n_ld is not None and n_ld <= 0:
            report.error(
                where,
                f"has n_ld={n_ld}; a non-positive equivalent length would give a "
                f"negative or zero fitting loss",
            )

        status = check_citation(report, where, row)
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

            status = check_citation(report, where, row)
            if status is not None:
                report.count(status)


def check_fluid_tables(report: Report, fluids_dir: Path) -> None:
    """Every fluid table the repository ships, checked unconditionally.

    Not driven by a spec's `data:` block, because the CLI and the property provider
    read these and no calc does, so there is no spec to hang the check off.
    """
    for path in sorted(fluids_dir.glob("*.csv")):
        rows = load_csv_rows(path)
        if not rows:
            report.error(display_path(path).as_posix(), "has a header but no data rows")
            continue
        check_fluids(report, {path.stem: rows})


def load_csv_rows(path: Path) -> list[dict[str, str]]:
    """A data file's rows, with the `#` comment banner stripped.

    The banner is not decoration - it carries the copyright statement and the
    meaning of `verify_status` - so the parser tolerates it rather than the file
    having to be machine-only.

    **Blank lines go too, and that is not tidiness.** `csv.DictReader` takes the
    first line it is handed as the header row, so a single blank line between the
    banner and the header makes every field name `None` and every value land in
    `restkey` - the whole file parses into rows that look like
    `{None: ['0', '999.8', ...]}`. The fluid tables have that blank line and the
    fittings registry does not, which is the only reason the fittings side was
    working: it was one blank line away from the same silent nothing.
    """
    with path.open(newline="", encoding="utf-8") as handle:
        content = [line for line in handle if line.strip() and not line.lstrip().startswith("#")]
    return [dict(row) for row in csv.DictReader(content)]


def load_fittings_csv(path: Path) -> dict[str, dict[str, str]]:
    """Read the fittings registry, skipping the leading comment block.

    The comment block is not decoration - it carries the copyright statement and
    the meaning of verify_status - so the parser has to tolerate it rather than
    the file having to be machine-only.
    """
    rows = load_csv_rows(path)
    if not rows:
        raise ValueError(f"{path} has no header row")
    return {row["fitting_id"]: row for row in rows}


def display_path(path: Path) -> Path:
    """A spec's path as a reader should see it.

    Specs normally live under the repository root. `--spec-dir` lets the test
    suite point this tool at a temporary tree, which is outside `ROOT`, so the
    relative form has to degrade rather than raise.
    """
    try:
        return path.relative_to(ROOT)
    except ValueError:
        return path


def check_schema(
    report: Report, validator: Draft202012Validator, spec_dir: Path
) -> list[tuple[Path, dict[str, Any]]]:
    """Validate each spec against the JSON Schema. Returns the ones that parsed."""
    parsed: list[tuple[Path, dict[str, Any]]] = []
    specs = sorted(spec_dir.rglob("*.toml"))
    if not specs:
        report.error("specs", f"no spec files found under {spec_dir}")
        return parsed

    for path in specs:
        rel = display_path(path)
        try:
            spec = tomllib.loads(path.read_text(encoding="utf-8"))
        except tomllib.TOMLDecodeError as exc:
            report.error(str(rel), f"TOML does not parse: {exc}")
            continue

        errors = sorted(validator.iter_errors(spec), key=lambda e: list(e.absolute_path))
        for err in errors:
            where = "/".join(str(p) for p in err.absolute_path) or "<root>"
            report.error(str(rel), f"schema violation at {where}: {err.message}")
        if not errors:
            parsed.append((rel, spec))
    return parsed


def check_identity(report: Report, rel: Path, spec: dict[str, Any]) -> None:
    """The id, the filename, and both implementation paths must agree.

    If they drift, the docs generator silently produces a page for a calc whose
    implementation cannot be found, and the signature-agreement test in Step 3 has
    nothing to compare against.
    """
    spec_id = spec["id"]
    stem = rel.stem
    if spec_id.split(".")[-1] != stem:
        report.error(
            str(rel),
            f"id '{spec_id}' does not end with the filename stem '{stem}'",
        )

    namespace, _, func = spec_id.rpartition(".")
    expected_python = f"azoth.{namespace}.{func}"
    if spec["implementations"]["python"] != expected_python:
        report.error(
            str(rel),
            f"implementations.python is '{spec['implementations']['python']}' "
            f"but the id implies '{expected_python}'",
        )
    expected_rust = f"azoth_{namespace}::{func}"
    if spec["implementations"]["rust"] != expected_rust:
        report.error(
            str(rel),
            f"implementations.rust is '{spec['implementations']['rust']}' "
            f"but the id implies '{expected_rust}'",
        )


def check_identifier_names(report: Report, rel: Path, spec: dict[str, Any]) -> None:
    """Names must be usable as parameters and result fields in both languages.

    Two different rules, because the two sections play different roles:

    * **Inputs** are function parameters, and they are deliberately the symbols
      from the published equation - `L`, `D`, `Re`, `rho`. Keeping them means a
      reader can check a signature against the paper it came from, and it is why
      the Rust implementations carry a scoped `allow(non_snake_case)`. The real
      requirement is that they are legal identifiers in both languages.

    * **Outputs** are attributes on a result object. Those must be snake_case,
      because `result.dP` is not an attribute either language would naturally
      produce - and a spec declaring `dP` while the code produces `dp` is a drift
      no numerical test can catch, since the values agree perfectly and only the
      attribute name differs. That happened: darcy_weisbach declared `dP`.

    Names in both sections must additionally not collide with language keywords.
    """
    import keyword

    legal = re.compile(r"^[A-Za-z_][A-Za-z0-9_]*$")
    snake = re.compile(r"^[a-z][a-z0-9_]*$")
    keywords = (
        set(keyword.kwlist)
        | set(keyword.softkwlist)
        | {
            # Rust keywords that would break a parameter or field name.
            "as",
            "break",
            "const",
            "continue",
            "crate",
            "else",
            "enum",
            "extern",
            "false",
            "fn",
            "for",
            "if",
            "impl",
            "in",
            "let",
            "loop",
            "match",
            "mod",
            "move",
            "mut",
            "pub",
            "ref",
            "return",
            "self",
            "static",
            "struct",
            "super",
            "trait",
            "true",
            "type",
            "unsafe",
            "use",
            "where",
            "while",
            "async",
            "await",
            "dyn",
            "abstract",
            "become",
            "box",
            "do",
            "final",
            "macro",
            "override",
            "priv",
            "try",
            "typeof",
            "unsized",
            "virtual",
            "yield",
        }
    )

    for section in ("inputs", "outputs"):
        for name in spec[section]:
            if not legal.match(name):
                report.error(
                    str(rel),
                    f"{section} key '{name}' is not a legal identifier, so it cannot be a "
                    f"parameter or field name in Rust or Python.",
                )
            elif name in keywords:
                report.error(
                    str(rel),
                    f"{section} key '{name}' is a reserved word in Rust or Python.",
                )
            elif section == "outputs" and not snake.match(name):
                suggestion = name.lower()
                report.error(
                    str(rel),
                    f"output '{name}' is not snake_case. Outputs become result fields, so "
                    f"they must be the name both languages produce - use '{suggestion}'. "
                    f"Keep '{name}' as the symbol in `equation` if that is what the source "
                    f"uses.",
                )


def check_model(report: Report, rel: Path, spec: dict[str, Any]) -> None:
    """The semantic checks a model spec needs beyond schema validation.

    `gen_models.py` validates against the schema before it generates anything, but a
    schema cannot say whether a bound is *evaluable* - which is the same gap
    `check_range_checks` fills for calcs, and the reason it is called here rather
    than left to the generator.

    The checks that do not apply are the calc ones: a model has no `equation`, no
    `worked_example` and no `tests` list, because what it pins down is a procedure.
    Its cases carry inputs and expected outputs like a calc's, so those are checked
    against the declared outputs here.
    """
    check_identity_model(report, rel, spec)
    check_range_checks(report, rel, spec)
    check_model_cases(report, rel, spec)
    check_prose_length(report, rel, spec)


def check_identity_model(report: Report, rel: Path, spec: dict[str, Any]) -> None:
    """The id, the filename and the declared implementations must agree.

    The same rule a calc follows, and for the same reason: the id is what
    `_models_gen.model()` is keyed by and what the two implementations are named
    after, so a mismatch is a model nothing can look up.
    """
    model_id = spec["id"]
    expected_stem = model_id.split(".")[-1]
    if rel.stem != expected_stem:
        report.error(
            str(rel),
            f"id '{model_id}' has final segment '{expected_stem}' but the file is "
            f"'{rel.stem}.toml'. The generated table is keyed by the id and the docs "
            f"page is named after the file, so the two have to agree.",
        )

    namespace = model_id.split(".")[0]
    expected_python = f"azoth.{namespace}.{expected_stem}"
    expected_rust = f"azoth_{namespace}::{expected_stem}"
    implementations = spec["implementations"]
    for language, expected in (("python", expected_python), ("rust", expected_rust)):
        actual = implementations[language]
        if actual != expected:
            report.error(
                str(rel),
                f"implementations.{language} is '{actual}' but the id implies '{expected}'.",
            )


#: The keys a case may state that the model does not take.
#:
#: A case names its fluid by component names and the runner resolves them; `associating`
#: is part of that resolution rather than of the model's signature, so it may appear in a
#: case's inputs and must not be handed to the function. Named here rather than matched by
#: a rule, because each one is a decision about the boundary.
BOUNDARY_KEYS = frozenset({"associating", "eos"})


def check_model_cases(report: Report, rel: Path, spec: dict[str, Any]) -> None:
    """Every case's inputs and expected outputs name declared quantities.

    A case naming an input the model does not take, or asserting an output it does
    not produce, is a test nothing can run - and the failure would otherwise be a
    KeyError in the generated table rather than a lint error here.

    **Except for a boundary key**, which is a property of the *runner* rather than
    of the model. A case states its fluid as component names, and the runner resolves
    them through the databank; whether that fluid runs the Wertheim association is
    part of how the runner builds it, because the model takes the mixture whole and
    would have nothing to say about a second, possibly contradictory, argument beside
    it. The set is closed and named here rather than matched by pattern, so a new one
    is a decision.
    """
    inputs = set(spec["inputs"])
    outputs = set(spec["outputs"])
    seen: set[str] = set()

    for case in spec["cases"]:
        case_id = case["id"]
        if case_id in seen:
            report.error(str(rel), f"duplicate case id '{case_id}'")
        seen.add(case_id)

        unknown_inputs = sorted(set(case["inputs"]) - inputs - BOUNDARY_KEYS)
        if unknown_inputs:
            report.error(
                str(rel),
                f"case '{case_id}' passes {unknown_inputs}, which the model does not "
                f"declare as inputs. Declared: {sorted(inputs)}",
            )
        unknown_expected = sorted(set(case["expected"]) - outputs)
        if unknown_expected:
            report.error(
                str(rel),
                f"case '{case_id}' asserts {unknown_expected}, which the model does "
                f"not declare as outputs. Declared: {sorted(outputs)}",
            )

        # The *declared* inputs the case provides: a boundary key is not an input, and
        # `eos` is a boundary key for one model and a declared input for another
        # (`ge_nrtl_flash` declares it, because its vapour is half of what it is).
        #
        # **An optional input may be left out, and leaving one out is not a gap.**
        # `mixer`'s `outlet_pressure` is the first the process layer has: a case that
        # omits it is exercising the branch where the outlet takes the lowest feed
        # pressure, and the model's own range-check warning says the bound went unrun.
        # Requiring it would make the two branches one case short of covering them.
        required = {
            name for name, declared in spec["inputs"].items() if not declared.get("optional")
        }
        supplied = set(case["inputs"]) & required
        if len(supplied) != len(required):
            missing = sorted(required - supplied)
            report.warn(
                str(rel),
                f"case '{case_id}' supplies {len(supplied)} of {len(required)} "
                f"required inputs (missing {missing}); a case that leaves one out is "
                f"not reproducing the model's whole signature.",
            )


def check_range_checks(report: Report, rel: Path, spec: dict[str, Any]) -> None:
    """Every range check must be evaluable, and must be honest about it.

    This is the check that stops the class of defect where a calc advertises a
    validated range it never actually evaluates - a range check that cannot run is
    worse than no check, because it reads as validation.
    """
    inputs = spec["inputs"]
    outputs = spec["outputs"]
    optional_inputs = {n for n, d in inputs.items() if d.get("optional", False)}

    for check in spec["valid_range"]:
        quantity = check["quantity"]

        # `enum` is a bound kind the schema permits and neither implementation can
        # evaluate. It is refused here rather than emitted and ignored, because a
        # check that never fires is the failure this whole file is written against -
        # and it was reachable: the schema's `anyOf` accepts `enum` as a bound, and
        # until `eos.rachford_rice_binary` no spec used `equals` either, so the whole
        # family of non-interval bounds sat there looking supported.
        if check.get("enum") is not None:
            report.error(
                str(rel),
                f"range check on '{quantity}' uses `enum`, which neither "
                f"implementation can evaluate. A range check resolves its quantity to "
                f"a float; an enum bound compares strings. Emitting it would produce "
                f"a check that silently never fires.",
            )
            continue

        known = quantity in inputs or quantity in outputs or quantity in STATE_KEYWORDS

        if not known:
            sources = check.get("computed_from")
            if not sources:
                report.error(
                    str(rel),
                    f"range check on '{quantity}' which is neither an input, an output, "
                    f"nor a state keyword, and has no computed_from. It cannot be evaluated.",
                )
                continue
            unknown = [s for s in sources if s not in inputs]
            if unknown:
                report.error(
                    str(rel),
                    f"range check on '{quantity}' says it is computed from {unknown}, "
                    f"which are not declared inputs",
                )
                continue
            if quantity not in outputs:
                report.warn(
                    str(rel),
                    f"range check on derived quantity '{quantity}' which is not an output, "
                    f"so callers never see the value that was checked",
                )

        # A check that depends on an optional input cannot always run. That is
        # acceptable, but the implementation must say so at runtime rather than
        # passing silently.
        deps = set(check.get("computed_from", []))
        if quantity in optional_inputs:
            deps.add(quantity)
        missing_optional = deps & optional_inputs
        if missing_optional:
            report.warn(
                str(rel),
                f"range check on '{quantity}' depends on optional input(s) "
                f"{sorted(missing_optional)}; the implementation must emit "
                f"RANGE_CHECK_SKIPPED when they are absent",
            )

        if check["severity"] == "error" and check.get("when") == "inside":
            report.warn(
                str(rel),
                f"range check on '{quantity}' is an 'inside' band with severity 'error'; "
                f"band conditions are usually warnings since the value is still computable",
            )


def check_worked_example(report: Report, rel: Path, spec: dict[str, Any]) -> None:
    """The worked example must be callable and must actually demonstrate the calc."""
    inputs = spec["inputs"]
    outputs = spec["outputs"]
    example = spec["worked_example"]

    required = {n for n, d in inputs.items() if not d.get("optional", False)}
    given = set(example["inputs"])
    missing = required - given
    if missing:
        report.error(
            str(rel),
            f"worked example omits required input(s) {sorted(missing)}; "
            f"the generated test could not call the implementation",
        )
    extra = given - set(inputs)
    if extra:
        report.error(
            str(rel),
            f"worked example passes {sorted(extra)} which are not declared inputs",
        )

    unknown_out = set(example["expected"]) - set(outputs)
    if unknown_out:
        report.error(
            str(rel),
            f"worked example expects {sorted(unknown_out)} which are not declared outputs",
        )

    # A worked example should demonstrate the calc's headline result, not a
    # peripheral one.
    demanded = {n for n, d in outputs.items() if not d.get("optional", False)}
    if not (set(example["expected"]) & demanded):
        report.error(
            str(rel),
            "worked example does not assert any of the calc's required outputs",
        )

    # No derivation check. One used to live here, behind `if False`, warning when a
    # worked example carried no `derivation` prose - and a rule that is switched off
    # is worse than a rule that was never written, because it reads as enforcement to
    # anyone skimming. The check is deleted rather than re-enabled: `derivation` is
    # optional prose about how an expected value was obtained, the example's own
    # inputs and expected values are what the tests compare, and nothing here can
    # tell whether arithmetic written in a sentence is true, so it is a reviewer's job.

    tol = example["tolerance"]
    if tol > MAX_SENSIBLE_TOLERANCE:
        report.warn(
            str(rel),
            f"worked example tolerance {tol} is looser than {MAX_SENSIBLE_TOLERANCE}; "
            f"at that width the test would pass for an implementation that was wrong",
        )


def check_tests(report: Report, rel: Path, spec: dict[str, Any]) -> None:
    """Tests must be uniquely named, and a calculation must ship a runnable example."""
    inputs = spec["inputs"]
    outputs = spec["outputs"]
    seen: set[str] = set()
    has_worked_example = False
    active_tests = 0

    for test in spec["tests"]:
        tid = test["id"]
        if tid in seen:
            report.error(str(rel), f"duplicate test id '{tid}'")
        seen.add(tid)

        if test["type"] == "worked_example":
            has_worked_example = True

        if test["status"] == "skipped":
            # Skipping is legitimate, but only when it says why - which the schema
            # enforces by requiring `skip_reason` on a skipped test.
            continue

        active_tests += 1

        if test["type"] == "reference":
            # Every key below is required by the schema for an active reference
            # test, so in a normal run they are present. They are still read
            # defensively: if the schema is ever relaxed and this tool is not,
            # reporting a clear error beats crashing with a KeyError, since a
            # traceback from a validator reads as a bug in the spec rather than in
            # the check.
            missing = [k for k in ("inputs", "expected", "tolerance") if k not in test]
            if missing:
                report.error(
                    str(rel),
                    f"active reference test '{tid}' is missing {sorted(missing)}; "
                    f"it cannot be run or checked",
                )
                continue

            unknown_in = set(test["inputs"]) - set(inputs)
            if unknown_in:
                report.error(
                    str(rel),
                    f"test '{tid}' passes {sorted(unknown_in)} which are not declared inputs",
                )
            unknown_out = set(test["expected"]) - set(outputs)
            if unknown_out:
                report.error(
                    str(rel),
                    f"test '{tid}' expects {sorted(unknown_out)} which are not declared outputs",
                )
            if test["tolerance"] > MAX_SENSIBLE_TOLERANCE:
                report.warn(
                    str(rel),
                    f"test '{tid}' tolerance {test['tolerance']} is looser than "
                    f"{MAX_SENSIBLE_TOLERANCE}",
                )

    if not has_worked_example:
        report.error(
            str(rel),
            "no worked_example test. Every calculation ships a worked example: it is "
            "what pins a number to something a reader can retrace, it costs a dozen "
            "lines, and it is the cheapest check in the registry. **Skipping it is "
            "allowed** - a calculation with nothing checkable may say so - but the "
            "test has to exist and its skip_reason has to say why.",
        )


def check_input_flags(report: Report, rel: Path, spec: dict[str, Any]) -> None:
    """The `interval` flag must be on an input that can actually have an interval.

    `interval: true` tells the Python boundary to refuse an offset unit (`degC`,
    `degF`) so an absolute temperature cannot be passed where a difference is meant.
    That only means anything for a unit that can express a difference, and there is
    exactly one: `K`. On, say, `m` the flag would silently do nothing - a check that
    looks like validation and is not, which is the failure mode this tool exists to
    catch.

    The second rule is the one with teeth. A temperature *difference* declared
    without the flag is the hazard itself: nothing then distinguishes `Q(30, "degC")`
    from `Q(30, "delta_degC")`, `pint` applies the 273.15 offset, and the calc returns
    a plausible number wrong by a factor the caller cannot see.
    """
    for name, declaration in spec["inputs"].items():
        unit = declaration.get("unit")
        marked = bool(declaration.get("interval", False))

        if marked and unit != "K":
            report.error(
                str(rel),
                f"input '{name}' is marked interval: true but its unit is "
                f"'{unit}'. Only 'K' can express a temperature difference; on any "
                f"other unit the flag does nothing.",
            )
            continue

        # A name is not proof, so this is a heuristic and is reported as such: it
        # catches the two names a difference is actually written under in this tree
        # rather than claiming to read the English in `description`.
        looks_like_a_difference = name.lower().startswith("dt") or name.lower().startswith("delta_")
        if looks_like_a_difference and unit == "K" and not marked:
            report.warn(
                str(rel),
                f"input '{name}' is a kelvin-dimensioned quantity whose name reads "
                f"as a difference but is not marked interval: true. If it is a "
                f"difference, an absolute 'degC' will be silently offset by 273.15. "
                f"If it is an absolute temperature, ignore this.",
            )


def check_data(report: Report, rel: Path, spec: dict[str, Any]) -> None:
    """Declared data files must exist and must satisfy the fitting_list contract."""
    fittings_used: set[str] = set()

    def collect(node: dict[str, Any]) -> None:
        value = node.get("fittings")
        if isinstance(value, list):
            fittings_used.update(value)

    collect(spec["worked_example"]["inputs"])
    for test in spec["tests"]:
        if test["status"] == "active" and test["type"] == "reference":
            collect(test.get("inputs", {}))

    needs_fittings = any(d.get("type") == "fitting_list" for d in spec["inputs"].values()) or bool(
        fittings_used
    )

    data = spec.get("data", {})
    if needs_fittings and "fittings" not in data:
        report.error(
            str(rel),
            "declares a fitting_list input but no `data.fittings` file to resolve it against",
        )
        return

    for key, raw_path in data.items():
        path = ROOT / raw_path
        if not path.exists():
            report.error(str(rel), f"data.{key} points at '{raw_path}' which does not exist")
            continue

        if key == "fittings":
            try:
                rows = load_fittings_csv(path)
            except (ValueError, KeyError) as exc:
                report.error(str(rel), f"could not read fittings data '{raw_path}': {exc}")
                continue

            unknown = fittings_used - set(rows)
            if unknown:
                report.error(
                    str(rel),
                    f"references fitting id(s) {sorted(unknown)} not present in {raw_path}",
                )

            for fid, row in rows.items():
                try:
                    n_ld = float(row["n_ld"])
                except ValueError:
                    report.error(str(rel), f"{raw_path}: fitting '{fid}' has non-numeric n_ld")
                    continue
                if n_ld <= 0:
                    report.error(str(rel), f"{raw_path}: fitting '{fid}' has non-positive n_ld")

                # One definition of the provenance rules, shared with
                # `check_user_data.py` rather than restated there. A user's row and
                # the repository's row are the same kind of claim.
                check_citation(report, f"{raw_path}: fitting '{fid}'", row)

                # Coefficients must be exact decimals, not exponent notation, so
                # that Python and Rust parse identical bit patterns.
                raw = row["n_ld"]
                if "e" in raw.lower():
                    report.warn(
                        str(rel),
                        f"{raw_path}: fitting '{fid}' n_ld '{raw}' uses exponent notation; "
                        f"prefer a plain decimal so both languages parse it identically",
                    )


#: The longest a string in a spec may be, and the longer limit a worked substitution
#: gets. A spec declares a value and a short label. A substitution is the one long
#: string the format keeps, because a number nobody can follow is a number somebody
#: typed - so its limit is set by the arithmetic, not by the reading.
MAX_PROSE = 300
MAX_SUBSTITUTION = 2000

#: The fields that hold a worked substitution rather than prose.
SUBSTITUTIONS: tuple[str, ...] = ("derivation", "source")

#: The fields that are the calculation rather than a description of it, and so carry no
#: limit: the equation as Python and the same equation as LaTeX.
UNBOUNDED: tuple[str, ...] = ("equation", "latex")


def check_prose_length(report: Report, rel: Path, spec: dict[str, Any]) -> None:
    """Every string in a spec fits in a sentence or two.

    The rule is `docs/src/architecture/spec-files.md`'s: a spec declares and does not
    explain. It was written down and then violated, because nothing measured it - a
    rationale of 1,469
    characters is a page in a file whose reader wants a number, and it is the kind of
    prose that goes on arguing after the value beside it has changed.

    Measured on the value with its whitespace collapsed, which is how a reader of the
    generated page sees it.
    """
    for where, value in strings(spec):
        field = where.rpartition("/")[2]
        if field in UNBOUNDED:
            continue
        limit = MAX_SUBSTITUTION if field in SUBSTITUTIONS else MAX_PROSE
        length = len(" ".join(value.split()))
        if length > limit:
            report.error(
                str(rel),
                f"{where} is {length} characters, over the {limit} a spec field allows. "
                f"A spec declares a value and a short label; an argument belongs in a "
                f"page of the book.",
            )


def strings(value: Any, at: str = "") -> Iterator[tuple[str, str]]:
    """Every string in a document, with the path it sits at."""
    if isinstance(value, str):
        yield at, value
    elif isinstance(value, dict):
        for key, item in value.items():
            yield from strings(item, f"{at}/{key}")
    elif isinstance(value, list):
        for index, item in enumerate(value):
            yield from strings(item, f"{at}/{index}")


def lint_spec(report: Report, rel: Path, spec: dict[str, Any]) -> None:
    check_identity(report, rel, spec)
    check_input_flags(report, rel, spec)
    check_identifier_names(report, rel, spec)
    check_range_checks(report, rel, spec)
    check_worked_example(report, rel, spec)
    check_tests(report, rel, spec)
    check_data(report, rel, spec)
    check_prose_length(report, rel, spec)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--quiet", action="store_true", help="only print problems")
    parser.add_argument(
        "--spec-dir",
        type=Path,
        default=SPEC_DIR,
        help=(
            "tree of specs to check (default: specs/calcs). Exists so the test "
            "suite can lint a synthetic spec that is not part of the registry."
        ),
    )
    args = parser.parse_args()

    # Every schema under specs/schema/ in one registry, keyed by `$id`: the calc
    # schema `$ref`s `unit.schema.json` for the unit enum, and the model and case
    # schemas `$ref` the calc schema in turn. See `schema_registry`.
    import schema_registry

    schema = schema_registry.load("calc.schema.json")
    Draft202012Validator.check_schema(schema)
    validator = Draft202012Validator(schema, registry=schema_registry.registry())

    report = Report()
    parsed = check_schema(report, validator, args.spec_dir)
    for rel, spec in parsed:
        lint_spec(report, rel, spec)

    # Models, against their own schema. The model schema `$ref`s the calc schema for
    # units and quantities, so the two share one definition of what a unit is - and
    # the reference has to be resolvable locally rather than fetched.
    if MODEL_SCHEMA_PATH.exists():
        model_schema = schema_registry.load("model.schema.json")
        Draft202012Validator.check_schema(model_schema)
        model_validator = Draft202012Validator(model_schema, registry=schema_registry.registry())
        # The instances come from their own files, against their own schema: a model
        # is a type, and the machines to run are not part of it. They are attached
        # here so the case rules below run unchanged on the same shape they always did.
        case_schema = schema_registry.load("case.schema.json")
        Draft202012Validator.check_schema(case_schema)
        case_validator = Draft202012Validator(case_schema, registry=schema_registry.registry())
        instances: dict[str, list[dict[str, Any]]] = {}
        for _rel, document in check_schema(report, case_validator, CASE_DIR):
            instances.setdefault(document["model"], []).extend(document["cases"])

        for rel, spec in check_schema(report, model_validator, MODEL_DIR):
            spec["cases"] = instances.get(spec["id"], [])
            if not spec["cases"]:
                report.error(
                    str(rel),
                    f"{spec['id']} has no instances. Every model ships at least one "
                    f"file under specs/cases/, or nothing verifies it.",
                )
            check_model(report, rel, spec)

    if not args.quiet:
        print(f"spec_lint: {len(parsed)} spec(s) checked against {SCHEMA_PATH.name}")

    # The fluid tables, unconditionally: they are not reached through a spec's
    # `data:` block, because the CLI and the property provider read them and no calc
    # does.
    check_fluid_tables(report, ROOT / "data" / "fluids")

    for warning in report.warnings:
        print(f"  warning  {warning}", file=sys.stderr)
    for error in report.errors:
        print(f"  ERROR    {error}", file=sys.stderr)

    if report.errors:
        print(
            f"\nspec_lint: FAILED with {len(report.errors)} error(s), "
            f"{len(report.warnings)} warning(s)",
            file=sys.stderr,
        )
        return 1

    if not args.quiet:
        print(f"spec_lint: OK ({len(report.warnings)} warning(s))")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
