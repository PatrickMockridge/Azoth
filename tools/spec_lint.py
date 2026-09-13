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

Exit status is non-zero if any error-level problem is found. Warnings do not fail
the run, but they are printed, because a warning nobody reads is the same as no
check at all.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

try:
    import yaml
except ImportError:  # pragma: no cover
    sys.exit("spec_lint requires PyYAML: pip install pyyaml")

try:
    from jsonschema import Draft202012Validator
except ImportError:  # pragma: no cover
    sys.exit("spec_lint requires jsonschema: pip install jsonschema")


ROOT = Path(__file__).resolve().parent.parent
SPEC_DIR = ROOT / "specs" / "calcs"
SCHEMA_PATH = ROOT / "specs" / "schema" / "calc.schema.json"

# Quantities that name a state rather than a number. These have no units and are
# not function parameters, so they are allowed in range checks without appearing
# in `inputs`.
STATE_KEYWORDS = {"phase"}

# The tolerance a worked example may declare before it stops being a meaningful
# check. Relative tolerance, so 1e-3 means "agrees to 0.1%".
MAX_SENSIBLE_TOLERANCE = 0.05

# Permitted provenance states for a data row.
VALID_VERIFY_STATUS = {"verified", "unverified", "estimated_dummy"}


@dataclass
class Report:
    """Collected problems, in the order they were found."""

    errors: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    # Placeholder data rows, keyed by file, counted so the run can shout about
    # them. Dummy coefficients cannot be caught by a test - there is nothing
    # correct to compare against - so the only defence is making them loud.
    estimated_rows: dict[str, int] = field(default_factory=dict)

    def error(self, spec: str, msg: str) -> None:
        self.errors.append(f"{spec}: {msg}")

    def warn(self, spec: str, msg: str) -> None:
        self.warnings.append(f"{spec}: {msg}")


def load_fittings_csv(path: Path) -> dict[str, dict[str, str]]:
    """Read the fittings registry, skipping the leading comment block.

    The comment block is not decoration - it carries the copyright statement and
    the meaning of verify_status - so the parser has to tolerate it rather than
    the file having to be machine-only.
    """
    rows: dict[str, dict[str, str]] = {}
    with path.open(newline="", encoding="utf-8") as handle:
        content = [line for line in handle if not line.lstrip().startswith("#")]
    reader = csv.DictReader(content)
    if reader.fieldnames is None:
        raise ValueError(f"{path} has no header row")
    for row in reader:
        rows[row["fitting_id"]] = dict(row)
    return rows


def check_schema(report: Report, validator: Draft202012Validator) -> list[tuple[Path, dict[str, Any]]]:
    """Validate each spec against the JSON Schema. Returns the ones that parsed."""
    parsed: list[tuple[Path, dict[str, Any]]] = []
    specs = sorted(SPEC_DIR.rglob("*.yaml"))
    if not specs:
        report.error("specs", f"no spec files found under {SPEC_DIR}")
        return parsed

    for path in specs:
        rel = path.relative_to(ROOT)
        try:
            spec = yaml.safe_load(path.read_text(encoding="utf-8"))
        except yaml.YAMLError as exc:
            report.error(str(rel), f"YAML does not parse: {exc}")
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
    expected_python = f"chemeng.{namespace}.{func}"
    if spec["implementations"]["python"] != expected_python:
        report.error(
            str(rel),
            f"implementations.python is '{spec['implementations']['python']}' "
            f"but the id implies '{expected_python}'",
        )
    expected_rust = f"chemeng_{namespace}::{func}"
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
    import re

    legal = re.compile(r"^[A-Za-z_][A-Za-z0-9_]*$")
    snake = re.compile(r"^[a-z][a-z0-9_]*$")
    keywords = set(keyword.kwlist) | set(keyword.softkwlist) | {
        # Rust keywords that would break a parameter or field name.
        "as", "break", "const", "continue", "crate", "else", "enum", "extern",
        "false", "fn", "for", "if", "impl", "in", "let", "loop", "match", "mod",
        "move", "mut", "pub", "ref", "return", "self", "static", "struct",
        "super", "trait", "true", "type", "unsafe", "use", "where", "while",
        "async", "await", "dyn", "abstract", "become", "box", "do", "final",
        "macro", "override", "priv", "try", "typeof", "unsized", "virtual", "yield",
    }

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

    if not example.get("derivation"):
        report.warn(
            str(rel),
            "worked example has no `derivation`; a reviewer cannot retrace where the "
            "expected value came from",
        )

    tol = example["tolerance"]
    if tol > MAX_SENSIBLE_TOLERANCE:
        report.warn(
            str(rel),
            f"worked example tolerance {tol} is looser than {MAX_SENSIBLE_TOLERANCE}; "
            f"at that width the test would pass for an implementation that was wrong",
        )


def check_tests(report: Report, rel: Path, spec: dict[str, Any]) -> None:
    """Tests must be complete, uniquely named, and consistent with verification status."""
    inputs = spec["inputs"]
    outputs = spec["outputs"]
    seen: set[str] = set()
    active_worked_example = False

    for test in spec["tests"]:
        tid = test["id"]
        if tid in seen:
            report.error(str(rel), f"duplicate test id '{tid}'")
        seen.add(tid)

        if test["status"] == "skipped":
            # Skipping is legitimate, but only when it says why. This is the rule
            # that keeps "TODO: source needed" from becoming a silent omission.
            continue

        if test["type"] == "worked_example":
            active_worked_example = True

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

    if not active_worked_example:
        report.error(
            str(rel),
            "no active worked_example test. Every calc must ship a runnable worked "
            "example; skipping it requires the calc to be marked source_needed and "
            "the reason recorded in the test's skip_reason.",
        )

    status = spec["verification"]["status"]
    if status == "source_needed":
        example_tests = [t for t in spec["tests"] if t["type"] == "worked_example"]
        if any(t["status"] == "active" for t in example_tests):
            report.error(
                str(rel),
                "verification.status is 'source_needed' but a worked_example test is "
                "active. An unverified calc must not claim a passing example - either "
                "find a source or skip the test.",
            )
    elif status == "unverified" and not spec["verification"].get("notes"):
        report.warn(
            str(rel),
            "verification.status is 'unverified' with no notes explaining what is "
            "unconfirmed",
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

    needs_fittings = any(
        d.get("type") == "fitting_list" for d in spec["inputs"].values()
    ) or bool(fittings_used)

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
                    report.error(
                        str(rel), f"{raw_path}: fitting '{fid}' has non-positive n_ld"
                    )

                status = row.get("verify_status")
                if status not in VALID_VERIFY_STATUS:
                    report.error(
                        str(rel),
                        f"{raw_path}: fitting '{fid}' has invalid verify_status '{status}'; "
                        f"expected one of {sorted(VALID_VERIFY_STATUS)}",
                    )
                elif status == "estimated_dummy":
                    report.estimated_rows[raw_path] = (
                        report.estimated_rows.get(raw_path, 0) + 1
                    )
                    # A row cannot simultaneously be a placeholder and verified.
                    # This is the copy-paste that would silently promote dummy
                    # data to trusted data.
                    citation = row.get("citation", "")
                    if "DUMMY" not in citation.upper():
                        report.error(
                            str(rel),
                            f"{raw_path}: fitting '{fid}' is marked estimated_dummy but its "
                            f"citation does not say so. The marker and the citation must agree.",
                        )
                elif status == "verified" and "DUMMY" in row.get("citation", "").upper():
                    report.error(
                        str(rel),
                        f"{raw_path}: fitting '{fid}' is marked verified but its citation "
                        f"still says it is a dummy value",
                    )
                # Coefficients must be exact decimals, not exponent notation, so
                # that Python and Rust parse identical bit patterns.
                raw = row["n_ld"]
                if "e" in raw.lower():
                    report.warn(
                        str(rel),
                        f"{raw_path}: fitting '{fid}' n_ld '{raw}' uses exponent notation; "
                        f"prefer a plain decimal so both languages parse it identically",
                    )


def check_numeric_literals(report: Report, rel: Path, spec: dict[str, Any]) -> None:
    """Catch numbers that YAML silently turned into strings.

    PyYAML implements YAML 1.1, whose float resolver requires an explicit sign on
    the exponent: `1.0e-3` is a float, but `1.0e12` is a STRING. That is a genuine
    trap, because the spec still looks correct and the schema violation surfaces
    far away from the typo. Worse, if the value lands somewhere unconstrained it
    would flow into the generated Rust registry as a quoted string.

    Only the subtrees where numbers are genuinely expected are checked, so quoted
    years and edition strings elsewhere are left alone.
    """
    targets: list[tuple[str, Any]] = [
        ("worked_example/tolerance", spec["worked_example"].get("tolerance")),
    ]
    for field_name in ("inputs", "expected"):
        for key, value in spec["worked_example"][field_name].items():
            targets.append((f"worked_example/{field_name}/{key}", value))

    for index, check in enumerate(spec["valid_range"]):
        for bound in ("min", "max"):
            if bound in check:
                targets.append((f"valid_range/{index}/{bound}", check[bound]))

    for key, value in spec.get("solver", {}).items():
        if key != "kind" and key != "convergence":
            targets.append((f"solver/{key}", value))

    for index, test in enumerate(spec["tests"]):
        if "tolerance" in test:
            targets.append((f"tests/{index}/tolerance", test["tolerance"]))
        for field_name in ("inputs", "expected"):
            for key, value in test.get(field_name, {}).items():
                targets.append((f"tests/{index}/{field_name}/{key}", value))

    for where, value in targets:
        if isinstance(value, str):
            report.error(
                str(rel),
                f"{where} is the string {value!r}, not a number. YAML 1.1 requires a "
                f"signed exponent (write 1.0e+12, not 1.0e12) - an unsigned exponent "
                f"silently parses as a string.",
            )


def lint_spec(report: Report, rel: Path, spec: dict[str, Any]) -> None:
    check_identity(report, rel, spec)
    check_identifier_names(report, rel, spec)
    check_numeric_literals(report, rel, spec)
    check_range_checks(report, rel, spec)
    check_worked_example(report, rel, spec)
    check_tests(report, rel, spec)
    check_data(report, rel, spec)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--quiet", action="store_true", help="only print problems")
    args = parser.parse_args()

    schema = json.loads(SCHEMA_PATH.read_text(encoding="utf-8"))
    Draft202012Validator.check_schema(schema)
    validator = Draft202012Validator(schema)

    report = Report()
    parsed = check_schema(report, validator)
    for rel, spec in parsed:
        lint_spec(report, rel, spec)

    if not args.quiet:
        print(f"spec_lint: {len(parsed)} spec(s) checked against {SCHEMA_PATH.name}")

    # Placeholder data is the one defect this tool cannot fail on, because there
    # is no correct value to compare against. So it is printed loudly instead,
    # on every run, where it cannot be missed.
    if report.estimated_rows:
        total = sum(report.estimated_rows.values())
        banner = "!" * 78
        print(f"\n{banner}", file=sys.stderr)
        print(
            f"!! NOT ENGINEERING DATA: {total} coefficient row(s) are ESTIMATED DUMMY VALUES.",
            file=sys.stderr,
        )
        for path, count in sorted(report.estimated_rows.items()):
            print(f"!!   {path}: {count} row(s)", file=sys.stderr)
        print(
            "!! These are placeholders for software testing. Do NOT use them to size\n"
            "!! equipment. Replace with values from the primary standard and set\n"
            "!! verify_status=verified before any real use.",
            file=sys.stderr,
        )
        print(f"{banner}\n", file=sys.stderr)

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
