#!/usr/bin/env python3
"""Generate language-level registries from the calc specs.

Reads every spec under specs/calcs/ and emits:

  crates/azoth-<namespace>/src/spec_gen.rs   one file per namespace, as Rust statics
  python/src/azoth/_registry_gen.py          every spec, as Python data

The Rust side is split by namespace because a crate is the unit of compilation: a
calc reads its bounds from a table in its own crate, and a crate must not have to
depend on another namespace's crate to see its own spec. The types those tables are
built from are shared, and live in `azoth_core::spec`.

The Python side is one file. It all ships as one package, and `azoth.hydraulics`
and `azoth.thermal` are subpackages of a single distribution, so there is nothing
for a split to buy.

Both outputs are committed and drift-checked in CI, so the specs are genuinely
the single source of truth rather than a document that is supposed to match the
code.

Why generate Rust rather than have Rust read the TOML at run time: the spec's
contents would be available only at run time, so a malformed bound would be a
runtime error instead of a build failure. Generating a `&'static` table moves
that to compile time.

Design note on the generated range checks: they are emitted as full struct
literals rather than through the convenience constructors in azoth-core. The
constructors assume inclusive bounds, so a spec with, say, an exclusive lower
and inclusive upper bound would be silently generated wrong. Struct literals
cannot be wrong that way, and a reviewer reading the generated file sees the
exact bound semantics rather than having to recall what a constructor does.

Usage:
    python tools/gen_registry.py            # write the files
    python tools/gen_registry.py --check    # fail if they are out of date
"""

from __future__ import annotations

import argparse
import re
import subprocess
import sys
import tomllib
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent.parent
SPEC_DIR = ROOT / "specs" / "calcs"
CRATES_DIR = ROOT / "crates"
PY_OUT = ROOT / "python" / "src" / "azoth" / "_registry_gen.py"

GENERATED_BANNER = "GENERATED FILE - DO NOT EDIT BY HAND."


def namespace_of(spec: dict[str, Any]) -> str:
    """The namespace a spec belongs to, which is its id's first segment."""
    namespace: str = spec["id"].split(".")[0]
    return namespace


def crate_for(namespace: str) -> str:
    """The crate that owns a namespace's generated tables.

    Derived rather than looked up in a hand-maintained map, so that adding a
    namespace stays "add a directory under specs/calcs/" - the same rule the rest
    of the tooling follows. A namespace whose crate does not exist is caught by
    `rust_out_for` rather than silently writing tables into some default crate.
    """
    return f"azoth-{namespace.replace('_', '-')}"


def rust_out_for(namespace: str) -> Path:
    """Where a namespace's generated Rust tables live.

    Fails rather than falling back. The fallback would be writing one namespace's
    specs into another namespace's crate, which compiles, passes the drift check
    against itself, and puts a calc's bounds in a crate that calc does not depend
    on - so the file would simply never be read.
    """
    crate = crate_for(namespace)
    source_dir = CRATES_DIR / crate / "src"
    if not source_dir.is_dir():
        sys.exit(
            f"gen_registry: specs/calcs/{namespace}/ exists but crates/{crate}/src/ "
            f"does not. A namespace needs a crate to generate into; either create "
            f"crates/{crate}/, or the namespace is misspelled."
        )
    return source_dir / "spec_gen.rs"


def rustfmt(source: str) -> str:
    """Format generated Rust with rustfmt.

    The output is committed and CI runs `cargo fmt --check` over the whole
    workspace, so generated code has to be formatted like hand-written code.
    Formatting here rather than emitting carefully aligned strings means the
    generator never has to know rustfmt's rules - and if rustfmt is unavailable
    the source is still valid, just unformatted, which `cargo fmt --check` will
    report rather than something silently drifting.
    """
    try:
        proc = subprocess.run(
            ["rustfmt", "--edition", "2024", "--emit", "stdout"],
            input=source,
            capture_output=True,
            text=True,
            check=False,
        )
    except FileNotFoundError:
        print("gen_registry: rustfmt not found; emitting unformatted Rust", file=sys.stderr)
        return source
    if proc.returncode != 0 or not proc.stdout.strip():
        print(
            f"gen_registry: rustfmt failed, emitting unformatted Rust\n{proc.stderr}",
            file=sys.stderr,
        )
        return source
    return proc.stdout


def rust_str(value: str) -> str:
    """A Rust string literal for a bounded piece of prose.

    Whitespace runs are collapsed to single spaces, so the message renders
    identically on both sides and the generated file stays readable.
    """
    flat = " ".join(value.split())
    escaped = flat.replace("\\", "\\\\").replace('"', '\\"')
    return f'"{escaped}"'


def rust_f64(value: float) -> str:
    """A Rust `f64` literal.

    A spec may write a whole number (`0`), whose Python repr produces an integer
    literal in Rust that will not compile where an `f64` is expected, so the type
    is forced first. `repr` round-trips exactly, so no precision is lost.
    """
    return repr(float(value))


def rust_opt_f64(value: float | None) -> str:
    return "None" if value is None else f"Some({rust_f64(value)})"


def rust_slice_str(items: list[str]) -> str:
    inner = ", ".join(rust_str(i) for i in items)
    return f"&[{inner}]"


def camel_variant(screaming_snake: str) -> str:
    """`TRANSITIONAL_FLOW` -> `TransitionalFlow`, matching the Rust enum."""
    return "".join(part.capitalize() for part in screaming_snake.split("_"))


#: Expected keys this generator does not carry and that a model's own test reads from
#: the raw spec instead.
#:
#: Named rather than allowed silently: a case that asserts something *nothing* reads is
#: worse than one that asserts nothing, because it reads as verified. An entry here is a
#: claim that somewhere does read it, and the test named beside each is where.
SELF_ASSERTED_EXPECTATIONS: frozenset[str] = frozenset(
    {
        # `stability_test`'s per-trial phase compositions, which are a matrix and have no
        # field in `TestCase`. `crates/azoth-eos/tests/stability_test.rs` asserts their
        # structure and `python/tests/models/test_stability_test.py` reads this same spec
        # and compares both backends against them.
        "w",
        # `reactive_phase_equilibrium`'s element matrix, for the same reason and with the
        # same shape. `crates/azoth-reactions/tests/reactive_phase_equilibrium.rs` reads it
        # through `TestCase::matrix` and compares it row by row against the matrix the
        # operation built, so the specification is what constrains it either way.
        "a_matrix",
        # `splitter`'s per-outlet compositions, which are a matrix because a `many` port
        # writes its `z` one row per outlet. The Python runner compares the values, as it
        # does for `a_matrix`; the two kernels are compared with each other by
        # `test_cross_impl`, which reads the field off both results, so a shape one side
        # got wrong is caught there rather than here.
        "products_z",
    }
)


def collect_numbers(mapping: dict[str, Any]) -> list[tuple[str, float]]:
    """Scalar numbers, **excluding flags**.

    `isinstance(True, int)` is True in Python, so a boolean collected as a number arrives
    in Rust as `1.0` and a reader of the generated table cannot tell a flag from a
    quantity. A case's `associating` is the first of them; it crosses as a boolean.
    """
    return [
        (k, v)
        for k, v in mapping.items()
        if isinstance(v, (int, float)) and not isinstance(v, bool)
    ]


def collect_flags(mapping: dict[str, Any]) -> list[tuple[str, bool]]:
    """Scalar booleans: a switch whose value is yes or no."""
    return [(k, v) for k, v in mapping.items() if isinstance(v, bool)]


def collect_lists(mapping: dict[str, Any]) -> list[tuple[str, list[str]]]:
    """Lists of strings: fitting ids, which are identifiers rather than numbers."""
    return [
        (k, v)
        for k, v in mapping.items()
        if isinstance(v, list) and v and all(isinstance(x, str) for x in v)
    ]


def collect_strings(mapping: dict[str, Any]) -> list[tuple[str, str]]:
    """Scalar strings: an enum or categorical input, whose value is a name."""
    return [(k, v) for k, v in mapping.items() if isinstance(v, str)]


def collect_vectors(mapping: dict[str, Any]) -> list[tuple[str, list[float]]]:
    """Lists of numbers: a composition, or one constant per component."""
    return [
        (k, v)
        for k, v in mapping.items()
        if isinstance(v, list) and v and all(isinstance(x, (int, float)) for x in v)
    ]


def collect_matrices(mapping: dict[str, Any]) -> list[tuple[str, list[float]]]:
    """Lists of lists of numbers, flattened row-major.

    The dimension is not emitted: every matrix this registry carries is square and
    the same size as the vectors beside it, so a consumer reshapes against their
    length rather than against a second copy of the same number.
    """
    out: list[tuple[str, list[float]]] = []
    for key, value in mapping.items():
        if not isinstance(value, list) or not value:
            continue
        if not all(isinstance(row, list) for row in value):
            continue
        out.append((key, [float(x) for row in value for x in row]))
    return out


def emit_range_check(check: dict[str, Any], spec_id: str) -> str:
    severity = "Severity::Error" if check["severity"] == "error" else "Severity::Warning"
    band = "Band::Inside" if check.get("when", "outside") == "inside" else "Band::Outside"
    # Warning codes are emitted by their Rust variant name. The schema restricts
    # `code` to the known set, and spec_lint checks the same strings against the
    # Python enum, so a typo here cannot reach the generated file.
    code = "WarningCode::" + camel_variant(check.get("code", "OUT_OF_VALID_RANGE"))
    # A bound's rationale is optional and is not checked here. It was required, on
    # the argument that a bound with no reason is a bound nobody dares change - which
    # is true, and is not this tool's business. A bound with no stated reason is a
    # bound to ask about in review.
    rationale = check.get("rationale", "")
    # `enum` is a bound kind the schema permits and neither implementation can
    # evaluate: a range check resolves a quantity to a float, and an enum bound
    # compares strings. Emitting it would produce a check that silently never fires,
    # so the generator refuses instead. `spec_lint` refuses it earlier, with the
    # same reasoning; this is the second gate, because a spec reaching the generator
    # is a spec that got past the linter.
    if check.get("enum") is not None:
        raise SystemExit(
            f"{spec_id}: range check on '{check['quantity']}' uses `enum`, which "
            f"neither implementation can evaluate. A range check resolves its "
            f"quantity to a float; an enum bound compares strings. Emitting it would "
            f"produce a check that never fires."
        )
    return (
        "            check: RangeCheck {\n"
        f"                quantity: {rust_str(check['quantity'])},\n"
        f"                min: {rust_opt_f64(check.get('min'))},\n"
        f"                min_inclusive: {str(check.get('min_inclusive', True)).lower()},\n"
        f"                max: {rust_opt_f64(check.get('max'))},\n"
        f"                max_inclusive: {str(check.get('max_inclusive', True)).lower()},\n"
        f"                equals: {rust_opt_f64(check.get('equals'))},\n"
        f"                band: {band},\n"
        f"                severity: {severity},\n"
        f"                code: {code},\n"
        f"                rationale: {rust_str(rationale)},\n"
        "            },\n"
    )


def emit_test_case(
    test_id: str,
    kind: str,
    status: str,
    tolerance: float,
    inputs: dict[str, Any],
    expected: dict[str, Any],
    indent: str,
    property_name: str | None = None,
    skip_reason: str | None = None,
) -> str:
    numbers = collect_numbers(inputs)
    lists = collect_lists(inputs)
    strings = collect_strings(inputs)
    flags = collect_flags(inputs)
    vectors = collect_vectors(inputs)
    matrices = collect_matrices(inputs)
    expected_numbers = collect_numbers(expected)
    expected_vectors = collect_vectors(expected)
    expected_strings = collect_strings(expected)
    # A case's expectation this generator cannot carry is refused rather than dropped:
    # an assertion nothing reads is worse than no assertion, because it reads as one.
    carried = (
        {k for k, _ in expected_numbers}
        | {k for k, _ in expected_vectors}
        | {k for k, _ in expected_strings}
    )
    if unreachable := sorted(set(expected) - carried - SELF_ASSERTED_EXPECTATIONS):
        raise ValueError(
            f"case '{test_id}' asserts {unreachable}, which this generator cannot "
            f"carry. Declare it as a number, a vector or a string, add it to "
            f"SELF_ASSERTED_EXPECTATIONS with the test that reads it, or assert it in "
            f"the model's own test"
        )

    def pairs(items: list[tuple[str, float]]) -> str:
        if not items:
            return "&[]"
        inner = ", ".join(f"({rust_str(k)}, {rust_f64(v)})" for k, v in items)
        return f"&[{inner}]"

    def flag_pairs(items: list[tuple[str, bool]]) -> str:
        if not items:
            return "&[]"
        inner = ", ".join(f"({rust_str(k)}, {'true' if v else 'false'})" for k, v in items)
        return f"&[{inner}]"

    def list_pairs(items: list[tuple[str, list[str]]]) -> str:
        if not items:
            return "&[]"
        inner = ", ".join(f"({rust_str(k)}, {rust_slice_str(v)})" for k, v in items)
        return f"&[{inner}]"

    def string_pairs(items: list[tuple[str, str]]) -> str:
        if not items:
            return "&[]"
        inner = ", ".join(f"({rust_str(k)}, {rust_str(v)})" for k, v in items)
        return f"&[{inner}]"

    def number_slice_pairs(items: list[tuple[str, list[float]]]) -> str:
        if not items:
            return "&[]"
        inner = ", ".join(
            f"({rust_str(k)}, &[{', '.join(rust_f64(x) for x in v)}])" for k, v in items
        )
        return f"&[{inner}]"

    prop = f"Some({rust_str(property_name)})" if property_name else "None"
    skip = f"Some({rust_str(skip_reason)})" if skip_reason else "None"

    return (
        f"{indent}TestCase {{\n"
        f"{indent}    id: {rust_str(test_id)},\n"
        f"{indent}    kind: {rust_str(kind)},\n"
        f"{indent}    property: {prop},\n"
        f"{indent}    status: {rust_str(status)},\n"
        f"{indent}    skip_reason: {skip},\n"
        f"{indent}    tolerance: {rust_f64(tolerance)},\n"
        f"{indent}    numbers: {pairs(numbers)},\n"
        f"{indent}    flags: {flag_pairs(flags)},\n"
        f"{indent}    lists: {list_pairs(lists)},\n"
        f"{indent}    strings: {string_pairs(strings)},\n"
        f"{indent}    vectors: {number_slice_pairs(vectors)},\n"
        f"{indent}    matrices: {number_slice_pairs(matrices)},\n"
        f"{indent}    expected: {pairs(expected_numbers)},\n"
        f"{indent}    expected_vectors: {number_slice_pairs(expected_vectors)},\n"
        f"{indent}    expected_strings: {string_pairs(expected_strings)},\n"
        f"{indent}}},\n"
    )


def emit_rust(specs: list[dict[str, Any]], source_files: list[str], namespace: str) -> str:
    header: list[str] = []
    body: list[str] = []
    header.append(
        f"""//! {GENERATED_BANNER}
//!
//! Generated by `tools/gen_registry.py` from:
//!
"""
    )
    for src in source_files:
        header.append(f"//!   - {src}\n")
    header.append(
        f"""//!
//! Tables for the `{namespace}` namespace. Every namespace has its own generated
//! file, because a crate is the unit of compilation and a calculation must be able
//! to read its own bounds without depending on another namespace's crate. The
//! types these tables are built from are shared - see `azoth_core::spec`.
//!
//! These tables are what make the specs authoritative at runtime rather than
//! merely descriptive. Each calc reads its own range checks from here, so a
//! bound changed in a spec file changes the code's behaviour with no second
//! edit - and `cargo test` fails if the two ever disagree.
"""
    )
    # Per-calc statics.
    for spec in specs:
        ident = spec["id"].split(".")[-1].upper()
        inputs = spec["inputs"]

        checks: list[str] = []
        for check in spec["valid_range"]:
            on_input = check["quantity"] in inputs
            checks.append("        SpecCheck {\n")
            checks.append(f"            on_input: {str(on_input).lower()},\n")
            checks.append(emit_range_check(check, spec["id"]))
            checks.append("        },\n")

        solver = spec.get("solver")
        solver_str = "None"
        if solver:
            # `initial_guess` is optional in the spec: `fixed_point` iterates from a
            # declared start and the schema requires one, while `cubic_roots` forms
            # its roots analytically and has no starting point to declare. Emitting
            # `None` rather than a placeholder keeps that absence visible in the
            # generated table instead of inventing a value nothing reads.
            guess = solver.get("initial_guess")
            guess_str = "None" if guess is None else f"Some({rust_f64(guess)})"
            solver_str = (
                "Some(SolverSpec {\n"
                f"            kind: {rust_str(solver['kind'])},\n"
                f"            tolerance: {rust_f64(solver['tolerance'])},\n"
                f"            max_iterations: {solver['max_iterations']},\n"
                f"            initial_guess: {guess_str},\n"
                f"            convergence: {rust_str(solver['convergence'])},\n"
                "        })"
            )

        example = spec["worked_example"]
        body.append(f"/// Registry entry for `{spec['id']}`.\n")
        body.append(f"static {ident}_CHECKS: &[SpecCheck] = &[\n{''.join(checks)}];\n\n")

        others = [t for t in spec["tests"] if t["type"] != "worked_example"]
        body.append(f"static {ident}_TESTS: &[TestCase] = &[\n")
        for test in others:
            body.append(
                emit_test_case(
                    test_id=test["id"],
                    kind=test["type"],
                    status=test["status"],
                    tolerance=test.get("tolerance", example["tolerance"]),
                    inputs=test.get("inputs", {}),
                    expected=test.get("expected", {}),
                    indent="    ",
                    property_name=test.get("property"),
                    skip_reason=test.get("skip_reason"),
                )
            )
        body.append("];\n\n")

        body.append(
            f"""/// Registered spec for `{spec["id"]}`.
///
/// Public and addressable directly, so a calc can hold `&{ident}_SPEC` with no
/// lookup and no failure path. A calc whose spec is missing is a build-time
/// invariant, not a runtime condition, and this shape makes it unrepresentable
/// rather than something to handle.
pub static {ident}_SPEC: CalcSpec = CalcSpec {{
    id: {rust_str(spec["id"])},
    checks: {ident}_CHECKS,
    solver: {solver_str},
    worked_example: {
                emit_test_case(
                    # The spec's own worked-example test id, not a synthesized one, so the two
                    # languages report the same name for the same case.
                    test_id=next(
                        (t["id"] for t in spec["tests"] if t["type"] == "worked_example"),
                        spec["id"].split(".")[-1] + "_worked_example",
                    ),
                    kind="worked_example",
                    status=next(
                        (t["status"] for t in spec["tests"] if t["type"] == "worked_example"),
                        "active",
                    ),
                    tolerance=example["tolerance"],
                    inputs=example["inputs"],
                    expected=example["expected"],
                    indent="        ",
                    # strip() leaves the trailing comma that emit_test_case adds for list
                    # elements; the template supplies its own, and two in a row will not
                    # compile.
                )
                .strip()
                .removesuffix(",")
            },
    tests: {ident}_TESTS,
}};

"""
        )

    names = ", ".join(f"&{s['id'].split('.')[-1].upper()}_SPEC" for s in specs)
    body.append(
        f"""/// Every calculation in the registry, sorted by id.
static ALL_SPECS: &[&CalcSpec] = &[{names}];

/// All specs, in a stable order.
#[must_use]
pub fn specs() -> &'static [&'static CalcSpec] {{
    ALL_SPECS
}}

/// Look up a spec by id.
#[must_use]
pub fn spec(id: &str) -> Option<&'static CalcSpec> {{
    ALL_SPECS.iter().copied().find(|s| s.id == id)
}}
"""
    )
    body_text = "".join(body)

    # Import exactly the names the body uses. A namespace whose calculations are all
    # explicit has no `SolverSpec` in its tables, and importing it anyway is an unused
    # import - which `-D warnings` turns into a build failure for that crate. Deriving
    # the list rather than maintaining one per namespace keeps the generated file
    # honest about what it actually references.
    candidates = [
        "Band",
        "CalcSpec",
        "RangeCheck",
        "Severity",
        "SpecCheck",
        "SolverSpec",
        "TestCase",
        "WarningCode",
    ]
    needed = [name for name in candidates if re.search(rf"\b{name}\b", body_text)]
    imports = "use azoth_core::{\n" + "".join(f"    {name},\n" for name in needed) + "};\n\n"

    return "".join(header) + imports + body_text


def emit_python(specs: list[dict[str, Any]], source_files: list[str]) -> str:
    out: list[str] = []
    out.append(
        f'''"""{GENERATED_BANNER}

Generated by `tools/gen_registry.py` from:
'''
    )
    for src in source_files:
        out.append(f"  - {src}\n")
    out.append(
        '''"""

from __future__ import annotations

from typing import Any, Final

#: Every calc spec, as plain data, in a stable order.
CALCS: Final[tuple[dict[str, Any], ...]] = (
'''
    )
    for spec in specs:
        out.append("    {\n")
        out.append(f'        "id": {spec["id"]!r},\n')
        out.append(f'        "name": {spec["name"]!r},\n')
        out.append(f'        "equation": {spec["equation"]!r},\n')
        out.append(f'        "source": {spec["source"]!r},\n')
        out.append(f'        "inputs": { {k: dict(v) for k, v in spec["inputs"].items()}!r},\n')
        out.append(f'        "outputs": { {k: dict(v) for k, v in spec["outputs"].items()}!r},\n')
        out.append(f'        "valid_range": {spec["valid_range"]!r},\n')
        out.append(f'        "assumptions": {spec["assumptions"]!r},\n')
        out.append(f'        "references": {spec["references"]!r},\n')
        out.append(f'        "solver": {spec.get("solver")!r},\n')
        out.append(f'        "data": {spec.get("data")!r},\n')
        out.append(f'        "worked_example": {spec["worked_example"]!r},\n')
        out.append(f'        "implementations": {spec["implementations"]!r},\n')
        out.append(f'        "tests": {spec["tests"]!r},\n')
        out.append("    },\n")
    out.append(
        """)


#: Lookup by calc id.
BY_ID: Final[dict[str, dict[str, Any]]] = {c["id"]: c for c in CALCS}


def spec(calc_id: str) -> dict[str, Any]:
    \"\"\"Fetch a spec by id, raising a clear error if it is not a calc we ship.\"\"\"
    try:
        return BY_ID[calc_id]
    except KeyError:
        raise KeyError(
            f"unknown calc {calc_id!r}; known ids: {sorted(BY_ID)}"
        ) from None
"""
    )
    return "".join(out)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--check",
        action="store_true",
        help="do not write; fail if the generated files are out of date",
    )
    args = parser.parse_args()

    paths = sorted(SPEC_DIR.rglob("*.toml"))
    if not paths:
        sys.exit(f"gen_registry: no specs found under {SPEC_DIR}")

    # Paired and sorted together. The specs used to be sorted by id while the
    # source paths kept directory order - which happened to agree, because a
    # namespace's directory name is also its id prefix, so the two orders always
    # matched. The lists are consumed in parallel and are now split per namespace,
    # so an accidental disagreement would attribute one namespace's spec files to
    # another's tables and be very hard to see.
    loaded: list[tuple[str, dict[str, Any]]] = [
        (str(p.relative_to(ROOT)), tomllib.loads(p.read_text(encoding="utf-8"))) for p in paths
    ]
    loaded.sort(key=lambda pair: pair[1]["id"])

    # Guard against two specs claiming the same id: the generated tables key on
    # it, and a duplicate would silently drop one calc from the registry.
    ids = [spec["id"] for _, spec in loaded]
    duplicates = {i for i in ids if ids.count(i) > 1}
    if duplicates:
        sys.exit(f"gen_registry: duplicate calc id(s): {sorted(duplicates)}")

    by_namespace: dict[str, list[tuple[str, dict[str, Any]]]] = {}
    for source_file, spec in loaded:
        by_namespace.setdefault(namespace_of(spec), []).append((source_file, spec))

    outputs: dict[Path, str] = {
        PY_OUT: emit_python([spec for _, spec in loaded], [src for src, _ in loaded])
    }
    for namespace in sorted(by_namespace):
        entries = by_namespace[namespace]
        outputs[rust_out_for(namespace)] = rustfmt(
            emit_rust([spec for _, spec in entries], [src for src, _ in entries], namespace)
        )

    stale: list[Path] = []
    for path, content in outputs.items():
        path.parent.mkdir(parents=True, exist_ok=True)
        if args.check:
            current = path.read_text(encoding="utf-8") if path.exists() else ""
            if current != content:
                stale.append(path)
        else:
            path.write_text(content, encoding="utf-8")

    if args.check:
        if stale:
            for path in stale:
                print(f"gen_registry: {path.relative_to(ROOT)} is out of date", file=sys.stderr)
            print("Run `python tools/gen_registry.py` to regenerate.", file=sys.stderr)
            return 1
        print(f"gen_registry: {len(loaded)} spec(s), generated files up to date")
        return 0

    for path in outputs:
        print(f"gen_registry: wrote {path.relative_to(ROOT)}")
    print(f"gen_registry: {len(loaded)} spec(s)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
