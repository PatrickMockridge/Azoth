#!/usr/bin/env python3
"""Generate the calculation catalogue a sheet editor reads.

`specs/calcs/` is the authoritative declaration of every standalone calculation, but the Rust
registry carries only the bounds and test cases. A form needs the name, the inputs with their
kinds and units, the outputs and a worked example - all of which are in the spec and nowhere in
the generated tables. This tool projects them into `crates/azoth-sheets/src/catalogue_gen.rs`.

    python tools/gen_sheets.py            # write
    python tools/gen_sheets.py --check    # fail if the tree is not what this emits

**The catalogue is a UI projection, not a second spec.** Two editorial rules are applied here and
stated rather than hidden: the four calculations with table-shaped inputs or outputs are carried
with `available = false` and a refusal, and prose that names the upstream implementation is
dropped. A summary cannot be dropped - it is the sheet's first line - so an empty or upstream
summary stops the generator rather than reaching the browser blank. The library, CLI and docs keep
the full record; a front end that showed the upstream Java class path was the defect that made
this projection sanitized.
"""

from __future__ import annotations

import json
import re
import subprocess
import sys
import tomllib
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import gen_python_transport as gpt
import gen_python_wrappers as gpw
import rust_index

# The generated Rust is one long template line per field; wrapping it would make the emitted
# code harder to read than the line it replaces.
# ruff: noqa: E501

ROOT = Path(__file__).resolve().parent.parent
SPECS = ROOT / "specs" / "calcs"
OUT = ROOT / "crates" / "azoth-sheets" / "src" / "catalogue_gen.rs"
DISPATCH_OUT = ROOT / "crates" / "azoth-sheets" / "src" / "dispatch_gen.rs"

#: The input and output kinds a sheet field can draw. A calculation that needs anything else is
#: listed and refused rather than silently given a field that cannot express it.
SIMPLE = {"quantity", "enum", "boolean", "string"}

#: The upstream implementation's name, in any case. A UI string that carries it is not shipped.
UPSTREAM = "neqsim"

REFUSAL = (
    "this calculation takes a table of rows as input, which the sheet editor does not yet provide"
)


def documents() -> list[tuple[Path, dict]]:
    return [
        (path, tomllib.loads(path.read_text(encoding="utf-8")))
        for path in sorted(SPECS.rglob("*.toml"))
    ]


def kind(declaration: dict) -> str:
    return str(declaration.get("type", "quantity"))


def complex_inputs(spec: dict) -> bool:
    return any(kind(declaration) not in SIMPLE for declaration in spec.get("inputs", {}).values())


def complex_outputs(spec: dict) -> bool:
    return any(
        kind(declaration) not in SIMPLE or "fields" in declaration
        for declaration in spec.get("outputs", {}).values()
    )


def clean(text: str) -> str:
    """A string the UI may show, or empty where it cannot be shown without naming the upstream."""
    return "" if UPSTREAM in text.lower() else text


def clean_list(values: list[str]) -> list[str]:
    return [value for value in values if UPSTREAM not in value.lower()]


def humanise(calc_id: str) -> str:
    """A name for a spec whose own name names the upstream implementation."""
    leaf = calc_id.split(".", 1)[-1].replace("_", " ")
    return leaf[:1].upper() + leaf[1:]


def name_of(spec: dict) -> str:
    name = str(spec["name"])
    return humanise(spec["id"]) if UPSTREAM in name.lower() else name


def summary_of(spec: dict) -> str:
    """The authored purpose line, which the UI must be able to show verbatim.

    **Not a fallback and not a clean-up.** A summary that is empty or names the upstream
    implementation is a spec defect: the browser would show a blank or drop the line, and the
    sheet would be left without the one explanation the reader needs before filling a field.
    """
    summary = " ".join(str(spec.get("summary", "")).split())
    if summary == "":
        raise SystemExit(f"gen_sheets: {spec['id']} has no summary") from None
    if UPSTREAM in summary.lower():
        raise SystemExit(
            f"gen_sheets: {spec['id']}'s summary names the upstream implementation"
        ) from None
    return summary


def dimension_of(unit: object) -> str | None:
    if unit is None:
        return None
    entry = rust_index._units().get(str(unit))
    return None if entry is None else str(entry["dimension"])


def ranges(spec: dict, name: str) -> list[dict]:
    out = []
    for check in spec.get("valid_range", []):
        if check.get("quantity") != name or "computed_from" in check:
            continue
        equals = check.get("equals")
        numeric_equals = (
            float(equals)
            if isinstance(equals, (int, float)) and not isinstance(equals, bool)
            else None
        )
        out.append(
            {
                "min": None if check.get("min") is None else float(check["min"]),
                "min_inclusive": bool(check.get("min_inclusive", True)),
                "max": None if check.get("max") is None else float(check["max"]),
                "max_inclusive": bool(check.get("max_inclusive", True)),
                "equals": numeric_equals,
                "band": str(check.get("when", "outside")),
                "severity": str(check["severity"]),
                "code": str(check.get("code", "OUT_OF_VALID_RANGE")),
                "rationale": clean(str(check.get("rationale", ""))),
            }
        )
    return out


def inputs(spec: dict) -> list[dict]:
    out = []
    for name, declaration in spec.get("inputs", {}).items():
        unit = declaration.get("unit")
        default = declaration.get("default")
        out.append(
            {
                "name": name,
                "kind": kind(declaration),
                "required": not bool(declaration.get("optional", False)),
                "unit": None if unit is None else str(unit),
                "dimension": dimension_of(unit),
                "values": [str(value) for value in declaration.get("values", [])],
                "description": clean(str(declaration.get("description", ""))),
                "range": ranges(spec, name),
                "default": (
                    None
                    if default is None
                    else {"number": float(default)}
                    if not isinstance(default, bool)
                    else {"flag": bool(default)}
                ),
            }
        )
    return out


def outputs(spec: dict) -> list[dict]:
    out = []
    for name, declaration in spec.get("outputs", {}).items():
        unit = declaration.get("unit")
        out.append(
            {
                "name": name,
                "kind": kind(declaration),
                "unit": None if unit is None else str(unit),
                "dimension": dimension_of(unit),
                "values": [str(value) for value in declaration.get("values", [])],
                "description": clean(str(declaration.get("description", ""))),
                "optional": bool(declaration.get("optional", False)),
            }
        )
    return out


def example(spec: dict, block: str = "inputs") -> list[dict]:
    declarations = spec.get("inputs" if block == "inputs" else "outputs", {})
    kinds = {name: kind(declaration) for name, declaration in declarations.items()}
    out = []
    for name, value in spec.get("worked_example", {}).get(block, {}).items():
        if isinstance(value, (list, dict)):
            # A table-shaped example belongs to a calculation this tranche refuses; the sheet
            # editor has no field for it.
            continue
        if kinds.get(name) == "boolean":
            out.append({"name": name, "value": {"flag": bool(value)}})
        elif kinds.get(name) in ("enum", "string"):
            out.append({"name": name, "value": {"name": str(value)}})
        else:
            out.append({"name": name, "value": {"number": float(value)}})
    return out


def rust_str(text: str) -> str:
    return json.dumps(text, ensure_ascii=False)


def opt_str(text: str | None) -> str:
    return "None" if text is None else f"Some({rust_str(text)})"


def f64(value: float) -> str:
    return repr(float(value))


def opt_f64(value: float | None) -> str:
    return "None" if value is None else f"Some({f64(value)})"


def strings(values: list[str]) -> str:
    return "&[" + ", ".join(rust_str(value) for value in values) + "]"


def example_value(value: dict) -> str:
    if "number" in value:
        return f"ExampleValue::Number({f64(value['number'])})"
    if "flag" in value:
        return f"ExampleValue::Flag({'true' if value['flag'] else 'false'})"
    return f"ExampleValue::Name({rust_str(value['name'])})"


def emit_input(entry: dict) -> str:
    range_rows = ", ".join(
        "CalcRange { "
        f"min: {opt_f64(row['min'])}, min_inclusive: {'true' if row['min_inclusive'] else 'false'}, "
        f"max: {opt_f64(row['max'])}, max_inclusive: {'true' if row['max_inclusive'] else 'false'}, "
        f"equals: {opt_f64(row['equals'])}, band: {rust_str(row['band'])}, "
        f"severity: {rust_str(row['severity'])}, code: {rust_str(row['code'])}, "
        f"rationale: {rust_str(row['rationale'])} "
        "}"
        for row in entry["range"]
    )
    default = "None" if entry["default"] is None else f"Some({example_value(entry['default'])})"
    return (
        "CalcInput { "
        f"name: {rust_str(entry['name'])}, kind: {rust_str(entry['kind'])}, "
        f"required: {'true' if entry['required'] else 'false'}, "
        f"unit: {opt_str(entry['unit'])}, dimension: {opt_str(entry['dimension'])}, "
        f"values: {strings(entry['values'])}, description: {rust_str(entry['description'])}, "
        f"range: &[{range_rows}], default: {default} "
        "}"
    )


def emit_output(entry: dict) -> str:
    return (
        "CalcOutput { "
        f"name: {rust_str(entry['name'])}, kind: {rust_str(entry['kind'])}, "
        f"unit: {opt_str(entry['unit'])}, dimension: {opt_str(entry['dimension'])}, "
        f"values: {strings(entry['values'])}, description: {rust_str(entry['description'])}, "
        f"optional: {'true' if entry['optional'] else 'false'} "
        "}"
    )


def emit(spec: dict) -> str:
    unavailable = complex_inputs(spec) or complex_outputs(spec)
    name = name_of(spec)
    summary = summary_of(spec)
    equation = spec.get("equation")
    equation = None if equation is None or UPSTREAM in str(equation).lower() else str(equation)
    assumptions = clean_list([str(value) for value in spec.get("assumptions", [])])
    references = clean_list([str(value) for value in spec.get("references", [])])
    ins = ", ".join(emit_input(entry) for entry in inputs(spec))
    outs = ", ".join(emit_output(entry) for entry in outputs(spec))
    examples = ", ".join(
        f"CalcExample {{ name: {rust_str(entry['name'])}, value: {example_value(entry['value'])} }}"
        for entry in example(spec)
    )
    expected = ", ".join(
        f"CalcExample {{ name: {rust_str(entry['name'])}, value: {example_value(entry['value'])} }}"
        for entry in example(spec, "expected")
    )
    tolerance = float(spec.get("worked_example", {}).get("tolerance", 1e-12))
    return (
        "    CalcForm { "
        f"id: {rust_str(spec['id'])}, name: {rust_str(name)}, "
        f"summary: {rust_str(summary)}, "
        f"namespace: {rust_str(spec['id'].split('.', 1)[0])}, "
        f"available: {'false' if unavailable else 'true'}, "
        f"refusal: {opt_str(REFUSAL if unavailable else None)}, "
        f"equation: {opt_str(equation)}, assumptions: {strings(assumptions)}, "
        f"references: {strings(references)}, inputs: &[{ins}], outputs: &[{outs}], "
        f"worked_example: &[{examples}], expected: &[{expected}], "
        f"tolerance: {f64(tolerance)} "
        "},"
    )


def render_catalogue() -> str:
    rows = "\n".join(emit(spec) for _, spec in documents())
    return f"""//! GENERATED FILE - DO NOT EDIT BY HAND.
//!
//! Generated by `tools/gen_sheets.py` from `specs/calcs/`.
//! Regenerate with `python tools/gen_sheets.py`; CI runs `--check` and fails on a difference.
//!
//! **A UI projection, not a second spec.** Calculations whose inputs or outputs are tables are
//! carried with `available: false`; prose that names the upstream implementation is dropped.

// A worked example is the spec's own number, and one of them is π/4 written out. Clippy's
// approximation lint is about code that meant the constant; this is a value a spec pinned.
#![allow(clippy::approx_constant)]

use crate::form::{{CalcExample, CalcForm, CalcInput, CalcOutput, CalcRange, ExampleValue}};

/// Every registered calculation, in id order.
pub static CALCS: &[CalcForm] = &[
{rows}
];
"""


def formatted(text: str) -> str:
    """The generated Rust, through `rustfmt`, so `cargo fmt --check` and `--check` agree.

    **A generated file is still a Rust file.** Without this, `cargo fmt` reformats it and the
    next `--check` calls the tree stale - a drift with no source change. `gen_palette` already
    shells out to `rustfmt` for the same reason.
    """
    return subprocess.run(
        ["rustfmt", "--edition", "2024", "--emit", "stdout"],
        input=text,
        capture_output=True,
        text=True,
        check=True,
    ).stdout


def main() -> None:
    written = {OUT: formatted(render_catalogue()), DISPATCH_OUT: formatted(render_dispatch())}
    if "--check" in sys.argv:
        stale = []
        for path, text in written.items():
            current = path.read_text(encoding="utf-8") if path.exists() else ""
            if current != text:
                stale.append(str(path.relative_to(ROOT)))
        if stale:
            print(
                f"gen_sheets: {', '.join(stale)} is out of date. Run `python tools/gen_sheets.py`.",
                file=sys.stderr,
            )
            raise SystemExit(1)
        print("gen_sheets: catalogue and dispatch up to date")
        return
    for path, text in written.items():
        path.write_text(text, encoding="utf-8")
        print(f"gen_sheets: wrote {path.relative_to(ROOT)}")


# ---------------------------------------------------------------------------------------------
# The dispatch half
# ---------------------------------------------------------------------------------------------


class CannotGenerate(Exception):
    """This calculation needs a shape the generic emitter does not have."""


#: Kernel spellings that name a declared input under a different word. Read from the seven
#: signatures the wrapper generator left hand-written rather than guessed at: three packing/tray
#: kernels call the packing `name`, and the reaction kernel calls the reaction `reaction_name`.
ALIASES = {
    "name": "packing",
    "reaction_name": "reaction",
    "temperature": "T",
}


def resolve(param: str, spec: dict) -> str:
    inputs = spec["inputs"]
    if param in inputs:
        return param
    lower = {name.lower(): name for name in inputs}
    if param.lower() in lower:
        return lower[param.lower()]
    alias = ALIASES.get(param)
    if alias is not None and alias in inputs:
        return alias
    raise CannotGenerate(f"{spec['id']}: kernel parameter `{param}` names no declared input")


def input_expr(
    kind: str, decl: dict, name: str, aliases: set[str], source: str, crate: str, file: Path
) -> tuple[list[str], str]:
    """The statements and argument for one scalar or struct parameter."""
    unit = decl.get("unit")
    if kind == "f64":
        return [], f"number(inputs, {rust_str(name)})?"
    if kind == "Option<f64>":
        return [], f"optional_number(inputs, {rust_str(name)})?"
    if kind == "bool":
        return [], f"flag(inputs, {rust_str(name)})?"
    if kind == "Option<bool>":
        return [], f"optional_flag(inputs, {rust_str(name)})?"
    if kind in ("usize", "u32", "u64", "i64"):
        return [], f"unsigned(inputs, {rust_str(name)})? as {kind}"
    if kind in ("String",):
        return [], f"text(inputs, {rust_str(name)})?.to_string()"
    if kind == "Option<String>":
        return [], f"optional_text(inputs, {rust_str(name)})?.map(str::to_string)"
    if kind == "&str":
        return [], f"text(inputs, {rust_str(name)})?"
    if kind == "Option<&str>":
        return [], f"optional_text(inputs, {rust_str(name)})?"
    if kind in ("&[f64]", "Vec<f64>"):
        return [], f"numbers(inputs, {rust_str(name)})?"
    if kind == "Option<&[f64]>":
        return [], f"optional_numbers(inputs, {rust_str(name)})?.as_deref()"
    inner = kind[len("Option<") : -1] if kind.startswith("Option<") and kind.endswith(">") else kind
    if inner in aliases:
        ctor = (
            gpw.INTERVAL_CONSTRUCTORS["TemperatureInterval"]
            if decl.get("interval")
            else gpw.base_ctor(unit)
        )
        if ctor is None:
            # A dimensionless quantity has no constructor; it is the number itself.
            if decl.get("interval"):
                raise CannotGenerate(f"{spec_name(decl)} has an interval with no constructor")
            return [], f"number(inputs, {rust_str(name)})?"
        if kind.startswith("Option<"):
            return [], f"optional_number(inputs, {rust_str(name)})?.map({ctor})"
        return [], f"{ctor}(number(inputs, {rust_str(name)})?)"
    if kind.startswith("&") and kind[1:] not in aliases and not kind.startswith(("&[", "&str")):
        raise CannotGenerate(f"reference to struct `{kind[1:]}`")
    fields = gpw.record_fields(kind)
    if fields is not None:
        raise CannotGenerate(f"by-value struct `{kind}`")
    # A bare capitalised name that is not an alias is an enum this boundary parses.
    resolved = gpw.type_path(inner, source, crate, file)
    error = gpw.from_str_error(resolved)
    if kind.startswith("Option<"):
        if error == "infallible":
            return (
                [
                    f"    let {name}_parsed: Option<{resolved}> = optional_text(inputs, {rust_str(name)})?",
                    "        .map(|value| value.parse().unwrap_or_default());",
                ],
                f"{name}_parsed",
            )
        lines = [
            f"    let {name}_parsed: Option<{resolved}> = match optional_text(inputs, {rust_str(name)})? {{",
            "        Some(value) => match value.parse() {",
            "            Ok(parsed) => Some(parsed),",
            f"            Err(error) => return Err(Outcome::Error {{ message: error.to_string(), field: Some({rust_str(name)}.to_string()) }}),",
            "        },",
            "        None => None,",
            "    };",
        ]
        return lines, f"{name}_parsed"
    if error == "infallible":
        lines = [
            f"    let {name}_parsed: {resolved} = text(inputs, {rust_str(name)})?.parse().unwrap_or_default();"
        ]
    else:
        lines = [
            f"    let {name}_parsed: {resolved} = match text(inputs, {rust_str(name)})?.parse() {{",
            "        Ok(parsed) => parsed,",
            f"        Err(error) => return Err(Outcome::Error {{ message: error.to_string(), field: Some({rust_str(name)}.to_string()) }}),",
            "    };",
        ]
    return lines, f"{name}_parsed"


def spec_name(decl: dict) -> str:
    return str(decl.get("description", ""))[:40]


def struct_expr(
    type_name: str,
    spec: dict,
    aliases: set[str],
    source: str,
    crate: str,
    file: Path,
) -> tuple[list[str], str]:
    """One kernel state struct, built field by field from the declared inputs."""
    fields = gpw.record_fields(type_name)
    if fields is None:
        raise CannotGenerate(f"struct `{type_name}` is not declared in a crate")
    statements: list[str] = []
    arms: list[str] = []
    for field, field_type in fields:
        public = resolve(field, spec)
        decl = spec["inputs"][public]
        if field_type == "String":
            expr = f"text(inputs, {rust_str(public)})?.to_string()"
        elif field_type == "Option<String>":
            expr = f"optional_text(inputs, {rust_str(public)})?.map(str::to_string)"
        else:
            more, expr = input_expr(field_type, decl, public, aliases, source, crate, file)
            statements.extend(more)
        arms.append(f"{field}: {expr}")
    path = gpw.type_path(type_name, source, crate, file)
    return statements, f"{path} {{ {', '.join(arms)} }}"


def output_expr(
    name: str,
    rust_field: str,
    rust_type: str,
    declaration: dict,
) -> str:
    """One output's JSON insertion, or the statements that insert it conditionally."""
    kind = str(declaration.get("type", "quantity"))
    unit = declaration.get("unit")
    optional, depth, element = gpt.parse_type(rust_type)
    if depth != 0:
        raise CannotGenerate(f"output {name}: {rust_type} is a vector the sheet cannot draw")
    if kind == "enum":
        accessor = rust_index.string_accessor(element)
        value = f"r.{rust_field}.{accessor}()"
        if optional:
            return (
                f"            if let Some(value) = r.{rust_field} {{ "
                f"outputs.insert({rust_str(name)}.to_string(), json!(value.{accessor}())); }}"
            )
        return f"            outputs.insert({rust_str(name)}.to_string(), json!({value}));"
    if kind == "string":
        if optional:
            return (
                f"            if let Some(value) = r.{rust_field}.as_ref() {{ "
                f"outputs.insert({rust_str(name)}.to_string(), json!(value)); }}"
            )
        return f"            outputs.insert({rust_str(name)}.to_string(), json!(r.{rust_field}));"
    if kind == "boolean":
        if optional:
            return (
                f"            if let Some(value) = r.{rust_field} {{ "
                f"outputs.insert({rust_str(name)}.to_string(), json!(value)); }}"
            )
        return f"            outputs.insert({rust_str(name)}.to_string(), json!(r.{rust_field}));"
    # A quantity: a number for a dimensionless field, a quantity object otherwise.
    quantity = gpt.dimensioned(unit)
    piece = "value" if element in ("f64",) or element in gpt.INTEGERS else "value.value"
    if quantity:
        payload = f'json!({{"magnitude_si": {piece}, "unit": {rust_str(str(unit))}}})'
    else:
        payload = f"json!({piece})"
    if optional:
        return (
            f"            if let Some(value) = r.{rust_field} {{ "
            f"outputs.insert({rust_str(name)}.to_string(), {payload}); }}"
        )
    source = (
        f"r.{rust_field}"
        if element in ("f64",) or element in gpt.INTEGERS
        else f"r.{rust_field}.value"
    )
    if quantity:
        payload = f'json!({{"magnitude_si": {source}, "unit": {rust_str(str(unit))}}})'
    else:
        payload = f"json!({source})"
    return f"            outputs.insert({rust_str(name)}.to_string(), {payload});"


def dispatch_function(index: int, spec: dict, kernels: dict, aliases: set[str]) -> str:
    function = spec["implementations"]["python"].split(".")[-1]
    if function not in kernels:
        raise CannotGenerate(f"{spec['id']}: no kernel signature for `{function}`")
    kernel, crate, file = kernels[function]
    source = rust_index._blank_comments(file.read_text(encoding="utf-8"))
    path = rust_index._item_path(file, function)
    statements: list[str] = []
    arguments: list[str] = []
    for param, kind in kernel:
        if kind.startswith("&") and kind[1:] not in aliases and not kind.startswith(("&[", "&str")):
            raise CannotGenerate(f"{spec['id']}: {kind} is a struct reference")
        fields = gpw.record_fields(kind) if kind not in aliases else None
        if fields is not None:
            more, expr = struct_expr(kind, spec, aliases, source, crate, file)
            statements.extend(more)
            arguments.append(expr)
            continue
        public = resolve(param, spec)
        decl = spec["inputs"][public]
        more, expr = input_expr(kind, decl, public, aliases, source, crate, file)
        statements.extend(more)
        arguments.append(expr)
    result = next((r for r in rust_index.result_types() if r.calc_id == spec["id"]), None)
    if result is None:
        raise CannotGenerate(f"{spec['id']}: no result type")
    pairs = dict(zip(result.public_fields, result.rust_fields, strict=True))
    types = dict(zip(result.rust_fields, result.rust_field_types, strict=True))
    inserts: list[str] = []
    for name, declaration in spec.get("outputs", {}).items():
        rust_field = pairs.get(name)
        if rust_field is None:
            raise CannotGenerate(f"{spec['id']}: output `{name}` is not in the result")
        inserts.append(output_expr(name, rust_field, types[rust_field], declaration))
    body = "\n".join(statements)
    inserts_text = "\n".join(inserts)
    return f"""fn calc_{index}(inputs: &Value) -> Result<Outcome, Outcome> {{
{body}
    Ok(match {path}({", ".join(arguments)}) {{
        Ok(r) => {{
            let mut outputs = BTreeMap::new();
{inserts_text}
            Outcome::Ok(CalcOutcome {{
                outputs,
                warnings: warnings(r.warnings()),
            }})
        }}
        Err(error) => refused(error),
    }})
}}
"""


def render_dispatch() -> str:
    kernels = gpw.kernel_signatures()
    aliases = gpw.uom_aliases()
    functions: list[str] = []
    arms: list[str] = []
    for index, (_, spec) in enumerate(documents()):
        if complex_inputs(spec) or complex_outputs(spec):
            arms.append(
                f"        {rust_str(spec['id'])} => Outcome::Error {{ message: {rust_str(REFUSAL)}.to_string(), field: None }},"
            )
            continue
        try:
            functions.append(dispatch_function(index, spec, kernels, aliases))
        except CannotGenerate as refusal:
            print(f"gen_sheets: cannot dispatch {refusal}", file=sys.stderr)
            raise SystemExit(1) from None
        arms.append(
            f"        {rust_str(spec['id'])} => calc_{index}(inputs).unwrap_or_else(|outcome| outcome),"
        )
    functions_text = "\n".join(functions)
    used = sorted(name for name in gpw._ALL_CTORS if re.search(rf"\b{name}\(", functions_text))
    ctors = ", ".join(used)
    return f"""//! GENERATED FILE - DO NOT EDIT BY HAND.
//!
//! Generated by `tools/gen_sheets.py` from `specs/calcs/`.
//! Regenerate with `python tools/gen_sheets.py`; CI runs `--check` and fails on a difference.
//!
//! **The boundary is SI.** A quantity input is an SI magnitude; a quantity output is
//! `{{magnitude_si, unit}}`. A kernel's `AzothError` is an outcome rather than an exception.

#![allow(clippy::approx_constant)]

use std::collections::BTreeMap;

use serde_json::{{Value, json}};

use azoth_core::CalcResult;
use azoth_core::units::{{{ctors}}};

use crate::dispatch::Outcome;
#[allow(unused_imports)]
use crate::dispatch_support::{{flag, number, numbers, optional_flag, optional_number, optional_numbers, optional_text, refused, text, unsigned, warnings}};
use crate::workspace::CalcOutcome;

/// Run one registered calculation.
#[must_use]
pub fn run(id: &str, inputs: &Value) -> Outcome {{
    match id {{
{chr(10).join(arms)}
        _ => Outcome::Error {{
            message: format!("`{{id}}` is not a calculation in this build"),
            field: None,
        }},
    }}
}}

{functions_text}
"""


if __name__ == "__main__":
    main()
