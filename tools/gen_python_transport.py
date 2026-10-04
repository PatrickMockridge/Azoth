#!/usr/bin/env python3
"""Generate the extension's transport structs, one per registered result.

`crates/azoth-python/src/results.rs` carried a `Py*` struct, a `#[pymethods] __repr__` and a
`From<&Kernel>` impl for every registered result - 194 of them, a hand-written block each. A new
calculation cost an edit there beside its spec, its two kernels and its test, and the block was
a *function of data already in the tree*: the field list and order come from `CalcResult::FIELDS`,
the Rust field names and types from the struct itself (which is the one thing no spec carries,
so `tools/rust_index.py` reads it), and the unit a dimensioned field transports in is the spec's
own `outputs.<name>.unit`.

**Two things are editorial rather than derived, and this file states them rather than hiding
them.** A field's doc comment would be a hand paraphrase; it is emitted as the spec's own
`description`, so the prose lives in one place. A `__repr__`'s choice of which fields to show
would be an editorial judgement per type; it is emitted uniformly - every field but `warnings`,
then the warning count - so the string is a function of the field list.

**What a field transports is decided by its declared unit.** A dimensionless output is a number
and a dimensioned one is a [`PyQty`], whatever the kernel happens to store: the kernels are typed
against `uom::si::f64`, where a bare `f64` and a `uom` alias are the same number in two spellings.
That rule is the *reason* this layer is generated rather than copied - the hand-written types
disagreed with it on 37 fields (`flooding_velocity` was a `float` in `m/s` while
`pressure_drop_per_meter` was a `Qty` in the same struct), and a rule that is applied by a
generator cannot drift field by field.

    python tools/gen_python_transport.py            # write
    python tools/gen_python_transport.py --check    # fail if the tree is not what this emits
"""

from __future__ import annotations

import subprocess
import sys
import tomllib
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import rust_index

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "crates" / "azoth-python" / "src" / "transport_gen.rs"

#: Rust integer widths a result field can carry. The width is the kernel's own choice - an
#: iteration counter is a `u32`, a tray index a `usize` - so it is read from the type rather
#: than chosen here.
INTEGERS = ("u8", "u16", "u32", "u64", "usize", "i8", "i16", "i32", "i64", "isize")

#: A unit that means "no unit". Anything else is dimensioned and transports as a `PyQty`.
DIMENSIONLESS = (None, "dimensionless")

PY_TYPE = {
    "bool": "bool",
    "string": "String",
    "enum": "String",
    "warning": "PyWarning",
    "kcomponent": "PyKComponent",
}


def specs() -> dict[str, dict]:
    """Every spec, by id - calcs and models both, since both register results."""
    out: dict[str, dict] = {}
    for namespace in ("calcs", "models"):
        for path in sorted((ROOT / "specs" / namespace).rglob("*.toml")):
            document = tomllib.loads(path.read_text(encoding="utf-8"))
            out[document["id"]] = document
    return out


def parse_type(ty: str) -> tuple[bool, int, str]:
    """A Rust field type as `(optional, vector depth, element type)`.

    Depth is `0` for a scalar, `1` for `Vec<T>` and `2` for `Vec<Vec<T>>`. `Option` is only
    ever outermost in this tree, so it is peeled once rather than tracked as a depth.
    """
    optional = False
    if ty.startswith("Option<") and ty.endswith(">"):
        optional = True
        ty = ty[len("Option<") : -1]
    depth = 0
    while ty.startswith("Vec<") and ty.endswith(">"):
        depth += 1
        ty = ty[len("Vec<") : -1]
    return optional, depth, ty


def element_kind(element: str, declared: str) -> str:
    """What one element is, from its Rust type and the spec's own `type`.

    A capitalised name that is not one of the framework's is either a `uom` alias or an enum,
    and the two are spelled the same way in Rust. The spec already claims one or the other -
    `type = "enum"` with its values - so the spec is asked rather than a name list kept here.
    """
    if element == "f64":
        return "number"
    if element in INTEGERS:
        return element
    if element == "bool":
        return "bool"
    if element == "String":
        return "string"
    if element == "Warning":
        return "warning"
    if element == "KComponent":
        return "kcomponent"
    return "enum" if declared == "enum" else "quantity"


def dimensioned(unit: object) -> bool:
    return unit not in DIMENSIONLESS


def py_element(kind: str, unit: object) -> str:
    """The Python-visible type of one element.

    A `uom` alias is a quantity only when its unit is dimensioned: `Ratio` is `uom`'s
    dimensionless type, and `azoth-core`'s own convention is that an unadorned ratio is a plain
    `f64` in the Python API. The declared unit decides for both spellings.
    """
    if kind in ("number", "quantity"):
        return "PyQty" if dimensioned(unit) else "f64"
    return PY_TYPE.get(kind, kind)


def element_expr(kind: str, unit: object, src: str, element: str = "") -> str:
    """One element's conversion, from a source expression bound to `{src}`."""
    if kind == "number":
        if dimensioned(unit):
            return f'PyQty {{ magnitude_si: {src}, unit: "{unit}".to_string() }}'
        return src
    if kind == "quantity":
        if dimensioned(unit):
            return f'PyQty {{ magnitude_si: {src}.value, unit: "{unit}".to_string() }}'
        return f"{src}.value"
    if kind == "enum":
        return f"{src}.{rust_index.string_accessor(element)}().to_string()"
    if kind == "string":
        return f"{src}.clone()"
    return src


def identity(kind: str, unit: object) -> bool:
    """Whether an element transports as itself, so a container can `clone()` instead of `map`."""
    return (
        (kind == "number" and not dimensioned(unit))
        or kind in INTEGERS
        or kind in ("bool", "string")
    )


def field_expr(
    field: str, kind: str, depth: int, optional: bool, unit: object, element: str
) -> str:
    """The right-hand side of one transport field's initialiser."""
    base = f"r.{field}"
    if kind == "warning":
        return f"transport(&{base})"
    if kind == "kcomponent":
        return (
            f"{base}.iter().map(|c| PyKComponent {{ "
            f"fitting_id: c.fitting_id.clone(), n_ld: c.n_ld, k: c.k }}).collect()"
        )
    if optional:
        if identity(kind, unit):
            return base
        return f"{base}.map(|v| {element_expr(kind, unit, 'v', element)})"
    if depth == 0:
        return element_expr(kind, unit, base, element)
    if identity(kind, unit):
        return f"{base}.clone()"
    # An element bound from a slice is a reference, so a plain `f64` needs its value taken;
    # a `uom` alias or an enum derefs through the method call and is left alone.
    src = "*v" if kind == "number" else "v"
    if depth == 1:
        return f"{base}.iter().map(|v| {element_expr(kind, unit, src, element)}).collect()"
    return (
        f"{base}.iter().map(|row| "
        f"row.iter().map(|v| {element_expr(kind, unit, src, element)}).collect()).collect()"
    )


def doc(text: str) -> str:
    """One `///` line from arbitrary prose, so a description's newlines cannot break the doc."""
    return " ".join(text.split())


def field_type(kind: str, depth: int, optional: bool, unit: object) -> str:
    inner = py_element(kind, unit)
    rendered = "Vec<" * depth + inner + ">" * depth
    return f"Option<{rendered}>" if optional else rendered


def struct_for(result: rust_index.ResultType, spec: dict) -> str:
    """One result's struct, its `__repr__` and its `From` impl.

    **The transported field is the public name, and the initialiser reads the Rust one.** They
    are the same string on 174 of the 194 results; on the `eos` flashes they are not, and the
    transport must expose `T` - what the spec, `FIELDS` and the Python dataclass all call it -
    while the kernel struct spells the same field `temperature`.
    """
    outputs = spec.get("outputs") or {}
    fields: list[tuple[str, str, str, str]] = []
    for (public, rust), rust_type in zip(result.fields, result.rust_field_types, strict=True):
        declared = outputs.get(public) or {}
        if public == "warnings":
            kind, depth, optional, unit = "warning", 1, False, None
            text = "Caveats, deduplicated."
        else:
            optional, depth, element = parse_type(rust_type)
            kind = element_kind(element, declared.get("type", "scalar"))
            unit = declared.get("unit")
            text = doc(declared.get("description", public))
        fields.append(
            (
                public,
                field_type(kind, depth, optional, unit),
                text,
                field_expr(rust, kind, depth, optional, unit, element),
            )
        )

    non_snake = any(public != public.lower() for public, *_ in fields)
    lines = [
        f"/// Result of `{result.calc_id}`, transported.",
        "#[pyclass(",
        "    frozen,",
        "    skip_from_py_object,",
        '    module = "azoth._core",',
        f'    name = "{result.rust_name}"',
        ")]",
        "#[derive(Debug, Clone, PartialEq)]",
    ]
    if non_snake:
        lines.append("#[allow(non_snake_case)] // the spec's own symbol")
    lines.append(f"pub struct Py{result.rust_name} {{")
    for public, ty, text, _ in fields:
        lines.append(f"    /// {text}")
        lines.append("    #[pyo3(get)]")
        lines.append(f"    pub {public}: {ty},")
    lines.append("}")
    lines.append("")
    lines.extend(pymethods(result, fields))
    lines.append("")
    lines.extend(from_impl(result, fields))
    return "\n".join(lines)


def pymethods(result: rust_index.ResultType, fields: list[tuple[str, str, str, str]]) -> list[str]:
    """The `__repr__`, uniform over the field list rather than chosen per type."""
    shown = [public for public, _, _, _ in fields if public != "warnings"]
    parts = ", ".join(f"{public}={{:?}}" for public in shown)
    args = ", ".join(f"self.{public}" for public in shown)
    body = (
        f"{result.rust_name}({parts}, {{}} warning(s))"
        if parts
        else f"{result.rust_name}({{}} warning(s))"
    )
    arguments = f"{args}, self.warnings.len()" if args else "self.warnings.len()"
    return [
        "#[pymethods]",
        f"impl Py{result.rust_name} {{",
        "    fn __repr__(&self) -> String {",
        "        format!(",
        f'            "{body}",',
        f"            {arguments}",
        "        )",
        "    }",
        "}",
    ]


def from_impl(result: rust_index.ResultType, fields: list[tuple[str, str, str, str]]) -> list[str]:
    """`From<&Kernel>` - one field per line, so a diff names the field that moved."""
    lines = [
        f"impl From<&{result.item_path}> for Py{result.rust_name} {{",
        f"    fn from(r: &{result.item_path}) -> Self {{",
        "        Self {",
    ]
    for public, _, _, expr in fields:
        lines.append(f"            {public}: {expr},")
    lines.extend(["        }", "    }", "}"])
    return lines


def emit() -> str:
    """The whole file, in id order so it is a stable diff."""
    rust_index.check()
    known = specs()
    blocks = []
    for result in rust_index.result_types():
        spec = known.get(result.calc_id)
        if spec is None:
            sys.exit(f"gen_python_transport: {result.calc_id} has no spec")
        blocks.append(struct_for(result, spec))
    header = (
        "//! GENERATED FILE - DO NOT EDIT BY HAND.\n"
        "//!\n"
        "//! Generated by `tools/gen_python_transport.py` from:\n"
        "//!   - specs/calcs/**/*.toml and specs/models/**/*.toml, for each output's unit,\n"
        "//!     declared kind and description\n"
        "//!   - the `impl CalcResult` blocks under `crates/*/src`, for each result type's\n"
        "//!     field names and Rust types (`tools/rust_index.py`)\n"
        "//!\n"
        "//! The transport model is `frozen` and carries a `#[pyo3(get)]` per field, so a result\n"
        "//! can be read from Python but not built there: every one of these is produced by a\n"
        "//! calculation, and a hand-built one would be a value with no provenance behind it.\n"
        "//!\n"
        "//! Regenerate with `python tools/gen_python_transport.py`; CI runs `--check` and fails\n"
        "//! on any difference.\n"
        "\n"
        "use crate::results::{PyKComponent, PyQty, PyWarning, transport};\n"
        "use pyo3::prelude::*;\n"
        "\n"
    )
    return header + "\n".join(blocks) + "\n"


def format_rust(text: str) -> str:
    """Run the emitted Rust through rustfmt, which the committed file is compared against."""
    try:
        proc = subprocess.run(
            ["rustfmt", "--edition", "2024"],
            input=text,
            capture_output=True,
            text=True,
            check=False,
        )
    except FileNotFoundError:
        sys.exit(
            "gen_python_transport: rustfmt is not on PATH, and the generated file is compared "
            "against rustfmt's own output - install the Rust toolchain"
        )
    if proc.returncode != 0:
        sys.exit(f"gen_python_transport: rustfmt refused the output:\n{proc.stderr}")
    return proc.stdout


def main() -> None:
    check = "--check" in sys.argv
    rendered = format_rust(emit())
    existing = OUT.read_text(encoding="utf-8") if OUT.exists() else None
    if existing == rendered:
        print(f"gen_python_transport: {OUT.relative_to(ROOT)} is current")
        return
    if check:
        sys.exit(
            f"gen_python_transport: {OUT.relative_to(ROOT)} is out of date; run "
            "tools/gen_python_transport.py"
        )
    OUT.write_text(rendered, encoding="utf-8")
    print(f"gen_python_transport: wrote {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
