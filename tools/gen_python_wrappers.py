#!/usr/bin/env python3
"""Generate the `#[pyfunction]` wrappers whose signature is the spec's declared input list.

`crates/azoth-python/src/{hydraulics,eos,thermal,reactions,standards}.rs` carried one
hand-written wrapper per calculation, each converting SI magnitudes into the kernel's own
quantities and mapping the result into its transport. The conversion is a function of two things
already in the tree - the kernel's parameter list, and the unit the spec declares for each input -
so the wrapper is emitted rather than written.

**The coverage rule is read from the kernel, not from a list.** An id is covered when the
kernel's parameter names are the spec's declared input names; that is the plan's own "the
signature *is* the declared input list". A mixture-expanded model takes a `Mixture` the spec does
not declare, so its kernel names do not match and this file leaves it alone - which is why the
`expands_to` declaration the plan wants for those is *not* a prerequisite here.

Two things are reused rather than re-derived: the optional-last ordering is `gen_stub`'s
`ordered_parameters` (Python has no syntax for a defaulted parameter before a required one), and
a unit's constructor is `rust_index.rust_ctor`, the vocabulary's own map.

    python tools/gen_python_wrappers.py            # write
    python tools/gen_python_wrappers.py --check    # fail if the tree is not what this emits
    python tools/gen_python_wrappers.py --survey   # what it covers, and what it leaves
"""

from __future__ import annotations

import re
import subprocess
import sys
import tomllib
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import gen_stub
import rust_index

ROOT = Path(__file__).resolve().parent.parent
CRATES = ROOT / "crates"
OUT = ROOT / "crates" / "azoth-python" / "src" / "wrappers_gen.rs"

#: Kernel parameter types that cross as themselves, with the type the binding should declare.
PASSTHROUGH = {
    "f64": "f64",
    "bool": "bool",
    "usize": "usize",
    "u32": "u32",
    "u64": "u64",
    "i64": "i64",
    "&str": "&str",
    "Option<f64>": "Option<f64>",
    "Option<bool>": "Option<bool>",
    "Option<usize>": "Option<usize>",
    "Option<u32>": "Option<u32>",
    "Option<&str>": "Option<&str>",
    "&[f64]": "Vec<f64>",
    "&[String]": "Vec<String>",
    "&[Vec<f64>]": "Vec<Vec<f64>>",
    "&[[f64; 3]]": "Vec<[f64; 3]>",
    "Option<&[f64]>": "Option<Vec<f64>>",
}

#: Kernel parameters that need a local binding before the call, and what it is.
NEEDS_LOCAL = {"&[&str]": ("Vec<String>", "Vec<&str>", "iter().map(String::as_str).collect()")}


def uom_aliases() -> set[str]:
    """The `uom` quantity aliases `azoth-core` re-exports, read from its own `pub use`.

    Read rather than listed, for the reason every other index here is: a type that is a `uom`
    quantity converts through a constructor, and one that is an enum is parsed from a string.
    Deciding which by convention would be a second place to change when the re-export does.
    """
    source = rust_index._blank_comments(
        (CRATES / "azoth-core" / "src" / "units.rs").read_text(encoding="utf-8")
    )
    match = re.search(r"pub use uom::si::f64::\{([^}]*)\}", source)
    if match is None:
        sys.exit("gen_python_wrappers: azoth-core/src/units.rs no longer re-exports uom aliases")
    return {name.strip() for name in match.group(1).split(",") if name.strip()}


def si_ctors() -> dict[str, str]:
    """Each dimension's SI base constructor, from the vocabulary's own order.

    **The binding takes SI magnitudes**, whatever unit the spec declares - `azoth.hydraulics`
    converts once in Python and hands over the base value. So the constructor is not the declared
    unit's (`millimeters` for an input declared in `mm`, which would scale it a second time) but
    the base unit's for the same dimension. The vocabulary lists the SI base first in every
    dimension, which is the property this reads rather than a name rule kept here.
    """
    table = rust_index._table()
    out: dict[str, str] = {}
    for unit in table["units"]:
        ctor = unit.get("rust_ctor")
        if ctor:
            out.setdefault(str(unit["dimension"]), str(ctor))
    return out


def kernel_signatures() -> dict[str, tuple[list[tuple[str, str]], str, str]]:
    """Every `pub fn` the kernels declare, by name, as `(parameters, crate, item path)`.

    The crate is carried so an enum parameter can be named in full where the binding uses it.

    By name rather than by path: a spec's `implementations.rust` names the item, and which module
    it lives in is the crate's business. A name two crates both use is dropped rather than picked,
    so nothing is emitted against a signature that might be the other one.

    **`azoth-python` is not scanned.** It is this binding, and every wrapper it carries has the
    kernel's own name - so including it would make each kernel look like a duplicate of itself.
    """
    out: dict[str, tuple[list[tuple[str, str]], str, str]] = {}
    seen: dict[str, int] = {}
    for path in sorted(CRATES.glob("*/src/**/*.rs")):
        parts = path.relative_to(CRATES).parts
        if parts[0] == "azoth-python":
            continue
        module = parts[0].replace("-", "_")
        text = rust_index._blank_comments(path.read_text(encoding="utf-8"))
        for name, params in signatures(text):
            if name in seen:
                seen[name] += 1
                out.pop(name, None)
                continue
            seen[name] = 1
            out[name] = (params, module, path)
    return {name: value for name, value in out.items() if seen[name] == 1}


def signatures(text: str) -> list[tuple[str, list[tuple[str, str]]]]:
    """Every top-level `pub fn name(...)` in one file, as `(name, parameters)`."""
    out = []
    for match in re.finditer(r"^pub fn (\w+)\(", text, re.M):
        depth, i = 0, match.end() - 1
        while i < len(text):
            if text[i] == "(":
                depth += 1
            elif text[i] == ")":
                depth -= 1
                if depth == 0:
                    break
            i += 1
        inner, parts, depth, start = text[match.end() : i], [], 0, 0
        for j, ch in enumerate(inner):
            if ch in "<([":
                depth += 1
            elif ch in ">)]":
                depth -= 1
            elif ch == "," and depth == 0:
                parts.append(inner[start:j])
                start = j + 1
        parts.append(inner[start:])
        params: list[tuple[str, str]] = []
        for part in parts:
            found = re.match(r"(?:mut\s+)?(\w+)\s*:\s*(.+)$", " ".join(part.split()))
            if found:
                params.append((found.group(1), found.group(2)))
        out.append((match.group(1), params))
    return out


def type_path(name: str, source: str, crate: str, file: Path) -> str:
    """Where a kernel names a type it takes, read from that file's own `use`.

    A parameter's type is spelled short in the kernel - `AntoineForm`, `Cubic` - and pyo3 needs
    the path. The kernel's `use` line is where the path is written, so it is read rather than
    guessed: `Cubic` comes from `azoth_eos` and `LiquidViscosityLadder` from a submodule, and a
    rule that prefixed the kernel's own crate would name the wrong one for both.
    """
    for match in re.finditer(r"^use\s+([^;]+);", source, re.M):
        path = " ".join(match.group(1).split())
        if path.endswith("}") and "{" in path:
            base, _, inner = path.partition("::{")
            if re.search(rf"\b{re.escape(name)}\b", inner):
                return f"{_root(base, crate)}::{name}"
            continue
        if path.rsplit("::", 1)[-1] == name:
            return f"{_root(path.rsplit('::', 1)[0], crate)}::{name}"
    if re.search(rf"^pub (?:enum|struct) {re.escape(name)}\b", source, re.M):
        # Defined in the kernel's own file: the path is that file's, not a `use`.
        return rust_index._item_path(file, name)
    return f"{crate}::{name}"


def _root(base: str, crate: str) -> str:
    """A `use` path's root, with `crate::` resolved to the package the binding names."""
    if base == "crate":
        return crate
    if base.startswith("crate::"):
        return base.replace("crate", crate, 1)
    if base.startswith(("self::", "super::")):
        return f"{crate}::{base.split('::', 1)[1]}"
    return base


def specs_by_id() -> dict[str, dict]:
    out: dict[str, dict] = {}
    for namespace in ("calcs", "models"):
        for path in sorted((ROOT / "specs" / namespace).rglob("*.toml")):
            document = tomllib.loads(path.read_text(encoding="utf-8"))
            out[document["id"]] = document
    return out


def base_ctor(unit: str | None) -> str | None:
    """The SI base unit's constructor for a declared unit's dimension, or `None`."""
    if not unit:
        return None
    dimension = rust_index._units()[unit]["dimension"] if unit in rust_index._units() else None
    return si_ctors().get(str(dimension))


class Refusal(Exception):
    """This file cannot derive the wrapper, so it stays hand-written."""


def convert(kind: str, unit: str | None, name: str, aliases: set[str]) -> tuple[str, str]:
    """`(declared type, argument)` for one kernel parameter, or a refusal.

    The kernel's own type decides the conversion: a `uom` alias goes through the vocabulary's
    constructor for the unit the spec declares, a bare `f64` crosses as it is, and an enum is
    parsed from the `&str` pyo3 hands over.
    """
    if kind in PASSTHROUGH:
        declared = PASSTHROUGH[kind]
        # A slice crosses as an owned `Vec`, and the kernel wants a borrow of it.
        return declared, (f"&{name}" if declared.startswith("Vec<") else name)
    if kind in NEEDS_LOCAL:
        declared, local, body = NEEDS_LOCAL[kind]
        return declared, f"&{name}_refs"
    if kind.startswith("Option<") and kind.endswith(">"):
        inner = kind[len("Option<") : -1]
        if inner in aliases:
            ctor = base_ctor(unit)
            if ctor is None:
                raise Refusal(f"Option<{inner}> has no constructor for unit {unit!r}")
            return "Option<f64>", f"{name}.map({ctor})"
        raise Refusal(f"Option<{inner}> is neither a uom alias nor a passthrough")
    if kind.startswith("&[") and kind.endswith("]"):
        inner = kind[2:-1]
        if inner not in aliases:
            raise Refusal(f"&[{inner}] is not a uom alias")
        ctor = base_ctor(unit)
        if ctor is None:
            raise Refusal(f"&[{inner}] has no constructor for unit {unit!r}")
        return "Vec<f64>", f"&{name}.iter().map(|v| {ctor}(*v)).collect::<Vec<_>>()"
    if kind in aliases:
        ctor = base_ctor(unit)
        if ctor is None:
            raise Refusal(f"{kind} has no base constructor for unit {unit!r}")
        if rust_index.uom_quantity(str(unit)) != kind:
            # A unit's constructor yields the quantity its *dimension* names, and a kernel that
            # wants a different one - `TemperatureInterval` for a `dT` declared in `K` - needs a
            # conversion this file has no rule for.
            raise Refusal(f"unit {unit!r} constructs {rust_index.uom_quantity(str(unit))}, not {kind}")
        return "f64", f"{ctor}({name})"
    # A capitalised name that is not a `uom` alias is an enum: parsed from the string pyo3 gives.
    return "&str", name


def locals_for(kind: str, name: str) -> list[str]:
    """The `let` statements a parameter needs, before the call."""
    if kind == "&[&str]":
        return [f"    let {name}_refs: Vec<&str> = {name}.iter().map(String::as_str).collect();"]
    return []


def enum_local(kind: str, name: str, resolved: str) -> tuple[list[str], str]:
    """An enum parameter, parsed into a local. `(statements, argument)`."""
    if kind.startswith("Option<"):
        inner = kind[len("Option<") : -1]
        return (
            [
                f"    let {name}_parsed: Option<{resolved}> = {name}",
                "        .map(|value| value.parse().map_err(pyo3::exceptions::PyValueError::new_err))",
                "        .transpose()?;",
            ],
            f"{name}_parsed",
        )
    return (
        [
            f"    let {name}_parsed: {resolved} = {name}",
            "        .parse()",
            "        .map_err(pyo3::exceptions::PyValueError::new_err)?;",
        ],
        f"{name}_parsed",
    )


def body(entry, kernel, crate: str, path: str, source: str, file: Path, aliases: set[str]) -> str:
    calc_id, module, function, spec = entry
    inputs = spec["inputs"]
    params = dict(kernel)

    # What each parameter is, and whether anything about it was refused.
    plan: list[tuple[str, str, str, bool]] = []  # (parameter, kind, argument, is_enum)
    for name, kind in kernel:
        is_enum = kind not in PASSTHROUGH and kind not in aliases and not kind.startswith(("&[", "Option<"))
        if kind.startswith("Option<") and kind.endswith(">"):
            is_enum = kind[len("Option<") : -1] not in aliases and kind not in PASSTHROUGH
        try:
            _declared, argument = convert(kind, inputs[name].get("unit"), name, aliases)
        except Refusal as refusal:
            raise Refusal(str(refusal)) from None
        plan.append((name, kind, argument if not is_enum else name, is_enum))

    statements: list[str] = []
    arguments: list[str] = []
    for name, kind, argument, is_enum in plan:
        statements.extend(locals_for(kind, name))
        if is_enum:
            inner = kind[len("Option<") : -1] if kind.startswith("Option<") else kind
            parsed, argument = enum_local(kind, name, type_path(inner, source, crate, file))
            statements.extend(parsed)
        arguments.append(argument)

    # The declared signature, with the optional ones last, which is `gen_stub`'s rule.
    declarations = []
    for name, decl in gen_stub.ordered_parameters(inputs):
        declared, _ = convert(dict(kernel)[name], decl.get("unit"), name, aliases)
        if decl.get("optional") and not declared.startswith("Option<") and declared != "&str" and not declared.startswith("Vec<"):
            declared = f"Option<{declared}>"
        declarations.append(f"{name}: {declared}")

    rust_name = next(r.rust_name for r in rust_index.result_types() if r.calc_id == calc_id)
    names = [name for name, _ in gen_stub.ordered_parameters(inputs)]
    signature_list = ", ".join(names)
    lines = [
        f"/// {spec.get('name', function)}.",
        "///",
        "/// All arguments are SI magnitudes. See the module documentation for why.",
        "#[pyfunction]",
        f"#[pyo3(signature = ({signature_list}))]",
        f'#[pyo3(text_signature = "({signature_list})")]',
    ]
    if any(name != name.lower() for name in names):
        lines.append("#[allow(non_snake_case)] // symbols from the published equation")
    if len(declarations) > 6:
        lines.append(
            "#[allow(clippy::too_many_arguments)] // The signature is the spec's declared inputs."
        )
    lines.append(f"pub fn {function}(")
    lines.append("    py: Python<'_>,")
    lines.extend(f"    {declaration}," for declaration in declarations)
    lines.append(f") -> PyResult<crate::transport_gen::Py{rust_name}> {{")
    lines.extend(statements)
    lines.append(f"    {path}({', '.join(arguments)})")
    lines.append(f"        .map(|r| crate::transport_gen::Py{rust_name}::from(&r))")
    lines.append("        .map_err(|e| to_pyerr(py, e))")
    lines.append("}")
    return "\n".join(lines)


def covered() -> list[tuple[tuple[str, str, str, dict], str]]:
    """Every id this file can emit, with its rendered wrapper."""
    kernels = kernel_signatures()
    aliases = uom_aliases()
    specs = specs_by_id()
    out = []
    for calc_id, module, function in rust_index.implementations():
        if module == "process":
            continue
        found = kernels.get(function)
        spec = specs[calc_id]
        entry = (calc_id, module, function, spec)
        if found is None:
            continue
        kernel, kernel_module, kernel_file = found
        kernel_path = rust_index._item_path(kernel_file, function)
        if set(n for n, _ in kernel) != set(spec["inputs"]):
            continue
        try:
            out.append((
                entry,
                body(entry, kernel, kernel_module, kernel_path,
                     rust_index._blank_comments(kernel_file.read_text(encoding="utf-8")),
                     kernel_file, aliases),
            ))
        except Refusal:
            continue
    return out


def ctors(text: str) -> list[str]:
    """The vocabulary constructors the emitted wrappers call."""
    return sorted(set(re.findall(r"\b([a-z_]+)\((?=[A-Za-z_])", text)) & set(_ALL_CTORS))


def _all_ctors() -> set[str]:
    table = rust_index._table()
    return {str(u["rust_ctor"]) for u in table["units"] if u.get("rust_ctor")}


_ALL_CTORS = _all_ctors()


def covered_ids() -> set[tuple[str, str]]:
    """The `(module, function)` pairs this file emits, which `gen_python_register` points at."""
    return {(entry[1], entry[2]) for entry, _ in covered()}


def emit() -> str:
    bodies = "\n\n".join(text for _, text in covered())
    used = ctors(bodies)
    imports = (
        "use crate::errors::to_pyerr;\n"
        + (f"use azoth_core::units::{{{', '.join(used)}}};\n" if used else "")
        + "use pyo3::prelude::*;\n"
    )
    header = (
        "//! GENERATED FILE - DO NOT EDIT BY HAND.\n"
        "//!\n"
        "//! Generated by `tools/gen_python_wrappers.py` from each spec's declared inputs and the\n"
        "//! kernel's own signature - the two agree by name on every id this file covers.\n"
        "//!\n"
        "//! Regenerate with `python tools/gen_python_wrappers.py`; CI runs `--check` and fails on\n"
        "//! any difference.\n"
        "\n"
    )
    return header + imports + "\n" + bodies + "\n"


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
            "gen_python_wrappers: rustfmt is not on PATH, and the generated file is compared "
            "against rustfmt's own output - install the Rust toolchain"
        )
    if proc.returncode != 0:
        sys.exit(f"gen_python_wrappers: rustfmt refused the output:\n{proc.stderr}")
    return proc.stdout


def main() -> None:
    if "--survey" in sys.argv:
        taken = covered()
        total = len(rust_index.implementations())
        print(f"covered: {len(taken)} of {total}")
        names = {cid for (cid, _, _, _), _ in taken}
        for calc_id, module, function in rust_index.implementations():
            if calc_id not in names:
                print(f"   left hand-written: {calc_id}")
        return
    check = "--check" in sys.argv
    rendered = format_rust(emit())
    existing = OUT.read_text(encoding="utf-8") if OUT.exists() else None
    if existing == rendered:
        print(f"gen_python_wrappers: {OUT.relative_to(ROOT)} is current")
        return
    if check:
        sys.exit(
            f"gen_python_wrappers: {OUT.relative_to(ROOT)} is out of date; run "
            "tools/gen_python_wrappers.py"
        )
    OUT.write_text(rendered, encoding="utf-8")
    print(f"gen_python_wrappers: wrote {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
