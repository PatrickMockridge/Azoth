#!/usr/bin/env python3
"""Generate the Rust backend's Python functions.

`python/src/azoth/_rust_bridge.py` carries one function per registered id, and each is the same
three statements: look the spec up, call the extension, rebuild the Python result dataclass from
the transport object. Nothing in it is arithmetic - the arithmetic is in Rust - so every line of it
is a function of declarations that already exist:

* **the signature** is the public wrapper's, because the wrapper is the caller: it does
  `resolve(...)(**kwargs)` and the bridge is what comes back. `python/src/azoth/{eos,hydraulics,
  process,reactions,standards,thermal}/__init__.py` states it once, and that is the contract.
* **the call's arguments** are `_core.pyi`'s parameter list, in its order, each converted
  according to the kind and unit the spec declares for it.
* **the result** is the transport struct's fields, in `CalcResult::FIELDS` order, each rebuilt
  according to its Rust transport type.

**Why this one matters more than the others.** A hand-written bridge is a second declaration of a
contract stated elsewhere, and nothing compared the two. Three of these functions were an argument
short for four commits - `absorption_column`, `stripping_column` and `packed_column` never received
`column_diameter` or `max_allowable_fs_factor`, so every call through `azoth.process.*` raised
`TypeError` on the Rust backend while the Python backend answered normally. A count does not catch
that: both parameters are optional in `_core`'s own signature, so the short call was accepted and
the caller's value silently dropped.

**Coverage is read, not listed.** An id is covered when the public wrapper's parameter names are
`_core`'s, in that order - which is the rule the whole file is built on, since the call is
positional - and when every parameter is an input the spec declares. A kernel whose boundary is a
mixture or a record the spec does not name fails both and is left hand-written.

    python tools/gen_python_bridge.py            # write
    python tools/gen_python_bridge.py --check    # fail if the tree is not what this emits
    python tools/gen_python_bridge.py --survey   # what it covers, and what it leaves
"""

from __future__ import annotations

import ast
import re
import sys
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parent))
import rust_index

ROOT = Path(__file__).resolve().parent.parent
PY_SRC = ROOT / "python" / "src" / "azoth"
OUT = PY_SRC / "_rust_bridge_gen.py"
TRANSPORT = ROOT / "crates" / "azoth-python" / "src" / "transport_gen.rs"
CORE_PYI = PY_SRC / "_core.pyi"

#: The packages whose `__init__.py` holds the public wrapper for a registered id.
NAMESPACES = ("eos", "hydraulics", "process", "reactions", "standards", "thermal")

#: The Python class an enum field is rebuilt as, where the Rust type's own name is not the name
#: the module imported. Four of the ten enums are spelled with a leading underscore and six are
#: not, and nothing derives the difference - they are names in a module.
ENUM_ALIASES = {
    "Phase": "_Phase",
    "StabilityVerdict": "_StabilityVerdict",
    "TpMultiflashSeed": "_TpMultiflashSeed",
    "HenryStatus": "_HenryStatus",
}

#: The one enum that crosses as its own spelling rather than as the enum: the spec declares
#: `dataset` an enum and `PitzerPhaseResult.dataset` is annotated `str`, so the Python field is a
#: string and a caller reads the name.
ENUMS_THAT_STAY_STRINGS = {"PitzerDataset"}


class Refusal(Exception):
    """This file cannot derive the function, so it stays hand-written."""


def _blank(text: str) -> str:
    return rust_index._blank_comments(text)


def public_signatures() -> dict[str, ast.FunctionDef]:
    """Every public wrapper, by the name it is called by.

    `hydraulics.tray_hydraulics` is `azoth.hydraulics.tray_hydraulics`, so the id's last segment
    is the key. Read as an AST node rather than as text, because the header is emitted verbatim
    and a default has to survive the trip.
    """
    out: dict[str, ast.FunctionDef] = {}
    for namespace in NAMESPACES:
        path = PY_SRC / namespace / "__init__.py"
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        for node in tree.body:
            if isinstance(node, ast.FunctionDef):
                out.setdefault(node.name, node)
    return out


def core_signatures() -> dict[str, tuple[list[str], str]]:
    """`_core.pyi`: each function's parameters, and the result class it answers with."""
    tree = ast.parse(CORE_PYI.read_text(encoding="utf-8"), filename=str(CORE_PYI))
    out: dict[str, tuple[list[str], str]] = {}
    for node in tree.body:
        if not isinstance(node, ast.FunctionDef):
            continue
        params = node.args.posonlyargs + node.args.args + node.args.kwonlyargs
        out[node.name] = ([p.arg for p in params], ast.unparse(node.returns))
    return out


def transport_fields() -> dict[str, dict[str, str]]:
    """`transport_gen.rs`: field name -> Rust transport type, by class name without the `Py`."""
    text = _blank(TRANSPORT.read_text(encoding="utf-8"))
    out: dict[str, dict[str, str]] = {}
    for block in re.finditer(r"^pub struct (Py\w+) \{(.*?)^\}", text, re.M | re.S):
        out[block.group(1)[2:]] = {
            field.group(1): field.group(2)
            for field in re.finditer(r"^\s*pub (\w+): (.+),$", block.group(2), re.M)
        }
    return out


def rust_result(calc_id: str) -> rust_index.ResultType:
    """The Rust result struct for an id. Carries its name, and its fields' types."""
    for result in rust_index.result_types():
        if result.calc_id == calc_id:
            return result
    raise Refusal(f"{calc_id} has no result type")


def rust_result_fields(calc_id: str) -> dict[str, str]:
    """The Rust result struct's own field types, by the name that field has in `FIELDS` order.

    Read for one thing: `String` in the transport is both an enum and a string, and the Rust type
    is what says which - and names the enum's Python class.
    """
    result = rust_result(calc_id)
    return {
        rust: kind for (_, rust), kind in zip(result.fields, result.rust_field_types, strict=False)
    }


_SPECS: dict[str, dict[str, Any]] = {}


def specs() -> dict[str, dict[str, Any]]:
    """Every spec, by id - calcs and models both."""
    if not _SPECS:
        from azoth._models_gen import MODELS
        from azoth._registry_gen import CALCS

        _SPECS.update({entry["id"]: entry for entry in [*CALCS, *MODELS]})
    return _SPECS


# --- the return: one expression per transport type --------------------------


def _enum_class(rust: str | None) -> str:
    """The Python class for a Rust enum field.

    The Rust type is the *struct's*, so an optional enum reads `Option<Phase>` and the wrapper has
    to be peeled before the name is a name.
    """
    if rust is None:
        raise Refusal("an enum field whose Rust type does not name it")
    for wrapper in ("Option<", "Vec<"):
        if rust.startswith(wrapper) and rust.endswith(">"):
            rust = rust[len(wrapper) : -1]
    return ENUM_ALIASES.get(rust, rust)


def return_expr(field: str, transport: str, rust: str | None, declared: str | None) -> str:
    """The right-hand side that rebuilds one result field from the transport object.

    **The rule is the transport's Rust type**, which is why this table is short and total: a
    quantity is a `PyQty` carrying its SI magnitude and the unit it is stated in, a plain number is
    the number, and a warning is a warning. `String` is the only case needing a second look - an
    enum crosses as its own spelling and a plain string as itself, and the spec's `type` says which.
    """
    src = f"result.{field}"
    if transport == "PyQty":
        return f"from_si({src}.magnitude_si, {src}.unit)"
    if transport == "Option<PyQty>":
        return f"None if {src} is None else from_si({src}.magnitude_si, {src}.unit)"
    if transport == "Vec<PyQty>":
        return f"tuple(from_si(q.magnitude_si, q.unit) for q in {src})"
    if transport == "Vec<Vec<PyQty>>":
        return f"tuple(tuple(from_si(v.magnitude_si, v.unit) for v in row) for row in {src})"
    if transport == "Vec<PyWarning>":
        return f"_warnings({src})"
    if transport == "Vec<PyKComponent>":
        return f"tuple(KComponent(fitting_id=c.fitting_id, n_ld=c.n_ld, k=c.k) for c in {src})"
    if transport in ("Vec<f64>", "Vec<String>", "Vec<u32>"):
        return f"tuple({src})"
    if transport == "Vec<Vec<f64>>":
        return f"tuple(tuple(row) for row in {src})"
    if transport in ("f64", "u32", "u64", "usize", "i64", "bool", "Option<f64>"):
        return src
    if transport == "String":
        if declared != "enum" or rust in ENUMS_THAT_STAY_STRINGS:
            return src
        return f"{_enum_class(rust)}({src})"
    if transport == "Option<String>":
        if declared == "enum":
            return f"None if {src} is None else {_enum_class(rust)}({src})"
        return src
    raise Refusal(f"{field}: the transport type {transport} has no rule")


# --- the arguments: one expression per declared input ----------------------


def argument_expr(name: str, declaration: dict[str, Any]) -> str:
    """The argument expression for one `_core` parameter, from the spec's declaration of it.

    A bare number and a dimensionless quantity are the same thing to the extension, so the
    dimensionless case crosses unchanged; a dimensioned one goes through `input_to_si`, which reads
    the unit *and* the `interval` flag from the spec rather than restating them - a call site that
    writes `to_si(dT, "K", "dT")` is a flag the spec holds and the code ignores.
    """
    kind = declaration.get("type", "quantity")
    unit = declaration.get("unit")
    if kind == "vector":
        if unit is None or unit == "dimensionless":
            return f"list({name})"
        raise Refusal(f"{name}: a unit-carrying vector has no rule yet")
    if kind == "matrix":
        raise Refusal(f"{name}: a matrix input has no rule yet")
    if kind == "components":
        # The extension takes an owned `Vec<String>`; everything else that is not a number is
        # already the thing the signature declares, an enum and a boolean included.
        return f"list({name})"
    if kind in ("string", "enum", "boolean"):
        return name
    if unit is None or unit == "dimensionless":
        return name
    if declaration.get("optional"):
        return f'None if {name} is None else input_to_si(spec, "{name}", {name})'
    return f'input_to_si(spec, "{name}", {name})'


def coerce(expression: str, annotation: str | None) -> str:
    """The expression as the annotation the extension declares for it.

    A spec says an input is a number; the kernel says whether that number is a count. `models`,
    `number_of_stages` and a solver's `max_iterations` are `int` in the signature pyo3 publishes
    and `f64` in the spec's vocabulary, and a float reaches Rust as a `float` and is refused as an
    integer - which is a `TypeError` at the boundary rather than a wrong answer.
    """
    if annotation == "int" and not expression.startswith("int("):
        return f"int({expression})"
    if annotation == "float" and not expression.startswith("float("):
        return f"float({expression})"
    if annotation == "bool":
        return f"bool({expression})"
    return expression


# --- the whole function -----------------------------------------------------


def emit_function(calc_id: str, used: set[str]) -> str:
    """One bridge function, as source, or a `Refusal` naming what stopped it.

    `used` collects the names the function needs imported, so the module's import block is
    the union over the functions rather than a second list kept beside them.
    """
    function = calc_id.rpartition(".")[2]
    public = public_signatures().get(function)
    if public is None:
        raise Refusal(f"no public wrapper named {function}")
    if public.args.vararg or public.args.kwarg:
        raise Refusal("the public wrapper takes *args")

    params = [p.arg for p in public.args.posonlyargs + public.args.args]
    params += [p.arg for p in public.args.kwonlyargs]
    core_params, result_class = core_signatures()[function]
    if params != core_params:
        raise Refusal(
            f"the public wrapper's parameters are not _core's: "
            f"+{[p for p in params if p not in core_params]} "
            f"-{[p for p in core_params if p not in params]}"
        )

    spec = specs()[calc_id]
    declared = spec.get("inputs", {})
    undeclared = [p for p in params if p not in declared]
    if undeclared:
        raise Refusal(f"{undeclared} is not an input the spec declares")

    annotations = {
        p.arg: (ast.unparse(p.annotation) if p.annotation else None)
        for p in public.args.posonlyargs + public.args.args + public.args.kwonlyargs
    }
    arguments = [coerce(argument_expr(p, declared[p]), annotations.get(p)) for p in params]
    needs_spec = any("spec," in a for a in arguments)
    lookup = (
        f'_spec_for("{calc_id}")'
        if calc_id in {e["id"] for e in _calcs()}
        else f'_models_gen.model("{calc_id}")'
    )

    # **The transport's class name is the Rust struct's, which is not always the Python one**:
    # `eos.tp_flash_saft` answers a `SaftFlashResult` in Rust and a `TpFlashSaftResult` here, so
    # the Python name comes from the stub and the transport is looked up by the Rust name.
    fields = transport_fields()[rust_result(calc_id).rust_name]
    rust = rust_result_fields(calc_id)
    outputs = spec.get("outputs", {})
    returned = [
        (field, return_expr(field, kind, rust.get(field), (outputs.get(field) or {}).get("type")))
        for field, kind in fields.items()
    ]
    used.add(result_class)
    # Read from the same decision `return_expr` made rather than from the emitted text: a scan for
    # `Name(` would import `Phase` for a `_Phase(` the module never declares.
    for field, kind in fields.items():
        if kind != "String" or (outputs.get(field) or {}).get("type") != "enum":
            continue
        if rust.get(field) not in ENUMS_THAT_STAY_STRINGS:
            used.add(_enum_class(rust.get(field)))
    if any("KComponent(" in expression for _, expression in returned):
        used.add("KComponent")

    # `ast.unparse` of an `arguments` node is the list *without* its parentheses, so they are
    # put back here - the defaults, the annotations and any `*` marker all survive the trip.
    lines = [
        f"def {function}({ast.unparse(public.args)}) -> {ast.unparse(public.returns)}:",
        f'    """``{calc_id}``, computed in Rust."""',
    ]
    if needs_spec:
        lines.append(f"    spec = {lookup}")
    lines.append(f"    result = _core.{function}(")
    lines += [f"        {a}," for a in arguments]
    lines.append("    )")
    lines.append(f"    return {result_class}(")
    lines += [f"        {field}={expr}," for field, expr in returned]
    lines.append("    )")
    return "\n".join(lines)


def _calcs() -> list[dict[str, Any]]:
    from azoth._registry_gen import CALCS

    return list(CALCS)


def covered() -> tuple[list[tuple[str, str, set[str]]], list[tuple[str, str]]]:
    """`(covered, refused)`, and each covered entry carries the names it needs imported."""
    out, left = [], []
    for calc_id, _, _ in rust_index.implementations():
        used: set[str] = set()
        try:
            source = emit_function(calc_id, used)
        except Refusal as refusal:
            left.append((calc_id, str(refusal)))
            continue
        out.append((calc_id, source, used))
    return out, left


def emit() -> str:
    """`_rust_bridge_gen.py`: every id this file can derive, in id order."""
    out, _ = covered()
    bodies = "\n\n\n".join(source for _, source, _ in out)

    used: set[str] = set().union(*(names for _, _, names in out)) if out else set()
    # The four aliased enums are imported under their other name, so the alias is what appears in
    # the block and `from ... import X as _X` is what binds it.
    plain = {alias: rust for rust, alias in ENUM_ALIASES.items() if alias in used}

    imports = [
        "from collections.abc import Sequence",
        "",
        "from azoth import _core, _models_gen",
        "from azoth._registry_gen import spec as _spec_for",
        "from azoth.core.bridge_support import _warnings",
        "from azoth.core.result import (",
    ]
    # isort compares names case-insensitively, so `HeaterResult` sorts before
    # `HeatOfVaporizationResult` and a plain `sorted()` would emit a block ruff rejects.
    imports += [f"    {name}," for name in sorted(used - set(plain), key=str.lower)]
    imports += [")"]
    # One block each, which is what isort settles on for a name imported under another name.
    for alias, rust in sorted(plain.items(), key=lambda item: item[1].lower()):
        imports += ["from azoth.core.result import (", f"    {rust} as {alias},", ")"]
    imports.append("from azoth.core.units import Q, from_si, input_to_si")

    header = (
        '"""The Rust backend\'s functions, emitted from the extension\'s own declarations.\n'
        "\n"
        "Generated by `tools/gen_python_bridge.py` - see its module documentation for the three\n"
        "declarations each one is a function of. **Do not edit by hand.**\n"
        "\n"
        "`python/src/azoth/_rust_bridge.py` re-exports these and keeps the ids this file cannot\n"
        "derive, whose boundary is a mixture or a record the spec does not name.\n"
        '"""\n'
        "\n"
        "from __future__ import annotations\n"
        "\n" + "\n".join(imports) + "\n"
        "\n"
        "__all__ = [\n"
    )
    names = sorted(calc_id.rpartition(".")[2] for calc_id, _, _ in out)
    header += "".join(f'    "{name}",\n' for name in names) + "]\n\n\n"
    return header + bodies + "\n"


def main() -> None:
    if "--survey" in sys.argv:
        out, left = covered()
        print(f"covered: {len(out)} of {len(out) + len(left)}")
        for calc_id, why in left:
            print(f"   left hand-written: {calc_id}  ({why})")
        return
    check = "--check" in sys.argv
    rendered = emit()
    existing = OUT.read_text(encoding="utf-8") if OUT.exists() else None
    if existing == rendered:
        print(f"gen_python_bridge: {OUT.relative_to(ROOT)} is current")
        return
    if check:
        sys.exit(f"gen_python_bridge: {OUT.relative_to(ROOT)} is out of date")
    OUT.write_text(rendered, encoding="utf-8")
    print(f"gen_python_bridge: wrote {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
