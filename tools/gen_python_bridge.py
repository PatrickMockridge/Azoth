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
BINDING = ROOT / "crates" / "azoth-python" / "src"
TRANSPORT = BINDING / "transport_gen.rs"
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

#: A parameter the caller may omit and the extension's own type does not accept as absent.
#:
#: The two are declarations of different things - "this argument is optional" and "this function
#: takes a string" - and the hand-written adapters bridged them with a per-site literal, because
#: the default is the class's own and not the first of the input's declared `values`. An id with
#: this shape and no entry here is refused rather than sent a `None` the extension will reject.
ARGUMENT_DEFAULTS: dict[tuple[str, str], str] = {
    ("reactions.reactive_ph_flash", "cubic"): '"srk"',
    ("reactions.reactive_tp_flash", "cubic"): '"srk"',
}

#: The parameters whose boundary conversion is the site's own rather than the spec's.
#:
#: **One entry, and it is the reason this is a table and not silent generation.** `tray_hydraulics`
#: declares `hole_diameter` in millimetres and `input_to_si` puts it in metres - which is what
#: every other input wants - but the two correlations that read the hole in NeqSim are written in
#: millimetres, so the site scales it back and the comment says so. Generating the argument without
#: it is a factor of a thousand on `downcommer_backup`, which is what a first run of this file
#: produced and what `test_cross_impl` refused. A second entry would need its own reason here.
ARGUMENT_SCALES: dict[tuple[str, str], tuple[str, str]] = {
    ("hydraulics.tray_hydraulics", "hole_diameter"): (
        " * 1000.0",
        "millimetres on both sides of the wire: the spec declares them and the correlations "
        "that read the hole are written in them",
    ),
}


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


def wrapper_types() -> dict[str, dict[str, str]]:
    """The extension's own parameters, by wrapper name, as the Rust types pyo3 publishes.

    **This is the authority for what the call must hand over**, and it is not the stub and not the
    public wrapper: a count is `f64` in the spec's vocabulary and `usize` in the `#[pyfunction]`, so
    a float reaches Rust as a float and is refused as an integer. `gen_stub.py` renders the stub
    from the spec, so it says `float` for the same parameter.
    """
    out: dict[str, dict[str, str]] = {}
    for name in (
        "wrappers_gen",
        "process",
        "eos",
        "hydraulics",
        "thermal",
        "reactions",
        "standards",
    ):
        path = BINDING / f"{name}.rs"
        if not path.exists():
            continue
        text = _blank(path.read_text(encoding="utf-8"))
        for match in re.finditer(
            # `\)` and not `\n\)`: a short wrapper is written on one line, and requiring the
            # newline made the non-greedy group run past it to the next function's close.
            r"pub fn (\w+)\(\s*(?:py: Python<'_>,\s*)?(.*?)\)\s*->",
            text,
            re.S,
        ):
            out[match.group(1)] = {
                field.group(1): field.group(2)
                for line in match.group(2).splitlines()
                if (field := re.match(r"\s*(\w+): (.+?),?$", line))
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


def argument_expr(name: str, declaration: dict[str, Any], annotation: str) -> str:
    """The argument expression for one `_core` parameter, from the spec's declaration of it.

    Three questions in this order: is it a container, does it carry a unit, and may it be absent.
    A bare number and a dimensionless quantity are the same thing to the extension, so the
    dimensionless case crosses unchanged; a dimensioned one goes through `input_to_si`, which reads
    the unit *and* the `interval` flag from the spec rather than restating them - a call site that
    writes `to_si(dT, "K", "dT")` is a flag the spec holds and the code ignores.
    """
    kind = declaration.get("type", "quantity")
    unit = declaration.get("unit")
    dimensioned = unit is not None and unit != "dimensionless"
    optional = bool(declaration.get("optional"))

    if kind == "matrix":
        body = _matrix(name, dimensioned, annotation)
    elif kind == "vector":
        body = _vector(name, dimensioned, annotation)
    elif kind == "components":
        # The extension takes an owned `Vec<String>`; everything else that is not a number is
        # already the thing the signature declares, an enum and a boolean included.
        body = f"list({name})"
    elif kind in ("string", "enum", "boolean") or not dimensioned:
        body = name
    else:
        body = f'input_to_si(spec, "{name}", {name})'
    if optional and kind not in ("string", "enum", "boolean"):
        return f"None if {name} is None else {body}"
    return body


def _vector(name: str, dimensioned: bool, annotation: str) -> str:
    """A vector crosses element by element, and the element's unit is the annotation's.

    `Sequence[Q]` is a vector of quantities and each entry converts; a vector of bare numbers
    is already SI magnitudes, which is what `_si` is - the tolerant half of the same function,
    which passes a number through and converts anything else.
    """
    if not dimensioned:
        return f"list({name})"
    if "Q" in annotation:
        return f'[input_to_si(spec, "{name}", v) for v in {name}]'
    return f'[_si(spec, "{name}", v) for v in {name}]'


def _matrix(name: str, dimensioned: bool, annotation: str) -> str:
    """A matrix is a `Vec<Vec<f64>>` either way; only a quantity's entries need converting."""
    if not dimensioned or "Q" not in annotation:
        return name
    return f'[[_si(spec, "{name}", v) for v in row] for row in {name}]'


#: The Rust parameter types that are counts, and so take an `int` on this side.
INTEGERS = ("usize", "u32", "u64", "u8", "i32", "i64", "u16")


def coerce(expression: str, rust_type: str | None, parameter: str) -> str:
    """The expression as the type pyo3 publishes for the parameter.

    A spec says an input is a number; the `#[pyfunction]` says whether that number is a *count*.
    `max_phases`, `number_of_stages` and a solver's `max_iterations` are `usize` in the signature
    the extension publishes and `f64` in the spec's vocabulary, and a float reaches Rust as a float
    and is refused as an integer - a `TypeError` at the boundary rather than a wrong answer.
    """
    optional = bool(rust_type and rust_type.startswith("Option<") and rust_type.endswith(">"))
    inner = rust_type[len("Option<") : -1] if optional else rust_type
    if inner not in INTEGERS or expression.startswith("int("):
        return expression
    # An optional count is absent or a count: `None if x is None else int(x)`. The guard is added
    # only where the expression *is* the parameter, because a dimensioned one already carries its
    # own - `None if x is None else input_to_si(...)` - and wrapping that in `int` would convert
    # the `None` branch.
    if optional:
        # The guard is already there - `argument_expr` writes one for anything that may be absent -
        # so the count is taken inside its else-branch rather than around the whole expression,
        # which would call `int(None)` on the branch that means "absent".
        guarded = re.fullmatch(rf"None if {re.escape(parameter)} is None else (.+)", expression)
        if guarded is None:
            raise Refusal(
                f"{parameter}: an optional count whose expression is `{expression}`, which this "
                f"file cannot put the conversion inside"
            )
        return f"None if {parameter} is None else int({guarded.group(1)})"
    return f"int({expression})"


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
    _, result_class = core_signatures()[function]
    # **The call's order is the extension's own, read from the `#[pyfunction]` itself.** The stub
    # describes that function and carries the *spec's* order, and for three ids the two disagree:
    # `reactions.equilibrium_constant` declares `(reaction, source, T)` in the spec and takes
    # `(source, reaction, T)`, so a positional call built from the stub swaps two strings.
    types = wrapper_types().get(function, {})
    core_params = list(types)
    if set(params) != set(core_params):
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
    # **The call is in `_core`'s order and the signature in the public wrapper's.** They are not
    # always the same: `absorption_column` declares its two component lists first and the public
    # wrapper states them where the caller expects them. Both orders are declarations, so neither
    # is copied into the other - the call follows the extension and the signature the caller.
    arguments: list[tuple[str, str | None]] = []
    for p in core_params:
        argument = coerce(argument_expr(p, declared[p], annotations.get(p) or ""), types.get(p), p)
        if (
            "None" in (annotations.get(p) or "")
            and types.get(p) is not None
            and not types[p].startswith("Option<")
        ):
            default = ARGUMENT_DEFAULTS.get((calc_id, p))
            if default is None:
                raise Refusal(
                    f"{p}: the signature allows it to be absent and the extension takes "
                    f"{types[p]}, with no default declared for this id"
                )
            argument = f"{argument} if {p} is not None else {default}"
        scale = ARGUMENT_SCALES.get((calc_id, p))
        # The reason travels with the argument, so a reader of the generated file sees why this
        # one is not the spec's conversion without opening this one.
        arguments.append((argument + scale[0], scale[1]) if scale else (argument, None))
    needs_spec = any("spec," in expression for expression, _ in arguments)
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
    if any("_si(" in expression for expression, _ in arguments):
        used.add("_si")
    if "warnings" in fields:
        used.add("_warnings")

    # `ast.unparse` of an `arguments` node is the list *without* its parentheses, so they are
    # put back here - the defaults, the annotations and any `*` marker all survive the trip.
    lines = [
        f"def {function}({ast.unparse(public.args)}) -> {ast.unparse(public.returns)}:",
        f'    """``{calc_id}``, computed in Rust."""',
    ]
    if needs_spec:
        lines.append(f"    spec = {lookup}")
    lines.append(f"    result = _core.{function}(")
    lines += [f"        {a},{f'  # {why}' if why else ''}" for a, why in arguments]
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

    # `_warnings` and `_si` come from the support module rather than from `azoth.core.result`,
    # and isort wants them named in the block they come from.
    support = sorted(name for name in used if name in ("_si", "_warnings"))
    imports = [
        "from collections.abc import Sequence",
        "",
        "from azoth import _core, _models_gen",
        "from azoth._registry_gen import spec as _spec_for",
        f"from azoth.core.bridge_support import {', '.join(support)}",
        "from azoth.core.result import (",
    ]
    # isort compares names case-insensitively, so `HeaterResult` sorts before
    # `HeatOfVaporizationResult` and a plain `sorted()` would emit a block ruff rejects.
    imports += [f"    {name}," for name in sorted(used - set(plain) - set(support), key=str.lower)]
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
