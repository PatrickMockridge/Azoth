"""Every calc id is wired into every hand-maintained list that has to know it.

Adding a calculation needs five new files - a spec, a kernel and a test in each
language - and then edits to the existing ones that still know an id by hand:
the Rust result type and its `lib.rs` line, the bridge adapter, the Python result
dataclass and the public package. None of those fails at generation time if
forgotten, and this file exists to make the omission fail loudly, naming the file
to go and edit. The registration lists the tranche has generated since - the
extension's id tables, its register call, its transport structs and its wrappers -
are not among them, because a generator that forgot one fails its own `--check`.

# What is checked elsewhere, and deliberately not repeated here

Three of the registration points already have an owner, and duplicating them would
be two places to fix when one changes:

* a spec with no entry in ``result_types()`` — ``test_registry_contract``
* a declared output that is not a field on the result dataclass — same file
* the namespace function existing — same file, from the spec's
  ``implementations.python``
* the Python and Rust field lists agreeing — ``test_cross_impl``

# What only this file checks

Five gaps, each of which was invisible until something called the calc:

* ``azoth._core.calc_ids()`` — the extension's own list. A calc missing from it is
  reported by no other test, because every check that iterates the registry
  iterates the *registry*, not the extension.
* the bridge function existing for a calc id. The bridge used to carry an id table
  whose absence was only reachable through ``resolve``; the table is gone and
  ``resolve`` derives the name, so what can now be missing is the *function*.
* that the bridge accepts the arguments the public API passes it. The two are
  separate hand-written declarations of one contract, so a parameter added to the
  public wrapper and its kernel but not to the bridge is a ``TypeError`` on the
  Rust backend alone — which is what it was, for ``column_diameter`` and
  ``max_allowable_fs_factor`` on the three columns, until this check was added.
* ``azoth._core.pyi`` — the stub. Nothing validated it at all.
* that the result class is a usable dataclass. This one is not hypothetical: the
  Haaland calc was added with its ``@dataclass`` decorator missing, every Rust test
  passed, and the Python path raised ``TypeError: HaalandResult() takes no
  arguments`` on first call. mypy does not catch it either, because the stub is a
  separate declaration the real class is never compared against.
"""

from __future__ import annotations

import ast
import dataclasses
import importlib
import inspect
import re
import sys
from pathlib import Path
from types import ModuleType
from typing import Any

import pytest

from azoth._dispatch import result_types
from azoth._registry_gen import CALCS

REPO_ROOT = Path(__file__).resolve().parents[2]
STUB_PATH = REPO_ROOT / "python" / "src" / "azoth" / "_core.pyi"

CALC_IDS = sorted(calc["id"] for calc in CALCS)
#: Calc id -> the function name it becomes in both languages.
FUNCTION_NAME = {calc_id: calc_id.rpartition(".")[2] for calc_id in CALC_IDS}


def _extension() -> ModuleType:
    """The compiled extension, imported by name.

    The same accessor `test_cross_impl.py` uses, and for the same reason: the
    module does not exist until the bindings are built, and an attribute mypy
    cannot resolve is a worse trade than a lookup that fails clearly.
    """
    try:
        return importlib.import_module("azoth._core")
    except ImportError as exc:  # pragma: no cover - the marker skips these
        raise AssertionError("azoth._core is not built; run `maturin develop`") from exc


def _stub_tree() -> ast.Module:
    """The parsed stub, so names can be read out of it rather than grepped."""
    return ast.parse(STUB_PATH.read_text(encoding="utf-8"), filename=str(STUB_PATH))


def _stub_function_names() -> set[str]:
    """Top-level function names declared in the stub."""
    return {node.name for node in _stub_tree().body if isinstance(node, ast.FunctionDef)}


def _stub_class_names() -> set[str]:
    """Top-level class names declared in the stub."""
    return {node.name for node in _stub_tree().body if isinstance(node, ast.ClassDef)}


def test_there_are_calcs_to_check() -> None:
    """Guard against every parametrized test below passing vacuously."""
    assert CALC_IDS, "no calcs in the registry"
    assert len(CALC_IDS) == len(CALCS)


# --- the result class is a real dataclass ---------------------------------


@pytest.mark.parametrize("calc_id", CALC_IDS)
def test_the_result_type_is_a_usable_dataclass(calc_id: str) -> None:
    """A result class without its decorator is a class you cannot construct.

    ``result_types()`` maps a calc id to a class, and every check that reads those
    classes reads them through ``dataclasses.fields`` - which raises rather than
    reporting a missing decorator. The class then looks registered and fails at the
    first call. This asserts the decorator, which nothing else did.
    """
    result_type: Any = result_types()[calc_id]
    # Read before the assert: `is_dataclass` narrows `result_type` to a dataclass type,
    # which has no `__name__` under mypy --strict even though a class always does.
    class_name: str = result_type.__name__
    assert dataclasses.is_dataclass(result_type), (
        f"{calc_id}: {class_name} in python/src/azoth/core/result.py is not "
        f"a dataclass - check for a missing @dataclass(frozen=True, slots=True, "
        f"eq=False) decorator"
    )
    fields = {f.name for f in dataclasses.fields(result_type)}
    assert "warnings" in fields, f"{calc_id}: {class_name} has no warnings field"


# --- the extension's own list --------------------------------------------


@pytest.mark.requires_rust
def test_the_extension_registers_exactly_the_registry() -> None:
    """``_core.calc_ids()`` is the same set as the generated registry.

    Asserted as an equality in both directions, so a calc dropped from the
    extension is caught alongside one added to the registry and forgotten. The
    generated registry is the source of truth here - it is what the rest of the
    machinery is built from - so this is the check that makes the generated
    ``crates/azoth-python/src/registry_tables_gen.rs`` self-policing.
    """
    extension_ids = set(_extension().calc_ids())
    registry_ids = set(CALC_IDS)
    assert extension_ids == registry_ids, (
        f"azoth._core.calc_ids() and the registry differ\n"
        f"  missing from the extension: {sorted(registry_ids - extension_ids)}\n"
        f"  not in the registry: {sorted(extension_ids - registry_ids)}\n"
        f"Add or remove the spec under specs/calcs/ and run "
        f"`python tools/gen_python_registry.py`, "
        f"register the function in crates/azoth-python/src/lib.rs (add_function), "
        f"and declare it in python/src/azoth/_core.pyi."
    )


@pytest.mark.requires_rust
@pytest.mark.parametrize("calc_id", CALC_IDS)
def test_the_extension_exposes_the_function(calc_id: str) -> None:
    """The extension has a callable under the calc's function name.

    ``calc_ids()`` naming a calc whose function was never registered would give a
    list that promises something the module cannot do.
    """
    function = FUNCTION_NAME[calc_id]
    assert hasattr(_extension(), function), (
        f"{calc_id}: azoth._core has no {function}. Register it in "
        f"crates/azoth-python/src/lib.rs with m.add_function, and declare it in "
        f"python/src/azoth/_core.pyi."
    )


# --- the bridge's id table ----------------------------------------------


def _registered_ids() -> list[str]:
    """Every id `resolve` accepts: the calcs *and* the models.

    `CALC_IDS` above is the calcs alone, because it is the extension's ``calc_ids()``
    table that has to match those. The bridge covers both halves, so the tests that ask
    what the bridge knows read this.
    """
    from azoth._models_gen import MODELS

    return sorted([calc["id"] for calc in CALCS] + [model["id"] for model in MODELS])


@pytest.mark.requires_rust
def test_the_bridge_covers_exactly_the_registry() -> None:
    """Every registered id resolves to a bridge function that exists.

    The bridge used to carry an id -> function table, and this test compared it with
    the registry. Both are gone: `_rust_bridge.resolve` derives the function from the
    id, because a calc's id is its address on the Rust side exactly as it is on the
    reference side. So what can go wrong now is a bridge function that was never
    *written*, and that is what this checks - the derivation is only as good as the
    names it looks up.
    """
    bridge = importlib.import_module("azoth._rust_bridge")
    missing = [
        calc_id for calc_id in _registered_ids() if not hasattr(bridge, calc_id.rpartition(".")[2])
    ]
    assert not missing, (
        f"{missing} have no bridge function. Add one named after the id's last "
        f"segment in python/src/azoth/_rust_bridge.py - there is no table to update."
    )
    with pytest.raises(KeyError):
        bridge.resolve("not.a_calc")


def _public_dispatch_sites() -> list[tuple[str, list[str], str]]:
    """Every ``resolve(...)(...)`` call site in the public packages.

    Returns ``(calc_id, keyword names, file)`` for each, read from the source rather
    than listed. The public wrapper is the caller `_rust_bridge` exists to serve: it
    resolves the id and calls whatever comes back with these keywords, so this *is* the
    signature the bridge has to accept - one hand-written declaration of a contract two
    files state separately.

    The id is written either as the module's own constant (``resolve(_ORIFICE_FLOW)``)
    or as a literal for the one site that has no constant, so both spellings are read.
    Every site is keyword-only today; a positional one would be a contract this cannot
    read, and it is refused rather than skipped.
    """
    sites: list[tuple[str, list[str], str]] = []
    root = REPO_ROOT / "python" / "src" / "azoth"
    for path in sorted(root.glob("*/__init__.py")):
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        constants = {
            node.targets[0].id: node.value.value
            for node in tree.body
            if isinstance(node, ast.Assign)
            and len(node.targets) == 1
            and isinstance(node.targets[0], ast.Name)
            and isinstance(node.value, ast.Constant)
            and isinstance(node.value.value, str)
        }
        for call in ast.walk(tree):
            if not isinstance(call, ast.Call) or not isinstance(call.func, ast.Call):
                continue
            inner = call.func
            if not (isinstance(inner.func, ast.Name) and inner.func.id == "resolve"):
                continue
            (argument,) = inner.args
            calc_id: str
            if isinstance(argument, ast.Constant):
                calc_id = str(argument.value)
            elif isinstance(argument, ast.Name) and argument.id in constants:
                calc_id = str(constants[argument.id])
            else:
                raise AssertionError(
                    f"{path.name}: `resolve({ast.unparse(argument)})` is not a literal or "
                    f"a module-level string constant, so the id it dispatches to cannot "
                    f"be read from this file"
                )
            assert not call.args, (
                f"{path.name}: the call to `resolve({calc_id!r})(...)` passes a positional "
                f"argument, which this check cannot compare against the bridge's own "
                f"parameter list"
            )
            keywords: list[str] = []
            for keyword in call.keywords:
                assert keyword.arg is not None, (
                    f"{path.name}: the call to `resolve({calc_id!r})(...)` unpacks keywords"
                )
                keywords.append(keyword.arg)
            sites.append((calc_id, keywords, path.name))
    return sites


@pytest.mark.requires_rust
def test_the_bridge_accepts_the_arguments_the_public_api_passes() -> None:
    """Every keyword the public wrapper sends, the bridge function takes.

    The two are one contract written twice, and nothing compared them: the bridge is
    called through ``resolve`` with exactly the keywords its caller names, so a
    parameter the wrapper forwards but the bridge does not declare is a ``TypeError``
    on the Rust backend and a working call on the Python one. That is what it was for
    ``column_diameter`` and ``max_allowable_fs_factor``: the capacity-limit port added
    both to the spec, the kernel, the reference and the public wrapper, and the bridge
    was never told, so ``azoth.process.absorption_column`` raised
    ``TypeError: absorption_column() got an unexpected keyword argument
    'column_diameter'`` on every call - reachable from the public API, and reported by
    no test, because no spec case sets either input.

    A count is not enough to catch it. Both are optional in ``_core``'s own signature,
    so the bridge's short call is accepted there and the caller's value is dropped;
    only the *names* separate the two.
    """
    bridge = importlib.import_module("azoth._rust_bridge")
    sites = _public_dispatch_sites()
    assert len(sites) == len(_registered_ids()), (
        f"{len(sites)} public resolve() call sites for {len(_registered_ids())} ids - "
        f"every registered calculation is reachable through the public API, so the two "
        f"counts must agree"
    )

    problems: list[str] = []
    for calc_id, keywords, source in sites:
        function = getattr(bridge, calc_id.rpartition(".")[2], None)
        if function is None:
            problems.append(f"{calc_id}: no bridge function ({source})")
            continue
        accepted = set(inspect.signature(function).parameters)
        missing = [keyword for keyword in keywords if keyword not in accepted]
        if missing:
            problems.append(f"{calc_id}: the bridge does not take {missing} ({source})")
    assert not problems, (
        "the bridge and the public API disagree about the arguments:\n  "
        + "\n  ".join(problems)
        + "\nAdd the parameter to python/src/azoth/_rust_bridge.py and forward it to "
        "azoth._core."
    )


# --- the type stub ------------------------------------------------------


def _rust_pyclasses() -> dict[str, tuple[str, ...]]:
    """Every `#[pyclass]` struct in the binding, by its Python-visible name.

    Read from the Rust *source* rather than from the compiled module, because a
    `#[pyclass]` with no `#[new]` cannot be constructed from Python and its fields
    therefore cannot be introspected - `dir()` on the class shows methods, not
    attributes. The source is the only place the fields are visible, and it is also the
    thing the stub is supposed to describe.

    The struct body is found by **counting braces** rather than by a regex. The first
    version of this used `[^}]*` and silently truncated any struct whose doc comments
    contained a `}` - and one does: `|x_k - x_{k-1}|`, in `PrZFactorResult`'s own
    comment. The fields after it were reported as missing from a struct that has them,
    which is a test that fails for a reason the reader cannot see.
    """
    source_dir = REPO_ROOT / "crates" / "azoth-python" / "src"
    # The whole attribute rather than the part up to the name: `get_all` is written after
    # it, and which of the two forms a struct uses is what decides where its fields are.
    attribute = re.compile(r"#\[pyclass\((?P<attributes>(?:[^()]|\([^()]*\))*?)\)\]", re.DOTALL)
    name = re.compile(r'name\s*=\s*"(?P<name>\w+)"')
    struct = re.compile(r"pub struct \w+ \{")

    found: dict[str, tuple[str, ...]] = {}
    for path in sorted(source_dir.glob("*.rs")):
        text = path.read_text(encoding="utf-8")
        for match in attribute.finditer(text):
            declared = name.search(match.group("attributes"))
            if declared is None:
                continue
            opening = struct.search(text, match.end())
            if opening is None:  # pragma: no cover - a pyclass with no struct
                continue
            depth = 0
            end = opening.end() - 1
            for index in range(end, len(text)):
                if text[index] == "{":
                    depth += 1
                elif text[index] == "}":
                    depth -= 1
                    if depth == 0:
                        end = index
                        break
            body = text[opening.end() : end]
            if re.search(r"\bget_all\b", match.group("attributes")):
                # `get_all` gives every field a getter, so every field is one - private
                # fields included, which is the one place this reads further than the
                # `#[pyo3(get)]` form does.
                fields = re.findall(r"^\s*(?:pub )?(\w+):", body, re.MULTILINE)
            else:
                fields = re.findall(r"#\[pyo3\(get\)\]\s*pub (\w+):", body)
            found[declared.group("name")] = tuple(fields)
    return found


def test_the_stub_matches_the_rust_transport_types() -> None:
    """Every field the stub declares for a transport type exists in the Rust struct.

    **This test exists because the drift it catches had already happened.** `FittingRow`
    declared `source_ref` and `source_locator` in the stub for as long as it took someone
    to read it, and both fields had been deleted from the Rust struct and from the
    project. The completeness test checked that the stub declares each *name*; nothing
    checked a field, so a stub describing attributes no object has type-checked clean and
    failed at the first attribute access.

    **Both directions.** A field the stub declares and the Rust does not have was the
    failure that shipped; a field the Rust carries with `#[pyo3(get)]` and the stub omits
    is the same defect one file over - a caller reaching for it type-checks clean and
    fails at the first attribute access. The parser reads only the `#[pyo3(get)]` fields,
    so a private Rust field is still allowed to be absent from the stub.
    """
    rust = _rust_pyclasses()
    assert rust, "the Rust source parser found no pyclasses, so it is broken not empty"

    stub = _stub_tree()
    checked = 0
    for node in stub.body:
        if not isinstance(node, ast.ClassDef) or node.name not in rust:
            continue
        declared = {
            statement.target.id
            for statement in node.body
            if isinstance(statement, ast.AnnAssign) and isinstance(statement.target, ast.Name)
        }
        missing = sorted(declared - set(rust[node.name]))
        assert not missing, (
            f"_core.pyi declares {node.name}.{missing}, which the Rust struct does not "
            f"have. The stub is describing attributes no object carries - run "
            f"`python tools/gen_stub.py` if the field was removed, or add it to the Rust "
            f"struct with `#[pyo3(get)]`."
        )
        undescribed = sorted(set(rust[node.name]) - declared)
        assert not undescribed, (
            f"_core.pyi does not declare {node.name}.{undescribed}, which the Rust struct "
            f"carries with `#[pyo3(get)]`. A caller reaching for one type-checks clean and "
            f"fails at the first attribute access - add it to `gen_stub.py`'s declaration "
            f"of {node.name} and run `python tools/gen_stub.py`."
        )
        checked += 1

    # A guard against the regex silently matching nothing, which would make every
    # assertion above vacuous and this test pass while checking nothing at all.
    assert checked >= 8, (
        f"only {checked} transport type(s) were compared against the Rust source, but "
        f"the stub declares more. The parser is matching less than it should."
    )


#: A Rust transport type that carries a quantity, and how the stub spells it.
#:
#: **Only these are compared, and the reason is the one this file keeps meeting.** A dimensioned
#: field is a `PyQty` in the transport and a `Qty` in the stub, and a *dimensionless* one is a
#: bare `f64` and a `float` - and in Rust the two are the same number, because the kernels are
#: typed against `uom::si::f64` where an alias and a float differ only in spelling. So a field
#: can be a quantity on one side of the wire and a number on the other without anything failing
#: to compile, which is exactly what happened: fourteen vector fields were `Vec<PyQty>` in the
#: transport and `list[float]` in the stub, the bridge read `.magnitude_si` off them, and mypy
#: was the only thing that noticed - sixty errors' worth, in code that was right.
_QUANTITY_STUB_TYPES = {
    "PyQty": "Qty",
    "Option<PyQty>": "Qty | None",
    "Vec<PyQty>": "list[Qty]",
    "Option<Vec<PyQty>>": "list[Qty] | None",
    "Vec<Vec<PyQty>>": "list[list[Qty]]",
}

TRANSPORT_GEN = REPO_ROOT / "crates" / "azoth-python" / "src" / "transport_gen.rs"


def _tools_module(name: str) -> ModuleType:
    """A module out of `tools/`, which is not a package and so is not importable by name."""
    sys.path.insert(0, str(REPO_ROOT / "tools"))
    try:
        return importlib.import_module(name)
    finally:
        sys.path.pop(0)


def _transport_field_types() -> dict[str, dict[str, str]]:
    """The generated transport's fields, by class name without the `Py` prefix.

    Read from `transport_gen.rs` rather than from `gen_stub.py`'s own rule, because the
    transport is what an attribute access actually reaches: `_core.<calc>(...).field` is a
    `PyQty` or a number according to this file, and the stub is a description of it.
    """
    rust_index = _tools_module("rust_index")
    text = rust_index._blank_comments(TRANSPORT_GEN.read_text(encoding="utf-8"))
    out: dict[str, dict[str, str]] = {}
    for block in re.finditer(r"^pub struct (Py\w+) \{(.*?)^\}", text, re.M | re.S):
        raw, body = block.group(1), block.group(2)
        out[raw.removeprefix("Py")] = {
            field.group(1): field.group(2)
            for field in re.finditer(r"^\s*pub (\w+): (.+),$", body, re.M)
        }
    return out


def test_the_stub_types_the_quantity_fields_the_transport_carries() -> None:
    """Every quantity the transport carries is a `Qty` in the stub, not a number.

    The sibling above compares the two files' *field names* and was written because a stub
    describing an attribute no object has type-checks clean. This is the same failure one
    column over: a stub describing an attribute as a `float` that is a quantity type-checks
    clean too, and then every read of `.magnitude_si` on it is an error.
    """
    transport = _transport_field_types()
    assert transport, f"{TRANSPORT_GEN.name} parsed to no structs"

    stubs = {
        node.name: {
            statement.target.id: ast.unparse(statement.annotation)
            for statement in node.body
            if isinstance(statement, ast.AnnAssign) and isinstance(statement.target, ast.Name)
        }
        for node in _stub_tree().body
        if isinstance(node, ast.ClassDef)
    }

    problems: list[str] = []
    checked = 0
    for class_name, fields in transport.items():
        declared = stubs.get(class_name)
        if declared is None:
            continue
        for field, rust_type in fields.items():
            expected = _QUANTITY_STUB_TYPES.get(rust_type)
            if expected is None:
                continue
            checked += 1
            if declared.get(field) != expected:
                problems.append(
                    f"{class_name}.{field}: the transport carries {rust_type}, so the stub "
                    f"must declare {expected}, and it declares {declared.get(field)!r}"
                )

    # The same guard the sibling carries: a mapping that stopped matching would make every
    # assertion below vacuous, and the count says whether it is still matching.
    assert checked >= 400, (
        f"only {checked} quantity field(s) were compared, but the transport carries far "
        f"more. The parser or the type mapping is matching less than it should."
    )
    assert not problems, (
        "\n  ".join(["_core.pyi and the transport disagree about a field's type:", *problems])
        + "\nRun `python tools/gen_stub.py` after correcting the annotation in "
        "python/src/azoth/core/result.py - the stub is emitted from it."
    )


@pytest.mark.parametrize("calc_id", CALC_IDS)
def test_the_stub_declares_the_calc_function(calc_id: str) -> None:
    """The stub declares the extension function for every calc.

    The stub is a separate declaration that nothing compares against the compiled
    module, so an omission is invisible: the code that depends on it type-checks
    against the extension's ``.so`` at runtime and against nothing at all in mypy.
    """
    function = FUNCTION_NAME[calc_id]
    assert function in _stub_function_names(), (
        f"{calc_id}: python/src/azoth/_core.pyi does not declare {function}(). The "
        f"stub has to list every function the extension exposes or the typed path "
        f"silently has no signature for it."
    )


@pytest.mark.parametrize("calc_id", CALC_IDS)
def test_the_stub_declares_the_result_class(calc_id: str) -> None:
    """The stub declares the result class for every calc.

    Same reasoning as the function check, and the same failure: an undeclared
    result class means the transport object a Rust path returns is untyped.
    """
    result_type: Any = result_types()[calc_id]
    class_name: str = result_type.__name__
    assert class_name in _stub_class_names(), (
        f"{calc_id}: python/src/azoth/_core.pyi does not declare "
        f"{class_name}. Add it alongside the other result classes."
    )


def test_the_stub_parses() -> None:
    """The stub is syntactically valid.

    Worth asserting separately because every check above reads the stub through
    ``ast.parse``: if it stopped parsing, each of those would raise inside ``ast``
    rather than reporting the real problem.
    """
    assert _stub_tree().body, "python/src/azoth/_core.pyi parsed to nothing"
