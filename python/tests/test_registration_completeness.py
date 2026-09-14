"""Every calc id is wired into every hand-maintained list that has to know it.

Adding a calculation needs five new files, and then a dozen edits to existing
ones - none of which fails at generation time if forgotten. This file exists to
make those omissions fail loudly, naming the file to go and edit.

# What is checked elsewhere, and deliberately not repeated here

Three of the registration points already have an owner, and duplicating them would
be two places to fix when one changes:

* a spec with no entry in ``result_types()`` — ``test_registry_contract``
* a declared output that is not a field on the result dataclass — same file
* the namespace function existing — same file, from the spec's
  ``implementations.python``
* the Python and Rust field lists agreeing — ``test_cross_impl``

# What only this file checks

Four gaps, each of which was invisible until something called the calc:

* ``azoth._core.calc_ids()`` — the extension's own list. A calc missing from it is
  reported by no other test, because every check that iterates the registry
  iterates the *registry*, not the extension.
* ``azoth._rust_bridge._IMPLEMENTATIONS`` — the bridge's id table. Its absence is
  only reachable through ``resolve``, so it surfaced only for a calc that had an
  active spec case and was asked for by id.
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
import re
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
    machinery is built from - so this is the check that makes the hand-maintained
    list in ``crates/azoth-python/src/results.rs`` self-policing.
    """
    extension_ids = set(_extension().calc_ids())
    registry_ids = set(CALC_IDS)
    assert extension_ids == registry_ids, (
        f"azoth._core.calc_ids() and the registry differ\n"
        f"  missing from the extension: {sorted(registry_ids - extension_ids)}\n"
        f"  not in the registry: {sorted(extension_ids - registry_ids)}\n"
        f"Add or remove the entry in crates/azoth-python/src/results.rs (calc_ids), "
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
    missing = sorted(calc_id for calc_id in CALC_IDS if not hasattr(bridge, FUNCTION_NAME[calc_id]))
    assert not missing, (
        f"{missing} have no bridge function. Add one named after the id's last "
        f"segment in python/src/azoth/_rust_bridge.py - there is no table to update."
    )
    with pytest.raises(KeyError):
        bridge.resolve("not.a_calc")


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
    attribute = re.compile(
        r'#\[pyclass\([^)]*?name\s*=\s*"(?P<name>\w+)"', re.DOTALL
    )
    struct = re.compile(r"pub struct \w+ \{")

    found: dict[str, tuple[str, ...]] = {}
    for path in sorted(source_dir.glob("*.rs")):
        text = path.read_text(encoding="utf-8")
        for match in attribute.finditer(text):
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
            found[match.group("name")] = tuple(
                re.findall(r"#\[pyo3\(get\)\]\s*pub (\w+):", body)
            )
    return found


def test_the_stub_matches_the_rust_transport_types() -> None:
    """Every field the stub declares for a transport type exists in the Rust struct.

    **This test exists because the drift it catches had already happened.** `FittingRow`
    declared `source_ref` and `source_locator` in the stub for as long as it took someone
    to read it, and both fields had been deleted from the Rust struct and from the
    project. The completeness test checked that the stub declares each *name*; nothing
    checked a field, so a stub describing attributes no object has type-checked clean and
    failed at the first attribute access.

    One-directional on purpose: the Rust struct may carry fields the stub omits only if
    they are not `#[pyo3(get)]`, which this parses for. A field the stub *declares* and
    the Rust does not have is the failure that shipped.
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
        checked += 1

    # A guard against the regex silently matching nothing, which would make every
    # assertion above vacuous and this test pass while checking nothing at all.
    assert checked >= 8, (
        f"only {checked} transport type(s) were compared against the Rust source, but "
        f"the stub declares more. The parser is matching less than it should."
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
