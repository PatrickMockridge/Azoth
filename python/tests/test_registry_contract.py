"""Contract tests between the specs, the generated registry, and the code.

These are the tests that make "docs cannot drift from code" true rather than
aspirational. They are all cheap and none of them is numerical: they compare
*names* and *shapes*, which is where the drift these tests exist to catch always
shows up.

The motivating failure is a calc whose spec declares an output ``dP`` while the
implementation produces ``dp``. Every numerical test passes, because the numbers
agree perfectly - only the attribute name is wrong, and the first person to hit it
is a user.
"""

from __future__ import annotations

import dataclasses
import importlib
import inspect
import json
from pathlib import Path
from types import ModuleType
from typing import Any

import pytest

from azoth._registry_gen import BY_ID, CALCS, spec
from azoth.core.result import RESULT_TYPES
from azoth.core.warnings import WarningCode

REPO_ROOT = Path(__file__).resolve().parents[2]
SCHEMA_PATH = REPO_ROOT / "specs" / "schema" / "calc.schema.json"


def namespace_module(calc_id: str) -> ModuleType:
    """The public package a calc lives in, e.g. ``azoth.thermal`` for ``thermal.x``.

    Derived from the id rather than imported by name, because the point of this file
    is to hold the spec/registry/code contract for *every* namespace. A hardcoded
    import silently excluded the second one: the tests passed, having checked only
    the namespace that happened to be named.
    """
    namespace, _, _ = calc_id.rpartition(".")
    return importlib.import_module(f"azoth.{namespace}")


def reference_module(calc_id: str, function_name: str) -> ModuleType:
    """The pure-Python reference module for a calc, in whichever namespace."""
    namespace, _, _ = calc_id.rpartition(".")
    return importlib.import_module(f"azoth.{namespace}.reference.{function_name}")


def test_every_spec_has_a_result_type() -> None:
    """A spec with no result type is a calc nothing can represent."""
    missing = sorted(c["id"] for c in CALCS if c["id"] not in RESULT_TYPES)
    assert not missing, f"specs with no registered result type: {missing}"
    extra = sorted(set(RESULT_TYPES) - set(BY_ID))
    assert not extra, f"result types with no spec: {extra}"


@pytest.mark.parametrize("calc", CALCS, ids=lambda c: c["id"])
def test_declared_outputs_are_result_fields(calc: dict[str, Any]) -> None:
    """Every declared output must exist as a field on the result dataclass.

    The contract is one-directional: a result may carry diagnostic fields the spec
    does not declare (``iterations``, ``components``), but a declared output that
    is not a field means the spec promises something the code cannot deliver.
    """
    # Any rather than a parameterised type: the registry maps ids to a
    # heterogeneous set of dataclasses, and `dataclasses.fields` is
    # deliberately runtime, so there is no useful static type to give it.
    result_type: Any = RESULT_TYPES[calc["id"]]
    fields = {f.name for f in dataclasses.fields(result_type)}
    declared = set(calc["outputs"])
    missing = declared - fields
    assert not missing, (
        f"{calc['id']}: spec declares outputs {sorted(missing)} that are not fields "
        f"on {result_type.__name__} (has {sorted(fields)})"
    )
    assert "warnings" in fields, f"{calc['id']}: result has no warnings field"


@pytest.mark.parametrize("calc", CALCS, ids=lambda c: c["id"])
def test_every_result_carries_the_calc_id(calc: dict[str, Any]) -> None:
    """The spec's implementation paths must name the real functions."""
    namespace, _, function_name = calc["id"].rpartition(".")
    assert calc["implementations"]["python"] == f"azoth.{namespace}.{function_name}"

    namespace = namespace_module(calc["id"])
    assert hasattr(namespace, function_name), (
        f"{calc['id']}: azoth.{calc['id'].rpartition('.')[0]}.{function_name} does not exist"
    )

    reference = getattr(reference_module(calc["id"], function_name), function_name)
    assert callable(reference)


@pytest.mark.parametrize("calc", CALCS, ids=lambda c: c["id"])
def test_spec_inputs_match_function_signatures(calc: dict[str, Any]) -> None:
    """Spec inputs, the dispatched function's parameters, and the reference
    implementation's parameters must all be the same set.

    This is the anti-drift check on the *input* side, and it is the one a reviewer
    would otherwise have to do by eye across three files in two languages.
    """
    namespace, _, function_name = calc["id"].rpartition(".")

    declared = set(calc["inputs"])
    public = namespace_module(calc["id"])
    dispatched = set(inspect.signature(getattr(public, function_name)).parameters)
    reference = set(
        inspect.signature(
            getattr(reference_module(calc["id"], function_name), function_name)
        ).parameters
    )

    assert dispatched == declared, (
        f"{calc['id']}: azoth.{namespace}.{function_name} takes {sorted(dispatched)} "
        f"but the spec declares {sorted(declared)}"
    )
    assert reference == declared, (
        f"{calc['id']}: the reference implementation takes {sorted(reference)} "
        f"but the spec declares {sorted(declared)}"
    )
    # The namespace above is derived from the id, so this asserts the derivation
    # agrees with what the spec's own implementation path claims - the check that
    # used to be a hardcoded equality against "hydraulics".
    assert calc["implementations"]["python"].startswith(f"azoth.{namespace}.")


@pytest.mark.parametrize("calc", CALCS, ids=lambda c: c["id"])
def test_optional_inputs_are_optional_in_the_signature(calc: dict[str, Any]) -> None:
    """A spec declaring an input optional must be matched by a default in code.

    Otherwise the "omit it and get RANGE_CHECK_SKIPPED" behaviour is unreachable
    and the optionality is a lie.
    """
    _, _, function_name = calc["id"].rpartition(".")
    module = reference_module(calc["id"], function_name)
    parameters = inspect.signature(getattr(module, function_name)).parameters

    for name, declaration in calc["inputs"].items():
        is_optional = bool(declaration.get("optional", False))
        has_default = parameters[name].default is not inspect.Parameter.empty
        assert is_optional == has_default, (
            f"{calc['id']}.{name}: spec says optional={is_optional} but the signature "
            f"{'has' if has_default else 'has no'} default"
        )


def test_warning_codes_match_the_schema() -> None:
    """The Python warning codes must be exactly the set the schema allows.

    The schema is the contract both languages are held to, so a code added on one
    side only is caught here rather than discovered when two results that should
    compare equal do not.
    """
    schema = json.loads(SCHEMA_PATH.read_text(encoding="utf-8"))
    allowed = set(schema["$defs"]["range_check"]["properties"]["code"]["enum"])
    implemented = {code.value for code in WarningCode}
    assert implemented == allowed, (
        f"warning codes differ from the schema\n"
        f"  only in Python: {sorted(implemented - allowed)}\n"
        f"  only in schema: {sorted(allowed - implemented)}"
    )


def test_range_check_codes_used_by_specs_exist() -> None:
    """Every code a spec names must be one the implementation can emit."""
    for calc in CALCS:
        for raw in calc["valid_range"]:
            if "code" in raw:
                assert raw["code"] in {c.value for c in WarningCode}, (
                    f"{calc['id']}: range check on {raw['quantity']!r} names unknown "
                    f"code {raw['code']!r}"
                )


def test_registry_matches_the_spec_files() -> None:
    """The generated registry must be current.

    ``tools/gen_registry.py --check`` is the CI gate; this is the same check
    reachable from the test suite, so a developer who edits a spec and runs pytest
    finds out immediately rather than in CI.
    """
    import subprocess
    import sys

    result = subprocess.run(
        [sys.executable, str(REPO_ROOT / "tools" / "gen_registry.py"), "--check"],
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, (
        f"generated registries are out of date - run `python tools/gen_registry.py`\n"
        f"{result.stdout}{result.stderr}"
    )


#: Files that hand-maintain a list of the implemented calculations.
#:
#: These are the two places a new calc has to be announced by hand, and until this
#: test existed nothing checked either of them. They are prose rather than data, so
#: no generator can own them; the best available defence is a test that notices
#: when the prose has fallen behind the registry.
CALC_LIST_FILES = ("README.md", "docs/src/index.md")


@pytest.mark.parametrize("relative_path", CALC_LIST_FILES)
def test_every_calc_is_announced_in_the_hand_written_lists(relative_path: str) -> None:
    """Every registered calc must be named in the hand-written lists.

    This is the drift that actually happens: a calc is added, the spec, both
    implementations and the tests all follow from the machinery, and the two places
    a *reader* looks are the two places nothing generates. The failure is silent -
    the library gains a calculation and the front page does not mention it.

    One-directional on purpose. The lists may name things the registry does not,
    and today they should: both still say orifice, control valve, relief valve and
    pump calculations are not implemented, which is true and will stop being true.
    Asserting the reverse would forbid a reader-facing note about work in progress.

    That prose is also why this test cannot be a complete guard. It can catch a
    calc missing from the list; it cannot tell whether the surrounding sentences
    about what is *not* implemented are still accurate. Whoever adds the first
    orifice calc has to read those paragraphs, and this test is what sends them
    there.
    """
    text = (REPO_ROOT / relative_path).read_text(encoding="utf-8")
    missing = sorted(calc["id"] for calc in CALCS if calc["id"] not in text)
    assert not missing, (
        f"{relative_path} does not mention {missing}. Every calc in the registry has "
        f"to appear in the hand-written lists, because that is where a reader finds "
        f"out it exists."
    )


def test_spec_lookup_rejects_unknown_ids() -> None:
    """A wrong calc id must fail loudly rather than returning something."""
    with pytest.raises(KeyError, match="unknown calc"):
        spec("hydraulics.no_such_calc")
