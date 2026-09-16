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
import re
import subprocess
import sys
from pathlib import Path
from types import ModuleType
from typing import Any

import pytest

import _helpers as h
from azoth import _models_gen
from azoth._dispatch import result_types
from azoth._registry_gen import BY_ID, CALCS, spec
from azoth.core.warnings import WarningCode

REPO_ROOT = Path(__file__).resolve().parents[2]
SCHEMA_PATH = REPO_ROOT / "specs" / "schema" / "calc.schema.json"
SPEC_PAGE = REPO_ROOT / "docs" / "src" / "architecture" / "specification.md"

#: `azoth has **30 ids** — 21 calculations and 9 models.`
_COUNT_CLAIM = re.compile(
    r"\*\*(?P<ids>\d+) ids\*\* — (?P<calcs>\d+) calculations and (?P<models>\d+) models"
)


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


def test_the_specification_states_the_real_number_of_ids() -> None:
    """`spec.md` says how many ids there are, and the registry is the authority.

    The page is normative and a reader arriving at it takes the count as a fact
    about the library. Nothing checked it, and it went stale: the unit-operation
    tier was deleted and the page went on saying **38 ids — 21 calculations and 17
    models** with 21 and 9 in the registry. A count is cheap to check and expensive
    to notice.

    The three numbers are compared separately rather than by their sum, because a
    change that moved one calc to the model tree would leave the total right and the
    page wrong.
    """
    page = SPEC_PAGE.read_text(encoding="utf-8")
    match = _COUNT_CLAIM.search(page)
    assert match, (
        f"{SPEC_PAGE.relative_to(REPO_ROOT)} no longer states an id count in the form "
        f"this test reads. Either the page changed and this regex did not, or the "
        f"count was removed - and a claim nothing checks is the one that goes stale."
    )
    assert int(match.group("calcs")) == len(CALCS), (
        f"spec.md says {match.group('calcs')} calculations; the registry has {len(CALCS)}"
    )
    assert int(match.group("models")) == len(_models_gen.MODELS), (
        f"spec.md says {match.group('models')} models; the registry has {len(_models_gen.MODELS)}"
    )
    assert int(match.group("ids")) == len(CALCS) + len(_models_gen.MODELS), (
        f"spec.md says {match.group('ids')} ids, which is not the "
        f"{len(CALCS) + len(_models_gen.MODELS)} the two numbers beside it give"
    )


@pytest.mark.requires_rust
def test_every_declared_range_check_runs() -> None:
    """No calc may declare a bound its own cases never evaluate.

    The calc-side half of the same assertion `test_model_contract.py` makes for models,
    and the reason it is behavioural rather than a name check: a declared quantity whose
    name resolves proves nothing about whether the implementation's resolver has an arm
    for it. The worked example was `process.mixer` — deleted since, along with the
    whole unit-operation tier — where `T` was a declared input and neither
    implementation's resolver produced it, so the bound was skipped on every call while
    reading as enforced.

    A calc's bounds are checked statically by `tools/spec_lint.py`, which is stricter here
    than it is for models: `computed_from` must name declared inputs, so a bound on a
    derived quantity has to say what it is derived from. What no static check can see is
    the other end — whether the code that resolves the name actually does.
    """
    from azoth._dispatch import resolve, use_backend

    exercised = 0
    for calc in CALCS:
        declared = [check["quantity"] for check in calc.get("valid_range", [])]
        if not declared:
            continue
        allowed = h.range_checks_that_may_skip(calc)
        for case in h.all_tests(calc):
            if case["status"] != "active" or case["type"] not in ("worked_example", "reference"):
                continue
            kwargs = h.kwargs_for(calc, case["inputs"])
            for backend in ("python", "rust"):
                with use_backend(backend):
                    result = resolve(calc["id"])(**kwargs)
                skipped = {
                    warning.field
                    for warning in result.warnings
                    if warning.code == WarningCode.RANGE_CHECK_SKIPPED
                }
                unexpected = sorted(skipped - allowed)
                assert not unexpected, (
                    f"{calc['id']}::{case['id']} on {backend}: the spec declares bounds on "
                    f"{unexpected}, and the implementation reports it could not evaluate "
                    f"{'them' if len(unexpected) > 1 else 'it'}. Declared: {declared}. "
                    f"Resolvable: {sorted(set(declared) - skipped)}"
                )
                exercised += len(declared) - len(skipped)
    assert exercised, "no range check was exercised - the loop above proved nothing"


def test_every_spec_has_a_result_type() -> None:
    """Every calc and model resolves to a result type, and nothing else does.

    The mapping is derived from each implementation's return annotation rather than
    listed, so `missing` here means an implementation that annotates nothing - which
    is a calc nothing can represent. `extra` means an id in neither registry, which
    would be a name a caller could ask for and get a shape for.
    """
    known = set(BY_ID) | {m["id"] for m in _models_gen.MODELS}
    missing = sorted(c["id"] for c in CALCS if c["id"] not in result_types())
    assert not missing, f"specs with no resolvable result type: {missing}"
    extra = sorted(set(result_types()) - known)
    assert not extra, f"result types with no spec: {extra}"


@pytest.mark.parametrize("calc", CALCS, ids=lambda c: c["id"])
def test_declared_outputs_are_result_fields(calc: dict[str, Any]) -> None:
    """Every declared output must exist as a field on the result dataclass.

    A declared output that is not a field is a spec promising something the code
    cannot deliver. The converse is checked by
    ``test_every_result_field_is_declared_in_the_spec`` below; both directions are
    required, so a field and a declared output are the same set.
    """
    # Any rather than a parameterised type: the registry maps ids to a
    # heterogeneous set of dataclasses, and `dataclasses.fields` is
    # deliberately runtime, so there is no useful static type to give it.
    result_type: Any = result_types()[calc["id"]]
    fields = {f.name for f in dataclasses.fields(result_type)}
    declared = set(calc["outputs"])
    missing = declared - fields
    assert not missing, (
        f"{calc['id']}: spec declares outputs {sorted(missing)} that are not fields "
        f"on {result_type.__name__} (has {sorted(fields)})"
    )
    assert "warnings" in fields, f"{calc['id']}: result has no warnings field"


#: Fields every result carries that no spec declares. The caveats channel is framework
#: rather than physics: `Warning` is its own cross-language contract, asserted
#: structurally in `python/tests/_helpers.py`, and restating it in 38 specs would be 38
#: copies of one fact.
FRAMEWORK_FIELDS = frozenset({"warnings"})


@pytest.mark.parametrize("calc", CALCS, ids=lambda c: c["id"])
def test_every_result_field_is_declared_in_the_spec(calc: dict[str, Any]) -> None:
    """The reverse of the test above: every field is a declared output.

    A field with no declared output is a value the spec does not name, so nothing
    holds it to anything - it is compared across the two languages without a spec
    saying it should exist, and a result struct emitted from ``outputs`` would drop
    it. Four calcs carried one before this test existed: ``iterations``,
    ``converged`` and ``residual`` on the two solvers, ``f`` and ``regime`` on
    darcy_weisbach, and crane_k_factors' per-fitting ``components``.
    """
    result_type: Any = result_types()[calc["id"]]
    fields = {f.name for f in dataclasses.fields(result_type)} - FRAMEWORK_FIELDS
    declared = set(calc["outputs"])
    undeclared = sorted(fields - declared)
    assert not undeclared, (
        f"{calc['id']}: {undeclared} appear on {result_type.__name__} but the spec "
        f"declares no such output, so nothing holds them to anything"
    )


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


#: Parameters a public wrapper may take that are not spec inputs.
#:
#: The keycard is not an input to a calculation - it is the authority the call reads
#: its data under - so it does not belong in a spec's `inputs:` and it does not reach
#: the kernels. Named rather than matched by a pattern, so adding a second one is a
#: deliberate edit to this set rather than a name that happens to look like plumbing.
_CAPABILITY_PARAMETERS = frozenset({"card"})


@pytest.mark.parametrize("calc", CALCS, ids=lambda c: c["id"])
def test_spec_inputs_match_function_signatures(calc: dict[str, Any]) -> None:
    """Spec inputs, the dispatched function's parameters, and the reference
    implementation's parameters must all be the same set.

    This is the anti-drift check on the *input* side, and it is the one a reviewer
    would otherwise have to do by eye across three files in two languages.

    The public wrapper is allowed the capability parameters above; the **reference
    implementation is not**, and that asymmetry is the point. A card is resolved at
    the boundary, so a kernel that took one would be reading the caller's authority
    itself rather than the values it was handed - which is the thing making the card
    a parameter is supposed to prevent.
    """
    namespace, _, function_name = calc["id"].rpartition(".")

    declared = set(calc["inputs"])
    public = namespace_module(calc["id"])
    signature = inspect.signature(getattr(public, function_name))
    dispatched = set(signature.parameters) - _CAPABILITY_PARAMETERS
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

    # A capability parameter must be keyword-only and defaulted, or it is a
    # positional argument a caller can shift by accident and a required one that
    # makes every existing call site fail.
    for name, parameter in signature.parameters.items():
        if name not in _CAPABILITY_PARAMETERS:
            continue
        assert parameter.kind is inspect.Parameter.KEYWORD_ONLY, (
            f"{calc['id']}: {name} is not keyword-only, so a caller can pass it "
            f"positionally and shift the inputs beside it"
        )
        assert parameter.default is not inspect.Parameter.empty, (
            f"{calc['id']}: {name} has no default, so every call must pass it"
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


#: The book's front door, which lists every calculation.
#:
#: It is hand-written prose - the account of how the pieces fit - with one
#: *generated* block spliced between markers by `tools/gen_docs.py`. The README is
#: not here: it links the book rather than carrying its own list, so this is the one
#: place a new calc has to be announced, and the drift is silent when it is not - the
#: library gains a calculation and the front page does not mention it.
CALC_LIST_FILES = ("docs/src/index.md",)


@pytest.mark.parametrize("relative_path", CALC_LIST_FILES)
def test_every_calc_is_announced_in_the_hand_written_lists(relative_path: str) -> None:
    """Every registered calc **and model** must be named in the hand-written lists.

    This is the drift that actually happens: a calc is added, the spec, both
    implementations and the tests all follow from the machinery, and the two places
    a *reader* looks are the two places nothing generates. The failure is silent -
    the library gains a calculation and the front page does not mention it.

    **Models were outside this check until they were not**, and the omission was
    caught by reading rather than by a failure: `eos.stability_test` arrived with a
    spec, two implementations, tests and a generated docs page, and was absent from
    both lists for the whole time. The test was one-directional in the wrong place -
    it named `CALCS`, so the newer half of the registry was never announced. That is
    the same shape as the batch API being missed from a count of registration points:
    a list that was complete when it was written, and a second thing that grew beside
    it.

    One-directional otherwise, on purpose. The lists may name things the registry
    does not, and today they should: both still say relief-valve sizing is not
    implemented, which is true and will stop being true. Asserting the reverse would
    forbid a reader-facing note about work in progress.

    That prose is also why this test cannot be a complete guard. It can catch an id
    missing from the list; it cannot tell whether the surrounding sentences about
    what is *not* implemented are still accurate. Whoever ships the first flowsheet
    has to read those paragraphs, and this test is what sends them there.
    """
    text = (REPO_ROOT / relative_path).read_text(encoding="utf-8")
    missing = sorted(
        entry["id"] for entry in [*CALCS, *_models_gen.MODELS] if entry["id"] not in text
    )
    assert not missing, (
        f"{relative_path} does not mention {missing}. Every registered id - calc or "
        f"model - has to appear in the list, because that is where a reader finds out "
        f"it exists. Run `python tools/gen_docs.py`."
    )


def test_the_announced_list_is_generated_rather_than_hand_edited() -> None:
    """The list is spliced in by the generator, so a hand-edit is caught here.

    The test above proves every id is *mentioned*; this one proves nobody typed it.
    Both are needed: a hand-edited list passes the first until the day someone adds a
    calc and forgets, and a generated list that has gone stale passes it too - right
    up until `gen_docs --check` is run, which is what this does.

    That check is also what the `docs-drift` CI job runs, so this is the same
    guarantee one step earlier - at the point where the person who forgot is still
    looking at the problem.
    """
    result = subprocess.run(
        [sys.executable, str(REPO_ROOT / "tools" / "gen_docs.py"), "--check"],
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, (
        f"a generated file is out of date - run `python tools/gen_docs.py`\n"
        f"{result.stdout}{result.stderr}"
    )


def test_spec_lookup_rejects_unknown_ids() -> None:
    """A wrong calc id must fail loudly rather than returning something."""
    with pytest.raises(KeyError, match="unknown calc"):
        spec("hydraulics.no_such_calc")
