"""Run every validation case under ``validation/``.

These are the *external* checks. The tests elsewhere in the suite are generated
from the specs and check that the code does what the spec says; these check
whether the spec was right, because they come from outside the registry and can
disagree with it. See ``validation/README.md``.

A mismatch fails the build. That is the whole point: a validation case nobody
runs is a claim in a JSON file.

# Models as well as calculations

``calc`` names an id in *either* registry. A model needs this at least as much as
a calculation does: what a model's spec pins down is a procedure, and a procedure
is not something a published source states - so an external case is the only kind
of check on it that does not come from the spec itself. It also means a case's
``expected`` may be a vector, an enum member, or absent (for a model whose answer
is ``None``), so the comparison handles all four shapes rather than only floats.
"""

from __future__ import annotations

import importlib
import json
from pathlib import Path
from typing import Any

import pytest

import _helpers as h
from azoth import _models_gen, _registry_gen

REPO_ROOT = Path(__file__).resolve().parents[3]
VALIDATION_DIR = REPO_ROOT / "validation"

REQUIRED_KEYS = {"id", "calc", "source", "inputs", "expected", "tolerance"}

#: Statuses that oblige the case to explain itself. A case whose arithmetic is
#: ours must say so, or a passing test reads as an external validation it is not.
MUST_BE_EXPLAINED = {"unverified", "source_needed"}


def _case_paths() -> list[Path]:
    return sorted(p for p in VALIDATION_DIR.rglob("*.json"))


def _load(path: Path) -> dict[str, Any]:
    # Annotated on its own line: `json.loads` returns Any, and returning it
    # directly is an implicit Any at the boundary of every test here.
    case: dict[str, Any] = json.loads(path.read_text(encoding="utf-8"))
    return case


def _spec(calc_id: str) -> dict[str, Any] | None:
    """The spec behind an id, from whichever registry declares it."""
    if calc_id in _registry_gen.BY_ID:
        return _registry_gen.BY_ID[calc_id]
    return _models_gen.model(calc_id) or None


def _mixture_kwargs(case: dict[str, Any]) -> dict[str, Any]:
    """Resolve a case's component names into the mixture a flash takes.

    The one place a case's JSON does not map straight onto a call. A mixture is an
    object, and JSON has no way to hold one - so `eos.pt_flash` is called with the
    composition vector and the *names* the spec declares, and `mixture_of` turns them
    into the object through the databank.

    The names in `inputs` are the ones `specs/models/eos/pt_flash.toml` declares, so
    the "inputs name declared quantities" check below still means something, and a
    reader of the JSON can look each one up in the model's documentation.

    Until this change the adapter built the components from three per-component
    vectors the case carried, which made the case's fluid a second fluid - described
    by a JSON file rather than by NeqSim's tables. A validation case exists to check
    the *model* against something outside it, and a fluid the case invented made that
    check partly a check of the case.
    """
    from azoth import ureg
    from azoth.eos.components import mixture_of

    q = ureg.Quantity
    inputs = case["inputs"]
    return {
        "mixture": mixture_of(list(inputs["components"]))[0],
        "T": q(inputs["T"], "K"),
        "P": q(inputs["P"], "Pa"),
        "z": list(inputs["z"]),
    }


def _enthalpy_kwargs(case: dict[str, Any]) -> dict[str, Any]:
    """Resolve a case's component names into the mixture and ideal-gas model.

    Same shape as :func:`_mixture_kwargs`, one step further: `mixture_of` returns
    both the `Mixture` and the `IdealGasModel` the enthalpy model takes, and the
    case supplies the cubic root (`compressibility`) it is evaluated at rather
    than solving for it - which root describes the phase is a choice the flash
    has already made.
    """
    from azoth import ureg
    from azoth.eos.components import mixture_of

    q = ureg.Quantity
    inputs = case["inputs"]
    mixture, ideal_gas = mixture_of(list(inputs["components"]))
    return {
        "mixture": mixture,
        "ideal_gas": ideal_gas,
        "T": q(inputs["T"], "K"),
        "P": q(inputs["P"], "Pa"),
        "z": list(inputs["z"]),
        "compressibility": inputs["compressibility"],
    }


#: Ids whose call arguments are not a straight copy of the case's `inputs`.
#:
#: One entry today. It grows by one line per model whose arguments are objects
#: rather than numbers, which is the point of it being a table rather than a branch
#: buried in `_call` - a reader can see the whole of the exception list at once.
ARGUMENT_BUILDERS = {
    "eos.pt_flash": _mixture_kwargs,
    "eos.molar_enthalpy_entropy": _enthalpy_kwargs,
}


def _call(case: dict[str, Any]) -> Any:
    """Run a case through the public API of its own namespace.

    The module is derived from the id rather than imported by hand, so a case in a
    namespace this file has never heard of still runs - which is what makes adding
    one a data change rather than a code change.
    """
    namespace = case["calc"].split(".")[0]
    function_name = case["calc"].rpartition(".")[2]
    module = importlib.import_module(f"azoth.{namespace}")
    builder = ARGUMENT_BUILDERS.get(case["calc"])
    if builder is not None:
        kwargs = builder(case)
    else:
        # A bare number becomes a quantity in the unit the spec declares, the same
        # conversion every other test in the suite uses. Doing it here rather than in
        # the case means a validation case states numbers in the spec's units and
        # needs no unit handling of its own.
        spec = _spec(case["calc"])
        assert spec is not None
        kwargs = h.kwargs_for(spec, case["inputs"])
    return getattr(module, function_name)(**kwargs)


def _compare(got: Any, want: Any, tolerance: float, where: str, name: str) -> None:
    """Compare one expected output against one actual one, whatever shape it is."""
    if want is None:
        assert got is None, f"{where} ({name}): expected no value, got {got!r}"
        return

    if isinstance(want, list):
        # A vector output - a composition, or one K-value per component.
        actual = [float(v) for v in got]
        assert len(actual) == len(want), (
            f"{where} ({name}): got {len(actual)} entries, expected {len(want)}"
        )
        for i, (a, w) in enumerate(zip(actual, want, strict=True)):
            h.assert_close(a, float(w), tolerance, f"{where} ({name}[{i}])")
        return

    if isinstance(want, str):
        # An enum output, compared by its spec spelling.
        value = getattr(got, "value", got)
        assert value == want, f"{where} ({name}): got {value!r}, expected {want!r}"
        return

    if isinstance(want, bool):
        assert bool(got) is want, f"{where} ({name}): got {got!r}"
        return

    value = getattr(got, "magnitude", got)
    h.assert_close(float(value), float(want), tolerance, f"{where} ({name})")


CASES = _case_paths()


def test_there_is_at_least_one_case() -> None:
    """Guard against a refactor that silently stops running any of them."""
    assert CASES, f"no validation cases found under {VALIDATION_DIR}"


@pytest.mark.parametrize("path", CASES, ids=lambda p: p.stem)
def test_case_is_well_formed(path: Path) -> None:
    """A case must carry what the runner needs, and name a real id.

    Checked separately from running it, so a malformed case reports as a
    malformed case rather than as an arithmetic failure.
    """
    case = _load(path)
    where = path.relative_to(REPO_ROOT)

    missing = REQUIRED_KEYS - set(case)
    assert not missing, f"{where} is missing {sorted(missing)}"
    assert _spec(case["calc"]) is not None, (
        f"{where} names {case['calc']!r}, which is in neither the calc nor the model registry"
    )
    assert case["tolerance"] > 0, f"{where} has a non-positive tolerance"

    source = case["source"]
    assert source.get("verification") in {"verified", "unverified", "source_needed"}, (
        f"{where} has an unrecognised source.verification"
    )


@pytest.mark.parametrize("path", CASES, ids=lambda p: p.stem)
def test_unconfirmed_sources_are_explained(path: Path) -> None:
    """A case that is not backed by a confirmed source has to say why.

    Without this, `source.verification: source_needed` is a field nobody reads,
    and a green test tick reads as "validated against the standard" when it
    means "the arithmetic we wrote matches the arithmetic we wrote".
    """
    case = _load(path)
    where = path.relative_to(REPO_ROOT)
    status = case["source"]["verification"]

    if status in MUST_BE_EXPLAINED:
        notes = str(case["source"].get("notes", "")).strip()
        assert notes, (
            f"{where}: source.verification is {status!r} but there are no notes "
            f"explaining what is unconfirmed. A reader has to be able to tell this "
            f"from a case validated against a source."
        )


@pytest.mark.parametrize("path", CASES, ids=lambda p: p.stem)
def test_the_inputs_name_declared_quantities(path: Path) -> None:
    """Every input a case supplies must be one the spec declares.

    A case passing an input the implementation does not take would otherwise fail
    as a `TypeError` from inside the call, which reads as a broken function rather
    than as a broken case.
    """
    case = _load(path)
    spec = _spec(case["calc"])
    assert spec is not None
    where = path.relative_to(REPO_ROOT)

    unknown = sorted(set(case["inputs"]) - set(spec["inputs"]))
    assert not unknown, (
        f"{where}: passes {unknown}, which {case['calc']} does not declare. "
        f"Declared: {sorted(spec['inputs'])}"
    )


@pytest.mark.parametrize("path", CASES, ids=lambda p: p.stem)
def test_case_matches_implementation(path: Path) -> None:
    """Run the case and compare, output by output.

    Compares each named output separately so a failure says which one disagreed,
    by how much, and against what - rather than reporting that two dictionaries
    differ. The source's verification status is deliberately *not* a skip
    condition: the arithmetic is worth checking either way, and a skip would
    report a case that ran as one that did not.
    """
    case = _load(path)
    where = path.relative_to(REPO_ROOT)
    result = _call(case)

    for name, want in case["expected"].items():
        assert hasattr(result, name), (
            f"{where}: expected output {name!r} is not a field of {type(result).__name__}"
        )
        _compare(getattr(result, name), want, float(case["tolerance"]), str(where), name)
