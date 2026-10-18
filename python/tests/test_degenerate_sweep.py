"""Every registered model and calc, on inputs no case states, on both backends.

The corpus is a set of states NeqSim was captured at, and every one of them carries
strictly positive compositions. So the two implementations agree everywhere the corpus
looks, and cannot be compared where it does not look - which is how a zero mole fraction
raised `math domain error` in Python while Rust answered, on a model landed and green for
months.

**The rule this holds them to is the coin landing:** for one input, either both backends
answer and the answers agree at the case's own tolerance, or both refuse with the same
`AzothError` class object and the same field. One answering while the other refuses is a
defect, and so is a raw exception - it is a missing check, not a refusal.

The variants are applied to a concrete case's inputs rather than synthesised from the
spec, because the registry's `length` and `shape` are prose. One numeric vector, one
matrix and the first two scalar quantities are degraded in turn; the rest of the case is
the state the case declares.
"""

from __future__ import annotations

import dataclasses
import math
from collections.abc import Iterator
from typing import Any

import pytest

import _helpers as h
from azoth import _models_gen, use_backend
from azoth._dispatch import Backend, resolve
from azoth._registry_gen import CALCS
from azoth.core.errors import AzothError

pytestmark = [pytest.mark.requires_rust, pytest.mark.slow_sweep]

Entry = tuple[str, dict[str, Any]]


def _entries() -> list[Entry]:
    return [("calc", calc) for calc in CALCS] + [("model", model) for model in _models_gen.MODELS]


def _active_cases(entry: dict[str, Any], kind: str) -> list[dict[str, Any]]:
    """The cases a diff would run, from the same helper the cross-impl test reads."""
    if kind == "calc":
        return [
            case
            for case in h.all_tests(entry)
            if case["status"] == "active" and case["type"] in ("worked_example", "reference")
        ]
    return list(entry["cases"])


def _leaves(
    entry: dict[str, Any], inputs: dict[str, Any]
) -> tuple[list[str], list[str], list[str]]:
    """One vector, the matrices, and the first two scalar quantities the case states."""
    vectors: list[str] = []
    matrices: list[str] = []
    scalars: list[str] = []
    for name, declaration in entry["inputs"].items():
        if name not in inputs:
            continue
        kind = declaration.get("type", "quantity")
        if kind == "vector":
            vectors.append(name)
        elif kind == "matrix":
            matrices.append(name)
        elif kind == "quantity":
            scalars.append(name)
    return vectors[:1], matrices[:1], scalars[:2]


def _variants(entry: dict[str, Any], case: dict[str, Any]) -> Iterator[tuple[str, dict[str, Any]]]:
    """The case's inputs, with one leaf degraded at a time."""
    inputs = dict(case["inputs"])
    vectors, matrices, scalars = _leaves(entry, inputs)

    for name in vectors:
        vector: list[float] = [float(item) for item in inputs[name]]
        if not vector:
            continue
        yield f"{name}=zeros", {**inputs, name: [0.0] * len(vector)}
        yield f"{name}[0]=0", {**inputs, name: [0.0, *vector[1:]]}
        yield f"{name}=nan", {**inputs, name: [math.nan] * len(vector)}
        yield f"{name}=inf", {**inputs, name: [math.inf] * len(vector)}
        yield f"{name}=negative", {**inputs, name: [-item for item in vector]}
        yield f"{name}=empty", {**inputs, name: []}

    for name in matrices:
        rows: list[list[float]] = [[float(item) for item in row] for row in inputs[name]]
        if not rows or not rows[0]:
            continue
        yield f"{name}=zeros", {**inputs, name: [[0.0] * len(row) for row in rows]}
        yield f"{name}[0]=zeros", {**inputs, name: [[0.0] * len(rows[0]), *rows[1:]]}

    for name in scalars:
        yield f"{name}=0", {**inputs, name: 0.0}
        yield f"{name}=-1", {**inputs, name: -1.0}
        yield f"{name}=nan", {**inputs, name: math.nan}
        yield f"{name}=inf", {**inputs, name: math.inf}


SWEEP = [
    (kind, entry, case, label, inputs)
    for kind, entry in _entries()
    for case in _active_cases(entry, kind)
    for label, inputs in _variants(entry, case)
]

SWEEP_IDS = [f"{entry['id']}::{case['id']}::{label}" for _, entry, case, label, _ in SWEEP]


def _kwargs(entry: dict[str, Any], kind: str, inputs: dict[str, Any]) -> dict[str, Any]:
    return h.kwargs_for(entry, inputs) if kind == "calc" else h.model_kwargs(entry, inputs)


def _outcome(entry: dict[str, Any], kind: str, inputs: dict[str, Any], backend: Backend) -> Any:
    kwargs = _kwargs(entry, kind, inputs)
    with use_backend(backend):
        try:
            return ("ok", resolve(entry["id"])(**kwargs))
        except Exception as error:
            return ("raised", error)


def _merge_flat_phases(result: Any) -> Any:
    """Combine phases with the same composition into one, summing their fractions.

    **A split between two phases of identical composition is a flat direction.** The
    specification says so where the model can reach it: the Newton solve converges on the
    compositions but the individual fractions of two identical phases are path-decided, so
    two kernels can report `0.02200/0.01164` and `0.02199/0.01165` for the same state.
    Their sum - the amount in that composition - is what the physics fixes, so it is what
    is compared.
    """
    fractions = getattr(result, "phase_fractions", None)
    compositions = getattr(result, "x", None)
    if not fractions or not compositions or len(fractions) != len(compositions):
        return result
    groups: list[list[Any]] = []
    for fraction, row in zip(fractions, compositions, strict=True):
        for group in groups:
            if all(abs(a - b) <= 1e-9 for a, b in zip(row, group[0], strict=True)):
                group[1] += fraction
                break
        else:
            groups.append([row, fraction])
    if len(groups) == len(fractions):
        return result
    return dataclasses.replace(
        result,
        phase_fractions=tuple(group[1] for group in groups),
        x=tuple(group[0] for group in groups),
    )


def _compare(
    left: Any,
    right: Any,
    tolerance: float,
    context: str,
    diagnostic_bound: float | None,
) -> None:
    """Both answer and agree, or both refuse the same way - and nothing else."""
    left_kind, left_value = left
    right_kind, right_value = right
    if left_kind == "ok" and right_kind == "ok":
        # **A non-converged result is a reported failure, not an answer.** Both sides saying
        # so is the agreement: the iterate each publishes is where its own path stopped, and
        # `reactions.chemical_equilibrium`'s `does_not_converge` case is built to reach one.
        # Holding two last iterates to the answer tolerance would be comparing the route.
        if (
            getattr(left_value, "converged", None) is False
            and getattr(right_value, "converged", None) is False
        ):
            return
        # **The solver's own bound is an absolute floor on its answer.** An iteration that
        # stopped at `1e-6` knows its output to `1e-6`, not to `1e-6` relative: a reflux
        # flash whose vapour fraction is `8.68e-7` is decided to that absolute precision,
        # and the two kernels land `1.2e-10` apart, which a relative comparison reads as
        # `1.4e-5` of a number neither promised.
        h.assert_results_equal(
            _merge_flat_phases(left_value),
            _merge_flat_phases(right_value),
            tolerance,
            context,
            diagnostic_bound,
            nan_agrees=True,
            absolute_floor=diagnostic_bound,
        )
        return
    if left_kind == "raised" and right_kind == "raised":
        if not isinstance(left_value, AzothError) or not isinstance(right_value, AzothError):
            raise AssertionError(
                f"{context}: a raw exception is not a refusal; "
                f"python={left_value!r}, rust={right_value!r}"
            )
        assert type(left_value) is type(right_value), (
            f"{context}: python raises {type(left_value).__name__}, "
            f"rust raises {type(right_value).__name__}"
        )
        assert left_value.field() == right_value.field(), (
            f"{context}: both raise {type(left_value).__name__} but on different fields: "
            f"python={left_value.field()!r}, rust={right_value.field()!r}"
        )
        return
    if left_kind == "ok":
        raise AssertionError(f"{context}: python answered and rust refused with {right_value!r}")
    raise AssertionError(f"{context}: rust answered and python refused with {left_value!r}")


@pytest.mark.parametrize(
    ("kind", "entry", "case", "label", "inputs"),
    SWEEP,
    ids=SWEEP_IDS,
)
def test_a_degenerate_input_lands_on_one_side(
    kind: str, entry: dict[str, Any], case: dict[str, Any], label: str, inputs: dict[str, Any]
) -> None:
    tolerance = (
        float(case.get("tolerance") or entry["worked_example"]["tolerance"])
        if kind == "calc"
        else float(case["tolerance"])
    )
    context = f"{entry['id']}::{case['id']}::{label}"
    python = _outcome(entry, kind, inputs, "python")
    rust = _outcome(entry, kind, inputs, "rust")
    _compare(python, rust, tolerance, context, h.convergence_tolerance(entry))
