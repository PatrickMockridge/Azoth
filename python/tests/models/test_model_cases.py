"""Every model's every case, against the value the case records.

**Why this is here.** ``test_cross_impl`` runs each model case through both kernels and
compares them *with each other* - it never reads ``case["expected"]``. So a model whose
case is not also run by a per-model test has recorded numbers that nothing checks, and
``eos.srk_cpa_phase`` and ``eos.pr_cpa_phase`` were exactly that: two registered models
whose recorded departures were unchecked, which is how a 28% error in the association's
enthalpy term survived in both.

A per-model test file is still the right place for a model's own physics. This is the
floor under all of them: the recorded value is the case's whole point, and a case whose
expectations are never compared is a comment with a schema.
"""

from __future__ import annotations

import dataclasses
from typing import Any

import pint
import pytest

import _helpers as h
from azoth import _models_gen

#: Models whose case records **per-phase vectors** - `beta`, `x`, `y`, `keq` - whose order
#: is the model's rather than the case's. `eos.tp_multiflash` returns its phases in its own
#: order and its own test in `crates/azoth-eos/tests/tp_multiflash.rs` matches them *by
#: composition* for exactly that reason; comparing them positionally here would report a
#: reordering as a divergence and say nothing about the physics.
#:
#: Named rather than detected, because "is this vector per-phase" is not a property of its
#: type. A model joining the list is a claim that its case cannot be compared positionally,
#: and that claim belongs next to the reason for it.
PER_PHASE_ORDERED = frozenset({"eos.tp_multiflash"})

#: Cases are run through the selected backend, so this is the Python kernel by default
#: and the Rust one under `AZOTH_REQUIRE_RUST`. The cross-implementation comparison is
#: `test_cross_impl`'s, and it is the other half of the same check.
MODEL_CASES = [
    (model, case)
    for model in _models_gen.MODELS
    if model["id"] not in PER_PHASE_ORDERED
    for case in model["cases"]
]


def _expected_result(result: Any, case: dict[str, Any]) -> Any:
    """A copy of ``result`` with every field the case records replaced by its value.

    Built this way rather than compared field by field here, so that the comparison is
    :func:`_helpers.assert_results_equal`'s - the same one the cross-implementation tests
    use, with the same three passes over scalars, quantities and vectors. A second
    comparison written for this file is a second place for the conventions to drift.
    """
    updates: dict[str, Any] = {}
    for spec_field in dataclasses.fields(result):
        name = spec_field.name
        if name not in case["expected"] or name in h._UNCOMPARED_FIELDS:
            continue
        declared = case["expected"][name]
        current = getattr(result, name)
        if isinstance(current, pint.Quantity):
            # A scalar output carrying a unit. The case records the magnitude in that
            # unit, which is what `spec_lint` holds it to.
            updates[name] = type(current)(declared, current.units)
        elif isinstance(current, tuple) and current and isinstance(current[0], pint.Quantity):
            updates[name] = tuple(type(current[0])(value, current[0].units) for value in declared)
        elif isinstance(declared, (tuple, list)):
            updates[name] = type(current)(declared)
        else:
            updates[name] = declared
    return dataclasses.replace(result, **updates)


@pytest.mark.parametrize(
    ("model", "case"),
    MODEL_CASES,
    ids=lambda x: x["id"] if isinstance(x, dict) else "",
)
def test_a_model_case_matches_what_it_records(model: dict[str, Any], case: dict[str, Any]) -> None:
    """Run one model case and compare every field it declares an expectation for."""
    from azoth._dispatch import resolve

    kwargs = h.model_kwargs(model, case["inputs"])
    result = resolve(model["id"])(**kwargs)
    h.assert_results_equal(
        result,
        _expected_result(result, case),
        float(case["tolerance"]),
        f"{model['id']}::{case['id']}",
        diagnostic_bound=h.convergence_tolerance(model),
    )
