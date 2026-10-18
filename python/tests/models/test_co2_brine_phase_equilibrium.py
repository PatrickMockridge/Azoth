"""``eos.co2_brine_phase_equilibrium``: the states its spec cases cannot carry.

A case's feed is a set of strictly positive mole numbers, because every captured NeqSim
state writes an absent species at a `1e-50`-style floor. **An exact zero is therefore a
state neither backend was ever run on**, and the two languages disagreed there: Python's
`math.log(0.0)` raises where Rust's `f64::ln(0.0)` is `-inf`. This file drives the feed the
reactive-brine coupling actually passes, four of whose nine species are at exact zero.
"""

from __future__ import annotations

from typing import cast

import pytest

import _helpers as h
from azoth import _models_gen, use_backend
from azoth._dispatch import Backend, resolve
from azoth.core.result import Co2BrinePhaseEquilibriumResult

MODEL_ID = "eos.co2_brine_phase_equilibrium"

SPEC = _models_gen.model(MODEL_ID)

#: `validation/neqsim/ReactiveBrineProbe.java`'s feed at 50 bara: the nine species
#: `chemicalReactionInit()` grows the five-component brine to, with the four products of
#: the reactions it has not run yet at exact zero.
COMPONENTS = ["CO2", "water", "Na+", "K+", "Cl-", "OH-", "H3O+", "HCO3-", "CO3--"]
MOLES = [
    10.0,
    55.508435061791985,
    0.540339675351135,
    0.28239120684612296,
    0.822730882197258,
    0.0,
    0.0,
    0.0,
    0.0,
]
ZERO = {5, 6, 7, 8}


INPUTS = {"components": COMPONENTS, "T": 283.15, "P": 5.0e6, "moles": MOLES}


def _result(backend: Backend) -> Co2BrinePhaseEquilibriumResult:
    kwargs = h.model_kwargs(SPEC, INPUTS)
    with use_backend(backend):
        return cast(Co2BrinePhaseEquilibriumResult, resolve(MODEL_ID)(**kwargs))


@pytest.mark.requires_rust
def test_an_exact_zero_composition_answers_the_same_on_both_backends() -> None:
    """The Python reference used to raise `math domain error` here; Rust answered."""
    python = _result("python")
    rust = _result("rust")
    h.assert_results_equal(
        python,
        rust,
        1.0e-8,
        f"{MODEL_ID}::exact-zero feed",
        diagnostic_bound=h.convergence_tolerance(SPEC),
    )


@pytest.mark.requires_rust
def test_a_zero_species_contributes_no_material() -> None:
    """The class's floor at `1e-50` is the most a species at exact zero may appear as."""
    result = _result("python")
    for row in result.x:
        for index in ZERO:
            assert 0.0 <= row[index] <= 1.0e-50
    assert sum(result.phase_fractions) == pytest.approx(1.0, abs=1.0e-10)
