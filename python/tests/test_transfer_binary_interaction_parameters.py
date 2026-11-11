"""`characterization.transfer_binary_interaction_parameters`: the recorded matrices.

The expected matrix is a `matrix` output, which a case's `expected` block carries but the
generated `TestCase` has no field for - so the Rust test checks the answer's shape, and the
values are pinned here from the model's raw spec. `test_cross_impl` compares the two kernels on
the same cases, which carries these values to the Rust side.
"""

from __future__ import annotations

import pytest

from azoth import _models_gen, characterization
from azoth.core.errors import InvalidInputError

MODEL_ID = "characterization.transfer_binary_interaction_parameters"
SPEC = _models_gen.model(MODEL_ID)


@pytest.mark.parametrize("case", SPEC["cases"], ids=lambda case: case["id"])
def test_the_answer_is_the_recorded_matrix(case: dict[str, object]) -> None:
    inputs = case["inputs"]
    assert isinstance(inputs, dict)
    result = characterization.transfer_binary_interaction_parameters(
        reference_components=inputs["reference_components"],
        target_components=inputs["target_components"],
        reference_kij=inputs["reference_kij"],
        target_kij=inputs["target_kij"],
    )
    expectations = case["expected"]
    assert isinstance(expectations, dict)
    expected = expectations["target_kij"]
    assert len(result.target_kij) == len(expected)
    for index, (row, want) in enumerate(zip(result.target_kij, expected, strict=True)):
        assert list(row) == pytest.approx(want, rel=1e-12, abs=0.0), f"row {index}"


def test_a_missing_second_base_name_is_refused() -> None:
    reference = [
        {"name": "methane", "pseudo_component": False},
        {"name": "CO2", "pseudo_component": False},
    ]
    target = [
        {"name": "methane", "pseudo_component": False},
        {"name": "C7_PC", "pseudo_component": True},
    ]
    with pytest.raises(InvalidInputError, match="CO2"):
        characterization.transfer_binary_interaction_parameters(
            reference_components=reference,
            target_components=target,
            reference_kij=[[0.0, 0.09], [0.09, 0.0]],
            target_kij=[[0.0, 0.0], [0.0, 0.0]],
        )


def test_a_pair_the_reference_states_as_zero_is_not_transferred() -> None:
    """The target's own value survives: a zero is the reference declining, not clearing it."""
    components = [
        {"name": "methane", "pseudo_component": False},
        {"name": "CO2", "pseudo_component": False},
    ]
    result = characterization.transfer_binary_interaction_parameters(
        reference_components=components,
        target_components=components,
        reference_kij=[[0.0, 0.0], [0.0, 0.0]],
        target_kij=[[0.0, 0.0973], [0.0973, 0.0]],
    )
    assert result.target_kij[0][1] == pytest.approx(0.0973, rel=1e-12)


def test_a_matrix_that_is_not_square_is_refused() -> None:
    components = [{"name": "methane", "pseudo_component": False}]
    with pytest.raises(InvalidInputError, match="reference_kij"):
        characterization.transfer_binary_interaction_parameters(
            reference_components=components,
            target_components=components,
            reference_kij=[[0.0, 0.0]],
            target_kij=[[0.0]],
        )
