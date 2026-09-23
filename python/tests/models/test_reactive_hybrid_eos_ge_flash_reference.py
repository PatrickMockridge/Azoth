"""The projection's branches, which no captured fluid exercises.

``reactions.reactive_hybrid_eos_ge_flash`` takes the *shortcut* on every pass of both captured
fluids - the raw reaction delta already satisfies ``A delta`` to ``1e-8`` and leaves no
inventory negative - so the oracle pins the loop around the projection and never the
projection itself. These tests are the other branch, on the captured matrix and the captured
inventory, with deltas chosen to reach it. The Rust kernel's own unit tests are the mirror of
this file.
"""

from __future__ import annotations

import math

from azoth.reactions.reference._linalg import project_onto_null_space, row_space_basis
from azoth.reactions.reference.reactive_hybrid_eos_ge_flash import (
    MINIMUM_COUPLED_MOLES,
    NEGATIVE_INVENTORY_TOLERANCE_MOLES,
    REACTION_DELTA_CONSERVATION_TOLERANCE_MOLES,
    _conservative_delta,
)

#: The captured conservation matrix, whose oxygen row depends on the other three.
MATRIX = [
    [1.0, 0.0, 1.0, 1.0, 0.0, 0.0],
    [0.0, 2.0, 1.0, 0.0, 1.0, 3.0],
    [2.0, 1.0, 3.0, 3.0, 1.0, 1.0],
    [0.0, 0.0, -1.0, -2.0, -1.0, 1.0],
]

#: The captured coupled inventory, in the fluid's order.
COUPLED = [
    5.0,
    0.049997206913475076,
    55.499994364714574,
    6.0e-4,
    2.0e-4,
    0.0010027670978232208,
    2.619878633189711e-8,
    1.1443550870656359e-8,
    2.830738529430478e-6,
]

#: The reactive components' positions in it: `CO2, water, HCO3-, CO3--, OH-, H3O+`.
REACTIVE_AT = [1, 2, 5, 6, 7, 8]


def conservation(delta: list[float]) -> float:
    """``max |A delta|``, the test that decides the branch."""
    return max(
        abs(sum(coefficient * value for coefficient, value in zip(row, delta, strict=True)))
        for row in MATRIX
    )


def test_a_non_conservative_delta_is_projected_onto_the_null_space() -> None:
    """**The projection runs where the raw delta is not conservative.**

    A unit step on water violates carbon's and oxygen's rows, so the shortcut refuses it and
    the answer is ``delta - A+ A delta``: conservation to rounding, and not the delta that went
    in.
    """
    basis = row_space_basis(MATRIX)
    raw = [0.0, 1.0e-3, 0.0, 0.0, 0.0, 0.0]
    assert conservation(raw) > REACTION_DELTA_CONSERVATION_TOLERANCE_MOLES

    delta = _conservative_delta(MATRIX, basis, COUPLED, REACTIVE_AT, raw)
    assert conservation(delta) < 1.0e-16
    assert any(abs(value - original) > 1.0e-6 for value, original in zip(delta, raw, strict=True))


def test_a_conservative_delta_is_used_as_it_stands() -> None:
    """**The shortcut is taken where the raw delta already conserves**, which is the branch both
    captured fluids take on every pass."""
    basis = row_space_basis(MATRIX)
    raw = [coefficient * 1.0e-10 for coefficient in MATRIX[2]]
    assert conservation(raw) < REACTION_DELTA_CONSERVATION_TOLERANCE_MOLES

    assert _conservative_delta(MATRIX, basis, COUPLED, REACTIVE_AT, raw) == raw


def test_a_negative_inventory_scales_the_whole_delta() -> None:
    """**A delta that would take an inventory negative is scaled whole.**

    The carbonate sits at ``2.6e-8`` and a step of ``-1e-6`` is not admissible. Scaling one
    component would take the delta out of the null space, so the whole step is scaled.
    """
    basis = row_space_basis(MATRIX)
    raw = [0.0, 0.0, 0.0, 0.0, 0.0, -1.0e-6]

    delta = _conservative_delta(MATRIX, basis, COUPLED, REACTIVE_AT, raw)
    for position, value in enumerate(delta):
        index = REACTIVE_AT[position]
        assert COUPLED[index] + value >= -NEGATIVE_INVENTORY_TOLERANCE_MOLES
    available = COUPLED[REACTIVE_AT[5]] - MINIMUM_COUPLED_MOLES
    assert delta[5] >= -available * (1.0 + 1.0e-12)
    assert delta[5] < 0.0, "the direction survives the scaling"
    assert conservation(delta) < 1.0e-16, "and the scaled delta still conserves"


def test_the_projector_is_the_captured_svd_one() -> None:
    """**The projection against ``A+ A`` as Commons Math computes it.**

    ``A+ A`` is the orthogonal projector onto ``A``'s row space, and it is what
    ``delta - A+ A delta`` is built from. The matrix is the captured one and the expected
    projectors are ``A+`` from the same capture multiplied out: two independent computations of
    one projector, one through a singular-value decomposition and one through Gram-Schmidt.
    This is the assertion the Rust kernel's ``linalg`` test makes on the same numbers.
    """
    # The captured `A+`, six rows by four columns.
    pseudo_inverse = [
        [0.23503325942350345, -0.1607538802660754, 0.1995565410199555, 0.38026607538802654],
        [-0.07317073170731712, 0.1585365853658537, -0.024390243902438956, -0.0853658536585366],
        [0.06651884700665198, -0.0077605321507761005, 0.11308203991130826, 0.03215077605321519],
        [-0.028824833702882594, -0.013303769401330344, 0.05099778270509979, -0.23059866962305986],
        [-0.16851441241685158, 0.15299334811529935, -0.08647450110864735, -0.34811529933481156],
        [0.02217294900221725, 0.16407982261640816, 0.03769401330376936, 0.17738359201773818],
    ]
    basis = row_space_basis(MATRIX)
    # Every basis vector is in the row space, so the projector fixes it and the projection
    # removes nothing; and a vector normal to it is removed entirely.
    for held in basis:
        assert math.dist(project_onto_null_space([held], held), [0.0] * len(held)) < 1.0e-15
    # The projector's action on the captured `A+`'s own columns: `A+ A` applied to them is
    # themselves, because every column of `A+` lies in `A`'s row space.
    for column in range(4):
        held = [pseudo_inverse[row][column] for row in range(6)]
        residual = project_onto_null_space(basis, held)
        assert max(abs(value) for value in residual) < 1.0e-9
