"""``characterization.transfer_binary_interaction_parameters`` - a table moved between slates.

Spec: ``specs/models/characterization/transfer_binary_interaction_parameters.toml``. Oracle:
``validation/neqsim/captures/characterization_fluid_rebuild_probe.tsv``.

Mirrors ``crates/azoth-characterization/src/transfer_binary_interaction_parameters.rs``.

# What moves, and what stays

The reference's **nonzero** pairs are written onto the target - base components by name,
pseudo-components by position, and base-to-pseudo pairs alongside base-to-base ones. Everything
the reference does not state is left at the target's value, so the answer starts from the
target's own matrix: a databank-built target has fitted defaults, and dropping them silently is
the failure this shape exists to prevent.

# A missing name is skipped once and refused once

The class looks the reference's first name up in the target and ``continue``s when it is not
there, but its setter throws when the *second* name is missing. That asymmetry is reproduced.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any

from azoth import _models_gen
from azoth.core.errors import InvalidInputError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import TransferBinaryInteractionParametersResult
from azoth.core.units import input_records_to_si
from azoth.core.warnings import Warning

MODEL_ID = "characterization.transfer_binary_interaction_parameters"


def _components(rows: Sequence[Mapping[str, Any]]) -> list[tuple[str, bool]]:
    """The rows as `(name, pseudo)` pairs, in the order they are given."""
    return [(str(row["name"]), bool(row["pseudo_component"])) for row in rows]


def _check_matrix(field: str, matrix: Sequence[Sequence[float]], count: int) -> list[list[float]]:
    """A matrix is `N x N` for its component list, or the refusal naming it."""
    if len(matrix) != count or any(len(row) != count for row in matrix):
        raise InvalidInputError(
            field,
            f"an interaction matrix is {count} x {count}, one row and column per component",
        )
    return [[float(value) for value in row] for row in matrix]


def _pseudo_indices(components: Sequence[tuple[str, bool]]) -> list[int]:
    return [index for index, (_, pseudo) in enumerate(components) if pseudo]


def _index_of(components: Sequence[tuple[str, bool]], name: str) -> int | None:
    """The index of a component by name, or ``None``; a repeated name resolves to the first."""
    for index, (known, _) in enumerate(components):
        if known == name:
            return index
    return None


def _write_pair(matrix: list[list[float]], first: int, second: int, kij: float) -> None:
    """One pair, both triangles, the diagonal never touched."""
    matrix[first][second] = kij
    matrix[second][first] = kij


def transfer_binary_interaction_parameters(
    reference_components: Sequence[Mapping[str, Any]],
    target_components: Sequence[Mapping[str, Any]],
    reference_kij: Sequence[Sequence[float]],
    target_kij: Sequence[Sequence[float]],
) -> TransferBinaryInteractionParametersResult:
    """The reference's nonzero interaction pairs, written onto the target's own matrix.

    ``reference_components`` and ``target_components`` are lists of records, one per component
    in each fluid's own order, each stating its ``name`` and whether it is a
    ``pseudo_component``. ``reference_kij`` and ``target_kij`` are the two fluids' matrices,
    one row and one column per component.

    Base components are matched by name and pseudo-components by position. A pair the reference
    states as zero is not transferred, and every target value the reference does not overwrite
    survives.

    Raises:
        InvalidInputError: if either matrix is not `N x N` for its component list, or if a pair
            the reference states names a base component the target does not carry.
    """
    spec = _models_gen.model(MODEL_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    reference_rows: list[dict[str, Any]] = input_records_to_si(
        spec, "reference_components", reference_components
    )
    target_rows: list[dict[str, Any]] = input_records_to_si(
        spec, "target_components", target_components
    )
    reference = _components(reference_rows)
    target = _components(target_rows)
    reference_matrix = _check_matrix("reference_kij", reference_kij, len(reference))
    transferred = _check_matrix("target_kij", target_kij, len(target))

    reference_pseudo = _pseudo_indices(reference)
    target_pseudo = _pseudo_indices(target)

    for first_index, (first_name, first_pseudo) in enumerate(reference):
        if first_pseudo:
            continue
        first_in_target = _index_of(target, first_name)
        if first_in_target is None:
            continue
        for second_index, (second_name, second_pseudo) in enumerate(reference):
            if first_name == second_name:
                continue
            kij = reference_matrix[first_index][second_index]
            if kij == 0.0:
                continue
            if second_pseudo:
                position = reference_pseudo.index(second_index)
                if position >= len(target_pseudo):
                    continue
                second_in_target = target_pseudo[position]
            else:
                found = _index_of(target, second_name)
                if found is None:
                    raise InvalidInputError(
                        "target_components",
                        f"no component with name: {second_name} in the target fluid, and the "
                        "reference states a nonzero interaction with it",
                    )
                second_in_target = found
            _write_pair(transferred, first_in_target, second_in_target, kij)

    common = min(len(reference_pseudo), len(target_pseudo))
    for first in range(common):
        for second in range(first + 1, common):
            kij = reference_matrix[reference_pseudo[first]][reference_pseudo[second]]
            if kij != 0.0:
                _write_pair(transferred, target_pseudo[first], target_pseudo[second], kij)

    apply_checks(checks.derived, lambda _name: None, warnings)

    return TransferBinaryInteractionParametersResult(
        target_kij=tuple(tuple(row) for row in transferred),
        warnings=tuple(warnings),
    )
