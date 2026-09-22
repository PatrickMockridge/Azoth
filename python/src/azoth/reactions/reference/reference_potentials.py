"""``reactions.reference_potentials`` - the independent basis and the potentials from it.

Spec: ``specs/models/reactions/reference_potentials.toml``

The Python twin of ``crates/azoth-reactions/src/reference_potentials.rs``, written to
mirror it line for line. NeqSim builds this in three passes and the port keeps them:
``removeJunkReactions``, ``removeDependentReactions`` and ``calcReferencePotentials``.

# The three places a naive port goes wrong

**The reaction set is not a parameter.** It is what the source carries, filtered to what
the fluid can run.

**The rank tests see only the stoichiometry.** The matrix is one column wider than the
coefficients - the extra column holds ``-R T ln K`` - and the rank calls are on the
narrower matrix, so what is ranked is integers. See :mod:`._linalg`.

**The row order is part of the answer**, because both the reduction and the column choice
are greedy over the order the rows arrive in.

# What is refused rather than reproduced

``calcReferencePotentials`` has a deadlock fallback that seeds an unreachable component
with its Gibbs energy of formation. That is a databank column this library carries as
``vendored``, and no captured state reaches the branch, so it is refused with a named
reason instead.
"""

from __future__ import annotations

import math
from typing import Any

from azoth.core.errors import InvalidInputError
from azoth.core.result import ReferencePotentialsResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.reactions.reference import _linalg, _tables
from azoth.reactions.reference.equilibrium_constant import GAS_CONSTANT

#: The magnitude below which a stoichiometric coefficient counts as absent, from
#: `calcReferencePotentials`.
_COEFFICIENT_FLOOR = 1e-10


def reference_potentials(components: list[str], source: str, T: Q) -> ReferencePotentialsResult:
    """The standard-state reference potentials of a fluid's reactive components.

    Args:
        components: the fluid's reactive components, by name, **in the order the
            potentials are wanted back**. The order is part of the input, not a
            presentation detail: the column indices the basis is chosen over are these
            positions.
        source: which of the three reaction tables supplies the candidates.
        T: absolute temperature the constants are evaluated at.

    Raises:
        InvalidInputError: if ``source`` is unknown, or the propagation cannot reach a
            component, or a rank test meets a non-integral coefficient.

    Example:
        >>> import azoth
        >>> q = azoth.ureg.Quantity
        >>> r = reference_potentials(
        ...     ["CO2", "water", "OH-", "H3O+", "HCO3-", "CO3--"],
        ...     "standard", q(298.15, "K"),
        ... )
        >>> r.rank
        3
    """
    t = input_to_si(_spec(), "T", T)

    # Refused before the table is read, so an unknown source is reported as one.
    _tables.source_path(source)

    width = len(components)
    position = {name: i for i, name in enumerate(components)}

    rows: list[list[float]] = []
    rhs: list[float] = []
    loaded_positions: list[int] = []
    loaded = 0
    for candidate_row in _tables.reactions(source):
        if not candidate_row.use_reaction:
            continue
        at = loaded
        loaded += 1

        coefficients = _tables.stoichiometry(candidate_row.name)
        # Every reactant present: reactants are the negative coefficients.
        if not all(component in position for component, nu in coefficients if nu < 0.0):
            continue

        coefficient_row = [0.0] * width
        for component, coefficient in coefficients:
            if component in position:
                coefficient_row[position[component]] = coefficient

        k1, k2, k3, k4 = candidate_row.coefficients
        ln_k = k1 + k2 / t + k3 * math.log(t) + k4 * t

        rows.append(coefficient_row)
        rhs.append(-GAS_CONSTANT * t * ln_k)
        loaded_positions.append(at)

    # `removeDependentReactions`: greedy, and order-dependent.
    independent_rows: list[list[float]] = []
    independent_rhs: list[float] = []
    survivor_mask = [0.0] * loaded
    for coefficients_row, value, at in zip(rows, rhs, loaded_positions, strict=True):
        candidate = [*independent_rows, coefficients_row]
        if _linalg.rank_of_integer_matrix(candidate) > len(independent_rows):
            independent_rows = candidate
            independent_rhs.append(value)
            survivor_mask[at] = 1.0

    if not independent_rows:
        return ReferencePotentialsResult(
            potentials=tuple(from_si(0.0, "J/mol") for _ in range(width)),
            independent=(0.0,) * width,
            survivors=tuple(survivor_mask),
            rank=0,
            warnings=(),
        )

    # The greedy column selection: a column is independent when adding it raises the
    # rank of the rows so far.
    n_rows = len(independent_rows)
    current: list[list[float]] = []
    independent_columns: list[int] = []
    dependent_columns: list[int] = []

    for column in range(width):
        candidate = [
            ([*current[i]] if current else []) + [row[column]]
            for i, row in enumerate(independent_rows)
        ]
        current_rank = _linalg.rank_of_integer_matrix(current) if current else 0
        if _linalg.rank_of_integer_matrix(candidate) > current_rank:
            current = candidate
            independent_columns.append(column)
            if len(independent_columns) == n_rows:
                dependent_columns.extend(range(column + 1, width))
                break
        else:
            dependent_columns.append(column)

    if len(independent_columns) < n_rows:
        raise InvalidInputError(
            "reactions",
            f"the reaction basis has rank {len(independent_columns)} against {n_rows} "
            f"reaction(s), so no reference potentials follow from it. NeqSim returns an "
            f"empty array here and `ChemicalReactionOperations.calcChemRefPot` falls back "
            f"to each component's Gibbs energy of formation, a databank column this "
            f"library does not carry",
        )

    solved = _linalg.solve_lu(current, [-value for value in independent_rhs])

    potentials = [0.0] * width
    computed = [False] * width
    for i, column in enumerate(independent_columns):
        potentials[column] = solved[i]
        computed[column] = True

    # The propagation: a dependent component is computable once a surviving reaction
    # exists in which every *other* component it names is known.
    for _ in range(len(dependent_columns) * 2 + 1):
        progress = False
        for column in dependent_columns:
            if computed[column]:
                continue
            for r, row in enumerate(independent_rows):
                nu = row[column]
                if abs(nu) < _COEFFICIENT_FLOOR:
                    continue
                if not all(
                    j == column or abs(row[j]) <= _COEFFICIENT_FLOOR or computed[j]
                    for j in range(width)
                ):
                    continue
                sum_others = sum(row[j] * potentials[j] for j in range(width) if j != column)
                potentials[column] = (independent_rhs[r] - sum_others) / nu
                computed[column] = True
                progress = True
                break
        if not progress:
            stuck = [components[column] for column in dependent_columns if not computed[column]]
            if not stuck:
                break
            raise InvalidInputError(
                "reactions",
                f"the propagation cannot reach {', '.join(stuck)} from the independent "
                f"set: no surviving reaction has exactly one unknown. NeqSim falls back "
                f"to the component's Gibbs energy of formation here, which is the "
                f"`GIBBSENERGYOFFORMATION` column azoth carries as `vendored`",
            )

    return ReferencePotentialsResult(
        potentials=tuple(from_si(value, "J/mol") for value in potentials),
        independent=tuple(1.0 if column in independent_columns else 0.0 for column in range(width)),
        survivors=tuple(survivor_mask),
        rank=len(independent_columns),
        warnings=(),
    )


def _spec() -> dict[str, Any]:
    """The generated spec, imported here so the module has no import-time cycle."""
    from azoth._models_gen import model

    return model("reactions.reference_potentials")
