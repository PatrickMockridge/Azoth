"""``reactions.reference_potentials`` - the independent basis and the potentials from it.

Spec: ``specs/models/reactions/reference_potentials.toml``

The Python twin of ``crates/azoth-reactions/src/reference_potentials.rs``, written to
mirror it line for line. NeqSim builds this in three passes and the port keeps them:
``removeJunkReactions``, ``removeDependentReactions`` and ``calcReferencePotentials``.

# The three places a naive port goes wrong

**The reaction set is not a parameter.** It is what the source carries, filtered to what
the fluid can run: a reaction is kept when every reactant it names is present **or** every
product is, which is ``reactantsContains``' fall-through. The second half was missing from
this port until ``PitzerStrictnessProbe`` measured a fluid it changes - ``MDEAprot`` names
``MDEA+`` as a reactant, which an ``MDEA``/``CO2``/water fluid does not carry, and ``MDEA``
and ``H3O+`` as products, which it does.

**The pitzer source refuses unvalidated active rows**, after both removals, so the refusal
names the survivors and not the table.

**The rank tests see only the stoichiometry.** The matrix is one column wider than the
coefficients - the extra column holds ``-R T ln K`` - and the rank calls are on the
narrower matrix, so what is ranked is integers. See :mod:`._linalg`.

**The row order is part of the answer**, because both the reduction and the column choice
are greedy over the order the rows arrive in.

# The two branches that answer from the databank

``calcReferencePotentials`` has two ways to answer from a component's Gibbs energy of
formation rather than from the reaction set, and both are reproduced:

* when the component rank falls below the reaction count it returns null, and
  ``ChemicalReactionOperations.calcChemRefPot`` reads that null as *every* component's
  Gibbs energy of formation (``ChemicalReactionOperations.java:513-519``);
* when the propagation deadlocks it seeds the first uncomputed dependent component with
  that same number and carries on, and a final sweep seeds whatever is left
  (``ChemicalReactionList.java:526-557``).

The number is :func:`_tables.formation_properties`'s, out of ``GIBBSENERGYOFFORMATION``,
and it is taken **raw and not negated**, so a seeded potential is not on the same footing
as a solved one - which is why ``independent`` is reported. **The fallback is not an edge
case**: of the 14,333 subsets of the species the three tables' loaded reactions name, 11,837
answer - 3,320 of them seeding at least one component - and 2,496 are refused by the pitzer
source's evidence gate, which the Rust twin's sweep test asserts.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
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
    kept_names: list[str] = []
    rhs: list[float] = []
    statuses: list[str | None] = []
    loaded_positions: list[int] = []
    loaded = 0
    for candidate_row in _tables.reactions(source):
        if not candidate_row.use_reaction:
            continue
        at = loaded
        loaded += 1

        coefficients = _tables.stoichiometry(candidate_row.name)

        # `reactantsContains`: **all reactants present, or all products present.**
        #
        # Not "every reactant", which is what this port had and what the spec said: the
        # class falls through to the *products* when a reactant is missing, so a reaction
        # whose products the fluid can hold is kept even though it cannot run forwards.
        # `MDEAprot` is the measured case - `MDEA+` is a reactant it lacks, `MDEA` and
        # `H3O+` are products it has, so NeqSim keeps it and the basis is one larger.
        #
        # Reactants are the negative coefficients and products everything else, including a
        # zero, which is how the class splits them. A side with no names does not satisfy
        # the rule, because the class's `test` starts false and only a completed scan sets
        # it.
        reactant_side = [name for name, nu in coefficients if nu < 0.0]
        product_side = [name for name, nu in coefficients if nu >= 0.0]
        reactants_present = bool(reactant_side) and all(name in position for name in reactant_side)
        products_present = bool(product_side) and all(name in position for name in product_side)
        if not reactants_present and not products_present:
            continue

        coefficient_row = [0.0] * width
        for component, coefficient in coefficients:
            if component in position:
                coefficient_row[position[component]] = coefficient

        k1, k2, k3, k4 = candidate_row.coefficients
        ln_k = k1 + k2 / t + k3 * math.log(t) + k4 * t

        rows.append(coefficient_row)
        kept_names.append(candidate_row.name)
        rhs.append(-GAS_CONSTANT * t * ln_k)
        statuses.append(candidate_row.validation_status)
        loaded_positions.append(at)

    # `removeDependentReactions`: greedy, and order-dependent.
    independent_rows: list[list[float]] = []
    independent_rhs: list[float] = []
    survivor_mask = [0.0] * loaded
    survivor_names: list[str] = []
    survivor_statuses: list[str | None] = []
    for coefficients_row, name, status, value, at in zip(
        rows, kept_names, statuses, rhs, loaded_positions, strict=True
    ):
        candidate = [*independent_rows, coefficients_row]
        if _linalg.rank_of_integer_matrix(candidate) > len(independent_rows):
            independent_rows = candidate
            independent_rhs.append(value)
            survivor_mask[at] = 1.0
            survivor_names.append(name)
            survivor_statuses.append(status)

    # `requireValidatedEvidenceForActiveReactions`: the pitzer source is the one that
    # requires evidence, and it requires it of the *survivors* - the check runs after both
    # removals upstream. Measured on `SystemPitzer`: 43 of the table's 47 rows are not
    # `VALIDATED`, a plain CO2/water fluid survives on the four that are, and an
    # `MDEA`/`CO2`/water fluid trips it on `MDEAprot` alone.
    if source == "pitzer":
        unvalidated = sorted(
            name
            for name, status in zip(survivor_names, survivor_statuses, strict=True)
            if (status or "").strip().upper() != "VALIDATED"
        )
        if unvalidated:
            raise InvalidInputError(
                "components",
                f"the `{source}` source rejects unvalidated active rows, and the reactions "
                f"without validated evidence are {', '.join(unvalidated)}",
            )

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

    # **A rank-deficient basis is an answer and not an error.** `calcReferencePotentials`
    # returns null when the component rank falls below the reaction count, and
    # `ChemicalReactionOperations.calcChemRefPot` reads that null as every component's
    # Gibbs energy of formation. Nothing is marked independent there, because nothing was
    # solved for; `rank` still reports the rank it fell short of.
    rank = len(independent_columns)
    potentials = [0.0] * width
    computed = [False] * width
    independent_mask = [0.0] * width

    if rank < n_rows:
        for column in range(width):
            potentials[column] = _formation_seed(components, column)
            computed[column] = True
    elif n_rows > 0:
        solved = _linalg.solve_lu(current, [-value for value in independent_rhs])
        for i, column in enumerate(independent_columns):
            potentials[column] = solved[i]
            computed[column] = True
            independent_mask[column] = 1.0

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
            if progress:
                continue
            # **The deadlock fallback, reproduced.** NeqSim seeds the first uncomputed
            # dependent component with its Gibbs energy of formation and carries on, one
            # per round that made no progress (`ChemicalReactionList.java:526-541`).
            seed = next((column for column in dependent_columns if not computed[column]), None)
            if seed is None:
                break
            potentials[seed] = _formation_seed(components, seed)
            computed[seed] = True

        # The final sweep (`:546-557`): whatever the rounds above left uncomputed takes
        # the same seed.
        for column in dependent_columns:
            if not computed[column]:
                potentials[column] = _formation_seed(components, column)

    return ReferencePotentialsResult(
        potentials=tuple(from_si(value, "J/mol") for value in potentials),
        independent=tuple(independent_mask),
        survivors=tuple(survivor_mask),
        rank=rank,
        warnings=(),
    )


def _formation_seed(components: Sequence[str], column: int) -> float:
    """NeqSim's seed for a component the propagation cannot reach.

    **The component's Gibbs energy of formation, raw and not negated**
    (``ChemicalReactionList.java:531``). A component with no databank row is refused
    rather than seeded with zero: NeqSim's ``Component`` field defaults to zero when the
    row is absent, and a component with no row cannot be built there at all, while here
    the name is a caller's string.
    """
    name = components[column]
    row = _tables.formation_properties(name)
    if row is None:
        raise InvalidInputError(
            "reactions",
            f"the reference potential of {name} cannot be reached by propagation and the "
            f"component databank has no row for it, so there is no Gibbs energy of "
            f"formation to seed it with as NeqSim's fallback does",
        )
    return row.gibbs_energy_of_formation


def _spec() -> dict[str, Any]:
    """The generated spec, imported here so the module has no import-time cycle."""
    from azoth._models_gen import model

    return model("reactions.reference_potentials")
