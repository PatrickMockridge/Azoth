"""``characterization.set_heavy_tbp_to_plus`` - the row a fluid's plus end begins at.

Spec: ``specs/models/characterization/set_heavy_tbp_to_plus.toml``. Oracle:
``validation/neqsim/captures/set_heavy_tbp_probe.tsv``.

Mirrors ``crates/azoth-characterization/src/set_heavy_tbp_to_plus.rs`` line for line.

# The carbon number is in the name

``PlusCharacterize.setHeavyTBPtoPlus`` picks the row whose **name** carries the largest carbon
number, and reads it out of the string with a positional slice::

    name.substring(3, 4).equals("_") ? name.substring(1, 3) : name.substring(1, 2)

So ``C7_PC`` reads 7, ``C10_PC`` reads 10, and a name whose fourth character is not the underscore
reads a single character that will not parse. The class catches that along with everything else in
its loop and **skips** the row; this refuses it, because a plus end chosen from whatever happened
to parse is a plausible answer and not the caller's.

# The row below is the row before

The class reads ``component(compNumber - 1)`` for ``densLastTBP``, so it is the winner's neighbour
by **index** and not the row with the next carbon number. A winner at index zero reads index -1.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any

from azoth import _models_gen
from azoth.core.errors import InvalidInputError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import SetHeavyTbpToPlusResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning

MODEL_ID = "characterization.set_heavy_tbp_to_plus"

#: Where the class slices for the underscore that says the carbon number is two digits.
SUFFIX_MARKER = 3

#: The class's ``0.03``, which it adds to a density in g/cm3. In the kg/m3 this declares, 30.
INTERCEPT_STEP = 30.0


def carbon_number(name: str) -> int | None:
    """``setHeavyTBPtoPlus``'s parse: the carbon number a component's name carries.

    ``None`` for a name the class's own slice would throw on or fail to parse - what its ``catch``
    swallows, and what this refuses.
    """
    # The class's first act is `name.substring(3, 4)`, which throws on a name of fewer than four
    # characters; the handler that catches a number that will not parse catches that too.
    if len(name) <= SUFFIX_MARKER:
        return None
    digits = name[1:3] if name[SUFFIX_MARKER] == "_" else name[1:2]
    if not digits.isdigit():
        return None
    return int(digits)


def _dimensionless(spec: Mapping[str, Any], name: str, value: float | Q) -> float:
    """A declared dimensionless input as its SI magnitude, quantity or not.

    A mole fraction carries no unit, so a bare number is one - there is no dimension to be wrong
    about. The dimensioned columns stay strict: a bare number where a density is expected is the
    mistake this library exists to make impossible.
    """
    if isinstance(value, int | float):
        return float(value)
    return input_to_si(spec, name, value)


def set_heavy_tbp_to_plus(
    names: Sequence[str],
    molar_mass: Sequence[Q],
    mole_fraction: Sequence[float],
    normal_liquid_density: Sequence[Q],
) -> SetHeavyTbpToPlusResult:
    """The TBP row a fluid's plus end begins at, chosen by the carbon number in its name.

    ``names`` are the fluid's TBP rows in the order it carries them, and the three vectors beside
    them are those rows' data - each aligned with ``names`` entry for entry.

    Raises:
        InvalidInputError: if the vectors are not all the same length, if the slate is empty, if a
            name carries no carbon number, or if the winner is the first row - the class reads the
            row below it, and there is none.
    """
    spec = _models_gen.model(MODEL_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []
    apply_checks(checks.on_input, lambda _name: None, warnings)

    count = len(names)
    if count == 0:
        raise InvalidInputError(
            "names", "a fluid of no TBP rows has no row for a plus end to begin at"
        )
    if (
        len(molar_mass) != count
        or len(mole_fraction) != count
        or len(normal_liquid_density) != count
    ):
        raise InvalidInputError(
            "names",
            f"a fluid's TBP slate is four aligned vectors: {count} names against "
            f"{len(molar_mass)} molar masses, {len(mole_fraction)} mole fractions and "
            f"{len(normal_liquid_density)} normal liquid densities",
        )

    masses = [input_to_si(spec, "molar_mass", value) for value in molar_mass]
    fractions = [_dimensionless(spec, "mole_fraction", value) for value in mole_fraction]
    densities = [
        input_to_si(spec, "normal_liquid_density", value) for value in normal_liquid_density
    ]

    # `plusCompNumber < firstPlusNumber`, strictly: a tie keeps the row that came first, which is
    # the class's own comparison against its running maximum.
    winner: int | None = None
    best = 0
    for index, name in enumerate(names):
        number = carbon_number(name)
        if number is None:
            raise InvalidInputError(
                f"names[{index}]",
                f"`{name}` carries no carbon number, so the walk has no place to put it",
            )
        if winner is None or best < number:
            best = number
            winner = index
    assert winner is not None  # the slate is not empty and every name parsed

    if winner == 0:
        raise InvalidInputError(
            "names",
            f"the plus end is the first row, `{names[0]}`, and the row below it is what the class "
            f"reads for `densLastTBP`",
        )
    below = winner - 1

    apply_checks(checks.derived, lambda _name: None, warnings)

    return SetHeavyTbpToPlusResult(
        plus_index=float(winner),
        plus_name=names[winner],
        plus_molar_mass=from_si(masses[winner], "kg/mol"),
        plus_mole_fraction=fractions[winner],
        plus_density=from_si(densities[winner], "kg/m**3"),
        last_tbp_density=from_si(densities[below], "kg/m**3"),
        density_intercept=from_si(densities[below] + INTERCEPT_STEP, "kg/m**3"),
        warnings=tuple(warnings),
    )
