"""The Python twin of ``reactions.reference_potentials``, on the two rules a case cannot pin.

The cases run through both kernels, but two things the model does are *refusals* or
*fallbacks* rather than numbers a case can state:

* the pitzer source's evidence gate, which throws for an active row that is not
  ``VALIDATED`` - NeqSim throws too, so there is no result to compare;
* ``reactantsContains``' product fallback, which the Rust sweep measures over 14,333
  subsets and which one case per kernel would only sample.

The oracle for both is ``validation/neqsim/PitzerStrictnessProbe.java`` and
``validation/neqsim/ReferenceFallbackProbe.java``, and the numbers below are theirs.
"""

from __future__ import annotations

import pytest

import azoth
from azoth.core.errors import InvalidInputError
from azoth.reactions.reference.reference_potentials import reference_potentials

#: The component list ``PitzerStrictnessProbe`` reports for an ``MDEA``/``CO2``/water fluid
#: *after* ``chemicalReactionInit`` has added the ions its chemistry needs. The gate fires on
#: that augmented pass, not on the feed, so this is the list that reproduces it.
AUGMENTED = ["MDEA", "water", "CO2", "OH-", "H3O+", "HCO3-"]

#: The CO2/water fluid the four ``VALIDATED`` rows carry.
CO2_WATER = ["CO2", "water", "OH-", "H3O+", "HCO3-", "CO3--"]


def test_the_pitzer_source_refuses_an_unvalidated_active_row() -> None:
    """``MDEAprot`` survives on this fluid and is not ``VALIDATED``, so the source refuses."""
    with pytest.raises(InvalidInputError, match="MDEAprot"):
        reference_potentials(AUGMENTED, "pitzer", azoth.ureg.Quantity(313.15, "K"))


def test_only_the_pitzer_source_requires_evidence() -> None:
    """The standard source carries no evidence column, so the same list is not refused."""
    result = reference_potentials(AUGMENTED, "standard", azoth.ureg.Quantity(298.15, "K"))
    # Four reactions, not three: `MDEAprot` is kept by the product fallback.
    assert sum(result.survivors) == 4
    assert result.rank == 4


def test_the_gate_does_not_fire_on_the_fluids_the_validated_rows_carry() -> None:
    """CO2/water survives on `CO2water`, `waterreac` and `carbonate` alone."""
    result = reference_potentials(CO2_WATER, "pitzer", azoth.ureg.Quantity(298.15, "K"))
    assert sum(result.survivors) == 3
    assert result.rank == 3
