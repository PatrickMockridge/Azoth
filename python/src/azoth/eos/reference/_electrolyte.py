"""The molality-scale surface every electrolyte phase needs.

The Python mirror of ``crates/azoth-eos/src/electrolyte.rs``, held to the same constants
by the tests on both sides.

An activity-coefficient phase works in mole fractions and every electrolyte model works
in **molality**: `PhasePitzer`'s ionic strength is ``1/2 sum m_i z_i^2``, its osmotic
coefficient is a sum over ``m_i``, and ``ComponentGePitzer.getMolality`` divides the mole
count by the solvent's mass. This module is that change of scale.

**What is not here is the osmotic coefficient itself.** Pitzer's, Deshmukh-Mather's and
Kent-Eisenberg's are each a different formula over the same molalities, so only the
relation from a coefficient to a water activity is shared.

The molality NeqSim means: ``Component.getMolality`` is ``getMolarity / (density / 1e3)``
- moles per kilogram of *solution*, with its own comment saying ``// return mol/kg``.
**The two electrolyte models that read a molality override it**, and both return
``n_i / getSolventWeight()``, which is moles per kilogram of *water*. This module
implements the override, because it is the one the models that need it use.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from typing import NamedTuple

from azoth.core.errors import InvalidInputError

#: The substance a phase's molality is measured against.
#:
#: NeqSim matches the component's *name* in ``PhasePitzer.getSolventWeight``. Matched
#: without regard to case here, which is the rule the rest of this library resolves names
#: by.
SOLVENT = "water"


class ElectrolyteComposition(NamedTuple):
    """One mole of mixture's worth of the molality-scale composition.

    Per mole of mixture rather than per phase, because that is what the mole fractions
    determine: the total mole count cancels out of ``n_i / m_solvent``, so a phase's
    molalities are a function of its composition and not of its size.
    """

    #: Each component's molality ``n_i / m_solvent``, in mol/kg. **The solvent's own entry
    #: is its reciprocal molar mass**, about ``55.5 mol/kg`` for water - a model that
    #: summed this column as "the solute" would include the solvent.
    molality: tuple[float, ...]
    #: The mass of solvent per mole of mixture, in kg.
    solvent_mass: float
    #: ``I = 1/2 sum m_i z_i^2``, in mol/kg.
    ionic_strength: float


def solvent_mass(names: Sequence[str], x: Sequence[float], molar_mass: Sequence[float]) -> float:
    """The mass of solvent per mole of mixture, in kg: ``PhasePitzer.getSolventWeight``.

    Summed over every component named :data:`SOLVENT`, which is NeqSim's rule - a mixture
    with two water-named components weights both.
    """
    return sum(
        xi * mass
        for name, xi, mass in zip(names, x, molar_mass, strict=False)
        if name.lower() == SOLVENT
    )


def ionic_strength(molality: Sequence[float], charge: Sequence[float]) -> float:
    """``I = 1/2 sum m_i z_i^2``, in mol/kg: ``PhasePitzer.getIonicStrength``.

    Every component contributes, and a neutral one contributes zero because its charge is
    zero - so this needs no test of which components are ions.
    """
    return 0.5 * sum(m * z * z for m, z in zip(molality, charge, strict=False))


def composition(
    names: Sequence[str],
    x: Sequence[float],
    molar_mass: Sequence[float],
    charge: Sequence[float],
) -> ElectrolyteComposition:
    """The molality-scale composition of a phase, from its mole fractions.

    Raises:
        InvalidInputError: if the vectors disagree in length, or if the mixture carries no
            solvent.

    **The refusal is a divergence from NeqSim, and a deliberate one.**
    ``ComponentGePitzer.getMolality`` returns ``0.0`` when the phase has no solvent
    weight, so every molality is zero, the ionic strength is zero, and Pitzer evaluates as
    a solution of nothing - a finite answer with no symptom. A phase with no water is not a
    phase this surface describes, so it is refused.
    """
    n = len(names)
    if len(x) != n or len(molar_mass) != n or len(charge) != n:
        raise InvalidInputError(
            "components",
            f"an electrolyte composition needs one entry per component: {n} names, "
            f"{len(x)} mole fractions, {len(molar_mass)} molar masses and {len(charge)} "
            f"charges",
        )

    mass = solvent_mass(names, x, molar_mass)
    if not mass > 0.0:
        raise InvalidInputError(
            "components",
            "the mixture carries no solvent, so every molality would be zero and the ionic "
            "strength with it. NeqSim returns 0.0 here and evaluates the model as a solution "
            "of nothing; a phase without water is not one this surface describes",
        )

    molality = tuple(xi / mass for xi in x)
    return ElectrolyteComposition(
        molality=molality,
        solvent_mass=mass,
        ionic_strength=ionic_strength(molality, charge),
    )


def ln_water_activity(
    osmotic_coefficient: float, sum_molalities: float, water_molar_mass: float
) -> float:
    """``ln(a_w) = -phi M_w sum_i m_i``, with ``M_w`` in kg/mol.

    NeqSim's own comment (``ComponentGePitzer``, above its ``lnaw``): ``// Water activity:
    ln(a_w) = -phi * M_w * sumM / 1000``, where the ``1000`` converts its ``18.015 g/mol``
    to kg/mol. The two agree - the databank's water molar mass is ``0.018015 kg/mol``
    exactly - so this takes the molar mass rather than restating the constant.

    ``sum_molalities`` is the sum the model's own ``phi`` was built over, which is **not**
    the solvent: ``PhasePitzer`` sums the charged components and the neutral solutes it has
    interactions for, and returns ``phi = 1`` when that sum is below ``1e-12``.
    """
    return -osmotic_coefficient * water_molar_mass * sum_molalities


def water_activity_coefficient(ln_a_w: float, x_water: float) -> float:
    """``gamma_w = a_w / x_w``: NeqSim's conversion from a water activity to the
    mole-fraction activity coefficient the fugacity kernel takes.

    Raises:
        InvalidInputError: if the water mole fraction is not positive, where the ratio is
            not a number. NeqSim guards the same case and returns ``gamma = 1``; refusing
            is the house rule for a state the expression is not defined at.
    """
    if not x_water > 0.0:
        raise InvalidInputError(
            "components",
            "the water mole fraction is not positive, so `a_w / x_w` is not a number. "
            "NeqSim returns `gamma = 1` here, which is the ideal-solution value rather than "
            "the one this state has",
        )
    return math.exp(ln_a_w) / x_water
