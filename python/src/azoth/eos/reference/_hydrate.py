"""The van der Waals-Platteeuw hydrate fugacity: the Python reference.

The mirror of ``crates/azoth-eos/src/hydrate.rs``, which states the physics and why the
route is the fitted one - NeqSim's ``ComponentHydratePVTsim``, whose guest constants are
``C = A/T exp(B/T)`` from the ``HydrateA*``/``B*`` columns, rather than the base class's
Kihara integral over the three ``*HYDRATE`` Lennard-Jones columns. It also states the two
places a port goes wrong silently: the ``C f`` product is taken in **bar**, and the
coefficient's ``f_w^ref`` is a *computed* fugacity rather than the pressure.
"""

from __future__ import annotations

import math

from azoth.core.errors import InvalidInputError, OutOfRangeError

#: NeqSim's ``ThermodynamicConstantsInterface.R``, which ``calcDeltaChemPot`` divides by.
R = 8.3144621

#: The reference temperature of the chemical-potential constant, K.
T_REFERENCE = 273.15

#: How many cavities of each type a water molecule belongs to: structure I then II, small
#: then large. NeqSim's ``ComponentHydrate`` constructor, and hardcoded rather than read.
CAVITIES_PER_WATER: tuple[tuple[float, float], tuple[float, float]] = (
    (1.0 / 23.0, 3.0 / 23.0),
    (2.0 / 17.0, 1.0 / 17.0),
)

#: The chemical-potential constants ``(dGf, dHf, Cp, dV)`` of ``calcDeltaChemPot``,
#: structure I then II, in J/mol, J/mol, J/(mol K) and m3/mol.
_CHEMICAL_POTENTIAL: tuple[tuple[float, float, float, float], ...] = (
    (1264.0, -4858.0, -39.16, 4.6e-6),
    (883.0, -5201.0, -39.16, 5.0e-6),
)

#: NeqSim works its fugacities in bara and this states its own in pascals.
_PA_PER_BAR = 1.0e5

#: The cavities of each type in one unit cell: structure I then II, small then large.
#: NeqSim's ``ComponentHydrate`` constructor, and the counts its ``46/54`` and ``136/160``
#: bounds come from - which are the *fully occupied* limits, not the composition at a state.
CAVITIES_PER_CELL: tuple[tuple[float, float], tuple[float, float]] = ((2.0, 6.0), (16.0, 8.0))

#: The water molecules in one unit cell.
WATER_PER_CELL: tuple[float, float] = (46.0, 136.0)


class HydrateGuest:
    """A substance's hydrate record: what the occupancy loops read, and no more."""

    __slots__ = ("former", "langmuir_a", "langmuir_b", "name")

    def __init__(
        self,
        name: str,
        langmuir_a: tuple[float, float, float, float],
        langmuir_b: tuple[float, float, float, float],
        former: bool,
    ) -> None:
        self.name = name
        self.langmuir_a = langmuir_a
        self.langmuir_b = langmuir_b
        self.former = former


class Hydration:
    """What a mixture needs to have a hydrate calculated for it."""

    __slots__ = ("guests", "water_index")

    def __init__(self, guests: list[HydrateGuest], water_index: int | None) -> None:
        self.guests = guests
        self.water_index = water_index


def langmuir(guest: HydrateGuest, structure: int, cavity: int, t: float) -> float:
    """``C = A/T exp(B/T)`` for one guest at one structure's cavity."""
    index = structure * 2 + cavity
    return guest.langmuir_a[index] / t * math.exp(guest.langmuir_b[index] / t)


def occupancy(
    guests: list[HydrateGuest],
    ref_fugacities: list[float],
    structure: int,
    cavity: int,
    t: float,
) -> list[float]:
    """One cavity type's occupancy by each guest, ``C f_i / (1 + sum_j C_j f_j)``.

    **The product ``C f`` is taken in bar**, whatever unit the caller states its fugacities
    in: the Langmuir constant is fitted against NeqSim's bar-valued ones and carries the
    reciprocal. Feeding it pascals saturates every cage.
    """
    denominator = 1.0
    for guest, fugacity in zip(guests, ref_fugacities, strict=True):
        if guest.former:
            denominator += langmuir(guest, structure, cavity, t) * fugacity / _PA_PER_BAR
    return [
        langmuir(guest, structure, cavity, t) * fugacity / _PA_PER_BAR / denominator
        if guest.former
        else 0.0
        for guest, fugacity in zip(guests, ref_fugacities, strict=True)
    ]


def _chemical_potential(structure: int, t: float, p: float) -> float:
    """``calcDeltaChemPot``: the change in chemical potential on forming the hydrate."""
    dgf, dhf, cp, dvolume = _CHEMICAL_POTENTIAL[structure]
    return (
        dgf / R / T_REFERENCE
        - (
            -dhf * (1.0 / R / t - 1.0 / R / T_REFERENCE)
            + cp / R * math.log(t / T_REFERENCE)
            + cp * T_REFERENCE / R * (1.0 / t - 1.0 / T_REFERENCE)
        )
        + dvolume / R / t * p
    )


def exponent(
    guests: list[HydrateGuest],
    ref_fugacities: list[float],
    t: float,
    p: float,
    structure: int,
) -> float:
    """The exponent of a structure's water fugacity coefficient: the part that is its own.

    ``f_w^hydrate/P = (f_w^ref/P) exp(sum_cav n_cav ln(1 - sum_j theta_j) + dMu)``.

    **The fluid's own water fugacity cancels out of this**, which is why it is not an
    argument. NeqSim writes the coefficient as
    ``(f_w^fluid/P) exp(sum + dMu + ln(f_w^ref/f_w^fluid))``, where the fluid's water fugacity
    is multiplied in and divided out again inside the logarithm - it is not what the
    hydrate's water fugacity depends on, and ``f_w^ref`` is. Written out that way the
    coefficient is ``0 * inf`` on a fluid with no water, which is a state a hydrate
    fraction's *bound* is, so the cancellation is taken here.

    Raises:
        OutOfRangeError: if a cavity type is fully occupied, where the sum's logarithm has no
            value. A refusal rather than a clamp: a saturation of one means the occupancy
            model has left the region it was fitted to.
    """
    val = 0.0
    for cavity, per_water in enumerate(CAVITIES_PER_WATER[structure]):
        occupied = sum(occupancy(guests, ref_fugacities, structure, cavity, t))
        if occupied >= 1.0:
            raise OutOfRangeError(
                "occupancy",
                occupied,
                f"cavity {cavity} of structure {structure + 1} is fully occupied at {t} K, "
                f"where the cavity sum's logarithm has no value",
            )
        val += per_water * math.log(1.0 - occupied)
    return val + _chemical_potential(structure, t, p)


def water_fugacity_coefficient(
    guests: list[HydrateGuest],
    ref_fugacities: list[float],
    t: float,
    p: float,
    structure: int,
    reference_water_fugacity: float,
) -> float:
    """The hydrate's water fugacity coefficient on one structure, ``f_w^hydrate/P``.

    ``f_w^ref`` is the reference water phase's fugacity and is **not** the pressure: NeqSim
    builds it from a one-component phase of the host's own class, so its water component is
    the host's and its fugacity is a real number. The caller states it, because only the
    caller knows which equation the fluid is.

    Raises:
        OutOfRangeError: from the structure's cavity sum.
    """
    return (
        math.exp(exponent(guests, ref_fugacities, t, p, structure)) * reference_water_fugacity / p
    )


def stable_structure(
    guests: list[HydrateGuest],
    ref_fugacities: list[float],
    t: float,
    p: float,
    reference_water_fugacity: float,
) -> tuple[int, float]:
    """The stable structure and its coefficient: the **lower** of the two.

    The reference term is the same for both - it is the *reference fluid's* fugacity over the
    pressure, and neither depends on the structure - so the comparison is the exponents'.
    """
    first = exponent(guests, ref_fugacities, t, p, 0)
    second = exponent(guests, ref_fugacities, t, p, 1)
    structure, chosen = (0, first) if first <= second else (1, second)
    return structure, math.exp(chosen) * reference_water_fugacity / p


def composition(
    guests: list[HydrateGuest],
    ref_fugacities: list[float],
    structure: int,
    t: float,
    water_index: int,
) -> list[float]:
    """The hydrate's mole fractions at a state, from **both** cavity types of a structure.

    A cell of structure I holds 46 waters with two small and six large cavities, and one of
    structure II 136 with sixteen and eight. A cavity of type ``cav`` holds guest ``i`` with
    probability ``theta_icav``, so the cell's guest count is the cavities' own weighted sum
    of the occupancies, and each guest's is the one weighted by its own.

    **The second cavity is the whole point.** NeqSim's ``updateHydrateComposition``
    distributes the guests by cavity 0 alone - the small cage - which for structure I is where
    a guest mostly is not: its ethane and propane fractions come out as the *feed's* own,
    because the loop writes a zero there and the field keeps what it had, and the phase's
    fractions then sum to ``1.1201`` rather than one.
    """
    counts = [0.0] * len(guests)
    total_guests = 0.0
    for cavity in range(2):
        per_cell = CAVITIES_PER_CELL[structure][cavity]
        for index, occupied in enumerate(occupancy(guests, ref_fugacities, structure, cavity, t)):
            counts[index] += per_cell * occupied
            total_guests += per_cell * occupied
    water = WATER_PER_CELL[structure]
    total = water + total_guests
    fractions = [count / total for count in counts]
    fractions[water_index] = water / total
    return fractions


def hydration_for(names: list[str]) -> Hydration:
    """The guest records for a set of names, refused if the fluid cannot form a hydrate.

    Raises:
        InvalidInputError: if a name is not in the databank, if there is no water, or if
            nothing in the mixture is a hydrate former.
    """
    from azoth.eos import components as databank

    guests = []
    for name in names:
        entry = databank.entry(name)
        guests.append(
            HydrateGuest(
                name=entry.name,
                langmuir_a=entry.hydrate_langmuir_a,
                langmuir_b=entry.hydrate_langmuir_b,
                former=entry.hydrate_former,
            )
        )

    water_index = next((i for i, g in enumerate(guests) if g.name == "water"), None)
    if water_index is None:
        raise InvalidInputError(
            "components",
            "no water: a hydrate is water's, so a fluid without it has no formation temperature",
        )
    if not any(guest.former for guest in guests):
        raise InvalidInputError(
            "components",
            "nothing in this mixture is a hydrate former, so no cage would have a guest and "
            "the hydrate's fugacity would be its empty one",
        )
    return Hydration(guests=guests, water_index=water_index)
