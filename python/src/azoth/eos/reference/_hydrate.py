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


def water_fugacity_coefficient(
    guests: list[HydrateGuest],
    ref_fugacities: list[float],
    t: float,
    p: float,
    water_index: int,
    structure: int,
    reference_water_fugacity: float,
) -> float:
    """The hydrate's water fugacity coefficient on one structure.

    ``f_w^ref`` is the reference water phase's fugacity and is **not** the pressure: NeqSim
    builds it from a one-component phase of the host's own class, so its water component is
    the host's and its fugacity is a real number. The caller states it, because only the
    caller knows which equation the fluid is.

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

    alpha_water = ref_fugacities[water_index]
    water_alpha_ref = math.log(reference_water_fugacity / alpha_water)
    return alpha_water * math.exp(val + _chemical_potential(structure, t, p) + water_alpha_ref) / p


def stable_structure(
    guests: list[HydrateGuest],
    ref_fugacities: list[float],
    t: float,
    p: float,
    water_index: int,
    reference_water_fugacity: float,
) -> tuple[int, float]:
    """The stable structure and its coefficient: the **lower** of the two."""
    first = water_fugacity_coefficient(
        guests, ref_fugacities, t, p, water_index, 0, reference_water_fugacity
    )
    second = water_fugacity_coefficient(
        guests, ref_fugacities, t, p, water_index, 1, reference_water_fugacity
    )
    return (0, first) if first <= second else (1, second)


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
