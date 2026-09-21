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

#: The reference pressure the empty-lattice vapour pressure is stated against, bar. NeqSim's
#: ``ThermodynamicConstantsInterface.referencePressure``.
_REFERENCE_PRESSURE_BAR = 1.01325

#: The empty structure's vapour pressure constants, structure I then II: Sloan (1990), as
#: NeqSim's ``ComponentHydrate.emptyHydrateVapourPressureConstant``.
_EMPTY_VAPOUR_PRESSURE = ((17.44, -6003.9), (17.332, -6017.6))

#: Avlonitis (1994)'s molar volume of the empty hydrate, ``(v0, k1, k2, k3)`` in cm3/mol and
#: per K, structure I then II. NeqSim's ``ComponentHydrate.getMolarVolumeHydrate``.
_HYDRATE_VOLUME = (
    (22.35, 3.1075e-4, 5.9537e-7, 1.3707e-10),
    (22.57, 1.9335e-4, 2.1768e-7, -1.4786e-10),
)

#: Which fitted hydrate model a calculation runs. NeqSim picks between them by the *system's
#: model name* rather than by anything about the fluid, so the name is a model-level choice.
PVTSIM = "pvtsim"
GUO_FINCH = "guo_finch"

#: The cavities of each type in one unit cell: structure I then II, small then large.
#: NeqSim's ``ComponentHydrate`` constructor, and the counts its ``46/54`` and ``136/160``
#: bounds come from - which are the *fully occupied* limits, not the composition at a state.
CAVITIES_PER_CELL: tuple[tuple[float, float], tuple[float, float]] = ((2.0, 6.0), (16.0, 8.0))

#: The water molecules in one unit cell.
WATER_PER_CELL: tuple[float, float] = (46.0, 136.0)


class HydrateGuest:
    """A substance's hydrate record: what the occupancy loops read, and no more."""

    __slots__ = ("former", "guo_finch_a", "guo_finch_b", "langmuir_a", "langmuir_b", "name")

    def __init__(
        self,
        name: str,
        langmuir_a: tuple[float, float, float, float],
        langmuir_b: tuple[float, float, float, float],
        guo_finch_a: tuple[float, float, float, float],
        guo_finch_b: tuple[float, float, float, float],
        former: bool,
    ) -> None:
        self.name = name
        self.langmuir_a = langmuir_a
        self.langmuir_b = langmuir_b
        self.guo_finch_a = guo_finch_a
        self.guo_finch_b = guo_finch_b
        self.former = former


class Hydration:
    """What a mixture needs to have a hydrate calculated for it."""

    __slots__ = ("guests", "water_index")

    def __init__(self, guests: list[HydrateGuest], water_index: int | None) -> None:
        self.guests = guests
        self.water_index = water_index


def langmuir(guest: HydrateGuest, model: str, structure: int, cavity: int, t: float) -> float:
    """``C = A/T exp(B/T)`` for one guest at one structure's cavity, on one model's pair."""
    index = structure * 2 + cavity
    a, b = (
        (guest.langmuir_a, guest.langmuir_b)
        if model == PVTSIM
        else (guest.guo_finch_a, guest.guo_finch_b)
    )
    return a[index] / t * math.exp(b[index] / t)


def empty_hydrate_vapour_pressure(structure: int, t: float) -> float:
    """The empty hydrate's vapour pressure on one structure, **in bar**.

    The Guo-Finch route's reference term, and the reason that route breaks the structure
    comparison: it is a function of the structure, so it does not cancel the way the PVTsim
    route's single reference does.
    """
    c0, c1 = _EMPTY_VAPOUR_PRESSURE[structure]
    return math.exp(c0 + c1 / t) * _REFERENCE_PRESSURE_BAR


def molar_volume_hydrate(structure: int, t: float) -> float:
    """The molar volume of the empty hydrate on one structure, m3/mol. Avlonitis (1994)."""
    v0, k1, k2, k3 = _HYDRATE_VOLUME[structure]
    dt = t - T_REFERENCE
    return v0 * (1.0 + k1 * dt + k2 * dt * dt + k3 * dt * dt * dt) / 1.0e6


def _reference_term(
    model: str, structure: int, t: float, p: float, reference_water_fugacity: float
) -> float:
    """What multiplies ``exp(exponent)`` to make the coefficient."""
    if model == PVTSIM:
        return reference_water_fugacity / p
    p_empty = empty_hydrate_vapour_pressure(structure, t) * _PA_PER_BAR
    return p_empty * math.exp(molar_volume_hydrate(structure, t) / (R * t) * (p - p_empty)) / p


def occupancy(
    guests: list[HydrateGuest],
    ref_fugacities: list[float],
    model: str,
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
            denominator += langmuir(guest, model, structure, cavity, t) * fugacity / _PA_PER_BAR
    return [
        langmuir(guest, model, structure, cavity, t) * fugacity / _PA_PER_BAR / denominator
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
    model: str,
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
        occupied = sum(occupancy(guests, ref_fugacities, model, structure, cavity, t))
        if occupied >= 1.0:
            raise OutOfRangeError(
                "occupancy",
                occupied,
                f"cavity {cavity} of structure {structure + 1} is fully occupied at {t} K, "
                f"where the cavity sum's logarithm has no value",
            )
        val += per_water * math.log(1.0 - occupied)
    # The Guo-Finch route's reference is the empty lattice's own vapour pressure, carried by
    # `_reference_term`, so it adds no chemical-potential constant.
    return val + (_chemical_potential(structure, t, p) if model == PVTSIM else 0.0)


def water_fugacity_coefficient(
    guests: list[HydrateGuest],
    ref_fugacities: list[float],
    model: str,
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
    return math.exp(exponent(guests, ref_fugacities, model, t, p, structure)) * _reference_term(
        model, structure, t, p, reference_water_fugacity
    )


def stable_structure(
    guests: list[HydrateGuest],
    ref_fugacities: list[float],
    model: str,
    t: float,
    p: float,
    reference_water_fugacity: float,
) -> tuple[int, float]:
    """The stable structure and its coefficient: the **lower** of the two.

    **The coefficients are compared, not the exponents.** The PVTsim route's reference is one
    number for both structures, so comparing exponents happens to agree there; the Guo-Finch
    route's is the empty lattice's vapour pressure on that structure, which does not cancel.
    """
    first = water_fugacity_coefficient(
        guests, ref_fugacities, model, t, p, 0, reference_water_fugacity
    )
    second = water_fugacity_coefficient(
        guests, ref_fugacities, model, t, p, 1, reference_water_fugacity
    )
    return (0, first) if first <= second else (1, second)


def composition(
    guests: list[HydrateGuest],
    ref_fugacities: list[float],
    model: str,
    structure: int,
    t: float,
    water_index: int,
) -> list[float]:
    """The hydrate's mole fractions at a state, from **both** cavity types of a structure.

    A cell of structure I holds 46 waters with two small and six large cavities, and one of
    structure II 136 with sixteen and eight. A cavity of type ``cav`` holds guest ``i`` with
    probability ``theta_icav``, so the cell's guest count is the cavities' own weighted sum
    of the occupancies, and each guest's is the one weighted by its own.

    **The second cavity is the whole point.** A guest occupies one cage type at a time and
    its occupancy in each is a different function of the state, so a composition built from one
    cavity type leaves the other's guests at zero. For structure I that is the large
    ``5^12 6^2`` cage, where most of the ethane and propane sit: the phase's fractions would not
    sum to one, and its bound would not move with the temperature.
    """
    counts = [0.0] * len(guests)
    total_guests = 0.0
    for cavity in range(2):
        per_cell = CAVITIES_PER_CELL[structure][cavity]
        for index, occupied in enumerate(
            occupancy(guests, ref_fugacities, model, structure, cavity, t)
        ):
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
                guo_finch_a=entry.hydrate_guo_finch_a,
                guo_finch_b=entry.hydrate_guo_finch_b,
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
