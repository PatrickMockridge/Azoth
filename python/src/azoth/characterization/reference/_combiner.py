"""``PseudoComponentCombiner``'s arithmetic, less the fluid it usually takes.

Mirrors ``crates/azoth-characterization/src/combiner.rs`` line for line. The class is 1841 lines
and every public method takes or returns a ``SystemInterface``; what this module holds is the part
that is a function of a table - extraction, the equal-mass cut points, the distribution of each
fluid's own rows onto them, and the blend of the per-fluid profiles.

**The arithmetic is the class's, down to the summation order.** A boundary is the sorting key of
the row whose running mass first reaches a target, not an interpolation; a row joins a group by a
strict ``key > boundary``; a group's molar mass is its accumulated mass over its accumulated
moles. The boiling point and the three critical properties are weighted by fluid *mass* fraction
times *mole* fraction within that fluid, where the other thirteen are mass-weighted.

# The row carries ``getCpD``'s number

``Component.getCpD`` returns the stored coefficient plus ``1e-10``. The class reads a source
through that getter and assigns through ``setCpD``, so a blend's coefficient is the mean of the
getter's values and the rebuilt component's getter adds the offset again. A row states the
getter's number, and the answer states the getter's number of the mean - reproduced, not corrected.
"""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, fields
from typing import Any, cast

from azoth.characterization.reference.delumping import sorting_key
from azoth.core.errors import InvalidInputError
from azoth.core.result import CharacterizationComponent
from azoth.core.units import Q, from_si

#: ``PseudoComponentCombiner.MASS_TOLERANCE``, on an accumulated mass or mole amount.
MASS_TOLERANCE = 1.0e-12

#: ``Component.getCpD``'s own offset: the getter returns the stored coefficient plus this.
CP_D_GETTER_OFFSET = 1.0e-10

#: The units the answer states each property in, as the spec declares them. A field absent here
#: is dimensionless and crosses as a bare number.
PROPERTY_UNITS = {
    "molar_mass": "kg/mol",
    "density": "kg/m**3",
    "normal_boiling_point": "K",
    "critical_temperature": "K",
    "critical_pressure": "Pa",
    "critical_volume": "m**3/mol",
    "critical_viscosity": "Pa*s",
    "triple_point_temperature": "K",
    "heat_of_fusion": "J/mol",
    "ideal_gas_enthalpy_of_formation": "J/mol",
    "cp_a": "J/(mol*K)",
    "cp_b": "J/(mol*K**2)",
    "cp_c": "J/(mol*K**3)",
    "cp_d": "J/(mol*K**4)",
}


@dataclass(frozen=True, slots=True)
class Properties:
    """The seventeen properties a contribution carries beyond its moles, molar mass and density.

    **An absent property is a NaN**, which is how NeqSim spells one in its row: every
    ``Double.isFinite`` guard in the class is a test for it.
    """

    normal_boiling_point: float
    critical_temperature: float
    critical_pressure: float
    acentric_factor: float
    critical_volume: float
    racket_z: float
    racket_z_cpa: float
    parachor: float
    critical_viscosity: float
    triple_point_temperature: float
    heat_of_fusion: float
    ideal_gas_enthalpy_of_formation: float
    cp_a: float
    cp_b: float
    cp_c: float
    cp_d: float
    attractive_m: float

    @classmethod
    def of(cls, row: Mapping[str, Any]) -> Properties:
        """The row's properties, with an absent field as a NaN."""
        return cls(
            **{
                field.name: math.nan if row.get(field.name) is None else float(row[field.name])
                for field in fields(cls)
            }
        )


PROPERTY_FIELDS: tuple[str, ...] = tuple(field.name for field in fields(Properties))


def _properties_from_values(values: Sequence[float]) -> Properties:
    """The properties an iterable in :data:`PROPERTY_FIELDS`' order states."""
    return Properties(**dict(zip(PROPERTY_FIELDS, values, strict=True)))


@dataclass(frozen=True, slots=True)
class Contribution:
    """One pseudo-component, as ``extractComponents`` lifts it."""

    moles: float
    molar_mass: float
    density: float
    properties: Properties
    mass: float

    @classmethod
    def of(cls, index: int, row: Mapping[str, Any]) -> Contribution:
        """The row as a contribution, or the refusal for a pseudo-component that cannot be one."""
        if row.get("molar_mass") is None:
            raise InvalidInputError(
                f"components[{index}].molar_mass",
                "a pseudo-component carries its own molar mass, and this row states none",
            )
        if row.get("density") is None:
            raise InvalidInputError(
                f"components[{index}].density",
                "a pseudo-component carries its own normal liquid density, and this row "
                "states none",
            )
        moles = float(row["moles"])
        molar_mass = float(row["molar_mass"])
        return cls(
            moles=moles,
            molar_mass=molar_mass,
            density=float(row["density"]),
            properties=Properties.of(row),
            mass=moles * molar_mass,
        )

    def sorting_key(self) -> float:
        """The boiling point, or the molar mass without a finite positive one."""
        return sorting_key(self.properties.normal_boiling_point, self.molar_mass)

    def scaled_by_mass(self, factor: float) -> Contribution:
        """The contribution with its mole amount - and hence its mass - scaled."""
        moles = self.moles * factor
        return Contribution(
            moles=moles,
            molar_mass=self.molar_mass,
            density=self.density,
            properties=self.properties,
            mass=moles * self.molar_mass,
        )


@dataclass(slots=True)
class Extraction:
    """One fluid's extraction: its base components by name, and its pseudo-components by key."""

    base: dict[str, float]
    pseudo: list[Contribution]


def extract(rows: Sequence[Mapping[str, Any]]) -> list[Extraction]:
    """``extractComponents``, for every fluid the rows name, in first-appearance order.

    A row whose mole amount is not positive is skipped - the class's own ``!(moles > 0)`` - so a
    zero is a row that is not on the grid rather than a row that contributes nothing.
    """
    fluids: list[float] = []
    for index, row in enumerate(rows):
        fluid = float(row["fluid"])
        if not (math.isfinite(fluid) and fluid >= 0.0 and fluid.is_integer()):
            raise InvalidInputError(
                f"components[{index}].fluid",
                "a fluid index is a whole number from zero, and this row's is not",
            )
        if fluid not in fluids:
            fluids.append(fluid)

    extractions = [Extraction(base={}, pseudo=[]) for _ in fluids]
    for index, row in enumerate(rows):
        if not float(row["moles"]) > 0.0:
            continue
        extraction = extractions[fluids.index(float(row["fluid"]))]
        if row["pseudo_component"]:
            extraction.pseudo.append(Contribution.of(index, row))
        else:
            name = str(row["name"])
            extraction.base[name] = extraction.base.get(name, 0.0) + float(row["moles"])

    for extraction in extractions:
        extraction.pseudo.sort(key=Contribution.sorting_key)
    return extractions


def merge_bases(extractions: Sequence[Extraction]) -> dict[str, float]:
    """The base components of every fluid, merged by name in the fluids' own order."""
    merged: dict[str, float] = {}
    for extraction in extractions:
        for name, moles in extraction.base.items():
            merged[name] = merged.get(name, 0.0) + moles
    return merged


def quantile_boundaries(contributions: Sequence[Contribution], target: int) -> list[float]:
    """``determineQuantileBoundaries``: the sorting key of the row each mass target falls in.

    The targets are equal **masses**, ``total * (i + 1) / target``, and the walk places a boundary
    at the row whose running mass first reaches one. A target nothing reaches is padded with an
    infinity, which :func:`distribute` then never crosses.
    """
    if target <= 1 or not contributions:
        return []
    ordered = sorted(contributions, key=Contribution.sorting_key)
    total_mass = sum(contribution.mass for contribution in ordered)
    if not total_mass > MASS_TOLERANCE:
        return []

    targets = [total_mass * (index + 1) / target for index in range(target - 1)]
    boundaries: list[float] = []
    cumulative = 0.0
    target_index = 0
    for contribution in ordered:
        next_cumulative = cumulative + contribution.mass
        while (
            target_index < len(targets)
            and next_cumulative >= targets[target_index] - MASS_TOLERANCE
        ):
            boundaries.append(contribution.sorting_key())
            target_index += 1
        cumulative = next_cumulative
    while len(boundaries) < len(targets):
        boundaries.append(math.inf)
    return boundaries


@dataclass(frozen=True, slots=True)
class Profile:
    """``PseudoComponentProfile``: one cut of a fluid, or of the blend."""

    moles: float
    molar_mass: float
    density: float
    properties: Properties

    @classmethod
    def empty(cls) -> Profile:
        """The class's empty profile: no moles, and every property absent."""
        return cls(
            moles=0.0,
            molar_mass=0.0,
            density=math.nan,
            properties=Properties(**dict.fromkeys(PROPERTY_FIELDS, math.nan)),
        )

    def mass(self) -> float:
        return self.moles * self.molar_mass

    def is_valid(self) -> bool:
        """``moles > 0 && molarMass > 0 && density > 0``: what makes a cut a component."""
        return self.moles > 0.0 and self.molar_mass > 0.0 and self.density > 0.0


class _Weighted:
    """A weighted mean that ignores a non-finite sample, as the class's accumulators do."""

    __slots__ = ("any", "sum", "weight")

    def __init__(self) -> None:
        self.sum = 0.0
        self.weight = 0.0
        self.any = False

    def add(self, weight: float, value: float) -> None:
        if math.isfinite(value):
            self.sum += weight * value
            self.weight += weight
            self.any = True

    def mean(self, mass: float) -> float:
        return self.sum / mass if self.any else math.nan

    def weighted_mean(self) -> float:
        """``sum(w X) / sum(w)``, or a NaN below the tolerance - ``WeightedMean.mean()``."""
        return self.sum / self.weight if self.weight > MASS_TOLERANCE else math.nan


class GroupBuilder:
    """``PseudoComponentGroupBuilder``: the running totals one group is built from."""

    __slots__ = ("density_mass", "mass", "moles", "sums", "volume")

    def __init__(self) -> None:
        self.mass = 0.0
        self.moles = 0.0
        self.volume = 0.0
        self.density_mass = 0.0
        self.sums = {name: _Weighted() for name in PROPERTY_FIELDS}

    def _add_properties(self, mass_portion: float, properties: Properties) -> None:
        for name in PROPERTY_FIELDS:
            self.sums[name].add(mass_portion, getattr(properties, name))

    def add_contribution(self, source: Contribution, mass_portion: float) -> None:
        """Add a row's mass and its mass-weighted properties."""
        if not mass_portion > MASS_TOLERANCE or not source.mass > 0.0:
            return
        ratio = mass_portion / source.mass
        self.mass += mass_portion
        self.moles += source.moles * ratio
        if source.density > 0.0:
            self.volume += mass_portion / source.density
            self.density_mass += mass_portion * source.density
        self._add_properties(mass_portion, source.properties)

    def add_profile(self, profile: Profile) -> None:
        """Add a whole profile, which is what blending the fluids does."""
        if not profile.is_valid():
            return
        mass_portion = profile.mass()
        self.mass += mass_portion
        self.moles += profile.moles
        if profile.density > 0.0:
            self.volume += mass_portion / profile.density
            self.density_mass += mass_portion * profile.density
        self._add_properties(mass_portion, profile.properties)

    def build(
        self,
        boiling_point: float | None = None,
        critical_temperature: float | None = None,
        critical_pressure: float | None = None,
        acentric_factor: float | None = None,
    ) -> Profile:
        """``buildWithOverrides``: the profile, or an empty one where the mass is at nothing."""
        if not self.mass > MASS_TOLERANCE or not self.moles > MASS_TOLERANCE:
            return Profile.empty()
        if self.volume > MASS_TOLERANCE:
            density = self.mass / self.volume
        elif self.density_mass > MASS_TOLERANCE:
            density = self.density_mass / self.mass
        else:
            density = math.nan
        overrides = {
            "normal_boiling_point": boiling_point,
            "critical_temperature": critical_temperature,
            "critical_pressure": critical_pressure,
            "acentric_factor": acentric_factor,
        }
        values = {}
        for name in PROPERTY_FIELDS:
            override = overrides.get(name)
            values[name] = override if override is not None else self.sums[name].mean(self.mass)
        return Profile(
            moles=self.moles,
            molar_mass=self.mass / self.moles,
            density=density,
            properties=Properties(**values),
        )


def distribute(
    contributions: Sequence[Contribution], boundaries: Sequence[float], target: int
) -> list[Profile]:
    """``distributeToProfiles``: each row into the group its key falls in.

    The walk is forward and the comparison is **strict**, so a key exactly on a boundary stays in
    the lower group. ``contributions`` must already be in ascending sorting key.
    """
    if target <= 0:
        return []
    if not contributions:
        return [Profile.empty() for _ in range(target)]
    builders = [GroupBuilder() for _ in range(target)]
    group = 0
    current = boundaries[0] if boundaries else math.inf
    for contribution in contributions:
        key = contribution.sorting_key()
        while group < target - 1 and key > current:
            group += 1
            current = boundaries[group] if group < len(boundaries) else math.inf
        builders[group].add_contribution(contribution, contribution.mass)
    return [builder.build() for builder in builders]


def combine_profiles(
    per_fluid_profiles: Sequence[Sequence[Profile]],
    fluid_mass_totals: Sequence[float],
    fluid_mole_totals: Sequence[float],
) -> list[Profile]:
    """``combineProfiles``: the per-fluid cuts blended into one slate.

    A fluid's weight is its share of the pooled pseudo-component **mass**, and within it a cut's
    share is its mole fraction of that fluid. The boiling point and the three critical properties
    are weighted by the product of the two; the remaining thirteen are mass-weighted by the
    builder.
    """
    if not per_fluid_profiles:
        return []
    total_fluid_mass = sum(fluid_mass_totals)
    combined: list[Profile] = []
    for group in range(len(per_fluid_profiles[0])):
        builder = GroupBuilder()
        numerators = {
            "normal_boiling_point": 0.0,
            "critical_temperature": 0.0,
            "critical_pressure": 0.0,
            "acentric_factor": 0.0,
        }
        denominators = dict.fromkeys(numerators, 0.0)
        for fluid, profiles in enumerate(per_fluid_profiles):
            profile = profiles[group]
            builder.add_profile(profile)
            fluid_fraction = (
                fluid_mass_totals[fluid] / total_fluid_mass
                if total_fluid_mass > MASS_TOLERANCE
                else 0.0
            )
            mole_fraction = (
                profile.moles / fluid_mole_totals[fluid]
                if fluid_mole_totals[fluid] > MASS_TOLERANCE
                else 0.0
            )
            weight = fluid_fraction * mole_fraction
            for name in numerators:
                value = getattr(profile.properties, name)
                if math.isfinite(value):
                    numerators[name] += weight * value
                    denominators[name] += weight
        overrides = {
            name: numerator / denominators[name] if denominators[name] > MASS_TOLERANCE else None
            for name, numerator in numerators.items()
        }
        combined.append(
            builder.build(
                boiling_point=overrides["normal_boiling_point"],
                critical_temperature=overrides["critical_temperature"],
                critical_pressure=overrides["critical_pressure"],
                acentric_factor=overrides["acentric_factor"],
            )
        )
    return combined


def compute_common_slate(
    per_fluid_profiles: Sequence[Sequence[Profile]],
    fluid_mole_totals: Sequence[float],
    weights: Sequence[float],
    target: int,
) -> list[Profile]:
    """``computeCommonSlate``: the shared slate, by Pedersen's Chapter 5.6 weighting.

    Each fluid's cut contributes weighted by the fluid's slate weight times the cut's mole
    fraction **of that fluid's total** - base components included, which is what makes a dilute
    fluid's pseudo-components count for less than its weight alone would say. The molar mass is
    the weighted mean of the cuts' own; the density is that mean over the weighted mean of their
    molar volumes, so a cut that states no positive density contributes to neither.
    """
    slate: list[Profile] = []
    for index in range(target):
        molar_mass = _Weighted()
        molar_volume = _Weighted()
        sums = {name: _Weighted() for name in PROPERTY_FIELDS}
        slate_moles = 0.0
        for fluid, profiles in enumerate(per_fluid_profiles):
            profile = profiles[index]
            if not profile.is_valid():
                continue
            total = fluid_mole_totals[fluid]
            mole_fraction = profile.moles / total if total > MASS_TOLERANCE else 0.0
            weight = weights[fluid] * mole_fraction
            if not weight > 0.0:
                continue
            slate_moles += profile.moles
            molar_mass.add(weight, profile.molar_mass)
            if profile.density > 0.0:
                molar_volume.add(weight, profile.molar_mass / profile.density)
            for name in PROPERTY_FIELDS:
                sums[name].add(weight, getattr(profile.properties, name))
        mean_molar_mass = molar_mass.weighted_mean()
        mean_molar_volume = molar_volume.weighted_mean()
        density = (
            mean_molar_mass / mean_molar_volume
            if math.isfinite(mean_molar_mass)
            and math.isfinite(mean_molar_volume)
            and mean_molar_volume > 0.0
            else math.nan
        )
        slate.append(
            Profile(
                moles=slate_moles,
                molar_mass=mean_molar_mass if math.isfinite(mean_molar_mass) else 0.0,
                density=density,
                properties=_properties_from_values(
                    [sums[name].weighted_mean() for name in PROPERTY_FIELDS]
                ),
            )
        )
    return slate


def implied_acentric_factor(m: float) -> float:
    """``AttractiveTermSrk.setm``'s side effect: the acentric factor the fitted ``m`` implies.

    **The setter rewrites the component's acentric factor**, solving
    ``-0.176 w**2 + 1.574 w + (0.48 - m) = 0`` with the class's own Newton-Halley iteration from a
    0.25 seed. ``applyTo`` assigns the averaged acentric factor first and ``m`` last, so where the
    row carries an ``m`` the answer's omega is this root and not the average.
    """
    eps = 1.0e-10

    def value(x: float) -> float:
        return -0.176 * x * x + 1.574 * x + (0.48 - m)

    def derivative(x: float) -> float:
        return -0.352 * x + 1.574

    def second(_x: float) -> float:
        return -0.352

    x = 0.25
    next_x = x
    alternate = x
    iterations = 0
    while True:
        iterations += 1
        x = next_x if abs(value(next_x)) < abs(value(alternate)) else alternate
        slope, curvature = derivative(x), second(x)
        radical = slope * slope - 2.0 * value(x) * curvature
        if radical > 0.0:
            next_x = x - slope / curvature + math.sqrt(radical) / curvature
            alternate = x - slope / curvature - math.sqrt(radical) / curvature
        else:
            next_x = x - value(x) / slope
            alternate = next_x
        if next_x < -0.99:
            next_x = -0.99
        if next_x > 1.5:
            next_x = 1.0
        if not abs(value(x)) > eps or iterations > 500:
            break
    return next_x


def _value(name: str, value: float) -> Q | float | None:
    """One property as the row states it: a quantity in its declared unit, or absent."""
    if not math.isfinite(value):
        return None
    unit = PROPERTY_UNITS.get(name)
    return value if unit is None else from_si(value, unit)


def profile_row(name: str, profile: Profile) -> CharacterizationComponent:
    """A profile as a result row, a finite property present and an absent one left out."""
    values: dict[str, Any] = {}
    for field in PROPERTY_FIELDS:
        value = getattr(profile.properties, field)
        if field == "cp_d":
            # The rebuilt component states what its getter states.
            value = value + CP_D_GETTER_OFFSET
        elif field == "acentric_factor" and math.isfinite(profile.properties.attractive_m):
            # `setm` runs after `setAcentricFactor` in `applyTo`.
            value = implied_acentric_factor(profile.properties.attractive_m)
        values[field] = _value(field, value)
    return CharacterizationComponent(
        name=name,
        pseudo_component=True,
        moles=from_si(profile.moles, "mol"),
        molar_mass=cast("Q | None", _value("molar_mass", profile.molar_mass)),
        density=cast("Q | None", _value("density", profile.density)),
        **values,
    )


def base_row(name: str, moles: float) -> CharacterizationComponent:
    """A base component as a result row: its name and amount, every property absent."""
    return CharacterizationComponent(
        name=name,
        pseudo_component=False,
        moles=from_si(moles, "mol"),
        molar_mass=None,
        density=None,
        **dict.fromkeys(PROPERTY_FIELDS, None),
    )
