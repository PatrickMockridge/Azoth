"""``characterization.oil_assay_derivation`` - an assay's cuts derived to molar masses and moles.

Spec: ``specs/models/characterization/oil_assay_derivation.toml``. Oracle:
``validation/neqsim/captures/oil_assay_derivation_probe.tsv``.

Mirrors ``crates/azoth-characterization/src/oil_assay_derivation.rs`` line for line.

# The arithmetic, and not the mutation

``OilAssayCharacterisation.apply`` resolves the cuts' mass fractions and then derives one component
per cut with a positive yield. In NeqSim that component goes into a ``SystemInterface``; here the
answer is the table - each component's name, molar mass and moles - and the fluid is the caller's
step.

# Three routes to one molar mass, in this order

A **standard component** is a substance the databank knows, and its molar mass comes out of it by
name. A **stated molar mass** is used as it is. Otherwise the closed form
``5.805e-5 * Tb^2.3776 / SG^0.9371`` is evaluated, in g/mol and g/cm3 - the class's own scale - on a
boiling point the cut states or one resolved from its Watson factor, and a gravity it states or one
resolved from its API gravity.

**A cut that states a molar mass still needs a gravity**, and not for the correlation: the class
reaches ``resolveDensity`` for every cut it derives, before it asks whether the cut has a molar mass
at all.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any

from azoth import _models_gen
from azoth.core.errors import InvalidInputError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import OilAssayDerivationResult
from azoth.core.units import Q, from_si, input_records_to_si, input_to_si, quantity
from azoth.core.warnings import Warning

MODEL_ID = "characterization.oil_assay_derivation"

#: ``apply``'s own gate, on the mass it rebuilt against the mass it started from.
MASS_BALANCE_TOLERANCE = 1.0e-10

#: The divisor's floor, the class's ``Math.max(totalAssayMass, 1.0e-30)``.
MASS_BALANCE_FLOOR = 1.0e-30

#: ``AssayCut.resolveDensity``'s API relation.
API_NUMERATOR = 141.5
API_OFFSET = 131.5

#: ``AssayCut.resolveMolarMass``'s coefficient and exponents.
MOLAR_MASS_COEFFICIENT = 5.805e-5
BOILING_POINT_EXPONENT = 2.3776
DENSITY_EXPONENT = 0.9371

#: ``AssayCut.resolveAverageBoilingPoint``'s Watson divisor.
WATSON_DIVISOR = 1.8


def specific_gravity(row: Mapping[str, Any]) -> float | None:
    """``AssayCut.resolveDensity``: the stated gravity, or the one its API gravity implies.

    g/cm3, which is the scale every correlation in the class is in. ``None`` where the cut states
    neither, which the derivation refuses - for every cut it derives and not only for a volume
    basis.
    """
    density = row["density"]
    if density is not None:
        return float(density) / 1000.0
    api_gravity = row["api_gravity"]
    if api_gravity is None:
        return None
    return API_NUMERATOR / (float(api_gravity) + API_OFFSET)


def average_boiling_point(row: Mapping[str, Any], gravity: float) -> float | None:
    """``AssayCut.resolveAverageBoilingPoint``: the stated boiling point, or the Watson one."""
    stated = row["average_boiling_point"]
    if stated is not None:
        return float(stated)
    watson = row["watson_characterization_factor"]
    if watson is None:
        return None
    return (float(watson) * gravity) ** 3 / WATSON_DIVISOR


def correlated_molar_mass(boiling_point: float, gravity: float) -> float:
    """``AssayCut.resolveMolarMass``'s closed form, in kg/mol from K and g/cm3."""
    return float(
        MOLAR_MASS_COEFFICIENT
        * boiling_point**BOILING_POINT_EXPONENT
        / gravity**DENSITY_EXPONENT
        / 1000.0
    )


def oil_assay_derivation(
    basis: str,
    cuts: Sequence[Mapping[str, Any]],
    total_assay_mass: Q,
) -> OilAssayDerivationResult:
    """An assay's cuts derived to one component each, by ``OilAssayCharacterisation.apply``.

    Raises:
        InvalidInputError: if two cuts share a name, if two standard components name one substance,
            if a standard component is on a volume basis or also states pseudo-component
            properties, if a cut reaches the correlation with neither a boiling point nor a gravity,
            if a standard component is not in the databank, or if the mass it rebuilt misses the
            total by more than ``1e-10`` relative.
    """
    # The databank lookup is deferred because `azoth.eos` and this namespace are both imported by
    # `azoth.characterization`, and one of the two has to reach for the other at call time.
    from azoth.characterization.reference.assay_mass_fractions import (
        assay_mass_fractions as resolve_mass_fractions,
    )
    from azoth.eos import components as databank

    spec = _models_gen.model(MODEL_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []
    apply_checks(checks.on_input, lambda _name: None, warnings)

    rows: list[dict[str, Any]] = input_records_to_si(spec, "cuts", cuts)
    total = input_to_si(spec, "total_assay_mass", total_assay_mass)

    # `validateUniqueCutNames`: by name, case-sensitively.
    names = [str(row["name"]) for row in rows]
    if len(set(names)) != len(names):
        raise InvalidInputError(
            "cuts", "two cuts carry one name, and a component cannot be named twice"
        )
    # `validateUniqueStandardComponentNames`: by lower-cased name.
    substances = [
        str(row["standard_component_name"]).lower()
        for row in rows
        if row["standard_component_name"] is not None
    ]
    if len(set(substances)) != len(substances):
        raise InvalidInputError(
            "cuts", "two cuts name the same standard component, which is one substance"
        )

    gravities = [specific_gravity(row) for row in rows]
    # The assay kernel's own resolution is over SI densities; a cut that states no gravity passes
    # none, which its mass branch never reads and its volume branch refuses.
    assay = resolve_mass_fractions(
        basis,
        [
            {
                "declared_fraction": float(row["declared_fraction"]),
                "density": None if gravity is None else quantity(gravity * 1000.0, "kg/m**3"),
            }
            for row, gravity in zip(rows, gravities, strict=True)
        ],
    )
    mass_fractions = assay.mass_fraction

    component_name: list[str] = []
    derived_fraction: list[float] = []
    molar_mass: list[float] = []
    moles: list[float] = []
    reconstructed = 0.0
    for index, row in enumerate(rows):
        mass_fraction = float(mass_fractions[index])
        if not mass_fraction > 0.0:
            continue
        substance = row["standard_component_name"]
        if substance is not None:
            if basis == "volume":
                raise InvalidInputError(
                    "cuts",
                    f"`{row['name']}` is the standard component `{substance}`, which needs a "
                    f"mass-basis assay fraction",
                )
            if any(
                row[field] is not None
                for field in (
                    "density",
                    "api_gravity",
                    "molar_mass",
                    "average_boiling_point",
                    "watson_characterization_factor",
                )
            ):
                raise InvalidInputError(
                    "cuts",
                    f"`{row['name']}` is the standard component `{substance}` and also states "
                    f"pseudo-component properties, which is one cut wanting two definitions",
                )
            try:
                entry = databank.entry(str(substance))
            except Exception as exc:
                raise InvalidInputError(
                    "cuts", f"`{substance}` is not a substance the databank carries"
                ) from exc
            if entry.molar_mass is None:
                raise InvalidInputError(
                    "cuts", f"the databank carries no molar mass for `{substance}`"
                )
            name = str(substance)
            cut_molar_mass = float(entry.molar_mass.to("kg/mol").magnitude)
        else:
            gravity = gravities[index]
            if gravity is None:
                raise InvalidInputError(
                    f"cuts[{index}].density",
                    f"`{row['name']}` states neither a density nor an API gravity, and the mass "
                    f"the moles are taken over needs one",
                )
            stated = row["molar_mass"]
            if stated is not None:
                cut_molar_mass = float(stated)
            else:
                boiling_point = average_boiling_point(row, gravity)
                if boiling_point is None:
                    raise InvalidInputError(
                        f"cuts[{index}].average_boiling_point",
                        f"`{row['name']}` states neither an average boiling point nor a Watson "
                        f"factor, and one of them is what the correlation is evaluated on",
                    )
                cut_molar_mass = correlated_molar_mass(boiling_point, gravity)
            name = f"{row['name']}_PC"
        cut_moles = total * mass_fraction / cut_molar_mass
        reconstructed += cut_moles * cut_molar_mass
        component_name.append(name)
        derived_fraction.append(mass_fraction)
        molar_mass.append(cut_molar_mass)
        moles.append(cut_moles)

    mass_balance_error = abs(reconstructed - total) / max(abs(total), MASS_BALANCE_FLOOR)
    if mass_balance_error != mass_balance_error or mass_balance_error > MASS_BALANCE_TOLERANCE:
        raise InvalidInputError(
            "total_assay_mass",
            f"the derived components come to {reconstructed} against a total of {total}, which is "
            f"{mass_balance_error} relative and past the class's {MASS_BALANCE_TOLERANCE}",
        )

    apply_checks(checks.derived, lambda _name: None, warnings)

    return OilAssayDerivationResult(
        component_name=tuple(component_name),
        mass_fraction=tuple(derived_fraction),
        molar_mass=tuple(from_si(value, "kg/mol") for value in molar_mass),
        moles=tuple(from_si(value, "mol") for value in moles),
        total_moles=from_si(sum(moles), "mol"),
        mass_balance_error=mass_balance_error,
        warnings=tuple(warnings),
    )
