"""``eos.pitzer_phase`` - the activity coefficients of a Pitzer electrolyte phase.

Spec: ``specs/models/eos/pitzer_phase.toml``. A *direct* model: no iteration, so no
algorithm block.

This is the pure-Python reference: a second, independent expression of the same physics
as the Rust kernel, held to it by ``test_cross_impl.py``. The arithmetic is
``_pitzer_phase``'s and the dataset rule is ``_pitzer_catalog``'s; what is here is the
composition the three branches are taken over, and **which branch each component takes**.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from typing import Any, NamedTuple

from azoth import _models_gen
from azoth.core.errors import InvalidInputError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import HenryStatus, PitzerPhaseResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning
from azoth.eos import components as _components
from azoth.eos.reference import _electrolyte as electrolyte
from azoth.eos.reference import _pitzer_catalog as catalog
from azoth.eos.reference import iapws_henry_law as _iapws
from azoth.eos.reference._henry import (
    INSOLUBLE_HENRY_COEFFICIENT,
    coefficient,
    effective_coefficient,
    is_capped,
)
from azoth.eos.reference._pitzer_catalog import find
from azoth.eos.reference._pitzer_phase import (
    ION_CHARGE,
    DatasetParameters,
    catalogue_interactions,
    debye_huckel_a_phi,
    ln_gamma,
    ln_gamma_neutral,
    ln_gamma_water,
    osmotic_coefficient,
    osmotic_neutral,
)
from azoth.eos.reference.antoine_vapor_pressure import saturation_pressure

MODEL_ID = "eos.pitzer_phase"

#: Pascals per bar. Every reference pressure NeqSim divides by is in bar.
_BAR_TO_PA = 1.0e5


class _Activity(NamedTuple):
    """The activity surface, as :func:`_activity_of` leaves it."""

    entries: list[Any]
    charge: list[float]
    gamma: list[float]
    ln_gamma: list[float]
    molality: list[float]
    ionic_strength: float
    osmotic: float
    solvent: int
    neutral_active: bool
    dataset: str


def _activity_of(components: Sequence[str], t_si: float, x: Sequence[float]) -> _Activity:
    """The activity-coefficient surface alone, which is what the reference phase is.

    Separate from :func:`pitzer_phase` because the reference phase must not compute a
    fugacity coefficient: NeqSim evaluates ``getActivityCoefficientInfDilWater`` through
    the excess-Gibbs energy, which is this surface, and a twin that took the whole model
    would recurse into itself without end.
    """
    return _activity_inner(components, t_si, x)


def _activity_inner(components: Sequence[str], t_si: float, x: Sequence[float]) -> _Activity:
    """The activity-coefficient arithmetic, with the input checks left to the caller."""
    n = len(components)

    n = len(components)
    if len(x) != n:
        raise InvalidInputError("x", f"a mixture of {n} components has {len(x)} mole fractions")
    bad = next((i for i, value in enumerate(x) if value < 0.0), None)
    if bad is not None:
        raise InvalidInputError("x", f"x[{bad}] is {x[bad]} but a mole fraction cannot be negative")
    total = sum(x)
    if abs(total - 1.0) > 1.0e-9:
        raise InvalidInputError(
            "x",
            f"the mole fractions sum to {total}, not to one. Renormalising them here "
            "would make a composition error invisible in every number downstream, so it "
            "is refused instead",
        )

    entries = [_components.entry(name) for name in components]
    for name, record in zip(components, entries, strict=True):
        if record.molar_mass is None:
            raise InvalidInputError(
                "components",
                f"`{name}` carries no molar mass, and a molality is a mole count over a "
                "solvent mass",
            )

    charge = [record.ionic_charge for record in entries]
    molar_mass = [
        record.molar_mass.to_base_units().magnitude  # type: ignore[union-attr]
        for record in entries
    ]
    names = [record.name for record in entries]
    composition = electrolyte.composition(names, x, molar_mass, charge)

    species = [
        catalog.Species(
            name=record.name,
            moles=moles,
            charge=record.ionic_charge,
            formula=record.formula or "",
            hydrocarbon=record.component_type == "hc",
        )
        for record, moles in zip(entries, x, strict=True)
    ]
    selection = catalog.select_dataset(species)
    audit = catalog.coverage(species, selection, composition.solvent_mass)
    catalog.require_complete(audit)
    catalogue = selection[0] == "phreeqc"
    parameters = DatasetParameters(names, catalogue)

    ions = [i for i in range(n) if abs(charge[i]) >= catalog.ACTIVE_CHARGE]
    neutrals = [
        i
        for i in range(n)
        if abs(charge[i]) < catalog.ACTIVE_CHARGE
        and names[i].lower() != "water"
        and not catalog.is_hydrocarbon(species[i])
    ]
    interactions = catalogue_interactions(names, charge, ions, neutrals) if catalogue else None
    neutral_active = interactions is not None and len(interactions) > 0
    neutral_osmotic = (
        osmotic_neutral(interactions, composition.molality, t_si) if interactions else 0.0
    )

    # `isNonTwoTwoBeta2Active` is set by `setBeta2`, which the catalogue loader calls for
    # every row it applies - so a non-2:2 pair with a `B2` is enough. The legacy loader
    # writes the array without the setter, so it never fires there.
    non_two_two_beta2 = catalogue and any(
        coefficient is not None and abs(coefficient[0]) > 1.0e-20
        for first in range(n)
        for second in range(first + 1, n)
        if not (abs(charge[first]) >= 1.5 and abs(charge[second]) >= 1.5)
        for coefficient in [find("B2", [names[first], names[second]])]
    )
    unequal_charge_same_sign = any(
        charge[first] * charge[second] > 0.0 and abs(charge[first] - charge[second]) >= 1.0e-12
        for first in range(n)
        for second in range(first + 1, n)
    )

    solvent = next((i for i, name in enumerate(names) if name.lower() == "water"), None)
    if solvent is None:
        raise InvalidInputError(
            "components",
            "the mixture names no water, so there is no solvent to measure a molality against",
        )

    # Two keyword sets rather than one: the ion branch takes the two flags the ion terms
    # are gated on, and only the water node reads `neutral_interactions_active`.
    ion_keyword = {
        "non_two_two_beta2": non_two_two_beta2,
        "unequal_charge_same_sign": unequal_charge_same_sign,
    }
    water_keyword = {**ion_keyword, "neutral_interactions_active": neutral_active}

    a_phi = debye_huckel_a_phi(t_si)
    osmotic = osmotic_coefficient(
        composition.molality,
        charge,
        t_si,
        a_phi,
        parameters,
        solvent,
        neutral_osmotic,
        **water_keyword,
    )

    # Named `logged` rather than `ln_gamma`: the module function is what fills it.
    logged: list[float] = []
    for i in range(n):
        if i == solvent and entries[i].reference_state == _components.SOLVENT:
            logged.append(
                ln_gamma_water(
                    composition.molality,
                    charge,
                    t_si,
                    a_phi,
                    parameters,
                    solvent,
                    x[solvent],
                    neutral_osmotic,
                    **water_keyword,
                )
            )
        elif abs(charge[i]) < ION_CHARGE:
            logged.append(
                ln_gamma_neutral(interactions, composition.molality, i, t_si)
                if interactions
                else 0.0
            )
        else:
            logged.append(
                ln_gamma(
                    composition.molality,
                    charge,
                    t_si,
                    a_phi,
                    parameters,
                    i,
                    common_ion_terms=catalogue,
                    **ion_keyword,
                )
            )

    return _Activity(
        entries=entries,
        charge=charge,
        gamma=[math.exp(value) for value in logged],
        ln_gamma=logged,
        molality=list(composition.molality),
        ionic_strength=composition.ionic_strength,
        osmotic=osmotic,
        solvent=solvent,
        neutral_active=neutral_active,
        dataset=selection[0],
    )


def pitzer_phase(
    components: Sequence[str],
    T: Q,
    P: Q,
    x: Sequence[float],
) -> PitzerPhaseResult:
    """The activity coefficients of a brine whose non-ideality is Pitzer's.

    ``components`` are names rather than a resolved record, because the two parameter
    datasets are keyed by them: each name is resolved against the databank here for its
    ionic charge, molar mass and ``REFERENCESTATETYPE``. ``T`` is the only state variable
    - a molality is ``n_i / m_water`` from the mole fractions alone, so the phase's size
    cancels - and there is no pressure argument, because none of the branches NeqSim's
    ``getGamma`` dispatches to reads the one it is given.

    The branch a component takes is its charge, except for water: water with the
    ``solvent`` reference state takes the osmotic route, any other neutral takes the
    sparse neutral layer, and an ion takes the extended Pitzer expression.

    Args:
        components: the substances the brine is made of, by name.
        T: absolute temperature.
        x: the liquid's mole fractions; non-negative and summing to one, with water
            among them.

    Returns:
        The activity coefficients, the fugacity coefficients, the molalities, the ionic
        strength, the osmotic coefficient, the water activity and which dataset answered.

    Raises:
        InvalidInputError: if a component is not in the databank, if it carries no molar
            mass, if ``x`` is not a composition, if the mixture carries no water, or if
            the dataset the selection rule chose does not cover the brine's topology.
        OutOfRangeError: if ``T`` is not positive.

    Example:
        >>> import azoth
        >>> q = azoth.ureg.Quantity
        >>> r = pitzer_phase(
        ...     ["water", "na+", "cl-"], q(298.15, "K"), q(1.0e5, "Pa"), [0.88, 0.06, 0.06]
        ... )
        >>> r.dataset
        'phreeqc'
        >>> round(r.osmotic_coefficient, 12)
        1.099026940549
    """
    spec = _models_gen.model(MODEL_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    t_si = input_to_si(spec, "T", T)
    p_si = input_to_si(spec, "P", P)
    apply_checks(checks.on_input, {"T": t_si, "P": p_si}.get, warnings)

    activity = _activity_of(components, t_si, x)
    bar = p_si / _BAR_TO_PA
    n = len(components)

    ln_phi: list[float] = []
    henry: list[Q] = []
    gamma_inf: list[float] = []

    for i in range(n):
        record = activity.entries[i]
        is_water = components[i].strip().lower() == "water"
        if abs(activity.charge[i]) < ION_CHARGE and not is_water:
            # `ComponentGePitzer.fugcoef`. The ratio is one where it cannot be taken,
            # which is NeqSim's own fallback for a zero mole fraction.
            ratio = activity.molality[i] / x[i] if x[i] > 0.0 else 0.0
            ratio = ratio if ratio > 0.0 and math.isfinite(ratio) else 1.0
            h = _neutral_henry(
                record,
                components[i],
                t_si,
                activity.neutral_active,
                _has_ions(activity.charge),
            )
            henry.append(from_si(h, "Pa"))
            gamma_inf.append(1.0)
            ln_phi.append(math.log(activity.gamma[i] * h * ratio / bar))
        elif record.reference_state == _components.SOLVENT and not _uses_iapws(components[i]):
            p0 = saturation_pressure(record, t_si, warnings)
            henry.append(from_si(0.0, "Pa"))
            gamma_inf.append(1.0)
            ln_phi.append(math.log(activity.gamma[i] * p0 / _BAR_TO_PA / bar))
        else:
            reference = _reference_phase_gamma(components, i, activity.solvent, t_si)
            raw = coefficient(record.henry, t_si)
            h = INSOLUBLE_HENRY_COEFFICIENT if is_capped(record, raw) else raw
            henry.append(from_si(h, "Pa"))
            gamma_inf.append(reference)
            ln_phi.append(math.log(activity.gamma[i] / reference * h / bar))

    return PitzerPhaseResult(
        gamma=tuple(activity.gamma),
        ln_gamma=tuple(activity.ln_gamma),
        ln_phi=tuple(ln_phi),
        henry=tuple(henry),
        gamma_inf=tuple(gamma_inf),
        molality=tuple(activity.molality),
        ionic_strength=activity.ionic_strength,
        osmotic_coefficient=activity.osmotic,
        water_activity=activity.gamma[activity.solvent] * x[activity.solvent],
        dataset=activity.dataset,
        warnings=tuple(warnings),
    )


def _neutral_henry(
    entry: Any,
    name: str,
    temperature_k: float,
    neutral_active: bool,
    has_ions: bool,
) -> float:
    """``ComponentGePitzer.getEffectiveHenryCoefficient(phase)``, in bar.

    The override selects the IAPWS pure-water table behind **three gates**: no ions in
    the phase, no active neutral interaction family, and a solute that is not CO2 or
    H2S. Where a gate closes, the database correlation answers instead.
    """
    gas = _iapws.gas_from_name(name)
    if gas is None:
        return effective_coefficient(entry, temperature_k)
    reactive = gas in ("co2", "h2s")
    if not neutral_active and (has_ions or reactive):
        return effective_coefficient(entry, temperature_k)
    table = _iapws.iapws_henry_law(name, _kelvin(temperature_k))
    if table.status is HenryStatus.GUIDELINE_EXTRAPOLATION:
        # `isUsable` fails outside the row's fitted window, and the guideline's own
        # consumer fails closed to the insoluble limit rather than extrapolating.
        return INSOLUBLE_HENRY_COEFFICIENT
    # The table is on the mole-fraction scale and the activity on the molality scale, so
    # the constant is converted by water's molar mass.
    value = table.henry.to_base_units().magnitude / _BAR_TO_PA * _iapws.WATER_MOLAR_MASS_KG_PER_MOL
    return INSOLUBLE_HENRY_COEFFICIENT if is_capped(entry, value) else value


def _kelvin(value: float) -> Q:
    """A bare kelvin quantity, for the table's own argument."""
    import azoth

    return azoth.ureg.Quantity(value, "K")


def _uses_iapws(name: str) -> bool:
    """``ComponentGE.usesIapwsAqueousReference``: the table overrides a solvent
    classification for a non-water component the table carries a row for.

    Vacuous on this databank - see the spec's assumptions - and kept because it is
    NeqSim's own predicate.
    """
    return name.strip().lower() != "water" and _iapws.gas_from_name(name) is not None


def _reference_phase_gamma(
    components: Sequence[str], solute: int, solvent: int, temperature_k: float
) -> float:
    """``PhaseGE.getActivityCoefficientInfDilWater``: the solute's gamma in a
    **two-component reference phase**.

    The reference phase is the solute at ``1e-10`` mol beside the solvent at ``10.0``
    mol (``Phase.initRefPhases``), so its composition is ``[1e-11, 1 - 1e-11]`` - and
    the arithmetic is this model's **activity** surface, which is why it is
    ``_activity_of`` and not ``pitzer_phase``: NeqSim evaluates the reference phase
    through the excess Gibbs energy, and a twin that took the whole model would recurse
    into itself without end.
    """
    dilute = [1.0e-11, 1.0 - 1.0e-11]
    reference = _activity_of([components[solute], components[solvent]], temperature_k, dilute)
    return reference.gamma[0]


def _has_ions(charge: Sequence[float]) -> bool:
    """Whether the phase carries any ion at all, which is one of the three gates."""
    return any(abs(value) >= ION_CHARGE for value in charge)
