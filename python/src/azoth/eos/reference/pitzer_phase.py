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

from azoth import _models_gen
from azoth.core.errors import InvalidInputError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import PitzerPhaseResult
from azoth.core.units import Q, input_to_si
from azoth.core.warnings import Warning
from azoth.eos import components as _components
from azoth.eos.reference import _electrolyte as electrolyte
from azoth.eos.reference import _pitzer_catalog as catalog
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

MODEL_ID = "eos.pitzer_phase"


def pitzer_phase(
    components: Sequence[str],
    T: Q,
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
        The activity coefficients, the molalities, the ionic strength, the osmotic
        coefficient, the water activity and which dataset answered.

    Raises:
        InvalidInputError: if a component is not in the databank, if it carries no molar
            mass, if ``x`` is not a composition, if the mixture carries no water, or if
            the dataset the selection rule chose does not cover the brine's topology.
        OutOfRangeError: if ``T`` is not positive.

    Example:
        >>> import azoth
        >>> q = azoth.ureg.Quantity
        >>> r = pitzer_phase(["water", "na+", "cl-"], q(298.15, "K"), [0.88, 0.06, 0.06])
        >>> r.dataset
        'phreeqc'
        >>> round(r.osmotic_coefficient, 12)
        1.099026940549
    """
    spec = _models_gen.model(MODEL_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []

    t_si = input_to_si(spec, "T", T)
    apply_checks(checks.on_input, {"T": t_si}.get, warnings)

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

    return PitzerPhaseResult(
        gamma=tuple(math.exp(value) for value in logged),
        ln_gamma=tuple(logged),
        molality=tuple(composition.molality),
        ionic_strength=composition.ionic_strength,
        osmotic_coefficient=osmotic,
        water_activity=math.exp(logged[solvent]) * x[solvent],
        dataset=selection[0],
        warnings=tuple(warnings),
    )
