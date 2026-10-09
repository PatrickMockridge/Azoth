"""``characterization.characterise_plus_fraction`` - a C7+ end characterised end to end.

Spec: ``specs/models/characterization/characterise_plus_fraction.toml``. Oracle:
``validation/neqsim/captures/plus_fraction_probe.tsv``.

Mirrors ``crates/azoth-characterization/src/characterise_plus_fraction.rs`` line for line.

# The front door, and nothing of its own

``Characterise.characterisePlusFraction`` is four statements: initialise, find the plus
fraction, **replace the model if the plus fraction is heavier than the model can take**, and run
the lumping only if the split returned true. Everything between is
:func:`azoth.characterization.reference.pedersen_plus_split`,
:func:`azoth.characterization.reference.whitson_gamma_split` and
:func:`azoth.characterization.reference.lumping`.

# The replacement is silent and ``selected_model`` is the only report of it

``getMPlus() > getMaxPlusMolarMass()`` is checked **before** the model is used, so a request for
Whitson Gamma above the threshold is swept away with the rest - and Whitson's own threshold is
``PedersenPlusModel``'s ``0.605``, inherited by its constructor, not the 2.10 its sibling name
suggests.

# A false split is refused rather than skipped

The class takes the boolean and, when it is false, leaves the fluid with its single plus row and
no lumps. An id that returns a table has no fluid to leave unchanged, so the split's own refusal
is propagated. A stated divergence.
"""

from __future__ import annotations

from azoth import _models_gen
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import CharacterisePlusFractionResult
from azoth.core.units import Q, input_to_si
from azoth.core.warnings import Warning
from azoth.characterization.reference.lumping import lumping
from azoth.characterization.reference.pedersen_plus_split import pedersen_plus_split
from azoth.characterization.reference.whitson_gamma_split import whitson_gamma_split

MODEL_ID = "characterization.characterise_plus_fraction"

#: The first carbon number the gamma model's cut table begins at, which its override never moves.
GAMMA_FIRST_CARBON_NUMBER = 1

PLUS_MODELS = ("pedersen", "pedersen_heavy_oil", "whitson_gamma")

#: ``getMaxPlusMolarMass()`` and the constructor's last carbon number, per model.
MAXIMUM_MOLAR_MASS = {
    "pedersen": 0.605,
    "pedersen_heavy_oil": 2.10,
    "whitson_gamma": 0.605,
}
LAST_CARBON_NUMBER = {
    "pedersen": 80,
    "pedersen_heavy_oil": 200,
    "whitson_gamma": 80,
}


def characterise_plus_fraction(
    molar_mass: Q,
    density: Q,
    mole_fraction: float,
    first_carbon_number: int,
    plus_model: str | None = None,
    number_of_lumps: int | None = None,
) -> CharacterisePlusFractionResult:
    """A C7+ end characterised end to end: model, split and lumps.

    ``selected_model`` on the result is **not always** ``plus_model``: a plus fraction heavier
    than the requested model's maximum is re-modelled to ``pedersen_heavy_oil``, silently, because
    the class replaces it before using it.

    Raises:
        InvalidInputError, OutOfRangeError: from whichever of the three ids refuses - including a
            plus fraction the split declines, which this raises rather than skipping.
    """
    spec = _models_gen.model(MODEL_ID)
    checks = checks_for(spec)
    warnings: list[Warning] = []
    apply_checks(checks.on_input, lambda _name: None, warnings)

    requested = plus_model or "pedersen"
    mass = input_to_si(spec, "molar_mass", molar_mass)
    # The class replaces the model before using it, so the request is discarded rather than
    # reported as an error.
    selected = (
        "pedersen_heavy_oil" if mass > MAXIMUM_MOLAR_MASS[requested] else requested
    )
    last = LAST_CARBON_NUMBER[selected]

    if selected == "whitson_gamma":
        # **The first carbon number is ignored here.** The gamma model's override never reads it,
        # so its cuts always begin at the outer object's own value of one.
        split = whitson_gamma_split(
            molar_mass=molar_mass,
            density=density,
            mole_fraction=mole_fraction,
            first_carbon_number=GAMMA_FIRST_CARBON_NUMBER,
            last_carbon_number=last,
        )
    else:
        split = pedersen_plus_split(
            molar_mass=molar_mass,
            density=density,
            mole_fraction=mole_fraction,
            first_carbon_number=first_carbon_number,
            last_carbon_number=last,
        )
    warnings.extend(split.warnings)

    # **The lumping is handed the plus fraction's own two numbers**, not the cut table's sums.
    grouped = lumping(
        molar_mass=molar_mass,
        mole_fraction=mole_fraction,
        cut_z=list(split.cut_z),
        cut_molar_mass=list(split.cut_molar_mass),
        cut_density=list(split.cut_density),
        number_of_lumps=number_of_lumps,
    )
    warnings.extend(grouped.warnings)

    apply_checks(checks.derived, lambda _name: None, warnings)

    return CharacterisePlusFractionResult(
        selected_model=selected,
        fraction_of_heavy_end=grouped.fraction_of_heavy_end,
        lump_molar_mass=grouped.lump_molar_mass,
        lump_density=grouped.lump_density,
        warnings=tuple(warnings),
    )
