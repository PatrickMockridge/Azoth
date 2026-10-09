"""Petroleum-fraction characterisation.

The front door for a fluid that has no databank row. A crude is an *assay* - a list of cuts,
each a molar mass and a density - and everything a cubic needs about a cut follows from those
two numbers by correlation. NeqSim's ``thermo/characterization/`` is that subsystem.

It is a namespace of its own rather than a family under :mod:`azoth.eos`, mirroring NeqSim's
own package split: what is here is not an equation of state and is not reached from one.

* :func:`assay_mass_fractions` - an oil assay's cuts resolved to a mass basis
* :func:`characterise_plus_fraction` - the whole chain, from a C7+ end to lumps
* :func:`characterize_to_reference` - a fluid re-cut onto another's cut slate
* :func:`tbp_cut_properties` - a cut's critical properties, by any of NeqSim's ten models
* :func:`tbp_closure` - a cut's molar mass, from its boiling point and gravity
* :func:`tbp_density` - the same pair inverted for gravity, which only one closure supports
* :func:`pedersen_plus_split` - a plus fraction divided into carbon-number cuts
* :func:`whitson_gamma_split` - the same, on Whitson's three-parameter gamma
* :func:`lumping` - a cut table grouped into equal-mass lumps
* :func:`tbp_grouping` - a phase's components grouped by boiling point

# Which implementation answers

As in the other namespaces, each function dispatches to the Rust extension when it is built and
to :mod:`azoth.characterization.reference` otherwise. Both are always reachable - see
:func:`azoth.backends` and :func:`azoth.use_backend`.
"""

from __future__ import annotations

from azoth._dispatch import resolve
from azoth.core.result import (
    AssayMassFractionsResult,
    CharacterisePlusFractionResult,
    CharacterizeToReferenceResult,
    LumpingResult,
    PedersenPlusSplitResult,
    TbpClosureResult,
    TbpCutPropertiesResult,
    TbpDensityResult,
    TbpGroupingResult,
    WhitsonGammaSplitResult,
)
from azoth.core.units import Q

__all__ = [
    "assay_mass_fractions",
    "characterise_plus_fraction",
    "characterize_to_reference",
    "lumping",
    "pedersen_plus_split",
    "tbp_closure",
    "tbp_cut_properties",
    "tbp_density",
    "tbp_grouping",
    "whitson_gamma_split",
]

_TBP_CUT_PROPERTIES = "characterization.tbp_cut_properties"
_TBP_CLOSURE = "characterization.tbp_closure"
_TBP_DENSITY = "characterization.tbp_density"
_TBP_GROUPING = "characterization.tbp_grouping"
_PEDERSEN_PLUS_SPLIT = "characterization.pedersen_plus_split"
_ASSAY_MASS_FRACTIONS = "characterization.assay_mass_fractions"
_CHARACTERISE_PLUS_FRACTION = "characterization.characterise_plus_fraction"
_CHARACTERIZE_TO_REFERENCE = "characterization.characterize_to_reference"
_LUMPING = "characterization.lumping"
_WHITSON_GAMMA_SPLIT = "characterization.whitson_gamma_split"


def assay_mass_fractions(
    basis: str,
    declared_fraction: list[float],
    density: list[Q] | None = None,
) -> AssayMassFractionsResult:
    """An oil assay's declared yields, resolved to a mass basis.

    An assay states each cut's yield as a mass fraction or as a liquid-volume fraction, never both
    and never neither, and the whole row shares one basis. The row must close on one within
    ``1e-3`` - a tenth of a percent is a **hard error** here, not a warning.

    ``density`` is one entry per cut. A volume basis needs it to convert; a mass basis needs it
    only for ``bulk_density``, which is absent without it.

    Raises:
        InvalidInputError: if the densities do not match the cuts, if a volume basis has none, or
            if the declared fractions do not close on one.

    See :func:`azoth.characterization.reference.assay_mass_fractions`.
    """
    return resolve(_ASSAY_MASS_FRACTIONS)(  # type: ignore[no-any-return]
        basis=basis,
        declared_fraction=declared_fraction,
        density=density,
    )


def characterize_to_reference(
    source_moles: list[float],
    source_molar_mass: list[Q],
    source_density: list[Q],
    source_boiling_point: list[Q],
    reference_molar_mass: list[Q],
    reference_boiling_point: list[Q],
) -> CharacterizeToReferenceResult:
    """A fluid's pseudo-components re-cut onto another fluid's cut slate.

    ``reference_*`` decides where the cuts are and nothing else; its densities and amounts are not
    read. The groups are **not** the reference rows - a group that comes out empty is dropped -
    so ``reference_index`` is what re-aligns the answer with the caller's own slate.

    Raises:
        InvalidInputError: if any table is not internally aligned.

    See :func:`azoth.characterization.reference.characterize_to_reference`.
    """
    return resolve(_CHARACTERIZE_TO_REFERENCE)(  # type: ignore[no-any-return]
        source_moles=source_moles,
        source_molar_mass=source_molar_mass,
        source_density=source_density,
        source_boiling_point=source_boiling_point,
        reference_molar_mass=reference_molar_mass,
        reference_boiling_point=reference_boiling_point,
    )


def characterise_plus_fraction(
    molar_mass: Q,
    density: Q,
    mole_fraction: float,
    first_carbon_number: int,
    plus_model: str | None = None,
    number_of_lumps: int | None = None,
) -> CharacterisePlusFractionResult:
    """A C7+ end characterised end to end: model, split and lumps.

    This is the chain the whole namespace exists for - one pseudo-component described by three
    numbers in, a table of lumps out. ``plus_model`` is ``pedersen``, ``pedersen_heavy_oil`` or
    ``whitson_gamma``.

    **``selected_model`` on the result is not always ``plus_model``.** A plus fraction heavier
    than the requested model's maximum molar mass is re-modelled to ``pedersen_heavy_oil``, and
    Whitson Gamma's maximum is ``0.605`` rather than the 2.10 its sibling's name suggests.

    Raises:
        InvalidInputError, OutOfRangeError: from the split or the lumping, including a plus
            fraction the split declines - which is raised here rather than skipped.

    See :func:`azoth.characterization.reference.characterise_plus_fraction`.
    """
    return resolve(_CHARACTERISE_PLUS_FRACTION)(  # type: ignore[no-any-return]
        molar_mass=molar_mass,
        density=density,
        mole_fraction=mole_fraction,
        first_carbon_number=first_carbon_number,
        plus_model=plus_model,
        number_of_lumps=number_of_lumps,
    )


def tbp_closure(
    closure: str,
    boiling_point: Q,
    density: Q,
    model: str | None = None,
) -> TbpClosureResult:
    """A cut's molar mass, from its normal boiling point and its specific gravity.

    ``closure`` is one of ``riazi_daubert_1980``, ``riazi_daubert_1987``, ``soreide`` or
    ``tbp_model``; the last needs ``model``.

    Raises:
        InvalidInputError: if a `tbp_model` closure has no `model`, or the boiling point is not
            attainable over the search bracket.

    See :func:`azoth.characterization.reference.tbp_closure`.
    """
    return resolve(_TBP_CLOSURE)(  # type: ignore[no-any-return]
        closure=closure, boiling_point=boiling_point, density=density, model=model
    )


def tbp_density(
    boiling_point: Q,
    molar_mass: Q,
    closure: str | None = None,
) -> TbpDensityResult:
    """A cut's normal liquid density, from its normal boiling point and its molar mass.

    Only ``riazi_daubert_1980`` supports this direction; the other three closures are refused,
    each naming the class that would close it.

    Raises:
        InvalidInputError: for the three closures that cannot be inverted.

    See :func:`azoth.characterization.reference.tbp_density`.
    """
    return resolve(_TBP_DENSITY)(  # type: ignore[no-any-return]
        boiling_point=boiling_point, molar_mass=molar_mass, closure=closure
    )


def tbp_cut_properties(
    molar_mass: Q,
    density: Q,
    model: str | None = None,
    boiling_point: Q | None = None,
) -> TbpCutPropertiesResult:
    """A TBP cut's critical properties, by any of NeqSim's ten models.

    ``model`` is one of ``pedersen_srk``, ``pedersen_srk_heavy_oil``, ``pedersen_pr``,
    ``pedersen_pr2``, ``pedersen_pr_heavy_oil``, ``riazi_daubert``, ``lee_kesler``, ``twu``,
    ``cavett`` or ``standing``; ``None`` takes NeqSim's own default, ``pedersen_srk``.

    ``attraction_exponent`` is ``None`` for the five models that have no exponent to offer.

    Raises:
        OutOfRangeError: if ``molar_mass`` or ``density`` is not positive, or if a supplied
            ``boiling_point`` is not positive.

    See :func:`azoth.characterization.reference.tbp_cut_properties`.
    """
    return resolve(_TBP_CUT_PROPERTIES)(  # type: ignore[no-any-return]
        model=model, molar_mass=molar_mass, density=density, boiling_point=boiling_point
    )


def tbp_grouping(
    boiling_point: list[Q],
    mole_fraction: list[float],
) -> TbpGroupingResult:
    """A phase's components grouped into boiling-point bins.

    Fourteen thresholds in degrees Celsius, the lowest at 69.2, each opening one of twenty bins -
    so bins 0 to 5 are always zero, and a component below the lowest threshold is in no bin at
    all rather than in the first. That is why the bins do not sum to one.

    ``boiling_point`` is the **stored kelvin**, the number the class subtracts 273.15 from before
    comparing; passing a Celsius value would shift every threshold by 273.

    Raises:
        InvalidInputError: if the two vectors are not the same length.

    See :func:`azoth.characterization.reference.tbp_grouping`.
    """
    return resolve(_TBP_GROUPING)(  # type: ignore[no-any-return]
        boiling_point=boiling_point,
        mole_fraction=mole_fraction,
    )


def pedersen_plus_split(
    molar_mass: Q,
    density: Q,
    mole_fraction: float,
    first_carbon_number: int,
    last_carbon_number: int,
) -> PedersenPlusSplitResult:
    """A plus fraction divided into carbon-number cuts, by Pedersen's two Newton solves.

    ``first_carbon_number`` opens the range and ``last_carbon_number`` is one past its last cut,
    which is NeqSim's own half-open convention: ``80`` is ``PedersenPlusModel``'s default range
    and ``200`` is the heavy-oil model's, and the two are the same two solves either way.

    Raises:
        InvalidInputError: if the plus fraction is lighter than the first cut the table gives it,
            or if the range is empty.
        OutOfRangeError: if a carbon number falls outside the tables the model carries.

    See :func:`azoth.characterization.reference.pedersen_plus_split`.
    """
    return resolve(_PEDERSEN_PLUS_SPLIT)(  # type: ignore[no-any-return]
        molar_mass=molar_mass,
        density=density,
        mole_fraction=mole_fraction,
        first_carbon_number=first_carbon_number,
        last_carbon_number=last_carbon_number,
    )


def whitson_gamma_split(
    molar_mass: Q,
    density: Q,
    mole_fraction: float,
    first_carbon_number: int,
    last_carbon_number: int,
    alpha: float | None = None,
    eta: Q | None = None,
    density_model: str | None = None,
    auto_estimate_shape: bool | None = None,
) -> WhitsonGammaSplitResult:
    """A plus fraction split into cuts by Whitson's three-parameter gamma distribution.

    The cuts are windows in molar mass that step from ``eta`` by 14 g/mol, so
    ``first_carbon_number`` labels the row and ``last_carbon_number`` is the only bound the
    windowing respects. ``density_model`` selects the gravity correlation, ``uop`` by default;
    ``auto_estimate_shape`` replaces ``alpha`` with a correlation on the plus fraction's Watson
    factor.

    Raises:
        InvalidInputError: if the plus fraction's molar mass is not above ``eta``.
        OutOfRangeError: if an input is outside its declared range.

    See :func:`azoth.characterization.reference.whitson_gamma_split`.
    """
    return resolve(_WHITSON_GAMMA_SPLIT)(  # type: ignore[no-any-return]
        molar_mass=molar_mass,
        density=density,
        mole_fraction=mole_fraction,
        first_carbon_number=first_carbon_number,
        last_carbon_number=last_carbon_number,
        alpha=alpha,
        eta=eta,
        density_model=density_model,
        auto_estimate_shape=auto_estimate_shape,
    )


def lumping(
    molar_mass: Q,
    mole_fraction: float,
    cut_z: list[float],
    cut_molar_mass: list[Q],
    cut_density: list[Q],
    number_of_lumps: int | None = None,
) -> LumpingResult:
    """A cut table grouped into equal-mass lumps.

    A split hands back sixty-odd cuts and no cubic wants sixty components, so the cuts are
    grouped into a handful of pseudo-components of equal *mass* - the partition runs on the
    accumulated ``cut_z * cut_molar_mass``, not on the cut count.

    ``molar_mass`` and ``mole_fraction`` are the **plus fraction's** own two numbers rather than
    the table's sums; the class reads them off the fluid and the difference moves a partition
    boundary. ``number_of_lumps`` is the class's ``numberOfPseudocomponents``, 7 by default.

    Raises:
        InvalidInputError: if the three vectors are not one table, or if more lumps are asked for
            than there are cuts.
        OutOfRangeError: if an input is outside its declared range.

    See :func:`azoth.characterization.reference.lumping`.
    """
    return resolve(_LUMPING)(  # type: ignore[no-any-return]
        molar_mass=molar_mass,
        mole_fraction=mole_fraction,
        cut_z=cut_z,
        cut_molar_mass=cut_molar_mass,
        cut_density=cut_density,
        number_of_lumps=number_of_lumps,
    )
