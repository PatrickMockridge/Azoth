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

from collections.abc import Mapping, Sequence
from typing import TYPE_CHECKING, Any

from azoth._dispatch import resolve
from azoth.core.errors import InvalidInputError
from azoth.core.result import (
    AssayMassFractionsResult,
    CharacterisePlusFractionResult,
    CharacterizeToReferenceResult,
    LumpingResult,
    OilAssayDerivationResult,
    PedersenPlusSplitResult,
    SetHeavyTbpToPlusResult,
    TbpClosureResult,
    TbpCutBoundariesResult,
    TbpCutPropertiesResult,
    TbpDensityResult,
    TbpGroupingResult,
    WhitsonGammaSplitResult,
)
from azoth.core.units import Q

if TYPE_CHECKING:
    from azoth import keycard


def card_from_cuts(
    names: Sequence[str],
    molar_mass: Sequence[Q],
    density: Sequence[Q],
    *,
    model: str | None = None,
) -> keycard.Keycard:
    """A keycard whose components are these cuts, each described by its own correlation.

    **This is not a calculation.** It is the step between the two things this namespace already
    answers and the one thing a cubic needs: a crude's cuts are a molar mass and a normal liquid
    density each, `tbp_cut_properties` turns one into the constants a cubic reads, and a card is
    how a substance the databank does not carry reaches a `Mixture`. So this composes ids that
    exist rather than adding arithmetic of its own, and every number it writes is one of theirs.

    Each cut's entry carries `Tc`, `Pc` and `omega` - what a cubic reads - beside the two numbers
    the cut *is*, and `alpha_param_1` wherever the model offers an attraction exponent. A cut whose
    model offers none is a component with Soave's `m`, which is what `alpha="srk"` reads.

    ``model`` is `tbp_cut_properties`' own vocabulary; ``None`` takes NeqSim's default.

    **The card is a value and this returns one.** It stores nothing and is in force nowhere: a
    caller holds it, passes it, and can write it out with :func:`azoth.keycard.to_toml`.

    Raises:
        InvalidInputError: if the three sequences are not the same length, or if a value is not a
            quantity in the unit its field declares.

    See :func:`azoth.eos.from_names` for the call that consumes it.
    """
    from azoth import keycard

    if not len(names) == len(molar_mass) == len(density):
        raise InvalidInputError(
            "names",
            f"a cut table is one name, one molar mass and one density per cut: {len(names)} names, "
            f"{len(molar_mass)} molar masses and {len(density)} densities",
        )
    components: dict[str, dict[str, dict[str, object]]] = {}
    for name, mass, rho in zip(names, molar_mass, density, strict=True):
        cut = tbp_cut_properties(mass, rho, model)
        stated: dict[str, dict[str, object]] = {
            "Tc": {"value": cut.tc.to("K").magnitude, "unit": "K"},
            "Pc": {"value": cut.pc.to("Pa").magnitude, "unit": "Pa"},
            "omega": {"value": cut.acentric_factor, "unit": "dimensionless"},
            "molar_mass": {"value": mass.to("kg/mol").magnitude, "unit": "kg/mol"},
            "normal_liquid_density": {
                "value": rho.to("kg/m**3").magnitude,
                "unit": "kg/m**3",
            },
        }
        if cut.attraction_exponent is not None:
            stated["alpha_param_1"] = {
                "value": cut.attraction_exponent,
                "unit": "dimensionless",
            }
        components[name] = stated
    return keycard.use({"schema_version": 2, "components": components})


__all__ = [
    "assay_mass_fractions",
    "card_from_cuts",
    "characterise_plus_fraction",
    "characterize_to_reference",
    "lumping",
    "oil_assay_derivation",
    "pedersen_plus_split",
    "set_heavy_tbp_to_plus",
    "tbp_closure",
    "tbp_cut_boundaries",
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
_SET_HEAVY_TBP_TO_PLUS = "characterization.set_heavy_tbp_to_plus"
_OIL_ASSAY_DERIVATION = "characterization.oil_assay_derivation"
_TBP_CUT_BOUNDARIES = "characterization.tbp_cut_boundaries"


def tbp_cut_boundaries(
    name_prefix: str,
    cumulative_volume_percent: Sequence[float],
    boiling_point: Sequence[Q],
    specific_gravity: Sequence[float],
) -> TbpCutBoundariesResult:
    """A cumulative TBP curve cut into an assay, one row per interval.

    ``cumulative_volume_percent`` and ``boiling_point`` are the curve's boundaries - N of each -
    and ``specific_gravity`` is one per interval, so one shorter than either. Each cut's share is
    the width of its interval, and its average boiling point is the midpoint of the pair that
    bounds it, which is kept beside the two boundaries themselves.

    The class's Celsius entry point is the same arithmetic with its temperatures offset, so a
    caller states them in whichever unit it holds them in.

    Raises:
        InvalidInputError: if the prefix is empty, if the vectors are not N, N and N-1 long, if the
            yields do not span zero to a hundred, if either the yields or the temperatures fail to
            increase strictly, or if a gravity is not finite and positive.

    See :func:`azoth.characterization.reference.tbp_cut_boundaries`.
    """
    return resolve(_TBP_CUT_BOUNDARIES)(  # type: ignore[no-any-return]
        name_prefix=name_prefix,
        cumulative_volume_percent=cumulative_volume_percent,
        boiling_point=boiling_point,
        specific_gravity=specific_gravity,
    )


def oil_assay_derivation(
    basis: str,
    cuts: Sequence[Mapping[str, Any]],
    total_assay_mass: Q,
) -> OilAssayDerivationResult:
    """An assay's cuts derived to one component each, with a molar mass and a mole amount.

    ``basis`` and the first two fields of every cut are ``assay_mass_fractions``' own; the rest are
    the routes to a molar mass, tried in this order:

    * ``standard_component_name`` - the name of a substance in the databank, whose molar mass is
      read out of it;
    * ``molar_mass`` - stated outright;
    * otherwise the closed form, on ``average_boiling_point`` or one resolved from
      ``watson_characterization_factor``, and on ``density`` or one resolved from ``api_gravity``.

    **A cut that states a molar mass still needs a gravity**, because the class reads the density
    for every cut it derives. A cut with no positive ``declared_fraction`` produces no component, so
    ``component_name`` is what re-aligns the answer with the caller's slate.

    Raises:
        InvalidInputError: if two cuts share a name, if two standard components name one substance,
            if a standard component is on a volume basis or also states pseudo-component
            properties, if a cut reaches the correlation with neither a boiling point nor a gravity,
            or if a standard component is not in the databank.

    See :func:`azoth.characterization.reference.oil_assay_derivation`.
    """
    return resolve(_OIL_ASSAY_DERIVATION)(  # type: ignore[no-any-return]
        basis=basis,
        cuts=cuts,
        total_assay_mass=total_assay_mass,
    )


def set_heavy_tbp_to_plus(
    names: Sequence[str],
    molar_mass: Sequence[Q],
    mole_fraction: Sequence[float],
    normal_liquid_density: Sequence[Q],
) -> SetHeavyTbpToPlusResult:
    """The TBP row a fluid's plus end begins at, chosen by the carbon number in its name.

    ``names`` are the fluid's TBP rows in the order it carries them - a NeqSim fluid names them
    ``C7_PC``, ``C10_PC`` and so on - and the three vectors beside them are those rows' data, each
    aligned with ``names`` entry for entry.

    **The carbon number is read out of the name** with the class's own positional slice, so
    ``C10_PC`` reads ten and a name whose fourth character is not the underscore reads one
    character. ``last_tbp_density`` is the density of the row **below** the winner in index order.

    Raises:
        InvalidInputError: if the columns are not all the same length, if the slate is empty, if a
            name carries no carbon number, or if the winner is the first row - the class reads the
            row below it, and there is none.

    See :func:`azoth.characterization.reference.set_heavy_tbp_to_plus`.
    """
    return resolve(_SET_HEAVY_TBP_TO_PLUS)(  # type: ignore[no-any-return]
        names=names,
        molar_mass=molar_mass,
        mole_fraction=mole_fraction,
        normal_liquid_density=normal_liquid_density,
    )


def assay_mass_fractions(
    basis: str,
    cuts: Sequence[Mapping[str, Any]],
) -> AssayMassFractionsResult:
    """An oil assay's declared yields, resolved to a mass basis.

    An assay states each cut's yield as a mass fraction or as a liquid-volume fraction, never both
    and never neither, and the whole row shares one basis. The row must close on one within
    ``1e-3`` - a tenth of a percent is a **hard error** here, not a warning.

    ``cuts`` is one mapping per cut, keyed by the field names the spec declares:

    * ``declared_fraction`` - the cut's yield on ``basis``, dimensionless, on every cut;
    * ``density`` - the cut's own normal liquid density, *optional per cut*, in kg/m3.

    A volume basis needs every cut's density and refuses the first that states none, naming it; a
    mass basis reads the densities for ``bulk_density`` alone, which is absent unless every cut
    states one.

    Raises:
        InvalidInputError: if the assay has no cuts, if a volume basis reaches a cut without a
            density, or if the declared fractions do not close on one.

    See :func:`azoth.characterization.reference.assay_mass_fractions`.
    """
    return resolve(_ASSAY_MASS_FRACTIONS)(  # type: ignore[no-any-return]
        basis=basis,
        cuts=cuts,
    )


def characterize_to_reference(
    source: Sequence[Mapping[str, Any]],
    reference: Sequence[Mapping[str, Any]],
    inherit_reference_properties: bool | None = None,
) -> CharacterizeToReferenceResult:
    """A fluid's pseudo-components re-cut onto another fluid's cut slate.

    ``source`` and ``reference`` are lists of records, one per pseudo-component, keyed by the field
    names the spec declares. A source row states ``moles`` and ``molar_mass``, and optionally
    ``density`` in kg/m3 and ``boiling_point`` in K; a reference row states ``molar_mass`` and
    optionally ``boiling_point`` and ``density``.

    ``reference`` decides where the cuts are. Set ``inherit_reference_properties`` and each group
    keeps its own moles and takes the **reference row's** molar mass and density instead of the ones
    its own rows average to - Pedersen's common-EoS slate, and the class's own variant of the same
    walk. Absent is the class's default entry point, where the reference's densities are not read.

    The groups are **not** the reference rows - a group that comes out empty is dropped - so
    ``reference_index`` is what re-aligns the answer with the caller's own slate.

    See :func:`azoth.characterization.reference.characterize_to_reference`.
    """
    return resolve(_CHARACTERIZE_TO_REFERENCE)(  # type: ignore[no-any-return]
        source=source,
        reference=reference,
        inherit_reference_properties=inherit_reference_properties,
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
    boiling_point: Q,
    density: Q,
    closure: str,
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
        boiling_point=boiling_point, density=density, closure=closure, model=model
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
