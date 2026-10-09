"""Petroleum-fraction characterisation.

The front door for a fluid that has no databank row. A crude is an *assay* - a list of cuts,
each a molar mass and a density - and everything a cubic needs about a cut follows from those
two numbers by correlation. NeqSim's ``thermo/characterization/`` is that subsystem.

It is a namespace of its own rather than a family under :mod:`azoth.eos`, mirroring NeqSim's
own package split: what is here is not an equation of state and is not reached from one.

* :func:`tbp_cut_properties` - a cut's critical properties, by any of NeqSim's ten models
* :func:`tbp_closure` - a cut's molar mass, from its boiling point and gravity
* :func:`tbp_density` - the same pair inverted for gravity, which only one closure supports
* :func:`pedersen_plus_split` - a plus fraction divided into carbon-number cuts
* :func:`whitson_gamma_split` - the same, on Whitson's three-parameter gamma
* :func:`lumping` - a cut table grouped into equal-mass lumps

# Which implementation answers

As in the other namespaces, each function dispatches to the Rust extension when it is built and
to :mod:`azoth.characterization.reference` otherwise. Both are always reachable - see
:func:`azoth.backends` and :func:`azoth.use_backend`.
"""

from __future__ import annotations

from azoth._dispatch import resolve
from azoth.core.result import (
    LumpingResult,
    PedersenPlusSplitResult,
    TbpClosureResult,
    TbpCutPropertiesResult,
    TbpDensityResult,
    WhitsonGammaSplitResult,
)
from azoth.core.units import Q

__all__ = [
    "lumping",
    "pedersen_plus_split",
    "tbp_closure",
    "tbp_cut_properties",
    "tbp_density",
    "whitson_gamma_split",
]

_TBP_CUT_PROPERTIES = "characterization.tbp_cut_properties"
_TBP_CLOSURE = "characterization.tbp_closure"
_TBP_DENSITY = "characterization.tbp_density"
_PEDERSEN_PLUS_SPLIT = "characterization.pedersen_plus_split"
_LUMPING = "characterization.lumping"
_WHITSON_GAMMA_SPLIT = "characterization.whitson_gamma_split"


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
