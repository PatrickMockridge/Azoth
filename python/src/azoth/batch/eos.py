"""Equations of state, over arrays.

Same arrangement as :mod:`azoth.batch.hydraulics` and :mod:`azoth.batch.thermal`:
plain numbers in the spec's canonical units in, SI base magnitudes out, and a loop
over the scalar reference rather than a second implementation of it.

The one difference here is that there is nothing to convert. Every quantity in this
namespace is dimensionless by construction - the reduced variables and
constitutive coefficients an equation of state is written in - so the batch path
and the scalar path carry exactly the same numbers, and this module has no unit
handling to get wrong.
"""

from __future__ import annotations

from array import array
from collections.abc import Sequence
from dataclasses import dataclass

from azoth.batch._core import run, sequence
from azoth.batch._result import BatchResult
from azoth.core.result import RootStructure
from azoth.core.warnings import Warning

__all__ = [
    "ChungConductivityBatch",
    "ChungViscosityBatch",
    "Co2WaterDiffusivityBatch",
    "CostaldMolarVolumeBatch",
    "HeatOfVaporizationBatch",
    "IdealGasCpBatch",
    "LiquidHeatCapacityBatch",
    "Matcop5PrumrAlphaBatch",
    "MatcopAlphaBatch",
    "MatcopPrAlphaBatch",
    "MatcopPrumrAlphaBatch",
    "MatcopPrumrNewAlphaBatch",
    "MollerupAlphaBatch",
    "NitricSulfuricAcidVaporPressureBatch",
    "ParachorSurfaceTensionBatch",
    "Pr78KappaBatch",
    "PrAlphaAbBatch",
    "PrDaneshAlphaBatch",
    "PrDelft1998AlphaBatch",
    "PrDepartureBatch",
    "PrGassem2001AlphaBatch",
    "PrKappaBatch",
    "PrLeeKeslerAlphaBatch",
    "PrMassDensityBatch",
    "PrMolarVolumeBatch",
    "PrPenelouxShiftBatch",
    "PrZFactorBatch",
    "PrsvKappaBatch",
    "RachfordRiceBinaryBatch",
    "RackettMolarVolumeBatch",
    "RkAlphaAbBatch",
    "RkDepartureBatch",
    "SchwartzentruberAlphaBatch",
    "SoreideWhitsonAlphaBatch",
    "SrkAlphaAbBatch",
    "SrkDepartureBatch",
    "SrkKappaBatch",
    "SrkPenelouxShiftBatch",
    "SrkZFactorBatch",
    "TwuKappaBatch",
    "TwucoonAlphaBatch",
    "TwucoonParamAlphaBatch",
    "TwucoonStatoilAlphaBatch",
    "TynCalusDiffusivityBatch",
    "UmrprAlphaBatch",
    "Vdw1fMixBinaryBatch",
    "WilkeChangDiffusivityBatch",
    "chung_conductivity",
    "chung_viscosity",
    "co2_water_diffusivity",
    "costald_molar_volume",
    "heat_of_vaporization",
    "ideal_gas_cp",
    "liquid_heat_capacity",
    "matcop5_prumr_alpha",
    "matcop_alpha",
    "matcop_pr_alpha",
    "matcop_prumr_alpha",
    "matcop_prumr_new_alpha",
    "mollerup_alpha",
    "nitric_sulfuric_acid_vapor_pressure",
    "parachor_surface_tension",
    "pr78_kappa",
    "pr_alpha_ab",
    "pr_danesh_alpha",
    "pr_delft1998_alpha",
    "pr_departure",
    "pr_gassem2001_alpha",
    "pr_kappa",
    "pr_lee_kesler_alpha",
    "pr_mass_density",
    "pr_molar_volume",
    "pr_peneloux_shift",
    "pr_z_factor",
    "prsv_kappa",
    "rachford_rice_binary",
    "rackett_molar_volume",
    "rk_alpha_ab",
    "rk_departure",
    "schwartzentruber_alpha",
    "soreide_whitson_alpha",
    "srk_alpha_ab",
    "srk_departure",
    "srk_kappa",
    "srk_peneloux_shift",
    "srk_z_factor",
    "twu_kappa",
    "twucoon_alpha",
    "twucoon_param_alpha",
    "twucoon_statoil_alpha",
    "tyn_calus_diffusivity",
    "umrpr_alpha",
    "vdw1f_mix_binary",
    "wilke_chang_diffusivity",
]

_PR_KAPPA = "eos.pr_kappa"
_MATCOP5_PRUMR_ALPHA = "eos.matcop5_prumr_alpha"
_MATCOP_ALPHA = "eos.matcop_alpha"
_MATCOP_PR_ALPHA = "eos.matcop_pr_alpha"
_MATCOP_PRUMR_ALPHA = "eos.matcop_prumr_alpha"
_MATCOP_PRUMR_NEW_ALPHA = "eos.matcop_prumr_new_alpha"
_MOLLERUP_ALPHA = "eos.mollerup_alpha"
_PR_ALPHA_AB = "eos.pr_alpha_ab"
_PR_DANESH_ALPHA = "eos.pr_danesh_alpha"
_PR_DELFT1998_ALPHA = "eos.pr_delft1998_alpha"
_PR_GASSEM2001_ALPHA = "eos.pr_gassem2001_alpha"
_PR_LEE_KESLER_ALPHA = "eos.pr_lee_kesler_alpha"
_PR_Z_FACTOR = "eos.pr_z_factor"
_PRSV_KAPPA = "eos.prsv_kappa"
_PR_DEPARTURE = "eos.pr_departure"
_VDW1F_MIX_BINARY = "eos.vdw1f_mix_binary"
_CHUNG_CONDUCTIVITY = "eos.chung_conductivity"
_CO2_WATER_DIFFUSIVITY = "eos.co2_water_diffusivity"
_PARACHOR_SURFACE_TENSION = "eos.parachor_surface_tension"
_CHUNG_VISCOSITY = "eos.chung_viscosity"
_COSTALD_MOLAR_VOLUME = "eos.costald_molar_volume"
_RACHFORD_RICE_BINARY = "eos.rachford_rice_binary"
_RACKETT_MOLAR_VOLUME = "eos.rackett_molar_volume"
_PR_MOLAR_VOLUME = "eos.pr_molar_volume"
_IDEAL_GAS_CP = "eos.ideal_gas_cp"
_PR_MASS_DENSITY = "eos.pr_mass_density"
_PR_PENELOUX_SHIFT = "eos.pr_peneloux_shift"
_SRK_PENELOUX_SHIFT = "eos.srk_peneloux_shift"
_HEAT_OF_VAPORIZATION = "eos.heat_of_vaporization"
_LIQUID_HEAT_CAPACITY = "eos.liquid_heat_capacity"
_SRK_KAPPA = "eos.srk_kappa"
_SCHWARTZENTRUBER_ALPHA = "eos.schwartzentruber_alpha"
_SOREIDE_WHITSON_ALPHA = "eos.soreide_whitson_alpha"
_SRK_ALPHA_AB = "eos.srk_alpha_ab"
_SRK_Z_FACTOR = "eos.srk_z_factor"
_SRK_DEPARTURE = "eos.srk_departure"
_RK_ALPHA_AB = "eos.rk_alpha_ab"
_RK_DEPARTURE = "eos.rk_departure"
_PR78_KAPPA = "eos.pr78_kappa"
_TWU_KAPPA = "eos.twu_kappa"
_TWUCOON_ALPHA = "eos.twucoon_alpha"
_TWUCOON_PARAM_ALPHA = "eos.twucoon_param_alpha"
_TWUCOON_STATOIL_ALPHA = "eos.twucoon_statoil_alpha"
_TYN_CALUS_DIFFUSIVITY = "eos.tyn_calus_diffusivity"
_NITRIC_SULFURIC_ACID_VAPOR_PRESSURE = "eos.nitric_sulfuric_acid_vapor_pressure"
_UMRPR_ALPHA = "eos.umrpr_alpha"
_WILKE_CHANG_DIFFUSIVITY = "eos.wilke_chang_diffusivity"


@dataclass(frozen=True, slots=True, eq=False, repr=False)
class PrKappaBatch(BatchResult):
    """Result of a batch :func:`azoth.eos.pr_kappa`."""

    #: The Peng-Robinson alpha-function coefficient per element. Dimensionless.
    kappa: array[float]


def _build(
    columns: dict[str, object],
    units: dict[str, str],
    warnings: tuple[tuple[Warning, ...], ...],
) -> PrKappaBatch:
    return PrKappaBatch(warnings=warnings, units=units, kappa=columns["kappa"])  # type: ignore[arg-type]


def pr_kappa(*, omega: Sequence[float]) -> PrKappaBatch:
    """The Peng-Robinson alpha-function coefficient, over an array of acentric factors.

    ``omega`` is dimensionless, so these are the same numbers the scalar API takes
    and there is no unit to be wrong about.

    A negative ``kappa`` is not an error and does not raise: an element whose
    acentric factor is below about -0.2334 produces a negative coefficient and
    carries an ``OUT_OF_VALID_RANGE`` warning in its own warning tuple, exactly as
    the scalar call does. Batch does not change the warning policy - see
    :func:`azoth.eos.pr_kappa` for the calculation itself.
    """
    result: PrKappaBatch = run(
        _PR_KAPPA,
        {"omega": sequence(omega, "omega")},
        _build,
    )
    return result


@dataclass(frozen=True, slots=True, eq=False, repr=False)
class NitricSulfuricAcidVaporPressureBatch(BatchResult):
    """Result of a batch :func:`azoth.eos.nitric_sulfuric_acid_vapor_pressure`."""

    #: The pure-component vapour pressure of water per element.
    p_water: array[float]
    #: The pure-component vapour pressure of nitric acid per element.
    p_nitric_acid: array[float]
    #: The pure-component vapour pressure of sulfuric acid per element.
    p_sulfuric_acid: array[float]


def _build_acid(
    columns: dict[str, object],
    units: dict[str, str],
    warnings: tuple[tuple[Warning, ...], ...],
) -> NitricSulfuricAcidVaporPressureBatch:
    return NitricSulfuricAcidVaporPressureBatch(
        warnings=warnings,
        units=units,
        p_water=columns["p_water"],  # type: ignore[arg-type]
        p_nitric_acid=columns["p_nitric_acid"],  # type: ignore[arg-type]
        p_sulfuric_acid=columns["p_sulfuric_acid"],  # type: ignore[arg-type]
    )


def nitric_sulfuric_acid_vapor_pressure(
    *, T: Sequence[float]
) -> NitricSulfuricAcidVaporPressureBatch:
    """The acid-system pure-component vapour pressures, over an array of temperatures.

    Outside the stated 190-298 K an element carries ``OUT_OF_VALID_RANGE`` in its own
    warning tuple rather than raising, exactly as the scalar call does - the arithmetic
    is defined there. An element at or below the nitric-acid form's 43 K pole does raise,
    because there it is a pole rather than a pressure.

    See :func:`azoth.eos.nitric_sulfuric_acid_vapor_pressure` for the calculation itself.
    """
    result: NitricSulfuricAcidVaporPressureBatch = run(
        _NITRIC_SULFURIC_ACID_VAPOR_PRESSURE,
        {"T": sequence(T, "T")},
        _build_acid,
    )
    return result


@dataclass(frozen=True, slots=True, eq=False, repr=False)
class MatcopAlphaBatch(BatchResult):
    """Result of a batch :func:`azoth.eos.matcop_alpha`."""

    #: The Mathias-Copeman alpha function per element. Dimensionless.
    alpha: array[float]


def _build_matcop(
    columns: dict[str, object],
    units: dict[str, str],
    warnings: tuple[tuple[Warning, ...], ...],
) -> MatcopAlphaBatch:
    return MatcopAlphaBatch(warnings=warnings, units=units, alpha=columns["alpha"])  # type: ignore[arg-type]


def matcop_alpha(
    *, mc1: Sequence[float], mc2: Sequence[float], mc3: Sequence[float], Tr: Sequence[float]
) -> MatcopAlphaBatch:
    """The Mathias-Copeman alpha function, over arrays."""
    result: MatcopAlphaBatch = run(
        _MATCOP_ALPHA,
        {
            "mc1": sequence(mc1, "mc1"),
            "mc2": sequence(mc2, "mc2"),
            "mc3": sequence(mc3, "mc3"),
            "Tr": sequence(Tr, "Tr"),
        },
        _build_matcop,
    )
    return result


@dataclass(frozen=True, slots=True, eq=False, repr=False)
class MatcopPrAlphaBatch(BatchResult):
    """Result of a batch :func:`azoth.eos.matcop_pr_alpha`."""

    #: The Mathias-Copeman alpha function per element. Dimensionless.
    alpha: array[float]


def _build_matcop_pr(
    columns: dict[str, object],
    units: dict[str, str],
    warnings: tuple[tuple[Warning, ...], ...],
) -> MatcopPrAlphaBatch:
    return MatcopPrAlphaBatch(warnings=warnings, units=units, alpha=columns["alpha"])  # type: ignore[arg-type]


def matcop_pr_alpha(
    *,
    omega: Sequence[float],
    mc1: Sequence[float],
    mc2: Sequence[float],
    mc3: Sequence[float],
    Tr: Sequence[float],
) -> MatcopPrAlphaBatch:
    """The Mathias-Copeman alpha with a Peng-Robinson fallback, over arrays."""
    result: MatcopPrAlphaBatch = run(
        _MATCOP_PR_ALPHA,
        {
            "omega": sequence(omega, "omega"),
            "mc1": sequence(mc1, "mc1"),
            "mc2": sequence(mc2, "mc2"),
            "mc3": sequence(mc3, "mc3"),
            "Tr": sequence(Tr, "Tr"),
        },
        _build_matcop_pr,
    )
    return result


@dataclass(frozen=True, slots=True, eq=False, repr=False)
class MatcopPrumrAlphaBatch(BatchResult):
    """Result of a batch :func:`azoth.eos.matcop_prumr_alpha`."""

    #: The Mathias-Copeman alpha function per element. Dimensionless.
    alpha: array[float]


def _build_matcop_prumr(
    columns: dict[str, object],
    units: dict[str, str],
    warnings: tuple[tuple[Warning, ...], ...],
) -> MatcopPrumrAlphaBatch:
    return MatcopPrumrAlphaBatch(warnings=warnings, units=units, alpha=columns["alpha"])  # type: ignore[arg-type]


def matcop_prumr_alpha(
    *,
    omega: Sequence[float],
    mc1: Sequence[float],
    mc2: Sequence[float],
    mc3: Sequence[float],
    Tr: Sequence[float],
) -> MatcopPrumrAlphaBatch:
    """The Mathias-Copeman alpha with the UMR-PR fallback, over arrays."""
    result: MatcopPrumrAlphaBatch = run(
        _MATCOP_PRUMR_ALPHA,
        {
            "omega": sequence(omega, "omega"),
            "mc1": sequence(mc1, "mc1"),
            "mc2": sequence(mc2, "mc2"),
            "mc3": sequence(mc3, "mc3"),
            "Tr": sequence(Tr, "Tr"),
        },
        _build_matcop_prumr,
    )
    return result


@dataclass(frozen=True, slots=True, eq=False, repr=False)
class MatcopPrumrNewAlphaBatch(BatchResult):
    """Result of a batch :func:`azoth.eos.matcop_prumr_new_alpha`."""

    #: The five-parameter Mathias-Copeman alpha function per element. Dimensionless.
    alpha: array[float]


def _build_matcop_prumr_new(
    columns: dict[str, object],
    units: dict[str, str],
    warnings: tuple[tuple[Warning, ...], ...],
) -> MatcopPrumrNewAlphaBatch:
    return MatcopPrumrNewAlphaBatch(warnings=warnings, units=units, alpha=columns["alpha"])  # type: ignore[arg-type]


def matcop_prumr_new_alpha(
    *,
    omega: Sequence[float],
    mc1: Sequence[float],
    mc2: Sequence[float],
    mc3: Sequence[float],
    mc4: Sequence[float],
    mc5: Sequence[float],
    Tr: Sequence[float],
) -> MatcopPrumrNewAlphaBatch:
    """The five-parameter Mathias-Copeman alpha, UMR-PR new variant, over arrays."""
    result: MatcopPrumrNewAlphaBatch = run(
        _MATCOP_PRUMR_NEW_ALPHA,
        {
            "omega": sequence(omega, "omega"),
            "mc1": sequence(mc1, "mc1"),
            "mc2": sequence(mc2, "mc2"),
            "mc3": sequence(mc3, "mc3"),
            "mc4": sequence(mc4, "mc4"),
            "mc5": sequence(mc5, "mc5"),
            "Tr": sequence(Tr, "Tr"),
        },
        _build_matcop_prumr_new,
    )
    return result


@dataclass(frozen=True, slots=True, eq=False, repr=False)
class Matcop5PrumrAlphaBatch(BatchResult):
    """Result of a batch :func:`azoth.eos.matcop5_prumr_alpha`."""

    #: The five-parameter Mathias-Copeman alpha function per element. Dimensionless.
    alpha: array[float]


def _build_matcop5_prumr(
    columns: dict[str, object],
    units: dict[str, str],
    warnings: tuple[tuple[Warning, ...], ...],
) -> Matcop5PrumrAlphaBatch:
    return Matcop5PrumrAlphaBatch(warnings=warnings, units=units, alpha=columns["alpha"])  # type: ignore[arg-type]


def matcop5_prumr_alpha(
    *,
    omega: Sequence[float],
    mc1: Sequence[float],
    mc2: Sequence[float],
    mc3: Sequence[float],
    mc4: Sequence[float],
    mc5: Sequence[float],
    Tr: Sequence[float],
) -> Matcop5PrumrAlphaBatch:
    """The five-parameter Mathias-Copeman alpha function, over arrays."""
    result: Matcop5PrumrAlphaBatch = run(
        _MATCOP5_PRUMR_ALPHA,
        {
            "omega": sequence(omega, "omega"),
            "mc1": sequence(mc1, "mc1"),
            "mc2": sequence(mc2, "mc2"),
            "mc3": sequence(mc3, "mc3"),
            "mc4": sequence(mc4, "mc4"),
            "mc5": sequence(mc5, "mc5"),
            "Tr": sequence(Tr, "Tr"),
        },
        _build_matcop5_prumr,
    )
    return result


@dataclass(frozen=True, slots=True, eq=False, repr=False)
class MollerupAlphaBatch(BatchResult):
    """Result of a batch :func:`azoth.eos.mollerup_alpha`."""

    #: The Mollerup alpha function per element. Dimensionless.
    alpha: array[float]


def _build_mollerup(
    columns: dict[str, object],
    units: dict[str, str],
    warnings: tuple[tuple[Warning, ...], ...],
) -> MollerupAlphaBatch:
    return MollerupAlphaBatch(warnings=warnings, units=units, alpha=columns["alpha"])  # type: ignore[arg-type]


def mollerup_alpha(
    *, p1: Sequence[float], p2: Sequence[float], p3: Sequence[float], Tr: Sequence[float]
) -> MollerupAlphaBatch:
    """The Mollerup alpha function, over arrays."""
    result: MollerupAlphaBatch = run(
        _MOLLERUP_ALPHA,
        {
            "p1": sequence(p1, "p1"),
            "p2": sequence(p2, "p2"),
            "p3": sequence(p3, "p3"),
            "Tr": sequence(Tr, "Tr"),
        },
        _build_mollerup,
    )
    return result


@dataclass(frozen=True, slots=True, eq=False, repr=False)
class PrDaneshAlphaBatch(BatchResult):
    """Result of a batch :func:`azoth.eos.pr_danesh_alpha`."""

    #: The Danesh alpha function per element. Dimensionless.
    alpha: array[float]


def _build_danesh(
    columns: dict[str, object],
    units: dict[str, str],
    warnings: tuple[tuple[Warning, ...], ...],
) -> PrDaneshAlphaBatch:
    return PrDaneshAlphaBatch(warnings=warnings, units=units, alpha=columns["alpha"])  # type: ignore[arg-type]


def pr_danesh_alpha(*, omega: Sequence[float], Tr: Sequence[float]) -> PrDaneshAlphaBatch:
    """The Danesh alpha function, over arrays."""
    result: PrDaneshAlphaBatch = run(
        _PR_DANESH_ALPHA,
        {"omega": sequence(omega, "omega"), "Tr": sequence(Tr, "Tr")},
        _build_danesh,
    )
    return result


@dataclass(frozen=True, slots=True, eq=False, repr=False)
class PrDelft1998AlphaBatch(BatchResult):
    """Result of a batch :func:`azoth.eos.pr_delft1998_alpha`."""

    #: The Peng-Robinson alpha, Delft (1998), per element. Dimensionless.
    alpha: array[float]


def _build_delft1998(
    columns: dict[str, object],
    units: dict[str, str],
    warnings: tuple[tuple[Warning, ...], ...],
) -> PrDelft1998AlphaBatch:
    return PrDelft1998AlphaBatch(warnings=warnings, units=units, alpha=columns["alpha"])  # type: ignore[arg-type]


def pr_delft1998_alpha(*, omega: Sequence[float], Tr: Sequence[float]) -> PrDelft1998AlphaBatch:
    """The Peng-Robinson alpha, Delft (1998), over arrays."""
    result: PrDelft1998AlphaBatch = run(
        _PR_DELFT1998_ALPHA,
        {"omega": sequence(omega, "omega"), "Tr": sequence(Tr, "Tr")},
        _build_delft1998,
    )
    return result


@dataclass(frozen=True, slots=True, eq=False, repr=False)
class PrGassem2001AlphaBatch(BatchResult):
    """Result of a batch :func:`azoth.eos.pr_gassem2001_alpha`."""

    #: The Gassem (2001) alpha function per element. Dimensionless.
    alpha: array[float]


def _build_gassem(
    columns: dict[str, object],
    units: dict[str, str],
    warnings: tuple[tuple[Warning, ...], ...],
) -> PrGassem2001AlphaBatch:
    return PrGassem2001AlphaBatch(
        warnings=warnings,
        units=units,
        alpha=columns["alpha"],  # type: ignore[arg-type]
    )


def pr_gassem2001_alpha(*, omega: Sequence[float], Tr: Sequence[float]) -> PrGassem2001AlphaBatch:
    """The Gassem (2001) alpha function, over arrays."""
    result: PrGassem2001AlphaBatch = run(
        _PR_GASSEM2001_ALPHA,
        {"omega": sequence(omega, "omega"), "Tr": sequence(Tr, "Tr")},
        _build_gassem,
    )
    return result


@dataclass(frozen=True, slots=True, eq=False, repr=False)
class PrLeeKeslerAlphaBatch(BatchResult):
    """Result of a batch :func:`azoth.eos.pr_lee_kesler_alpha`."""

    #: The Peng-Robinson alpha with a Soave-form m-factor, per element. Dimensionless.
    alpha: array[float]


def _build_lee_kesler(
    columns: dict[str, object],
    units: dict[str, str],
    warnings: tuple[tuple[Warning, ...], ...],
) -> PrLeeKeslerAlphaBatch:
    return PrLeeKeslerAlphaBatch(warnings=warnings, units=units, alpha=columns["alpha"])  # type: ignore[arg-type]


def pr_lee_kesler_alpha(*, omega: Sequence[float], Tr: Sequence[float]) -> PrLeeKeslerAlphaBatch:
    """The Peng-Robinson alpha with a Soave-form m-factor, over arrays."""
    result: PrLeeKeslerAlphaBatch = run(
        _PR_LEE_KESLER_ALPHA,
        {"omega": sequence(omega, "omega"), "Tr": sequence(Tr, "Tr")},
        _build_lee_kesler,
    )
    return result


@dataclass(frozen=True, slots=True, eq=False, repr=False)
class PrAlphaAbBatch(BatchResult):
    """Result of a batch :func:`azoth.eos.pr_alpha_ab`."""

    #: The alpha function per element. Dimensionless.
    alpha: array[float]
    #: The cubic's ``A`` per element. Dimensionless.
    a_reduced: array[float]
    #: The cubic's ``B`` per element. Dimensionless.
    b_reduced: array[float]


def _build_alpha_ab(
    columns: dict[str, object],
    units: dict[str, str],
    warnings: tuple[tuple[Warning, ...], ...],
) -> PrAlphaAbBatch:
    return PrAlphaAbBatch(
        warnings=warnings,
        units=units,
        alpha=columns["alpha"],  # type: ignore[arg-type]
        a_reduced=columns["a_reduced"],  # type: ignore[arg-type]
        b_reduced=columns["b_reduced"],  # type: ignore[arg-type]
    )


def pr_alpha_ab(
    *,
    kappa: Sequence[float],
    Tr: Sequence[float],
    Pr: Sequence[float],
) -> PrAlphaAbBatch:
    """The alpha function and reduced parameters, over arrays.

    All six quantities are dimensionless, so ``Tr`` and ``Pr`` here are the bare
    reduced values - not temperatures or pressures in any unit. See
    :func:`azoth.eos.pr_alpha_ab` for the calculation itself.
    """
    result: PrAlphaAbBatch = run(
        _PR_ALPHA_AB,
        {
            "kappa": sequence(kappa, "kappa"),
            "Tr": sequence(Tr, "Tr"),
            "Pr": sequence(Pr, "Pr"),
        },
        _build_alpha_ab,
    )
    return result


@dataclass(frozen=True, slots=True, eq=False, repr=False)
class PrZFactorBatch(BatchResult):
    """Result of a batch :func:`azoth.eos.pr_z_factor`."""

    #: Smallest admissible root per element. Dimensionless.
    z_min: array[float]
    #: Largest admissible root per element. Dimensionless.
    z_max: array[float]
    #: How many admissible roots each element had. A label column: an enum has no
    #: numeric form, and an index into a table would make a caller look the mapping
    #: up to read a value - the same arrangement `reynolds_number`'s regime uses.
    root_structure: tuple[RootStructure | None, ...]
    #: Newton steps the polish took per element.
    iterations: array[float]
    #: Whether each element's polish met its tolerance, as ``1.0`` or ``0.0``.
    converged: array[float]
    #: Final change between iterates per element.
    residual: array[float]


def _build_z_factor(
    columns: dict[str, object],
    units: dict[str, str],
    warnings: tuple[tuple[Warning, ...], ...],
) -> PrZFactorBatch:
    return PrZFactorBatch(
        warnings=warnings,
        units=units,
        z_min=columns["z_min"],  # type: ignore[arg-type]
        z_max=columns["z_max"],  # type: ignore[arg-type]
        root_structure=_structures(columns["root_structure"]),  # type: ignore[arg-type]
        iterations=columns["iterations"],  # type: ignore[arg-type]
        converged=columns["converged"],  # type: ignore[arg-type]
        residual=columns["residual"],  # type: ignore[arg-type]
    )


def _structures(
    labels: Sequence[str | None],
) -> tuple[RootStructure | None, ...]:
    """Rebuild the enum from the label column, `None` where the value was absent."""
    return tuple(None if label is None else RootStructure(label) for label in labels)


def pr_z_factor(*, a_reduced: Sequence[float], b_reduced: Sequence[float]) -> PrZFactorBatch:
    """The Peng-Robinson compressibility factor, over arrays.

    Both arguments are the cubic's dimensionless parameters, so there is no unit to
    be wrong about. The solver report travels with the answer, as it does in the
    scalar API, because the roots are what the polish produced and without
    ``converged`` a caller cannot tell a solution from a failure to converge.

    See :func:`azoth.eos.pr_z_factor` for the calculation itself.
    """
    result: PrZFactorBatch = run(
        _PR_Z_FACTOR,
        {
            "a_reduced": sequence(a_reduced, "a_reduced"),
            "b_reduced": sequence(b_reduced, "b_reduced"),
        },
        _build_z_factor,
    )
    return result


@dataclass(frozen=True, slots=True, eq=False, repr=False)
class PrsvKappaBatch(BatchResult):
    """Result of a batch :func:`azoth.eos.prsv_kappa`."""

    #: The PRSV alpha-function coefficient per element. Dimensionless.
    kappa: array[float]


def _build_prsv(
    columns: dict[str, object],
    units: dict[str, str],
    warnings: tuple[tuple[Warning, ...], ...],
) -> PrsvKappaBatch:
    return PrsvKappaBatch(warnings=warnings, units=units, kappa=columns["kappa"])  # type: ignore[arg-type]


def prsv_kappa(
    *,
    omega: Sequence[float],
    Tr: Sequence[float],
    kappa1: Sequence[float],
) -> PrsvKappaBatch:
    """The PRSV alpha-function coefficient, over arrays.

    ``kappa1`` is the per-substance parameter this library ships no values for, so
    a batch call needs one array of them from the caller - there is nothing to
    broadcast from. See :func:`azoth.eos.prsv_kappa` for the calculation itself.
    """
    result: PrsvKappaBatch = run(
        _PRSV_KAPPA,
        {
            "omega": sequence(omega, "omega"),
            "Tr": sequence(Tr, "Tr"),
            "kappa1": sequence(kappa1, "kappa1"),
        },
        _build_prsv,
    )
    return result


@dataclass(frozen=True, slots=True, eq=False, repr=False)
class PrDepartureBatch(BatchResult):
    """Result of a batch :func:`azoth.eos.pr_departure`."""

    #: Logarithm of the fugacity coefficient per element. Dimensionless.
    ln_phi: array[float]
    #: Departure enthalpy over ``R*T`` per element. Dimensionless.
    h_dep_rt: array[float]
    #: Departure entropy over ``R`` per element. Dimensionless.
    s_dep_r: array[float]
    #: Departure heat capacity over ``R`` per element. Dimensionless.
    cp_dep_r: array[float]


def _build_departure(
    columns: dict[str, object],
    units: dict[str, str],
    warnings: tuple[tuple[Warning, ...], ...],
) -> PrDepartureBatch:
    return PrDepartureBatch(
        warnings=warnings,
        units=units,
        ln_phi=columns["ln_phi"],  # type: ignore[arg-type]
        h_dep_rt=columns["h_dep_rt"],  # type: ignore[arg-type]
        s_dep_r=columns["s_dep_r"],  # type: ignore[arg-type]
        cp_dep_r=columns["cp_dep_r"],  # type: ignore[arg-type]
    )


def pr_departure(
    *,
    a_reduced: Sequence[float],
    b_reduced: Sequence[float],
    z: Sequence[float],
    kappa: Sequence[float],
    Tr: Sequence[float],
) -> PrDepartureBatch:
    """The Peng-Robinson fugacity coefficient and departures, over arrays.

    All five arguments are dimensionless, and so are the three outputs - the
    multiplication by ``R`` and ``T`` happens in the model layer. See
    :func:`azoth.eos.pr_departure` for the calculation itself.
    """
    result: PrDepartureBatch = run(
        _PR_DEPARTURE,
        {
            "a_reduced": sequence(a_reduced, "a_reduced"),
            "b_reduced": sequence(b_reduced, "b_reduced"),
            "z": sequence(z, "z"),
            "kappa": sequence(kappa, "kappa"),
            "Tr": sequence(Tr, "Tr"),
        },
        _build_departure,
    )
    return result


@dataclass(frozen=True, slots=True, eq=False, repr=False)
class Vdw1fMixBinaryBatch(BatchResult):
    """Result of a batch :func:`azoth.eos.vdw1f_mix_binary`."""

    #: The mixture's attraction parameter per element. Dimensionless.
    a_mix: array[float]
    #: The mixture's repulsion parameter per element. Dimensionless.
    b_mix: array[float]


def _build_vdw1f(
    columns: dict[str, object],
    units: dict[str, str],
    warnings: tuple[tuple[Warning, ...], ...],
) -> Vdw1fMixBinaryBatch:
    return Vdw1fMixBinaryBatch(
        warnings=warnings,
        units=units,
        a_mix=columns["a_mix"],  # type: ignore[arg-type]
        b_mix=columns["b_mix"],  # type: ignore[arg-type]
    )


def vdw1f_mix_binary(
    *,
    z1: Sequence[float],
    a1: Sequence[float],
    a2: Sequence[float],
    b1: Sequence[float],
    b2: Sequence[float],
    k12: Sequence[float],
) -> Vdw1fMixBinaryBatch:
    """Van der Waals one-fluid mixing, over arrays.

    Six dimensionless arrays in, two out. ``k12`` is the per-pair parameter this
    library ships no values for, so there is nothing to broadcast from. See
    :func:`azoth.eos.vdw1f_mix_binary` for the calculation itself.
    """
    result: Vdw1fMixBinaryBatch = run(
        _VDW1F_MIX_BINARY,
        {
            "z1": sequence(z1, "z1"),
            "a1": sequence(a1, "a1"),
            "a2": sequence(a2, "a2"),
            "b1": sequence(b1, "b1"),
            "b2": sequence(b2, "b2"),
            "k12": sequence(k12, "k12"),
        },
        _build_vdw1f,
    )
    return result


@dataclass(frozen=True, slots=True, eq=False, repr=False)
class RachfordRiceBinaryBatch(BatchResult):
    """Result of a batch :func:`azoth.eos.rachford_rice_binary`."""

    #: The vapour fraction per element. Dimensionless, and outside ``[0, 1]`` for
    #: elements whose feed is single phase - each of which carries its own warning.
    beta: array[float]


def _build_rachford_rice(
    columns: dict[str, object],
    units: dict[str, str],
    warnings: tuple[tuple[Warning, ...], ...],
) -> RachfordRiceBinaryBatch:
    return RachfordRiceBinaryBatch(
        warnings=warnings,
        units=units,
        beta=columns["beta"],  # type: ignore[arg-type]
    )


def rachford_rice_binary(
    *, z1: Sequence[float], K1: Sequence[float], K2: Sequence[float]
) -> RachfordRiceBinaryBatch:
    """The binary Rachford-Rice vapour fraction, over arrays.

    See :func:`azoth.eos.rachford_rice_binary` for the calculation itself, including
    what a ``beta`` outside ``[0, 1]`` means and why it is a warning.
    """
    result: RachfordRiceBinaryBatch = run(
        _RACHFORD_RICE_BINARY,
        {
            "z1": sequence(z1, "z1"),
            "K1": sequence(K1, "K1"),
            "K2": sequence(K2, "K2"),
        },
        _build_rachford_rice,
    )
    return result


@dataclass(frozen=True, slots=True, eq=False, repr=False)
class PrMolarVolumeBatch(BatchResult):
    """Result of a batch :func:`azoth.eos.pr_molar_volume`."""

    #: Molar volume per element.
    v: array[float]


def _build_molar_volume(
    columns: dict[str, object],
    units: dict[str, str],
    warnings: tuple[tuple[Warning, ...], ...],
) -> PrMolarVolumeBatch:
    return PrMolarVolumeBatch(warnings=warnings, units=units, v=columns["v"])  # type: ignore[arg-type]


def pr_molar_volume(
    *, z: Sequence[float], T: Sequence[float], P: Sequence[float]
) -> PrMolarVolumeBatch:
    """Molar volume, over arrays.

    Plain numbers in the spec's canonical units - ``T`` in K and ``P`` in Pa - and SI
    base magnitudes out, like every other batch call. See
    :func:`azoth.eos.pr_molar_volume` for the calculation itself.
    """
    result: PrMolarVolumeBatch = run(
        _PR_MOLAR_VOLUME,
        {
            "z": sequence(z, "z"),
            "T": sequence(T, "T"),
            "P": sequence(P, "P"),
        },
        _build_molar_volume,
    )
    return result


@dataclass(frozen=True, slots=True, eq=False, repr=False)
class PrMassDensityBatch(BatchResult):
    """Result of a batch :func:`azoth.eos.pr_mass_density`."""

    #: Mass density per element.
    rho: array[float]


def _build_mass_density(
    columns: dict[str, object],
    units: dict[str, str],
    warnings: tuple[tuple[Warning, ...], ...],
) -> PrMassDensityBatch:
    return PrMassDensityBatch(warnings=warnings, units=units, rho=columns["rho"])  # type: ignore[arg-type]


def pr_mass_density(*, M: Sequence[float], v: Sequence[float]) -> PrMassDensityBatch:
    """Mass density, over arrays.

    ``M`` is in kg/mol in the spec's canonical unit, as the scalar API declares - a
    batch call cannot check it, which is the one thing the batch API gives up and
    and the batch API cannot check it. See :func:`azoth.eos.pr_mass_density`.
    """
    result: PrMassDensityBatch = run(
        _PR_MASS_DENSITY,
        {"M": sequence(M, "M"), "v": sequence(v, "v")},
        _build_mass_density,
    )
    return result


@dataclass(frozen=True, slots=True, eq=False, repr=False)
class PrPenelouxShiftBatch(BatchResult):
    """Result of a batch :func:`azoth.eos.pr_peneloux_shift`."""

    #: Volume-translation parameter per element, in m**3/mol.
    c: array[float]


def _build_pr_peneloux_shift(
    columns: dict[str, object],
    units: dict[str, str],
    warnings: tuple[tuple[Warning, ...], ...],
) -> PrPenelouxShiftBatch:
    return PrPenelouxShiftBatch(warnings=warnings, units=units, c=columns["c"])  # type: ignore[arg-type]


def pr_peneloux_shift(
    *, omega: Sequence[float], Tc: Sequence[float], Pc: Sequence[float]
) -> PrPenelouxShiftBatch:
    """The Peng-Robinson Peneloux volume-translation parameter, over arrays."""
    result: PrPenelouxShiftBatch = run(
        _PR_PENELOUX_SHIFT,
        {
            "omega": sequence(omega, "omega"),
            "Tc": sequence(Tc, "Tc"),
            "Pc": sequence(Pc, "Pc"),
        },
        _build_pr_peneloux_shift,
    )
    return result


@dataclass(frozen=True, slots=True, eq=False, repr=False)
class SrkPenelouxShiftBatch(BatchResult):
    """Result of a batch :func:`azoth.eos.srk_peneloux_shift`."""

    #: Volume-translation parameter per element, in m**3/mol.
    c: array[float]


def _build_srk_peneloux_shift(
    columns: dict[str, object],
    units: dict[str, str],
    warnings: tuple[tuple[Warning, ...], ...],
) -> SrkPenelouxShiftBatch:
    return SrkPenelouxShiftBatch(warnings=warnings, units=units, c=columns["c"])  # type: ignore[arg-type]


def srk_peneloux_shift(
    *, omega: Sequence[float], Tc: Sequence[float], Pc: Sequence[float]
) -> SrkPenelouxShiftBatch:
    """The Soave-Redlich-Kwong Peneloux volume-translation parameter, over arrays."""
    result: SrkPenelouxShiftBatch = run(
        _SRK_PENELOUX_SHIFT,
        {
            "omega": sequence(omega, "omega"),
            "Tc": sequence(Tc, "Tc"),
            "Pc": sequence(Pc, "Pc"),
        },
        _build_srk_peneloux_shift,
    )
    return result


@dataclass(frozen=True, slots=True, eq=False, repr=False)
class Co2WaterDiffusivityBatch(BatchResult):
    """Result of a batch :func:`azoth.eos.co2_water_diffusivity`."""

    #: Binary diffusion coefficient per element, in m**2/s.
    d: array[float]


def _build_co2_water_diffusivity(
    columns: dict[str, object],
    units: dict[str, str],
    warnings: tuple[tuple[Warning, ...], ...],
) -> Co2WaterDiffusivityBatch:
    return Co2WaterDiffusivityBatch(warnings=warnings, units=units, d=columns["d"])  # type: ignore[arg-type]


def co2_water_diffusivity(*, T: Sequence[float]) -> Co2WaterDiffusivityBatch:
    """The CO2-in-water binary diffusivity, over arrays."""
    result: Co2WaterDiffusivityBatch = run(
        _CO2_WATER_DIFFUSIVITY,
        {"T": sequence(T, "T")},
        _build_co2_water_diffusivity,
    )
    return result


@dataclass(frozen=True, slots=True, eq=False, repr=False)
class ParachorSurfaceTensionBatch(BatchResult):
    """Result of a batch :func:`azoth.eos.parachor_surface_tension`."""

    #: Surface tension per element, in N/m.
    sigma: array[float]


def _build_parachor_surface_tension(
    columns: dict[str, object],
    units: dict[str, str],
    warnings: tuple[tuple[Warning, ...], ...],
) -> ParachorSurfaceTensionBatch:
    return ParachorSurfaceTensionBatch(
        warnings=warnings,
        units=units,
        sigma=columns["sigma"],  # type: ignore[arg-type]
    )


def parachor_surface_tension(
    *,
    parachor: Sequence[float],
    rho_l: Sequence[float],
    rho_v: Sequence[float],
    M: Sequence[float],
) -> ParachorSurfaceTensionBatch:
    """The surface tension from the parachor correlation, over arrays."""
    result: ParachorSurfaceTensionBatch = run(
        _PARACHOR_SURFACE_TENSION,
        {
            "parachor": sequence(parachor, "parachor"),
            "rho_l": sequence(rho_l, "rho_l"),
            "rho_v": sequence(rho_v, "rho_v"),
            "M": sequence(M, "M"),
        },
        _build_parachor_surface_tension,
    )
    return result


@dataclass(frozen=True, slots=True, eq=False, repr=False)
class ChungConductivityBatch(BatchResult):
    """Result of a batch :func:`azoth.eos.chung_conductivity`."""

    #: Gas thermal conductivity per element, in W/(m*K).
    k: array[float]


def _build_chung_conductivity(
    columns: dict[str, object],
    units: dict[str, str],
    warnings: tuple[tuple[Warning, ...], ...],
) -> ChungConductivityBatch:
    return ChungConductivityBatch(warnings=warnings, units=units, k=columns["k"])  # type: ignore[arg-type]


def chung_conductivity(
    *,
    Cv0: Sequence[float],
    M: Sequence[float],
    omega: Sequence[float],
    Tc: Sequence[float],
    Vc: Sequence[float],
    dipole: Sequence[float],
    kappa: Sequence[float],
    T: Sequence[float],
) -> ChungConductivityBatch:
    """The gas thermal conductivity, over arrays."""
    result: ChungConductivityBatch = run(
        _CHUNG_CONDUCTIVITY,
        {
            "Cv0": sequence(Cv0, "Cv0"),
            "M": sequence(M, "M"),
            "omega": sequence(omega, "omega"),
            "Tc": sequence(Tc, "Tc"),
            "Vc": sequence(Vc, "Vc"),
            "dipole": sequence(dipole, "dipole"),
            "kappa": sequence(kappa, "kappa"),
            "T": sequence(T, "T"),
        },
        _build_chung_conductivity,
    )
    return result


@dataclass(frozen=True, slots=True, eq=False, repr=False)
class ChungViscosityBatch(BatchResult):
    """Result of a batch :func:`azoth.eos.chung_viscosity`."""

    #: Gas dynamic viscosity per element, in Pa*s.
    mu: array[float]


def _build_chung_viscosity(
    columns: dict[str, object],
    units: dict[str, str],
    warnings: tuple[tuple[Warning, ...], ...],
) -> ChungViscosityBatch:
    return ChungViscosityBatch(warnings=warnings, units=units, mu=columns["mu"])  # type: ignore[arg-type]


def chung_viscosity(
    *,
    omega: Sequence[float],
    Tc: Sequence[float],
    Vc: Sequence[float],
    M: Sequence[float],
    dipole: Sequence[float],
    kappa: Sequence[float],
    T: Sequence[float],
    V: Sequence[float],
) -> ChungViscosityBatch:
    """The gas dynamic viscosity, over arrays."""
    result: ChungViscosityBatch = run(
        _CHUNG_VISCOSITY,
        {
            "omega": sequence(omega, "omega"),
            "Tc": sequence(Tc, "Tc"),
            "Vc": sequence(Vc, "Vc"),
            "M": sequence(M, "M"),
            "dipole": sequence(dipole, "dipole"),
            "kappa": sequence(kappa, "kappa"),
            "T": sequence(T, "T"),
            "V": sequence(V, "V"),
        },
        _build_chung_viscosity,
    )
    return result


@dataclass(frozen=True, slots=True, eq=False, repr=False)
class CostaldMolarVolumeBatch(BatchResult):
    """Result of a batch :func:`azoth.eos.costald_molar_volume`."""

    #: Saturated liquid molar volume per element, in m**3/mol.
    v: array[float]


def _build_costald_molar_volume(
    columns: dict[str, object],
    units: dict[str, str],
    warnings: tuple[tuple[Warning, ...], ...],
) -> CostaldMolarVolumeBatch:
    return CostaldMolarVolumeBatch(warnings=warnings, units=units, v=columns["v"])  # type: ignore[arg-type]


def costald_molar_volume(
    *,
    omega: Sequence[float],
    Tc: Sequence[float],
    Vc: Sequence[float],
    M: Sequence[float],
    rho_normal: Sequence[float],
    T: Sequence[float],
) -> CostaldMolarVolumeBatch:
    """The saturated liquid molar volume, over arrays."""
    result: CostaldMolarVolumeBatch = run(
        _COSTALD_MOLAR_VOLUME,
        {
            "omega": sequence(omega, "omega"),
            "Tc": sequence(Tc, "Tc"),
            "Vc": sequence(Vc, "Vc"),
            "M": sequence(M, "M"),
            "rho_normal": sequence(rho_normal, "rho_normal"),
            "T": sequence(T, "T"),
        },
        _build_costald_molar_volume,
    )
    return result


@dataclass(frozen=True, slots=True, eq=False, repr=False)
class RackettMolarVolumeBatch(BatchResult):
    """Result of a batch :func:`azoth.eos.rackett_molar_volume`."""

    #: Saturated liquid molar volume per element, in m**3/mol.
    v: array[float]


def _build_rackett_molar_volume(
    columns: dict[str, object],
    units: dict[str, str],
    warnings: tuple[tuple[Warning, ...], ...],
) -> RackettMolarVolumeBatch:
    return RackettMolarVolumeBatch(warnings=warnings, units=units, v=columns["v"])  # type: ignore[arg-type]


def rackett_molar_volume(
    *, omega: Sequence[float], Tc: Sequence[float], Pc: Sequence[float], T: Sequence[float]
) -> RackettMolarVolumeBatch:
    """The saturated liquid molar volume, over arrays."""
    result: RackettMolarVolumeBatch = run(
        _RACKETT_MOLAR_VOLUME,
        {
            "omega": sequence(omega, "omega"),
            "Tc": sequence(Tc, "Tc"),
            "Pc": sequence(Pc, "Pc"),
            "T": sequence(T, "T"),
        },
        _build_rackett_molar_volume,
    )
    return result


@dataclass(frozen=True, slots=True, eq=False, repr=False)
class HeatOfVaporizationBatch(BatchResult):
    """Result of a batch :func:`azoth.eos.heat_of_vaporization`."""

    #: Heat of vaporisation per element, in J/mol.
    hov: array[float]


def _build_heat_of_vaporization(
    columns: dict[str, object],
    units: dict[str, str],
    warnings: tuple[tuple[Warning, ...], ...],
) -> HeatOfVaporizationBatch:
    return HeatOfVaporizationBatch(warnings=warnings, units=units, hov=columns["hov"])  # type: ignore[arg-type]


def heat_of_vaporization(
    *,
    c0: Sequence[float],
    c1: Sequence[float],
    c2: Sequence[float],
    c3: Sequence[float],
    Tc: Sequence[float],
    T: Sequence[float],
) -> HeatOfVaporizationBatch:
    """The pure-component heat of vaporisation, over arrays."""
    result: HeatOfVaporizationBatch = run(
        _HEAT_OF_VAPORIZATION,
        {
            "c0": sequence(c0, "c0"),
            "c1": sequence(c1, "c1"),
            "c2": sequence(c2, "c2"),
            "c3": sequence(c3, "c3"),
            "Tc": sequence(Tc, "Tc"),
            "T": sequence(T, "T"),
        },
        _build_heat_of_vaporization,
    )
    return result


@dataclass(frozen=True, slots=True, eq=False, repr=False)
class LiquidHeatCapacityBatch(BatchResult):
    """Result of a batch :func:`azoth.eos.liquid_heat_capacity`."""

    #: Liquid heat capacity per element, in J/(mol*K).
    cp: array[float]


def _build_liquid_heat_capacity(
    columns: dict[str, object],
    units: dict[str, str],
    warnings: tuple[tuple[Warning, ...], ...],
) -> LiquidHeatCapacityBatch:
    return LiquidHeatCapacityBatch(warnings=warnings, units=units, cp=columns["cp"])  # type: ignore[arg-type]


def liquid_heat_capacity(
    *,
    c0: Sequence[float],
    c1: Sequence[float],
    c2: Sequence[float],
    c3: Sequence[float],
    c4: Sequence[float],
    T: Sequence[float],
) -> LiquidHeatCapacityBatch:
    """The pure-component liquid heat capacity, over arrays."""
    result: LiquidHeatCapacityBatch = run(
        _LIQUID_HEAT_CAPACITY,
        {
            "c0": sequence(c0, "c0"),
            "c1": sequence(c1, "c1"),
            "c2": sequence(c2, "c2"),
            "c3": sequence(c3, "c3"),
            "c4": sequence(c4, "c4"),
            "T": sequence(T, "T"),
        },
        _build_liquid_heat_capacity,
    )
    return result


@dataclass(frozen=True, slots=True, eq=False, repr=False)
class IdealGasCpBatch(BatchResult):
    """Result of a batch :func:`azoth.eos.ideal_gas_cp`."""

    #: Ideal-gas heat capacity per element, in J/(mol*K).
    cp: array[float]


def _build_ideal_gas_cp(
    columns: dict[str, object],
    units: dict[str, str],
    warnings: tuple[tuple[Warning, ...], ...],
) -> IdealGasCpBatch:
    return IdealGasCpBatch(
        warnings=warnings,
        units=units,
        cp=columns["cp"],  # type: ignore[arg-type]
    )


def ideal_gas_cp(
    *,
    cp_a: Sequence[float],
    cp_b: Sequence[float],
    cp_c: Sequence[float],
    cp_d: Sequence[float],
    cp_e: Sequence[float],
    T: Sequence[float],
) -> IdealGasCpBatch:
    """Ideal-gas heat capacity, over arrays.

    The five coefficients are in J/(mol*K) and one per kelvin per degree; ``T`` is in
    kelvin. ``cp`` comes back in ``J/(mol*K)``, which the result records.

    See :func:`azoth.eos.ideal_gas_cp`.
    """
    result: IdealGasCpBatch = run(
        _IDEAL_GAS_CP,
        {
            "cp_a": sequence(cp_a, "cp_a"),
            "cp_b": sequence(cp_b, "cp_b"),
            "cp_c": sequence(cp_c, "cp_c"),
            "cp_d": sequence(cp_d, "cp_d"),
            "cp_e": sequence(cp_e, "cp_e"),
            "T": sequence(T, "T"),
        },
        _build_ideal_gas_cp,
    )
    return result


@dataclass(frozen=True, slots=True, eq=False, repr=False)
class SrkKappaBatch(BatchResult):
    """Result of a batch :func:`azoth.eos.srk_kappa`."""

    #: The Soave-Redlich-Kwong alpha-function coefficient per element. Dimensionless.
    kappa: array[float]


def _build_srk_kappa(
    columns: dict[str, object],
    units: dict[str, str],
    warnings: tuple[tuple[Warning, ...], ...],
) -> SrkKappaBatch:
    return SrkKappaBatch(warnings=warnings, units=units, kappa=columns["kappa"])  # type: ignore[arg-type]


def srk_kappa(*, omega: Sequence[float]) -> SrkKappaBatch:
    """The Soave-Redlich-Kwong alpha-function coefficient, over an array.

    See :func:`azoth.eos.srk_kappa` for the calculation itself.
    """
    result: SrkKappaBatch = run(
        _SRK_KAPPA,
        {"omega": sequence(omega, "omega")},
        _build_srk_kappa,
    )
    return result


@dataclass(frozen=True, slots=True, eq=False, repr=False)
class SchwartzentruberAlphaBatch(BatchResult):
    """Result of a batch :func:`azoth.eos.schwartzentruber_alpha`."""

    #: The Schwartzentruber-Renon alpha function per element. Dimensionless.
    alpha: array[float]


def _build_schwartzentruber(
    columns: dict[str, object],
    units: dict[str, str],
    warnings: tuple[tuple[Warning, ...], ...],
) -> SchwartzentruberAlphaBatch:
    return SchwartzentruberAlphaBatch(warnings=warnings, units=units, alpha=columns["alpha"])  # type: ignore[arg-type]


def schwartzentruber_alpha(
    *,
    omega: Sequence[float],
    p1: Sequence[float],
    p2: Sequence[float],
    p3: Sequence[float],
    Tr: Sequence[float],
) -> SchwartzentruberAlphaBatch:
    """The Schwartzentruber-Renon alpha function, over arrays."""
    result: SchwartzentruberAlphaBatch = run(
        _SCHWARTZENTRUBER_ALPHA,
        {
            "omega": sequence(omega, "omega"),
            "p1": sequence(p1, "p1"),
            "p2": sequence(p2, "p2"),
            "p3": sequence(p3, "p3"),
            "Tr": sequence(Tr, "Tr"),
        },
        _build_schwartzentruber,
    )
    return result


@dataclass(frozen=True, slots=True, eq=False, repr=False)
class SoreideWhitsonAlphaBatch(BatchResult):
    """Result of a batch :func:`azoth.eos.soreide_whitson_alpha`."""

    #: The Soreide-Whitson alpha function per element. Dimensionless.
    alpha: array[float]


def _build_soreide_whitson(
    columns: dict[str, object],
    units: dict[str, str],
    warnings: tuple[tuple[Warning, ...], ...],
) -> SoreideWhitsonAlphaBatch:
    return SoreideWhitsonAlphaBatch(warnings=warnings, units=units, alpha=columns["alpha"])  # type: ignore[arg-type]


def soreide_whitson_alpha(
    *, salinity: Sequence[float], Tr: Sequence[float]
) -> SoreideWhitsonAlphaBatch:
    """The Soreide-Whitson alpha function, over arrays."""
    result: SoreideWhitsonAlphaBatch = run(
        _SOREIDE_WHITSON_ALPHA,
        {"salinity": sequence(salinity, "salinity"), "Tr": sequence(Tr, "Tr")},
        _build_soreide_whitson,
    )
    return result


@dataclass(frozen=True, slots=True, eq=False, repr=False)
class SrkAlphaAbBatch(BatchResult):
    """Result of a batch :func:`azoth.eos.srk_alpha_ab`."""

    #: The alpha function per element. Dimensionless.
    alpha: array[float]
    #: The cubic's ``A`` per element. Dimensionless.
    a_reduced: array[float]
    #: The cubic's ``B`` per element. Dimensionless.
    b_reduced: array[float]


def _build_srk_alpha_ab(
    columns: dict[str, object],
    units: dict[str, str],
    warnings: tuple[tuple[Warning, ...], ...],
) -> SrkAlphaAbBatch:
    return SrkAlphaAbBatch(
        warnings=warnings,
        units=units,
        alpha=columns["alpha"],  # type: ignore[arg-type]
        a_reduced=columns["a_reduced"],  # type: ignore[arg-type]
        b_reduced=columns["b_reduced"],  # type: ignore[arg-type]
    )


def srk_alpha_ab(
    *,
    kappa: Sequence[float],
    Tr: Sequence[float],
    Pr: Sequence[float],
) -> SrkAlphaAbBatch:
    """The Soave-Redlich-Kwong alpha function and reduced parameters, over arrays.

    See :func:`azoth.eos.srk_alpha_ab` for the calculation itself.
    """
    result: SrkAlphaAbBatch = run(
        _SRK_ALPHA_AB,
        {
            "kappa": sequence(kappa, "kappa"),
            "Tr": sequence(Tr, "Tr"),
            "Pr": sequence(Pr, "Pr"),
        },
        _build_srk_alpha_ab,
    )
    return result


@dataclass(frozen=True, slots=True, eq=False, repr=False)
class SrkZFactorBatch(BatchResult):
    """Result of a batch :func:`azoth.eos.srk_z_factor`."""

    #: Smallest admissible root per element. Dimensionless.
    z_min: array[float]
    #: Largest admissible root per element. Dimensionless.
    z_max: array[float]
    #: How many admissible roots each element had.
    root_structure: tuple[RootStructure | None, ...]
    #: Newton steps the polish took per element.
    iterations: array[float]
    #: Whether each element's polish met its tolerance, as ``1.0`` or ``0.0``.
    converged: array[float]
    #: Final change between iterates per element.
    residual: array[float]


def _build_srk_z_factor(
    columns: dict[str, object],
    units: dict[str, str],
    warnings: tuple[tuple[Warning, ...], ...],
) -> SrkZFactorBatch:
    return SrkZFactorBatch(
        warnings=warnings,
        units=units,
        z_min=columns["z_min"],  # type: ignore[arg-type]
        z_max=columns["z_max"],  # type: ignore[arg-type]
        root_structure=_structures(columns["root_structure"]),  # type: ignore[arg-type]
        iterations=columns["iterations"],  # type: ignore[arg-type]
        converged=columns["converged"],  # type: ignore[arg-type]
        residual=columns["residual"],  # type: ignore[arg-type]
    )


def srk_z_factor(*, a_reduced: Sequence[float], b_reduced: Sequence[float]) -> SrkZFactorBatch:
    """The Soave-Redlich-Kwong compressibility factor, over arrays.

    See :func:`azoth.eos.srk_z_factor` for the calculation itself.
    """
    result: SrkZFactorBatch = run(
        _SRK_Z_FACTOR,
        {
            "a_reduced": sequence(a_reduced, "a_reduced"),
            "b_reduced": sequence(b_reduced, "b_reduced"),
        },
        _build_srk_z_factor,
    )
    return result


@dataclass(frozen=True, slots=True, eq=False, repr=False)
class SrkDepartureBatch(BatchResult):
    """Result of a batch :func:`azoth.eos.srk_departure`."""

    #: Logarithm of the fugacity coefficient per element. Dimensionless.
    ln_phi: array[float]
    #: Departure enthalpy over ``R*T`` per element. Dimensionless.
    h_dep_rt: array[float]
    #: Departure entropy over ``R`` per element. Dimensionless.
    s_dep_r: array[float]
    #: Departure heat capacity over ``R`` per element. Dimensionless.
    cp_dep_r: array[float]


def _build_srk_departure(
    columns: dict[str, object],
    units: dict[str, str],
    warnings: tuple[tuple[Warning, ...], ...],
) -> SrkDepartureBatch:
    return SrkDepartureBatch(
        warnings=warnings,
        units=units,
        ln_phi=columns["ln_phi"],  # type: ignore[arg-type]
        h_dep_rt=columns["h_dep_rt"],  # type: ignore[arg-type]
        s_dep_r=columns["s_dep_r"],  # type: ignore[arg-type]
        cp_dep_r=columns["cp_dep_r"],  # type: ignore[arg-type]
    )


def srk_departure(
    *,
    a_reduced: Sequence[float],
    b_reduced: Sequence[float],
    z: Sequence[float],
    kappa: Sequence[float],
    Tr: Sequence[float],
) -> SrkDepartureBatch:
    """The Soave-Redlich-Kwong fugacity coefficient and departures, over arrays.

    See :func:`azoth.eos.srk_departure` for the calculation itself.
    """
    result: SrkDepartureBatch = run(
        _SRK_DEPARTURE,
        {
            "a_reduced": sequence(a_reduced, "a_reduced"),
            "b_reduced": sequence(b_reduced, "b_reduced"),
            "z": sequence(z, "z"),
            "kappa": sequence(kappa, "kappa"),
            "Tr": sequence(Tr, "Tr"),
        },
        _build_srk_departure,
    )
    return result


@dataclass(frozen=True, slots=True, eq=False, repr=False)
class RkAlphaAbBatch(BatchResult):
    """Result of a batch :func:`azoth.eos.rk_alpha_ab`."""

    #: The alpha function per element. Dimensionless.
    alpha: array[float]
    #: The cubic's ``A`` per element. Dimensionless.
    a_reduced: array[float]
    #: The cubic's ``B`` per element. Dimensionless.
    b_reduced: array[float]


def _build_rk_alpha_ab(
    columns: dict[str, object],
    units: dict[str, str],
    warnings: tuple[tuple[Warning, ...], ...],
) -> RkAlphaAbBatch:
    return RkAlphaAbBatch(
        warnings=warnings,
        units=units,
        alpha=columns["alpha"],  # type: ignore[arg-type]
        a_reduced=columns["a_reduced"],  # type: ignore[arg-type]
        b_reduced=columns["b_reduced"],  # type: ignore[arg-type]
    )


def rk_alpha_ab(*, Tr: Sequence[float], Pr: Sequence[float]) -> RkAlphaAbBatch:
    """The Redlich-Kwong alpha function and reduced parameters, over arrays.

    See :func:`azoth.eos.rk_alpha_ab` for the calculation itself.
    """
    result: RkAlphaAbBatch = run(
        _RK_ALPHA_AB,
        {
            "Tr": sequence(Tr, "Tr"),
            "Pr": sequence(Pr, "Pr"),
        },
        _build_rk_alpha_ab,
    )
    return result


@dataclass(frozen=True, slots=True, eq=False, repr=False)
class RkDepartureBatch(BatchResult):
    """Result of a batch :func:`azoth.eos.rk_departure`."""

    #: Logarithm of the fugacity coefficient per element. Dimensionless.
    ln_phi: array[float]
    #: Departure enthalpy over ``R*T`` per element. Dimensionless.
    h_dep_rt: array[float]
    #: Departure entropy over ``R`` per element. Dimensionless.
    s_dep_r: array[float]
    #: Departure heat capacity over ``R`` per element. Dimensionless.
    cp_dep_r: array[float]


def _build_rk_departure(
    columns: dict[str, object],
    units: dict[str, str],
    warnings: tuple[tuple[Warning, ...], ...],
) -> RkDepartureBatch:
    return RkDepartureBatch(
        warnings=warnings,
        units=units,
        ln_phi=columns["ln_phi"],  # type: ignore[arg-type]
        h_dep_rt=columns["h_dep_rt"],  # type: ignore[arg-type]
        s_dep_r=columns["s_dep_r"],  # type: ignore[arg-type]
        cp_dep_r=columns["cp_dep_r"],  # type: ignore[arg-type]
    )


def rk_departure(
    *,
    a_reduced: Sequence[float],
    b_reduced: Sequence[float],
    z: Sequence[float],
) -> RkDepartureBatch:
    """The Redlich-Kwong fugacity coefficient and departures, over arrays.

    See :func:`azoth.eos.rk_departure` for the calculation itself.
    """
    result: RkDepartureBatch = run(
        _RK_DEPARTURE,
        {
            "a_reduced": sequence(a_reduced, "a_reduced"),
            "b_reduced": sequence(b_reduced, "b_reduced"),
            "z": sequence(z, "z"),
        },
        _build_rk_departure,
    )
    return result


@dataclass(frozen=True, slots=True, eq=False, repr=False)
class Pr78KappaBatch(BatchResult):
    """Result of a batch :func:`azoth.eos.pr78_kappa`."""

    #: The 1978 Peng-Robinson alpha-function coefficient per element. Dimensionless.
    kappa: array[float]


def _build_pr78_kappa(
    columns: dict[str, object],
    units: dict[str, str],
    warnings: tuple[tuple[Warning, ...], ...],
) -> Pr78KappaBatch:
    return Pr78KappaBatch(warnings=warnings, units=units, kappa=columns["kappa"])  # type: ignore[arg-type]


def pr78_kappa(*, omega: Sequence[float]) -> Pr78KappaBatch:
    """The 1978 Peng-Robinson alpha-function coefficient, over an array.

    See :func:`azoth.eos.pr78_kappa` for the calculation itself.
    """
    result: Pr78KappaBatch = run(
        _PR78_KAPPA,
        {"omega": sequence(omega, "omega")},
        _build_pr78_kappa,
    )
    return result


@dataclass(frozen=True, slots=True, eq=False, repr=False)
class TwuKappaBatch(BatchResult):
    """Result of a batch :func:`azoth.eos.twu_kappa`."""

    #: Twu's alpha-function coefficient per element. Dimensionless.
    kappa: array[float]


def _build_twu_kappa(
    columns: dict[str, object],
    units: dict[str, str],
    warnings: tuple[tuple[Warning, ...], ...],
) -> TwuKappaBatch:
    return TwuKappaBatch(warnings=warnings, units=units, kappa=columns["kappa"])  # type: ignore[arg-type]


def twu_kappa(*, omega: Sequence[float]) -> TwuKappaBatch:
    """Twu's alpha-function coefficient, over an array.

    See :func:`azoth.eos.twu_kappa` for the calculation itself.
    """
    result: TwuKappaBatch = run(
        _TWU_KAPPA,
        {"omega": sequence(omega, "omega")},
        _build_twu_kappa,
    )
    return result


@dataclass(frozen=True, slots=True, eq=False, repr=False)
class TwucoonAlphaBatch(BatchResult):
    """Result of a batch :func:`azoth.eos.twucoon_alpha`."""

    #: The Twu-Coon alpha function per element. Dimensionless.
    alpha: array[float]


def _build_twucoon(
    columns: dict[str, object],
    units: dict[str, str],
    warnings: tuple[tuple[Warning, ...], ...],
) -> TwucoonAlphaBatch:
    return TwucoonAlphaBatch(warnings=warnings, units=units, alpha=columns["alpha"])  # type: ignore[arg-type]


def twucoon_alpha(*, omega: Sequence[float], Tr: Sequence[float]) -> TwucoonAlphaBatch:
    """The Twu-Coon alpha function, over arrays."""
    result: TwucoonAlphaBatch = run(
        _TWUCOON_ALPHA,
        {"omega": sequence(omega, "omega"), "Tr": sequence(Tr, "Tr")},
        _build_twucoon,
    )
    return result


@dataclass(frozen=True, slots=True, eq=False, repr=False)
class TwucoonParamAlphaBatch(BatchResult):
    """Result of a batch :func:`azoth.eos.twucoon_param_alpha`."""

    #: The Twu-Coon parameter alpha function per element. Dimensionless.
    alpha: array[float]


def _build_twucoon_param(
    columns: dict[str, object],
    units: dict[str, str],
    warnings: tuple[tuple[Warning, ...], ...],
) -> TwucoonParamAlphaBatch:
    return TwucoonParamAlphaBatch(warnings=warnings, units=units, alpha=columns["alpha"])  # type: ignore[arg-type]


def twucoon_param_alpha(
    *, a: Sequence[float], b: Sequence[float], c: Sequence[float], Tr: Sequence[float]
) -> TwucoonParamAlphaBatch:
    """The Twu-Coon parameter alpha function, over arrays."""
    result: TwucoonParamAlphaBatch = run(
        _TWUCOON_PARAM_ALPHA,
        {
            "a": sequence(a, "a"),
            "b": sequence(b, "b"),
            "c": sequence(c, "c"),
            "Tr": sequence(Tr, "Tr"),
        },
        _build_twucoon_param,
    )
    return result


@dataclass(frozen=True, slots=True, eq=False, repr=False)
class TwucoonStatoilAlphaBatch(BatchResult):
    """Result of a batch :func:`azoth.eos.twucoon_statoil_alpha`."""

    #: The Twu-Coon Statoil alpha function per element. Dimensionless.
    alpha: array[float]


def _build_twucoon_statoil(
    columns: dict[str, object],
    units: dict[str, str],
    warnings: tuple[tuple[Warning, ...], ...],
) -> TwucoonStatoilAlphaBatch:
    return TwucoonStatoilAlphaBatch(warnings=warnings, units=units, alpha=columns["alpha"])  # type: ignore[arg-type]


def twucoon_statoil_alpha(
    *, a: Sequence[float], b: Sequence[float], c: Sequence[float], Tr: Sequence[float]
) -> TwucoonStatoilAlphaBatch:
    """The Twu-Coon Statoil alpha function, over arrays."""
    result: TwucoonStatoilAlphaBatch = run(
        _TWUCOON_STATOIL_ALPHA,
        {
            "a": sequence(a, "a"),
            "b": sequence(b, "b"),
            "c": sequence(c, "c"),
            "Tr": sequence(Tr, "Tr"),
        },
        _build_twucoon_statoil,
    )
    return result


@dataclass(frozen=True, slots=True, eq=False, repr=False)
class TynCalusDiffusivityBatch(BatchResult):
    """Result of a batch :func:`azoth.eos.tyn_calus_diffusivity`."""

    #: Binary diffusion coefficient per element, in m**2/s.
    d: array[float]


def _build_tyn_calus_diffusivity(
    columns: dict[str, object],
    units: dict[str, str],
    warnings: tuple[tuple[Warning, ...], ...],
) -> TynCalusDiffusivityBatch:
    return TynCalusDiffusivityBatch(warnings=warnings, units=units, d=columns["d"])  # type: ignore[arg-type]


def tyn_calus_diffusivity(
    *, VA: Sequence[float], VB: Sequence[float], T: Sequence[float], eta: Sequence[float]
) -> TynCalusDiffusivityBatch:
    """The liquid binary diffusivity, over arrays."""
    result: TynCalusDiffusivityBatch = run(
        _TYN_CALUS_DIFFUSIVITY,
        {
            "VA": sequence(VA, "VA"),
            "VB": sequence(VB, "VB"),
            "T": sequence(T, "T"),
            "eta": sequence(eta, "eta"),
        },
        _build_tyn_calus_diffusivity,
    )
    return result


@dataclass(frozen=True, slots=True, eq=False, repr=False)
class UmrprAlphaBatch(BatchResult):
    """Result of a batch :func:`azoth.eos.umrpr_alpha`."""

    #: The UMR-PR alpha function per element. Dimensionless.
    alpha: array[float]


def _build_umrpr(
    columns: dict[str, object],
    units: dict[str, str],
    warnings: tuple[tuple[Warning, ...], ...],
) -> UmrprAlphaBatch:
    return UmrprAlphaBatch(warnings=warnings, units=units, alpha=columns["alpha"])  # type: ignore[arg-type]


def umrpr_alpha(*, omega: Sequence[float], Tr: Sequence[float]) -> UmrprAlphaBatch:
    """The UMR-PR alpha function, over arrays."""
    result: UmrprAlphaBatch = run(
        _UMRPR_ALPHA,
        {"omega": sequence(omega, "omega"), "Tr": sequence(Tr, "Tr")},
        _build_umrpr,
    )
    return result


@dataclass(frozen=True, slots=True, eq=False, repr=False)
class WilkeChangDiffusivityBatch(BatchResult):
    """Result of a batch :func:`azoth.eos.wilke_chang_diffusivity`."""

    #: Binary diffusion coefficient per element, in m**2/s.
    d: array[float]


def _build_wilke_chang_diffusivity(
    columns: dict[str, object],
    units: dict[str, str],
    warnings: tuple[tuple[Warning, ...], ...],
) -> WilkeChangDiffusivityBatch:
    return WilkeChangDiffusivityBatch(warnings=warnings, units=units, d=columns["d"])  # type: ignore[arg-type]


def wilke_chang_diffusivity(
    *,
    phi: Sequence[float],
    M: Sequence[float],
    T: Sequence[float],
    eta: Sequence[float],
    VA: Sequence[float],
) -> WilkeChangDiffusivityBatch:
    """The liquid binary diffusivity, over arrays."""
    result: WilkeChangDiffusivityBatch = run(
        _WILKE_CHANG_DIFFUSIVITY,
        {
            "phi": sequence(phi, "phi"),
            "M": sequence(M, "M"),
            "T": sequence(T, "T"),
            "eta": sequence(eta, "eta"),
            "VA": sequence(VA, "VA"),
        },
        _build_wilke_chang_diffusivity,
    )
    return result
