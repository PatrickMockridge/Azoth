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
    "HeatOfVaporizationBatch",
    "IdealGasCpBatch",
    "LiquidHeatCapacityBatch",
    "Pr78KappaBatch",
    "PrAlphaAbBatch",
    "PrDepartureBatch",
    "PrKappaBatch",
    "PrMassDensityBatch",
    "PrMolarVolumeBatch",
    "PrPenelouxShiftBatch",
    "PrZFactorBatch",
    "PrsvKappaBatch",
    "RachfordRiceBinaryBatch",
    "RackettMolarVolumeBatch",
    "RkAlphaAbBatch",
    "RkDepartureBatch",
    "SrkAlphaAbBatch",
    "SrkDepartureBatch",
    "SrkKappaBatch",
    "SrkPenelouxShiftBatch",
    "SrkZFactorBatch",
    "TwuKappaBatch",
    "Vdw1fMixBinaryBatch",
    "heat_of_vaporization",
    "ideal_gas_cp",
    "liquid_heat_capacity",
    "pr78_kappa",
    "pr_alpha_ab",
    "pr_departure",
    "pr_kappa",
    "pr_mass_density",
    "pr_molar_volume",
    "pr_peneloux_shift",
    "pr_z_factor",
    "prsv_kappa",
    "rachford_rice_binary",
    "rackett_molar_volume",
    "rk_alpha_ab",
    "rk_departure",
    "srk_alpha_ab",
    "srk_departure",
    "srk_kappa",
    "srk_peneloux_shift",
    "srk_z_factor",
    "twu_kappa",
    "vdw1f_mix_binary",
]

_PR_KAPPA = "eos.pr_kappa"
_PR_ALPHA_AB = "eos.pr_alpha_ab"
_PR_Z_FACTOR = "eos.pr_z_factor"
_PRSV_KAPPA = "eos.prsv_kappa"
_PR_DEPARTURE = "eos.pr_departure"
_VDW1F_MIX_BINARY = "eos.vdw1f_mix_binary"
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
_SRK_ALPHA_AB = "eos.srk_alpha_ab"
_SRK_Z_FACTOR = "eos.srk_z_factor"
_SRK_DEPARTURE = "eos.srk_departure"
_RK_ALPHA_AB = "eos.rk_alpha_ab"
_RK_DEPARTURE = "eos.rk_departure"
_PR78_KAPPA = "eos.pr78_kappa"
_TWU_KAPPA = "eos.twu_kappa"


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
