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
    "IdealGasCpBatch",
    "PrAlphaAbBatch",
    "PrDepartureBatch",
    "PrKappaBatch",
    "PrMassDensityBatch",
    "PrMolarVolumeBatch",
    "PrZFactorBatch",
    "PrsvKappaBatch",
    "RachfordRiceBinaryBatch",
    "Vdw1fMixBinaryBatch",
    "ideal_gas_cp",
    "pr_alpha_ab",
    "pr_departure",
    "pr_kappa",
    "pr_mass_density",
    "pr_molar_volume",
    "pr_z_factor",
    "prsv_kappa",
    "rachford_rice_binary",
    "vdw1f_mix_binary",
]

_PR_KAPPA = "eos.pr_kappa"
_PR_ALPHA_AB = "eos.pr_alpha_ab"
_PR_Z_FACTOR = "eos.pr_z_factor"
_PRSV_KAPPA = "eos.prsv_kappa"
_PR_DEPARTURE = "eos.pr_departure"
_VDW1F_MIX_BINARY = "eos.vdw1f_mix_binary"
_RACHFORD_RICE_BINARY = "eos.rachford_rice_binary"
_PR_MOLAR_VOLUME = "eos.pr_molar_volume"
_IDEAL_GAS_CP = "eos.ideal_gas_cp"
_PR_MASS_DENSITY = "eos.pr_mass_density"


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
    ``docs/src/batch.md`` says so. See :func:`azoth.eos.pr_mass_density`.
    """
    result: PrMassDensityBatch = run(
        _PR_MASS_DENSITY,
        {"M": sequence(M, "M"), "v": sequence(v, "v")},
        _build_mass_density,
    )
    return result


@dataclass(frozen=True, slots=True, eq=False, repr=False)
class IdealGasCpBatch(BatchResult):
    """Result of a batch :func:`azoth.eos.ideal_gas_cp`."""

    #: The polynomial's dimensionless value per element.
    cp_over_r: array[float]
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
        cp_over_r=columns["cp_over_r"],  # type: ignore[arg-type]
        cp=columns["cp"],  # type: ignore[arg-type]
    )


def ideal_gas_cp(
    *,
    a: Sequence[float],
    b: Sequence[float],
    c: Sequence[float],
    d: Sequence[float],
    T: Sequence[float],
) -> IdealGasCpBatch:
    """Ideal-gas heat capacity, over arrays.

    ``T`` is in kelvin. This is the first calc in this namespace whose unit is not
    dimensionless, so it is the first batch arm here that carries a unit back out -
    ``cp`` is in ``J/(mol*K)``, which the result records.

    See :func:`azoth.eos.ideal_gas_cp`.
    """
    result: IdealGasCpBatch = run(
        _IDEAL_GAS_CP,
        {
            "a": sequence(a, "a"),
            "b": sequence(b, "b"),
            "c": sequence(c, "c"),
            "d": sequence(d, "d"),
            "T": sequence(T, "T"),
        },
        _build_ideal_gas_cp,
    )
    return result
