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
    "PrAlphaAbBatch",
    "PrDepartureBatch",
    "PrKappaBatch",
    "PrZFactorBatch",
    "PrsvKappaBatch",
    "pr_alpha_ab",
    "pr_departure",
    "pr_kappa",
    "pr_z_factor",
    "prsv_kappa",
]

_PR_KAPPA = "eos.pr_kappa"
_PR_ALPHA_AB = "eos.pr_alpha_ab"
_PR_Z_FACTOR = "eos.pr_z_factor"
_PRSV_KAPPA = "eos.prsv_kappa"
_PR_DEPARTURE = "eos.pr_departure"


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
