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
from azoth.core.warnings import Warning

__all__ = ["PrAlphaAbBatch", "PrKappaBatch", "pr_alpha_ab", "pr_kappa"]

_PR_KAPPA = "eos.pr_kappa"
_PR_ALPHA_AB = "eos.pr_alpha_ab"


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
