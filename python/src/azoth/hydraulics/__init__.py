"""Hydraulics calculations.

The slice implemented here runs from the Reynolds number through to the
Darcy-Weisbach pressure drop for a straight pipe, with fitting losses available
separately:

* :func:`reynolds_number` - the flow regime, and the input every other calc needs
* :func:`friction_factor_colebrook` - the implicit, accurate friction factor
* :func:`friction_factor_swamee_jain` - the explicit approximation to it
* :func:`friction_factor_haaland` - a second explicit approximation, fitted differently
* :func:`crane_k_factors` - fitting losses by the equivalent-length method
* :func:`darcy_weisbach` - pressure drop over a straight pipe

Pipe *with* fittings is a composition of the last two, performed by the
``azoth pipe`` CLI rather than by a calc of its own, because the two losses are
computed by different methods and adding them is a modelling decision the caller
should be able to see.

# Which implementation answers

Each function below dispatches to the Rust extension when it is built, and to
:mod:`azoth.hydraulics.reference` otherwise. Both are always reachable - see
:func:`azoth.backends` and :func:`azoth.use_backend`.

# Warning before use

:func:`crane_k_factors` reads coefficients from
``data/fittings/crane_k_factors.csv``, where every row is currently an **estimated
dummy value** - a placeholder for software testing, not engineering data. Results
built from it carry an ``ESTIMATED_DATA`` warning. Nothing in this package should
be used for design work until that file is populated from a primary standard.
"""

from __future__ import annotations

from collections.abc import Sequence

from azoth._dispatch import resolve
from azoth.core.result import (
    ColebrookResult,
    DarcyWeisbachResult,
    HaalandResult,
    KFactorsResult,
    ReynoldsNumberResult,
    SwameeJainResult,
)
from azoth.core.units import Q

__all__ = [
    "crane_k_factors",
    "darcy_weisbach",
    "friction_factor_colebrook",
    "friction_factor_haaland",
    "friction_factor_swamee_jain",
    "reynolds_number",
]

_REYNOLDS_NUMBER = "hydraulics.reynolds_number"
_COLEBROOK = "hydraulics.friction_factor_colebrook"
_SWAMEE_JAIN = "hydraulics.friction_factor_swamee_jain"
_HAALAND = "hydraulics.friction_factor_haaland"
_CRANE_K = "hydraulics.crane_k_factors"
_DARCY_WEISBACH = "hydraulics.darcy_weisbach"


def reynolds_number(rho: Q, v: Q, D: Q, mu: Q) -> ReynoldsNumberResult:
    """Reynolds number for flow in a circular pipe.

    Returns the Reynolds number and its flow regime. A transitional result carries
    a ``TRANSITIONAL_FLOW`` warning, because the friction factor is indeterminate
    in that band.

    Raises:
        OutOfRangeError: if density, diameter or viscosity is not positive, or if
            velocity is negative.

    See :func:`azoth.hydraulics.reference.reynolds_number` for the full
    description; this wrapper only chooses which implementation runs.
    """
    return resolve(_REYNOLDS_NUMBER)(rho=rho, v=v, D=D, mu=mu)  # type: ignore[no-any-return]


def friction_factor_colebrook(re: float, relative_roughness: float) -> ColebrookResult:
    """Solve the Colebrook-White equation for the Darcy friction factor.

    Raises:
        OutOfRangeError: if ``re <= 0`` or ``relative_roughness < 0``.
        SolverNotConvergedError: if the iteration hits its cap.

    See :func:`azoth.hydraulics.reference.friction_factor_colebrook`.
    """
    return resolve(_COLEBROOK)(re=re, relative_roughness=relative_roughness)  # type: ignore[no-any-return]


def friction_factor_swamee_jain(re: float, relative_roughness: float) -> SwameeJainResult:
    """Explicit Swamee-Jain approximation to the Colebrook friction factor.

    Raises:
        OutOfRangeError: if ``re <= 0`` or ``relative_roughness < 0``.

    See :func:`azoth.hydraulics.reference.friction_factor_swamee_jain`.
    """
    return resolve(_SWAMEE_JAIN)(re=re, relative_roughness=relative_roughness)  # type: ignore[no-any-return]


def friction_factor_haaland(re: float, relative_roughness: float) -> HaalandResult:
    """Explicit Haaland approximation to the Colebrook friction factor.

    The second explicit approximation in this package, alongside Swamee-Jain. They
    were fitted differently and are accurate to about 1% and 2% respectively, so
    evaluating both at the same inputs is a cheap way to see how much the choice of
    explicit form matters.

    Raises:
        OutOfRangeError: if ``re <= 0`` or ``relative_roughness < 0``.

    See :func:`azoth.hydraulics.reference.friction_factor_haaland`.
    """
    return resolve(_HAALAND)(re=re, relative_roughness=relative_roughness)  # type: ignore[no-any-return]


def crane_k_factors(fittings: Sequence[str], f_t: float) -> KFactorsResult:
    """Total resistance coefficient for a list of fittings.

    Raises:
        UnknownFittingError: for an id not in the registry.
        OutOfRangeError: if ``f_t`` is not positive, or the list is empty.

    See :func:`azoth.hydraulics.reference.crane_k_factors`.
    """
    return resolve(_CRANE_K)(fittings=list(fittings), f_t=f_t)  # type: ignore[no-any-return]


def darcy_weisbach(
    f: float,
    L: Q,
    D: Q,
    rho: Q,
    v: Q,
    mu: Q | None = None,
) -> DarcyWeisbachResult:
    """Pressure drop over a length of straight pipe.

    ``mu`` is optional and is used only to check the flow regime. Omit it and the
    result carries a ``RANGE_CHECK_SKIPPED`` warning saying the regime went
    unchecked, rather than silently reporting a value that was never validated.

    Raises:
        OutOfRangeError: on a hard bound violation in ``f``, ``L``, ``D``, ``rho``
            or ``v``.

    See :func:`azoth.hydraulics.reference.darcy_weisbach`.
    """
    return resolve(_DARCY_WEISBACH)(  # type: ignore[no-any-return]
        f=f, L=L, D=D, rho=rho, v=v, mu=mu
    )
