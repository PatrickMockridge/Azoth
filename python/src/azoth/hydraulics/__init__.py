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
* :func:`pump_power` - shaft power from flow, head and efficiency
* :func:`orifice_flow` - flow through an orifice from its pressure difference
* :func:`control_valve_cv` - liquid flow through a control valve
* :func:`choked_flow_area` - the throat area a choked gas flow needs

Pipe *with* fittings is a composition of the last two, performed by the
``azoth pipe`` CLI rather than by a calc of its own, because the two losses are
computed by different methods and adding them is a modelling decision the caller
should be able to see.

# Which implementation answers

Each function below dispatches to the Rust extension when it is built, and to
:mod:`azoth.hydraulics.reference` otherwise. Both are always reachable - see
:func:`azoth.backends` and :func:`azoth.use_backend`.

# Not for design work yet

:func:`crane_k_factors` reads coefficients from
``data/fittings/crane_k_factors.csv``, where every row is currently an **estimated
dummy value** - a placeholder for software testing, not engineering data. Nothing
in this package should be used for design work until that file is populated from
a primary standard. The placeholder status is recorded in the data file and in
the calc's spec, not announced by a warning on the result.
"""

from __future__ import annotations

from collections.abc import Sequence

from azoth._dispatch import resolve
from azoth.core.result import (
    ChokedFlowAreaResult,
    ColebrookResult,
    ControlValveCvResult,
    DarcyWeisbachResult,
    HaalandResult,
    KFactorsResult,
    OrificeFlowResult,
    PumpPowerResult,
    ReynoldsNumberResult,
    SwameeJainResult,
)
from azoth.core.units import Q

__all__ = [
    "choked_flow_area",
    "control_valve_cv",
    "crane_k_factors",
    "darcy_weisbach",
    "friction_factor_colebrook",
    "friction_factor_haaland",
    "friction_factor_swamee_jain",
    "orifice_flow",
    "pump_power",
    "reynolds_number",
]

_REYNOLDS_NUMBER = "hydraulics.reynolds_number"
_COLEBROOK = "hydraulics.friction_factor_colebrook"
_SWAMEE_JAIN = "hydraulics.friction_factor_swamee_jain"
_HAALAND = "hydraulics.friction_factor_haaland"
_CRANE_K = "hydraulics.crane_k_factors"
_DARCY_WEISBACH = "hydraulics.darcy_weisbach"
_PUMP_POWER = "hydraulics.pump_power"
_ORIFICE_FLOW = "hydraulics.orifice_flow"
_CONTROL_VALVE_CV = "hydraulics.control_valve_cv"
_CHOKED_FLOW_AREA = "hydraulics.choked_flow_area"


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


def pump_power(rho: Q, q: Q, H: Q, eta: float) -> PumpPowerResult:
    """Shaft power a pump must be supplied with.

    ``eta`` is the pump's own efficiency, dimensionless and in ``(0, 1]``. ``H`` is
    the head *delivered*, in metres of the pumped fluid: the pump's internal losses
    are what ``eta`` accounts for, so they do not belong in ``H`` as well.

    Raises:
        OutOfRangeError: if ``rho`` is not positive, if ``q`` or ``H`` is negative,
            or if ``eta`` is outside ``(0, 1]``.

    See :func:`azoth.hydraulics.reference.pump_power`.
    """
    return resolve(_PUMP_POWER)(rho=rho, q=q, H=H, eta=eta)  # type: ignore[no-any-return]


def orifice_flow(d: Q, dP: Q, rho: Q, Cd: float) -> OrificeFlowResult:
    """Volumetric flow through an orifice.

    ``Cd`` is the discharge coefficient and is supplied rather than computed: the
    standard's coefficient equations are fitted expressions built on tables of
    experimental constants, and this library does not reproduce those. Supply the
    coefficient for ``q = Cd * A * sqrt(2*dP/rho)`` as written, with no
    velocity-of-approach factor added - see the reference implementation.

    ``d`` is the bore and is declared in millimetres, as bores are quoted; any length
    is accepted. ``dP`` is a magnitude and may not be negative.

    Raises:
        OutOfRangeError: if ``d`` or ``rho`` is not positive, if ``dP`` is negative,
            or if ``Cd`` is outside ``(0, 1]``.

    See :func:`azoth.hydraulics.reference.orifice_flow`.
    """
    return resolve(_ORIFICE_FLOW)(d=d, dP=dP, rho=rho, Cd=Cd)  # type: ignore[no-any-return]


def control_valve_cv(Cv: float, dP: Q, SG: float) -> ControlValveCvResult:
    """Liquid flow through a control valve.

    ``Cv`` is the valve flow coefficient in the **US** convention - gallons per
    minute of water at one psi - and is supplied rather than looked up: this library
    does not reproduce IEC 60534's coefficient tables. ``Cv`` is not dimensionless
    in any physical sense; see the reference implementation for why the schema has
    to declare it so and what the calculation does about it.

    ``SG`` is the liquid's specific gravity relative to water at 15.6 C.

    Raises:
        OutOfRangeError: if ``Cv`` or ``SG`` is not positive, or if ``dP`` is
            negative.

    See :func:`azoth.hydraulics.reference.control_valve_cv`.
    """
    return resolve(_CONTROL_VALVE_CV)(Cv=Cv, dP=dP, SG=SG)  # type: ignore[no-any-return]


def choked_flow_area(m_dot: Q, P0: Q, rho0: Q, k: float) -> ChokedFlowAreaResult:
    """Throat area required for a choked gas flow.

    This is the isentropic critical-flow relation, which is the physical basis of
    relief valve sizing. It is **not** relief valve sizing to a standard: the
    de-rating coefficients API 520 requires - discharge, back pressure, combination -
    are the caller's to apply, and their values are tabulated in the standard rather
    than reproduced here.

    ``k`` is the isentropic exponent, dimensionless and above 1 for every real gas.
    The flow must be choked for this area to be the right one, and that cannot be
    checked here; see the reference implementation.

    Raises:
        OutOfRangeError: if ``P0`` or ``rho0`` is not positive, if ``m_dot`` is
            negative, or if ``k`` is not greater than 1.

    See :func:`azoth.hydraulics.reference.choked_flow_area`.
    """
    return resolve(_CHOKED_FLOW_AREA)(m_dot=m_dot, P0=P0, rho0=rho0, k=k)  # type: ignore[no-any-return]


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
