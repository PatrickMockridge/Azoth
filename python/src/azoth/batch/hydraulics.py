"""Hydraulics, over arrays.

One call, N evaluations, arrays in and arrays out. The scalar functions in
:mod:`azoth.hydraulics` are unchanged and remain what this is checked against: the
batch layer is a loop over them, not a second implementation of them.

# Units, which are the one thing this gives up

Inputs are **plain numbers in the spec's canonical unit** and outputs are **SI base
magnitudes** with a unit map. So ``rho=[998.0]`` is kg/m**3 and ``result.dp`` is an
array of pascals.

The scalar API takes and returns `pint` quantities, and that is the difference worth
knowing about. Converting a sequence of quantities costs one `pint` operation per
element - the same order as the per-element boundary crossing this exists to remove - so
the conversion is one factor applied to a whole array instead. What is lost is the
scalar API's `UnitMismatchError`: passing metres where the spec says millimetres gets a
thousand-fold error with nothing to catch it. See `docs/src/batch.md`.
"""

from __future__ import annotations

from array import array
from collections.abc import Sequence
from dataclasses import dataclass

from azoth.batch._core import run, sequence
from azoth.batch._result import BatchResult
from azoth.core.result import FlowRegime
from azoth.core.warnings import Warning

__all__ = [
    "ChokedFlowAreaBatch",
    "ColebrookBatch",
    "ControlValveCvBatch",
    "DarcyWeisbachBatch",
    "HaalandBatch",
    "OrificeFlowBatch",
    "PumpPowerBatch",
    "ReynoldsNumberBatch",
    "SwameeJainBatch",
    "choked_flow_area",
    "control_valve_cv",
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
_DARCY_WEISBACH = "hydraulics.darcy_weisbach"
_PUMP_POWER = "hydraulics.pump_power"
_ORIFICE_FLOW = "hydraulics.orifice_flow"
_CONTROL_VALVE_CV = "hydraulics.control_valve_cv"
_CHOKED_FLOW_AREA = "hydraulics.choked_flow_area"


def _regimes(labels: Sequence[str | None]) -> tuple[FlowRegime | None, ...]:
    """Enum labels back into the enum, because that is what the scalar API returns.

    `None` survives as `None`: it means the calc could not form a regime - viscosity was
    omitted - which is not the same as any of the three real answers.
    """
    return tuple(FlowRegime(label) if label is not None else None for label in labels)


@dataclass(frozen=True, slots=True, eq=False, repr=False)
class ReynoldsNumberBatch(BatchResult):
    """Result of a batch :func:`azoth.hydraulics.reynolds_number`.

    Every field is a sequence of length ``len(result)``, in input order. See
    :class:`~azoth.batch._result.BatchResult` for the shape and the unit rule.
    """

    #: Reynolds number per element. Dimensionless.
    re: array[float]
    #: Flow regime per element, `None` where it could not be formed.
    regime: tuple[FlowRegime | None, ...]


@dataclass(frozen=True, slots=True, eq=False, repr=False)
class DarcyWeisbachBatch(BatchResult):
    """Result of a batch :func:`azoth.hydraulics.darcy_weisbach`.

    `re` and `regime` are `NaN` and `None` respectively when viscosity was omitted for
    the batch, and every element then carries a `RANGE_CHECK_SKIPPED` warning saying the
    regime went unchecked. `NaN` rather than zero because a Reynolds number of zero is a
    physical claim, and the input that would let this calc make one was withheld.
    """

    #: Pressure drop per element, in pascals.
    dp: array[float]
    #: Reynolds number per element, `NaN` where viscosity was omitted.
    re: array[float]
    #: Flow regime per element, `None` where viscosity was omitted.
    regime: tuple[FlowRegime | None, ...]


@dataclass(frozen=True, slots=True, eq=False, repr=False)
class ColebrookBatch(BatchResult):
    """Result of a batch :func:`azoth.hydraulics.friction_factor_colebrook`.

    The solver report travels with the answer, as it does in the scalar API, because
    ``f`` is the last iterate: without ``converged`` a caller cannot tell a solution from
    a failure to converge.
    """

    #: Darcy friction factor per element. Dimensionless.
    f: array[float]
    #: Iterations performed per element.
    iterations: array[float]
    #: Whether each element's iteration met its tolerance, as ``1.0`` or ``0.0``.
    #: A batch column is an array of numbers, so the flag is one; ``bool(...)`` recovers
    #: it, and the conversion is exact rather than merely close.
    converged: array[float]
    #: Final change between iterates per element.
    residual: array[float]


@dataclass(frozen=True, slots=True, eq=False, repr=False)
class SwameeJainBatch(BatchResult):
    """Result of a batch :func:`azoth.hydraulics.friction_factor_swamee_jain`."""

    #: Darcy friction factor per element. Dimensionless.
    f: array[float]


@dataclass(frozen=True, slots=True, eq=False, repr=False)
class HaalandBatch(BatchResult):
    """Result of a batch :func:`azoth.hydraulics.friction_factor_haaland`."""

    #: Darcy friction factor per element. Dimensionless.
    f: array[float]


@dataclass(frozen=True, slots=True, eq=False, repr=False)
class PumpPowerBatch(BatchResult):
    """Result of a batch :func:`azoth.hydraulics.pump_power`."""

    #: Shaft power per element, in watts.
    power: array[float]


@dataclass(frozen=True, slots=True, eq=False, repr=False)
class OrificeFlowBatch(BatchResult):
    """Result of a batch :func:`azoth.hydraulics.orifice_flow`."""

    #: Volumetric flow rate per element, in m**3/s.
    q: array[float]


@dataclass(frozen=True, slots=True, eq=False, repr=False)
class ControlValveCvBatch(BatchResult):
    """Result of a batch :func:`azoth.hydraulics.control_valve_cv`."""

    #: Volumetric flow rate per element, in m**3/s.
    q: array[float]


@dataclass(frozen=True, slots=True, eq=False, repr=False)
class ChokedFlowAreaBatch(BatchResult):
    """Result of a batch :func:`azoth.hydraulics.choked_flow_area`."""

    #: Throat area per element, in m**2.
    a: array[float]


def _build_colebrook(
    columns: dict[str, object],
    units: dict[str, str],
    warnings: tuple[tuple[Warning, ...], ...],
) -> ColebrookBatch:
    return ColebrookBatch(
        warnings=warnings,
        units=units,
        f=columns["f"],  # type: ignore[arg-type]
        iterations=columns["iterations"],  # type: ignore[arg-type]
        converged=columns["converged"],  # type: ignore[arg-type]
        residual=columns["residual"],  # type: ignore[arg-type]
    )


def _build_swamee_jain(
    columns: dict[str, object],
    units: dict[str, str],
    warnings: tuple[tuple[Warning, ...], ...],
) -> SwameeJainBatch:
    return SwameeJainBatch(warnings=warnings, units=units, f=columns["f"])  # type: ignore[arg-type]


def _build_haaland(
    columns: dict[str, object],
    units: dict[str, str],
    warnings: tuple[tuple[Warning, ...], ...],
) -> HaalandBatch:
    return HaalandBatch(warnings=warnings, units=units, f=columns["f"])  # type: ignore[arg-type]


def _build_pump_power(
    columns: dict[str, object],
    units: dict[str, str],
    warnings: tuple[tuple[Warning, ...], ...],
) -> PumpPowerBatch:
    return PumpPowerBatch(warnings=warnings, units=units, power=columns["power"])  # type: ignore[arg-type]


def _build_orifice(
    columns: dict[str, object],
    units: dict[str, str],
    warnings: tuple[tuple[Warning, ...], ...],
) -> OrificeFlowBatch:
    return OrificeFlowBatch(warnings=warnings, units=units, q=columns["q"])  # type: ignore[arg-type]


def _build_control_valve(
    columns: dict[str, object],
    units: dict[str, str],
    warnings: tuple[tuple[Warning, ...], ...],
) -> ControlValveCvBatch:
    return ControlValveCvBatch(warnings=warnings, units=units, q=columns["q"])  # type: ignore[arg-type]


def _build_choked_flow(
    columns: dict[str, object],
    units: dict[str, str],
    warnings: tuple[tuple[Warning, ...], ...],
) -> ChokedFlowAreaBatch:
    return ChokedFlowAreaBatch(warnings=warnings, units=units, a=columns["a"])  # type: ignore[arg-type]


def _build_reynolds(
    columns: dict[str, object],
    units: dict[str, str],
    warnings: tuple[tuple[Warning, ...], ...],
) -> ReynoldsNumberBatch:
    return ReynoldsNumberBatch(
        warnings=warnings,
        units=units,
        re=columns["re"],  # type: ignore[arg-type]
        regime=_regimes(columns["regime"]),  # type: ignore[arg-type]
    )


def _build_darcy(
    columns: dict[str, object],
    units: dict[str, str],
    warnings: tuple[tuple[Warning, ...], ...],
) -> DarcyWeisbachBatch:
    return DarcyWeisbachBatch(
        warnings=warnings,
        units=units,
        dp=columns["dp"],  # type: ignore[arg-type]
        re=columns["re"],  # type: ignore[arg-type]
        regime=_regimes(columns["regime"]),  # type: ignore[arg-type]
    )


def reynolds_number(
    *,
    rho: Sequence[float],
    v: Sequence[float],
    D: Sequence[float],
    mu: Sequence[float],
) -> ReynoldsNumberBatch:
    """Reynolds number and flow regime, over arrays. Inputs in kg/m**3, m/s, m, Pa*s.

    Keyword-only, and every argument is a sequence: a scalar here is an
    `InvalidInputError` rather than a value to broadcast, because a length-one array
    where a length-N array was meant produces plausible answers rather than an error.

    Raises:
        InvalidInputError: for a missing, unknown or mismatched input, or a scalar where
            a sequence belongs.
        OutOfRangeError: if any element violates a hard bound. The whole call raises;
            no partial result is returned.

    See :func:`azoth.hydraulics.reynolds_number` for the calculation itself.
    """
    result: ReynoldsNumberBatch = run(
        _REYNOLDS_NUMBER,
        {
            "rho": sequence(rho, "rho"),
            "v": sequence(v, "v"),
            "D": sequence(D, "D"),
            "mu": sequence(mu, "mu"),
        },
        _build_reynolds,
    )
    return result


def darcy_weisbach(
    *,
    f: Sequence[float],
    L: Sequence[float],
    D: Sequence[float],
    rho: Sequence[float],
    v: Sequence[float],
    mu: Sequence[float] | None = None,
) -> DarcyWeisbachBatch:
    """Darcy-Weisbach pressure drop, over arrays. Inputs in m, m, kg/m**3, m/s, Pa*s.

    ``mu`` is optional **per batch, not per element**: present for every element or for
    none. A partly-supplied optional input has no meaning here that a caller has asked
    for, and inventing one would need a per-element absent-value encoding. Omit it and
    every element carries a `RANGE_CHECK_SKIPPED` warning.

    Raises:
        InvalidInputError: for a missing, unknown or mismatched input.
        OutOfRangeError: if any element violates a hard bound.

    See :func:`azoth.hydraulics.darcy_weisbach` for the calculation itself.
    """
    arrays: dict[str, object] = {
        "f": sequence(f, "f"),
        "L": sequence(L, "L"),
        "D": sequence(D, "D"),
        "rho": sequence(rho, "rho"),
        "v": sequence(v, "v"),
    }
    if mu is not None:
        arrays["mu"] = sequence(mu, "mu")
    result: DarcyWeisbachBatch = run(_DARCY_WEISBACH, arrays, _build_darcy)
    return result


def friction_factor_colebrook(
    *, re: Sequence[float], relative_roughness: Sequence[float]
) -> ColebrookBatch:
    """Colebrook friction factor, over arrays. Both inputs dimensionless.

    See :func:`azoth.hydraulics.friction_factor_colebrook`. The solver report comes back
    with the answer, so ``converged`` is a column rather than a footnote.
    """
    result: ColebrookBatch = run(
        _COLEBROOK,
        {
            "re": sequence(re, "re"),
            "relative_roughness": sequence(relative_roughness, "relative_roughness"),
        },
        _build_colebrook,
    )
    return result


def friction_factor_swamee_jain(
    *, re: Sequence[float], relative_roughness: Sequence[float]
) -> SwameeJainBatch:
    """Swamee-Jain explicit friction factor, over arrays. Both inputs dimensionless.

    See :func:`azoth.hydraulics.friction_factor_swamee_jain`.
    """
    result: SwameeJainBatch = run(
        _SWAMEE_JAIN,
        {
            "re": sequence(re, "re"),
            "relative_roughness": sequence(relative_roughness, "relative_roughness"),
        },
        _build_swamee_jain,
    )
    return result


def friction_factor_haaland(
    *, re: Sequence[float], relative_roughness: Sequence[float]
) -> HaalandBatch:
    """Haaland explicit friction factor, over arrays. Both inputs dimensionless.

    See :func:`azoth.hydraulics.friction_factor_haaland`.
    """
    result: HaalandBatch = run(
        _HAALAND,
        {
            "re": sequence(re, "re"),
            "relative_roughness": sequence(relative_roughness, "relative_roughness"),
        },
        _build_haaland,
    )
    return result


def pump_power(
    *,
    rho: Sequence[float],
    q: Sequence[float],
    H: Sequence[float],
    eta: Sequence[float],
) -> PumpPowerBatch:
    """Pump shaft power, over arrays. Inputs in kg/m**3, m**3/s, m, and dimensionless.

    See :func:`azoth.hydraulics.pump_power`. ``power`` is the power the pump must be
    *supplied* with, which is why ``eta`` divides rather than multiplies.
    """
    result: PumpPowerBatch = run(
        _PUMP_POWER,
        {
            "rho": sequence(rho, "rho"),
            "q": sequence(q, "q"),
            "H": sequence(H, "H"),
            "eta": sequence(eta, "eta"),
        },
        _build_pump_power,
    )
    return result


def orifice_flow(
    *,
    d: Sequence[float],
    dP: Sequence[float],
    rho: Sequence[float],
    Cd: Sequence[float],
) -> OrificeFlowBatch:
    """Orifice flow rate, over arrays. ``d`` in **millimetres**, the rest m**3/s.

    ``d`` is millimetres because that is what the spec declares, and this is the one
    input in the tree where the declared unit is not the SI base - so it is the one place
    a caller who assumes SI gets a quiet factor of a thousand. Stated here and on the
    batch page rather than left to the spec.

    See :func:`azoth.hydraulics.orifice_flow`.
    """
    result: OrificeFlowBatch = run(
        _ORIFICE_FLOW,
        {
            "d": sequence(d, "d"),
            "dP": sequence(dP, "dP"),
            "rho": sequence(rho, "rho"),
            "Cd": sequence(Cd, "Cd"),
        },
        _build_orifice,
    )
    return result


def control_valve_cv(
    *, Cv: Sequence[float], dP: Sequence[float], SG: Sequence[float]
) -> ControlValveCvBatch:
    """Control-valve flow rate, over arrays. ``dP`` in Pa, the rest dimensionless.

    See :func:`azoth.hydraulics.control_valve_cv`, which is where the ``Cv``/``Kv`` unit
    convention is set out. ``dP`` is pascals here, not the psi the convention is quoted
    in: the convention fixes the *conversion*, and the spec's declared unit fixes the
    argument.
    """
    result: ControlValveCvBatch = run(
        _CONTROL_VALVE_CV,
        {
            "Cv": sequence(Cv, "Cv"),
            "dP": sequence(dP, "dP"),
            "SG": sequence(SG, "SG"),
        },
        _build_control_valve,
    )
    return result


def choked_flow_area(
    *,
    m_dot: Sequence[float],
    P0: Sequence[float],
    rho0: Sequence[float],
    k: Sequence[float],
) -> ChokedFlowAreaBatch:
    """Choked-flow throat area, over arrays. Inputs in kg/s, Pa, kg/m**3, dimensionless.

    See :func:`azoth.hydraulics.choked_flow_area`.
    """
    result: ChokedFlowAreaBatch = run(
        _CHOKED_FLOW_AREA,
        {
            "m_dot": sequence(m_dot, "m_dot"),
            "P0": sequence(P0, "P0"),
            "rho0": sequence(rho0, "rho0"),
            "k": sequence(k, "k"),
        },
        _build_choked_flow,
    )
    return result
