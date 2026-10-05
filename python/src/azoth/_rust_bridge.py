"""The Rust backend, presented with the same interface as the reference.

Each function here takes and returns exactly what
:mod:`azoth.hydraulics.reference` does, so the dispatcher can swap one for the
other and a caller cannot tell which answered. Two jobs sit between the two
representations:

* **Units in.** Quantities are converted to SI magnitudes before crossing into
  Rust. That conversion lives here, once, rather than in the extension - one
  conversion site used by the reference and the Rust path alike cannot disagree
  with itself about what a number is in.
* **Results out.** The extension's transport objects become the same frozen
  dataclasses the reference returns, so ``result.dp`` is a pint quantity either
  way and ``isinstance`` checks behave identically.

If this file is missing or out of step with the extension, the dispatcher raises
rather than falling back - a silent fallback would mean the cross-language
guarantee was being verified by nothing.
"""

from __future__ import annotations

from collections.abc import Callable, Mapping
from typing import Any

from azoth import _core, _models_gen

# **The generated half.** Every id whose signature is the public wrapper's and whose every
# boundary parameter is an input the spec declares, a field of the `params` record, or a boundary
# object's derived expansion is emitted by `tools/gen_python_bridge.py`; what is left in this file
# is the ids that generator refuses, each with its reason in `--survey`. `resolve` derives the
# id's function name and looks it up in this module's globals, which the star import fills.
from azoth._rust_bridge_gen import *  # noqa: F403

# The helpers the generated adapters and the hand-written ones both call. Imported rather than
# defined here because `azoth._rust_bridge_gen` imports them too, and a module cannot import the
# file that re-exports it. `from ... import _warnings` keeps the name reachable as
# `azoth._rust_bridge._warnings`, which is where `azoth.batch._core` reads it.
from azoth.core.bridge_support import _association_spec as _association_spec
from azoth.core.bridge_support import _si as _si
from azoth.core.bridge_support import _warnings as _warnings
from azoth.core.result import (
    HydrogenPhaseResult,
    PureSaturationResult,
)
from azoth.core.units import Q, from_si, input_to_si


def pure_saturation(Tc: Q, Pc: Q, omega: float, T: Q) -> PureSaturationResult:
    """The saturation pressure of a pure component, computed in Rust.

    A model rather than a calculation - its spec fixes a procedure and the Rust side
    composes the same kernels this bridge's scalar functions call, so the two
    implementations run the same search over the same arithmetic.
    """
    spec = _models_gen.model("eos.pure_saturation")
    result = _core.pure_saturation(
        # Not `input_to_si`: that reads the unit out of the spec's declaration, and the
        # spec no longer declares `Tc` or `Pc` - a component's constants come from the
        # databank by name, and a caller holding them holds quantities already. `T` is
        # still declared and still goes through the spec.
        Tc.to("K").magnitude,
        Pc.to("Pa").magnitude,
        # A plain float: `omega` is genuinely dimensionless, so it crosses as
        # the number the caller used - the same rule the other eos calcs follow.
        omega,
        input_to_si(spec, "T", T),
    )
    return PureSaturationResult(
        p_sat=from_si(result.p_sat.magnitude_si, result.p_sat.unit),
        ln_phi=result.ln_phi,
        iterations=result.iterations,
        residual=result.residual,
        warnings=_warnings(result.warnings),
    )


#: The same table for *models*, kept separate because the calc table is asserted to
#: be exactly the calc registry - ``test_registration_completeness`` compares it by
#: equality in both directions, so a model id in it would break that contract rather
#: than extend it. Models have their own id list for the same reason, and
#: ``test_model_contract.py`` holds this table to it.


def _boundary_result(raw: Any, result_type: Any, *, liquid_first: bool) -> Any:
    """Shared unpacking for the two phase-boundary models.

    The two differ only in whether the held phase is the liquid, so the `z_liquid`
    and `z_vapour` fields land in the opposite order. Doing that here rather than
    twice is the same argument the Rust side makes for one shared iteration.
    """
    held, incipient = (raw.z_liquid, raw.z_vapour) if liquid_first else (raw.z_vapour, raw.z_liquid)
    return result_type(
        pressure=from_si(raw.pressure.magnitude_si, raw.pressure.unit),
        incipient=tuple(raw.incipient),
        k=tuple(raw.k),
        z_liquid=held if liquid_first else incipient,
        z_vapour=incipient if liquid_first else held,
        min_t_over_tc=raw.min_t_over_tc,
        iterations=raw.iterations,
        residual=raw.residual,
        warnings=_warnings(raw.warnings),
    )


def _boundary_temperature_result(raw: Any, result_type: Any) -> Any:
    """Unpack a temperature boundary, whose scalar is a temperature rather than a
    pressure. The `z_liquid` and `z_vapour` fields already name the phases they say
    they do, so there is no swap here."""
    return result_type(
        temperature=from_si(raw.temperature.magnitude_si, raw.temperature.unit),
        incipient=tuple(raw.incipient),
        k=tuple(raw.k),
        z_liquid=raw.z_liquid,
        z_vapour=raw.z_vapour,
        min_t_over_tc=raw.min_t_over_tc,
        iterations=raw.iterations,
        residual=raw.residual,
        warnings=_warnings(raw.warnings),
    )


def hydrogen_phase(
    T: Q, P: Q, hydrogen_type: str = "normal", compressed_phase: str = "vapour"
) -> HydrogenPhaseResult:
    """The Leachman hydrogen phase state, computed in Rust."""
    spec = _models_gen.model("eos.hydrogen_phase")
    result = _core.hydrogen_phase(
        input_to_si(spec, "T", T), input_to_si(spec, "P", P), hydrogen_type, compressed_phase
    )
    return HydrogenPhaseResult(
        z_factor=result.z_factor,
        u=from_si(result.u.magnitude_si, result.u.unit),
        h=from_si(result.h.magnitude_si, result.h.unit),
        s=from_si(result.s.magnitude_si, result.s.unit),
        cv=from_si(result.cv.magnitude_si, result.cv.unit),
        cp=from_si(result.cp.magnitude_si, result.cp.unit),
        g=from_si(result.g.magnitude_si, result.g.unit),
        warnings=_warnings(result.warnings),
    )


def overlay_from(card: Any) -> Any:
    """A keycard as the extension reads it.

    The card's ``pint`` quantities become SI magnitudes here, the way every other
    dimensioned value crosses - the extension takes numbers, and one conversion site is
    what keeps the two languages working in the same units.

    **Nothing is resolved here.** This hands the card over as data; what a name resolves
    to is `azoth_eos::databank`'s answer, and `python/tests/test_data_agreement.py`
    compares it against `azoth.eos.components`' answer for the same card. Two
    implementations of one merge rule, held to each other the way the two kernels are.
    """
    components: list[_core.ComponentArguments] = [
        {
            "name": name,
            "tc": _si_of(stated.parameters.get("Tc")),
            "pc": _si_of(stated.parameters.get("Pc")),
            "omega": _plain(stated.parameters.get("omega")),
            "ion": stated.is_ion,
            "ionic_charge": _plain(stated.parameters.get("ionic_charge")),
            # Metres, which is the card's unit and the merge's input. The crossing to the
            # databank's ångström happens in one place, on the Rust side, so there is no
            # second conversion here to disagree with it.
            "deshmukh_mather_diameter": _si_of(stated.parameters.get("deshmukh_mather_diameter")),
            "dielectric": _dielectric(stated.parameters),
        }
        for name, stated in sorted(card.components.items())
    ]
    kij = [(a, b, value) for (a, b), value in sorted(card.kij.items())]
    return _core.overlay(components, kij)


def _si_of(value: Q | None) -> float | None:
    """A dimensioned card value as an SI magnitude, or ``None`` if the card omits it."""
    return None if value is None else float(value.to_base_units().magnitude)


def _dielectric(parameters: Mapping[str, Q]) -> list[float] | None:
    """The five coefficients the card states, each in its own unit, or ``None``."""
    names = ("dielectric_1", "dielectric_2", "dielectric_3", "dielectric_4", "dielectric_5")
    if not all(name in parameters for name in names):
        return None
    return [
        float(parameters["dielectric_1"].to("dimensionless").magnitude),
        float(parameters["dielectric_2"].to("K").magnitude),
        float(parameters["dielectric_3"].to("1/K").magnitude),
        float(parameters["dielectric_4"].to("1/K**2").magnitude),
        float(parameters["dielectric_5"].to("1/K**3").magnitude),
    ]


def _plain(value: Q | None) -> float | None:
    """A dimensionless card value as a bare float, or ``None`` if the card omits it."""
    return None if value is None else float(value.to("dimensionless").magnitude)


def resolve(calc_id: str) -> Callable[..., Any]:
    """The bridge function for a calc id, derived from the id rather than listed.

    A calculation's id is its address here exactly as it is on the reference side:
    `hydraulics.orifice_flow` names `orifice_flow` in this module. So there is no
    table, and adding a calculation adds no entry.

    The table this replaces was hand-maintained and asserted equal to the registry by
    a test - which is a list kept in step with a list, so the test could only ever
    tell you that you had forgotten something. Deriving it means there is nothing to
    forget, and a calc whose bridge function is missing fails the same completeness
    test one step earlier.

    Raises:
        KeyError: if the id is not in the registry at all. That is a programming
            error rather than a runtime condition - the reference and the extension
            are supposed to cover the same calcs, and a test asserts they do.
    """
    from azoth._models_gen import MODELS
    from azoth._registry_gen import CALCS

    if calc_id not in {entry["id"] for entry in [*CALCS, *MODELS]}:
        raise KeyError(
            f"no Rust implementation for {calc_id!r}: it is not in the registry. "
            f"The extension covers {len(CALCS)} calc(s) and {len(MODELS)} model(s)."
        ) from None

    resolved: Callable[..., Any] = globals()[calc_id.rpartition(".")[2]]
    return resolved


def _optional_si(quantity: object) -> Q | None:
    """A quantity the run may not have reached, as the dataclass spells it.

    The Rust transport writes `None` for a number a run did not reach - an ejector handed no
    flow, an exchanger's rating where one outlet was pinned - so the bridge's job here is to
    keep that absence rather than to invent a zero for it.
    """
    if quantity is None:
        return None
    magnitude_si = quantity.magnitude_si  # type: ignore[attr-defined]
    unit = quantity.unit  # type: ignore[attr-defined]
    return from_si(magnitude_si, unit)
