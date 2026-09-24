"""``process.filter`` - a fixed pressure drop, flashed at the feed's temperature.

Spec: ``specs/models/process/filter.toml``

The Python twin of ``crates/azoth-process/src/models/filter.rs`` over
``crates/azoth-process/src/kernels/filter.rs``, written to mirror them.

# It holds the temperature, which is what separates it from a throttling valve

``Filter.run`` reduces the pressure and runs a ``TPflash``, so the outlet is at the feed's
temperature and its enthalpy moves with the pressure. A valve is a ``PHflash`` and holds the
enthalpy instead. Measured on the fluid these cases share, one bar of drop moves the enthalpy
``24.8`` J/mol - the same pressure dependence the no-specification branch of
``process.heater`` prices at ``49.6`` over two.

# A drop past the inlet pressure is clamped rather than refused

``run`` applies ``min(max(0, dP), max(0, P_in - 1e-6 bar))`` and logs a warning, so the
outlet lands a millionth of a bar above vacuum. ``applied_drop`` carries what was applied and
the model warns when the two differ, because a result that carried the clamped value silently
would be indistinguishable from one that was asked for it.
"""

from __future__ import annotations

from dataclasses import dataclass

from azoth.core.errors import InvalidInputError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import FilterResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning, WarningCode
from azoth.eos import components as _components
from azoth.eos.reference.ph_flash import enthalpy_at


def filter(
    components: list[str],
    inlet_n: Q,
    inlet_z: list[float],
    inlet_p: Q,
    inlet_t: Q,
    pressure_drop: Q,
) -> FilterResult:
    """Drop a stream's pressure by a fixed amount at a constant temperature.

    Args:
        components: the fluid's substances, by name.
        inlet_n: the inlet's molar flow.
        inlet_z: the inlet composition.
        inlet_p: the inlet pressure.
        inlet_t: the inlet temperature, which is also the outlet's.
        pressure_drop: the pressure loss across the filter.

    Returns:
        The outlet's record and the drop that was applied.

    Raises:
        InvalidInputError: where the shapes disagree.

    Example:
        >>> import azoth
        >>> q = azoth.ureg.Quantity
        >>> r = filter(
        ...     ["methane", "n-butane"],
        ...     1.0,
        ...     q([0.9, 0.1]),
        ...     q(30.0, "bar"),
        ...     q(320.0, "K"),
        ...     q(1.0, "bar"),
        ... )
        >>> round(r.outlet_t.to("K").magnitude, 1)
        320.0
    """
    spec = _spec()
    checks = checks_for(spec)
    warnings: list[Warning] = []

    n = input_to_si(spec, "inlet_n", inlet_n)
    t = input_to_si(spec, "inlet_t", inlet_t)
    p = input_to_si(spec, "inlet_p", inlet_p)
    drop = input_to_si(spec, "pressure_drop", pressure_drop)

    apply_checks(checks.on_input, {"inlet_t": t, "pressure_drop": drop}.get, warnings)

    states = _route(components, n, inlet_z, p, t, drop)

    if states.applied_drop != drop:
        warnings.append(
            Warning(
                code=WarningCode.OUT_OF_VALID_RANGE,
                field="pressure_drop",
                message=(
                    f"a drop of {drop} Pa was asked for and {states.applied_drop} Pa was "
                    f"applied: the outlet is held a millionth of a bar above vacuum, which "
                    f"is the clamp `Filter.run` applies"
                ),
            )
        )

    return FilterResult(
        outlet_n=from_si(n, "mol/s"),
        outlet_z=tuple(inlet_z),
        outlet_p=from_si(states.p_out, "Pa"),
        outlet_t=from_si(t, "K"),
        outlet_h=from_si(states.h_out, "J/mol"),
        applied_drop=from_si(states.applied_drop, "Pa"),
        warnings=tuple(warnings),
    )


@dataclass(frozen=True, slots=True)
class FilterStates:
    """The filter's intermediates, in the order it computes them, in SI."""

    #: The drop applied, which is the requested one unless the clamp moved it.
    applied_drop: float
    #: The outlet pressure.
    p_out: float
    #: The outlet molar enthalpy, which is not the inlet's: the drop is isothermal.
    h_out: float


def _route(
    components: list[str],
    inlet_n: float,
    inlet_z: list[float],
    inlet_p: float,
    inlet_t: float,
    pressure_drop: float,
) -> FilterStates:
    """The filter's intermediates, in SI, in the order it computes them."""
    del inlet_n
    mixture, ideal_gas = _components.mixture_of(components, eos="pr")

    if not inlet_p > 0.0:
        raise InvalidInputError(
            "inlet_p", f"a filter needs an inlet pressure to drop from, and {inlet_p} is not one"
        )
    # `Math.max(0.0, inletPressure - 1.0e-6)` in the class's bara, which is 0.1 Pa.
    ceiling = max(0.0, inlet_p - 0.1)
    applied = min(max(0.0, pressure_drop), ceiling)
    p_out = inlet_p - applied

    h_out, _ = enthalpy_at(mixture, ideal_gas, inlet_t, p_out, inlet_z)

    return FilterStates(applied_drop=applied, p_out=p_out, h_out=h_out)


def _spec() -> dict[str, object]:
    """The generated spec, imported here so the module has no import-time cycle."""
    from azoth._models_gen import model

    return model("process.filter")
