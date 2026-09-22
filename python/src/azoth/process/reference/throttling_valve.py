"""``process.throttling_valve`` - the valve's kernel.

Spec: ``specs/models/process/throttling_valve.toml``

The Python twin of ``crates/azoth-process/src/models/throttling_valve.rs`` over
``crates/azoth-process/src/kernels/throttling_valve.rs``, written to mirror them.

# The kernel needed no change to become a port

``ThrottlingValve.run`` flashes a ``PHflash`` at the inlet's total enthalpy, and
``Stream::from_ph`` is the same operation, so the outlet's temperature is the one the
mixture reaches at the outlet pressure. And ``acceptNegativeDP`` defaults to ``true``, so
a stated pressure *above* the inlet is honoured rather than clamped - which is what the
capture's third row shows, 30 bara in and 40 bara out at the same enthalpy. The palette
entry was wrong about its parameters, not about its arithmetic.
"""

from __future__ import annotations

from azoth.core.range import apply_checks, checks_for
from azoth.core.result import ThrottlingValveResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning
from azoth.eos import components as _components
from azoth.eos.reference.ph_flash import enthalpy_at
from azoth.eos.reference.ph_flash import ph_flash as ph_flash_solve


def throttling_valve(
    components: list[str],
    inlet_n: Q,
    inlet_z: list[float],
    inlet_p: Q,
    inlet_t: Q,
    outlet_pressure: Q,
) -> ThrottlingValveResult:
    """Drop a stream to a lower pressure without heat or work.

    Args:
        components: the fluid's substances, by name.
        inlet_n: the inlet's molar flow.
        inlet_z: the inlet composition.
        inlet_p: the inlet pressure.
        inlet_t: the inlet temperature.
        outlet_pressure: the pressure the valve drops the stream to.

    Returns:
        The outlet's record: the inlet's flow, composition and molar enthalpy, the stated
        pressure, and the temperature the mixture reaches at it.

    Raises:
        OutOfRangeError: for a temperature or pressure the spec bounds.

    Example:
        >>> import azoth
        >>> q = azoth.ureg.Quantity
        >>> r = throttling_valve(
        ...     ["methane", "n-butane"],
        ...     q(1.0, "mol/s"),
        ...     [0.9, 0.1],
        ...     q(30.0, "bar"),
        ...     q(320.0, "K"),
        ...     q(10.0, "bar"),
        ... )
        >>> round(r.outlet_t.to("K").magnitude, 4)
        308.9296
    """
    spec = _spec()
    checks = checks_for(spec)
    warnings: list[Warning] = []

    n = input_to_si(spec, "inlet_n", inlet_n)
    t = input_to_si(spec, "inlet_t", inlet_t)
    p = input_to_si(spec, "inlet_p", inlet_p)

    # Every declared bound's quantity is fed here: a closure that supplies one and not
    # the other produces a `RANGE_CHECK_SKIPPED` warning for a check that could have run,
    # which reads as a defect in the caller's state rather than in this line.
    apply_checks(
        checks.on_input,
        {
            "inlet_t": t,
            "outlet_pressure": input_to_si(spec, "outlet_pressure", outlet_pressure),
        }.get,
        warnings,
    )

    mixture, ideal_gas = _components.mixture_of(components, eos="pr")
    # The inlet's own state, which the drop carries unchanged.
    h_in, _ = enthalpy_at(mixture, ideal_gas, t, p, inlet_z)
    solved = ph_flash_solve(mixture, ideal_gas, outlet_pressure, from_si(h_in, "J/mol"), inlet_z)

    return ThrottlingValveResult(
        outlet_n=from_si(n, "mol/s"),
        outlet_z=tuple(inlet_z),
        outlet_p=outlet_pressure,
        outlet_t=solved.T,
        outlet_h=from_si(h_in, "J/mol"),
        warnings=tuple(warnings),
    )


def _spec() -> dict[str, object]:
    """The generated spec, imported here so the module has no import-time cycle."""
    from azoth._models_gen import model

    return model("process.throttling_valve")
