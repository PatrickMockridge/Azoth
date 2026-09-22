"""``process.pump`` - the pump's kernel.

Spec: ``specs/models/process/pump.toml``

The Python twin of ``crates/azoth-process/src/models/pump.rs`` over
``crates/azoth-process/src/kernels/pump.rs``, written to mirror them.

# The work is an isentropic head, not ``v (P_out - P_in)``

NeqSim's ``Pump.run`` clones the fluid to the outlet pressure and runs a ``PSflash`` at the
inlet entropy, because ``calculateAsCompressor`` defaults to ``true`` and selects that
branch. So the work is ``(H(P_out, s_in) - H(P_in)) / eta`` - the enthalpy a real fluid
needs along an isentrope - and the incompressible form is its linearisation. Measured on the
captured fluid the two differ by ``1.2e-3``, which is why this is the isentropic route.

# The stream record

An inlet carries ``n, z, P, T`` and an outlet carries ``n, z, P, T, h``. ``h`` is a state
function of ``(T, P, z)`` - the same flash decides it - so accepting one as well would let a
case hand over a state that does not exist.
"""

from __future__ import annotations

from azoth.core.errors import InvalidInputError
from azoth.core.range import apply_checks, checks_for
from azoth.core.result import PumpResult
from azoth.core.units import Q, from_si, input_to_si
from azoth.core.warnings import Warning
from azoth.eos import components as _components

# **The modules' functions, imported by name.** `azoth.eos.reference` re-exports `ph_flash`
# and `ps_flash` as *functions*, so reaching the modules through the package gives the
# function where a module is meant. The `from ... import` form resolves the module path
# itself and is unambiguous.
from azoth.eos.reference.ph_flash import enthalpy_at
from azoth.eos.reference.ph_flash import ph_flash as ph_flash_solve
from azoth.eos.reference.ps_flash import entropy_at
from azoth.eos.reference.ps_flash import ps_flash as ps_flash_solve


def pump(
    components: list[str],
    inlet_n: Q,
    inlet_z: list[float],
    inlet_p: Q,
    inlet_t: Q,
    outlet_pressure: Q,
    isentropic_efficiency: float,
) -> PumpResult:
    """Raise a stream's pressure, adding the pump's work as enthalpy.

    Args:
        components: the fluid's substances, by name.
        inlet_n: the inlet's molar flow.
        inlet_z: the inlet composition.
        inlet_p: the inlet pressure.
        inlet_t: the inlet temperature.
        outlet_pressure: the pressure the pump raises the stream to.
        isentropic_efficiency: the pump's isentropic efficiency, in ``(0, 1]``.

    Returns:
        The outlet's record: its flow and composition (the inlet's), its pressure, and the
        temperature and molar enthalpy the isentropic head implies.

    Raises:
        InvalidInputError: where the shapes disagree or the efficiency is outside ``(0, 1]``.
        OutOfRangeError: for a temperature the spec bounds.

    Example:
        >>> import azoth
        >>> q = azoth.ureg.Quantity
        >>> r = pump(
        ...     ["n-butane"],
        ...     1.0,
        ...     [1.0],
        ...     q(5.0, "bar"),
        ...     q(250.0, "K"),
        ...     q(20.0, "bar"),
        ...     0.75,
        ... )
        >>> round(r.outlet_t.to("K").magnitude, 4)
        250.7621
    """
    spec = _spec()
    checks = checks_for(spec)
    warnings: list[Warning] = []

    n = input_to_si(spec, "inlet_n", inlet_n)
    t = input_to_si(spec, "inlet_t", inlet_t)
    p = input_to_si(spec, "inlet_p", inlet_p)
    p_out = input_to_si(spec, "outlet_pressure", outlet_pressure)

    apply_checks(
        checks.on_input,
        {"inlet_t": t, "isentropic_efficiency": isentropic_efficiency}.get,
        warnings,
    )

    if not 0.0 < isentropic_efficiency <= 1.0:
        raise InvalidInputError(
            "isentropic_efficiency",
            f"the pump's isentropic efficiency must be in (0, 1], and {isentropic_efficiency} "
            f"is not",
        )

    mixture, ideal_gas = _components.mixture_of(components, eos="pr")

    # The inlet's own state, which is where the entropy and the enthalpy come from.
    h_in, _ = enthalpy_at(mixture, ideal_gas, t, p, inlet_z)
    entropy, _ = entropy_at(mixture, ideal_gas, t, p, inlet_z)

    isentropic = ps_flash_solve(
        mixture, ideal_gas, outlet_pressure, from_si(entropy, "J/(mol*K)"), inlet_z
    )
    h_isentropic, _ = enthalpy_at(
        mixture, ideal_gas, isentropic.T.to("K").magnitude, p_out, inlet_z
    )

    h_out = h_in + (h_isentropic - h_in) / isentropic_efficiency
    outlet_t = ph_flash_solve(mixture, ideal_gas, outlet_pressure, from_si(h_out, "J/mol"), inlet_z)

    return PumpResult(
        outlet_n=n,
        outlet_z=tuple(inlet_z),
        outlet_p=outlet_pressure,
        outlet_t=outlet_t.T,
        outlet_h=from_si(h_out, "J/mol"),
        warnings=tuple(warnings),
    )


def _spec() -> dict[str, object]:
    """The generated spec, imported here so the module has no import-time cycle."""
    from azoth._models_gen import model

    return model("process.pump")
